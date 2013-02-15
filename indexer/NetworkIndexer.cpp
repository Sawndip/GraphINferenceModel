//
// BuildNetworkQuery
//
// Generates a network of SNOMED concepts (according to SNOMED concept relationships)
// from a Indri index.
//
// Created: 10-05-2012 by Bevan.Koopman@csiro.au

#include <math.h>
#include "glog/logging.h"

// lemon
#include <lemon/list_graph.h>
#include <lemon/lgf_writer.h>

// lemur
#include "common_headers.hpp"
#include "IndexManager.hpp"

using namespace lemon;
using namespace std;
using namespace lemur::api;
using namespace indri::index;

// SQLite access to concept, concept relationships
#include "Utils.cpp" 
#include "ConceptRelationships.cpp" 

// typedefs
typedef ListDigraph::Node Node;
typedef ListDigraph::Arc Arc;
typedef string CUI;
typedef int nodeid;

// Lemon graph data structures
lemur::api::Index *idx;
ListDigraph *g;
ListDigraph::NodeMap<CUI> *nodeToCUIMap;
ListDigraph::NodeMap<float> *nodeToTopoScoreMap;
ListDigraph::ArcMap<double> *simArcs;
ListDigraph::ArcMap<int> *reltypeArcs;
map<CUI, nodeid> cuiToNodeIdMap;

// Similarity measures
#include "doccosine.cpp" // calcaulate semantic similarity between concepts (as doccosine angle)
#include "cer.cpp" // CER semantic similarity calcaulatation
#include "pmi.cpp"

///////////////////////////////////////////////////////////////////////////////
// Diffusion factor functions
///////////////////////////////////////////////////////////////////////////////

map<DOCID_T, float> cer_background_model;

float idf(CUI cui)
{
    TERMID_T termId = idx->term(cui);
    int termDocCount = termId > 0 ? idx->docCount() > 0 : 1;    
    return log(idx->docCount() / (double) termDocCount); 
}

// doccoise between two term document vectors
float doccosine(TERMID_T t1, TERMID_T t2, string reltype) 
{
    return cosine(t1, t2);
}

// CER language model comparison
float cer(TERMID_T t1, TERMID_T t2, string reltype) 
{
    float lambda = 0.7;
    if(cer_background_model.size() == 0)
        cer_background_model = background_model();
    
    map<DOCID_T, float> theta1 = conceptLM(t1, lambda, cer_background_model);
    map<DOCID_T, float> theta2 = conceptLM(t2, lambda, cer_background_model);
    
    return dCER(theta1, cer_background_model, theta2);
}

// pmi similarity measure
float pmi(TERMID_T t1, TERMID_T t2, string reltype) 
{
    float thePmi = PositivePMI(t1, t2);
    return thePmi;
}
 
// edge weight based on semantic similarity
float semantic_similarity(Node u, Node v, string reltype) 
{
    float similarity = 0.0;
    if(idx->term((*nodeToCUIMap)[u]) > 0 && idx->term((*nodeToCUIMap)[v]) > 0)
    {
        similarity = doccosine(idx->term((*nodeToCUIMap)[u]), idx->term((*nodeToCUIMap)[v]), reltype);
    }
    return similarity;
}

// edge weight based on relationship type
float reltype_weight_umls(Node source, Node target, string reltype) {
    float weight = 0.5;
    if (reltype == "RB") {
        weight = 0.7;
    } else if (reltype == "RN") {
        weight = 0.3;
    }
    return weight;
}

// edge weight based on relationship type
float reltype_weight_snomed(Node source, Node target, string reltype) {
    return snomedRelationshipTypeWeight(reltype);
}

float reltype_weight(Node source, Node target, string reltype) {
    float the_reltype_weight;
    if((*nodeToCUIMap)[source].substr(0,1) == "C") // UMLS
    {
        the_reltype_weight = reltype_weight_umls(source, target, reltype);
    }
    else // SNOMED
    {
        the_reltype_weight = reltype_weight_snomed(source, target, reltype);
        // CURRENTLY APPLIED AT RETRIEVAL TIME
        //the_reltype_weight *= semanticType_weight_snomed((*nodeToCUIMap)[source], (*nodeToCUIMap)[target]);
    }
    VLOG(2) << "the_reltype_weight: " << reltype << "=" << the_reltype_weight;
    return the_reltype_weight;
}



// general diffusion function that calcaulate the weight between two nodes
double diffusionFactor(Node u, Node v, string reltype)
{
    float alpha = 0.5;
    return alpha * reltype_weight(u,v,reltype) + (1-alpha) * semantic_similarity(u,v,reltype);
    
//    stringstream ss;
//    ss << semantic_similarity(u,v,reltype) << " " << reltype;
//    VLOG(2) << "diffusionFactor=" << ss.str();
//    return ss.str();
}

///////////////////////////////////////////////////////////////////////////////
// Graph creation and serialisation functions
///////////////////////////////////////////////////////////////////////////////

// print progress bar % for jobs
void progress(string msg, int count, int total, int increments)
{
    if(total > increments && count % (total/100) == 0) {
        cout << msg << ": " << (int)(((float)count/total)*100) << "%\r";
        cout.flush();
    }
}

// gets all the related nodes and creates edges between related nodes based on relationship
// defined in DB
// Note: Recursive if new, non-index nodes are found
void addEdgesToRelatedNodes(Node &node, bool inEdge)
{
    VLOG(2) << "looking for edges for " << (*nodeToCUIMap)[node];
    
    // for each related concept (from UMLS)
    vector< vector<string> > relatedConcepts = findRelatedConceptsSQL((*nodeToCUIMap)[node], inEdge);
    VLOG(2) << "found " << relatedConcepts.size() << " edges ";
    for (vector< vector<string> >::iterator relIter = relatedConcepts.begin();
         relIter != relatedConcepts.end(); relIter++)
    {
        string relatedCUI = relIter->at(0);
        string reltype = relIter->at(1);
        
        VLOG(2) << "add edge " << (*nodeToCUIMap)[node] << " -> " << relatedCUI << endl;
        
        std::transform(relatedCUI.begin(), relatedCUI.end(), relatedCUI.begin(), ::tolower);
        
        ListDigraph::Node relatedNode;
        
        // check whether there is aleady a node for the related concept, if not create one
        if(cuiToNodeIdMap.count(relatedCUI) == 0)
        {
            VLOG(2) << "\t" << "creating new node for: " << relatedCUI << endl;
            relatedNode = g->addNode();
            cuiToNodeIdMap[relatedCUI] = g->id(relatedNode);
            (*nodeToCUIMap)[relatedNode] = relatedCUI;
            (*nodeToTopoScoreMap)[relatedNode] = idf((*nodeToCUIMap)[relatedNode]);
        } else
        {
            VLOG(2) << "\t" << "getting related node: " << relatedCUI << endl;
            relatedNode = g->nodeFromId(cuiToNodeIdMap[relatedCUI]);
        }
        
        VLOG(2) << "calc df for " << (*nodeToCUIMap)[node] << " -> " << (*nodeToCUIMap)[relatedNode] << ", " << reltype;
        // calculate the weight between concept and assing to arc weight
        double sim_weight = semantic_similarity(node,relatedNode,reltype);
        VLOG(2) << sim_weight << endl;

        Arc arc;
        if(inEdge) {
            arc = g->addArc(relatedNode, node);
        } else {
            arc = g->addArc(node, relatedNode);
        }

        (*simArcs)[arc] = sim_weight;
        (*reltypeArcs)[arc] = atoi(reltype.c_str());
    }
    
    VLOG(2) << "finished adding edges to " << (*nodeToCUIMap)[node] << endl;
}

// create a node for each term in the index
void createGraphNodesIndex()
{
    // cycle through each unique concept in the index
    int total = idx->termCountUnique();
    int  count = 0;
    for (int termId = 1; termId <= total; termId++)
    {
        cout << "Building network nodes:\t" << (int)((++count / (float)total)*100) << "%\r" << flush;
        
        //cout << termId << " (" << idx->term(termId) << ") " << idx->docCount(termId) << " docs" << endl ;
        
        // create the concept node
        ListDigraph::Node node = g->addNode();
        cuiToNodeIdMap[idx->term(termId)] = g->id(node);
        (*nodeToCUIMap)[node] = idx->term(termId);
        (*nodeToTopoScoreMap)[node] = idf((*nodeToCUIMap)[node]);

         
    }

    cout << endl << flush;
}

// create a node for each entry in the DB of related nodes
void createGraphNodesAll()
{
    int count = 0;
    vector<string> concepts = getAllConcepts();
    int total = concepts.size();

    for(vector<string>::iterator it = concepts.begin(); it != concepts.end(); it++)
    {
        cout << "Building network nodes:\t" << (int)((++count / (float)total)*100) << "%\r" << flush;
        
        ListDigraph::Node node = g->addNode();
        cuiToNodeIdMap[*it] = g->id(node);
        (*nodeToCUIMap)[node] = *it;
        (*nodeToTopoScoreMap)[node] = 1.0;
        
        VLOG(2) << count << " / " << total << endl;
    }
    cout << endl << flush;
}

// for each node create an edge to its other node based on terminology (UMLS or SNOMED) relationships
void createGraphEdges()
{
//    Node n = g->nodeFromId(cuiToNodeIdMap["116680003"]);
//    addEdgesToRelatedNodes(n);
    
    int nodeCount = 0;
    int nodeTotal = countNodes(*g);

    
    // for each node in the graph
    for (ListDigraph::NodeIt node(*g); node != INVALID; ++node)
    {
        cout << "Building network edges:\t" << (int)((++nodeCount / (float)nodeTotal)*100) << "%\r" << flush;
        VLOG(2) << "edge " << nodeCount << " / " << nodeTotal << endl;
        addEdgesToRelatedNodes(node, true);  // rn -> n
        //addEdgesToRelatedNodes(node, false); // n -> rn
    }
    cout << endl;

}

// iterates over all nodes and updates a score based on topology, based on
// Blanco & Lioma's Graph-IR paper
//
//                __          score(ni)
// score(n) =     \          -------------
//                /_         |outedge(ni)|
//            ni E inedge(n)
//
void calcTopoScore()
{
    int total = countNodes(*g);
    int count = 0;
    
    int NUM_ITERATIONS = 1;
    for (int i = 0; i < NUM_ITERATIONS; i++) {
        for (ListDigraph::NodeIt node(*g); node != INVALID; ++node)
        {
            float score = (*nodeToTopoScoreMap)[node];
            for (ListDigraph::InArcIt arc(*g, node); arc != INVALID; ++arc)
            {
                Node relatedNode = g->source(arc);
                int outEdges = 0;
                for (ListDigraph::OutArcIt outArc(*g, relatedNode); outArc != INVALID; ++outArc)
                    outEdges++;
                score += (*nodeToTopoScoreMap)[relatedNode] / outEdges;
            }
            VLOG(2) << "topo score for " << (*nodeToCUIMap)[node] << ": " << score << endl;

            (*nodeToTopoScoreMap)[node] = score;
            
            cout << "Calculating topo scores:\t" << (int)((++count / (float)total)*100) << "%\r" << flush;
        }
    }
    cout << endl;
}

// write the lemon graph to file as "graph.lgf" in the lemir index dir
void serialiseGraph(string indexLoc) 
{
    ofstream outfile;
    outfile.open((indexLoc+"/graph-active.lgf").c_str());
    digraphWriter(*g, outfile).nodeMap("Concepts", *nodeToCUIMap).nodeMap("TopoScore", *nodeToTopoScoreMap).arcMap("SimArcWeights", *simArcs).arcMap("RelTypeArc", *reltypeArcs).run();
    outfile.close();
}

///////////////////////////////////////////////////////////////////////////////
// Main functions
///////////////////////////////////////////////////////////////////////////////

// method that does all the work
void NetworkIndexer(string indexLoc, string conceptRelationshipDB)
{
    // load index and graph maps
    idx = lemur::api::IndexManager::openIndex(indexLoc);
    ListDigraph theGraph;
    g = &theGraph;
    ListDigraph::NodeMap<CUI> theNodeToCUIMap(*g);
    nodeToCUIMap = &theNodeToCUIMap;
    ListDigraph::NodeMap<float> theNodeToTopoScoreMap(*g);
    nodeToTopoScoreMap = &theNodeToTopoScoreMap;
    ListDigraph::ArcMap<double> theSimArcs(*g);
    simArcs = &theSimArcs;
    ListDigraph::ArcMap<int> theRelTypeArcs(*g);
    reltypeArcs = &theRelTypeArcs;
    
    prepareDB(conceptRelationshipDB);
        
    createGraphNodesIndex();
    createGraphEdges();
    calcTopoScore();
    
    cout << "Graph complete: " << countNodes(*g) << " nodes; " << countArcs(*g) << " edges." << endl;
    
    cout << "Serialising graph..." << flush;
    serialiseGraph(indexLoc);
    cout << "complete" << endl;
    
    delete idx;
}

// MAIN
int main(int argc, char* argv[])
{
    google::InitGoogleLogging(argv[0]);
    
    LOG_IF(FATAL, argc < 2) << "Usage: retrieval_inference_using_umls <IndriIndex> [<ConceptRelationshipDB>]" << endl;

    string indexLoc = argv[1];
    
    string conceptRelationshipDB = "../snomed_relations/snomed_rel.db";
    if (argc > 2) {
        conceptRelationshipDB = argv[2];
    }
    
    NetworkIndexer(indexLoc, conceptRelationshipDB);
}
