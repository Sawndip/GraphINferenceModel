//
// NetworkRetrieval
//
// Reads a lemon graph and indri index and performs retrieval using network retrieval model
//
//
// Created: 10-05-2012 by Bevan.Koopman@csiro.au

#include <math.h>
#include "glog/logging.h"


// tinyxml2
#include <tinyxml2.cpp>

// lemon
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/maps.h>
#include <lemon/dijkstra.h>
#include <lemon/bfs.h>

// lemur
#include "common_headers.hpp"
#include "IndexManager.hpp"
#include "ScoreAccumulator.hpp"
#include "IndexedReal.hpp"
#include "ResultFile.hpp"
#include "BasicDocStream.hpp"
#include "Param.hpp"

#include "../indexer/Utils.cpp"
// SQLite access to concept, concept relationships
#include "../indexer/ConceptRelationships.cpp" 
#include "report_to_visit_mapper.cpp"
#include "ConceptLookup.cpp"

using namespace lemon;
using namespace std;
using namespace lemur::api;
using namespace lemur::retrieval;
using namespace indri::index;

typedef ListDigraph::Node Node;
typedef ListDigraph::Arc Arc;
typedef int PHRASE_ID;

enum loglevel {ERROR, INFO, DEBUG, TRACE};
loglevel level = DEBUG;

// index and graphs
lemur::api::Index *idx;
ListDigraph *g;
ListDigraph::NodeMap<string> *nodeToCUIMap;
ListDigraph::NodeMap<float> *nodeToTopoScoreMap;
map<string, int> cuiToNodeId;
ListDigraph::ArcMap<double> *simArcs;
ListDigraph::ArcMap<string> *reltypeArcs;
map<Node, map<DOCID_T, double> > diffusionFactors;

#include "../indexer/doccosine.cpp"

// tuning params
float related_dampener = 0.6;
float semantic_type_weight_booster = 2;
double s = 0;
double relweightsweep = -1.0;

// debug maps
map<DOCID_T, double> lvl0Scores;
map<DOCID_T, double> lvl1Scores;

// rerank list
ResultFile lvl0Results(true);
IndexedRealVector *rerankResults;

// stopword list
set<string> stopwordList;

// empty document scores used for language model
map<int, map<TERMID_T, double> > emptyDocScores;

// GraphViz visualisation
ofstream *graphViz;
set<DOCID_T> relevantLvl0Docs;

// keep track how many nodes have contributed to a documents score
map<DOCID_T, int> docNodeCount;

///////////////////////////////////////////////////////////////////////////////
// Retrieve parameters
///////////////////////////////////////////////////////////////////////////////

namespace LocalParameter {
    std::string indexPath; // the index of the documents
    std::string queryFile; // the file of query stream
    std::string resultFile; // the name of the result file
    int resultCount; // the number of top ranked documents to return for each query
    std::string stopwordFile; // file containing list of stopword, one per line
    
    std::string weightScheme; // the weighting scheme
    double LMjm_lambda; // value of lambda in LM Jelinek-Mercer Smoothing
    int LMdirichlet_mu; // the value of mu in LM Dirichlet Smoothing
    double indri_tfidf_k1; // indri's tfidf term weighting param
    double indri_tfidf_b; // indri's tfidf doc length weighting param
    
    int depth; // how many edges to traverse away from the query
    double diffusionMix; // mix of semantic_similary and rel_type
    int diffusionOn; // use a diffusion factor (i.e. on) or just set it to 1 (i.e. off)
    double minDiffusionFactor; // alternative to depth setting were we keep traversing while DF is less than this value
    int maxNodesVisited; // keep traversing until visit this number of nodes
    double dfBackgroundSmoothing; // for calculating the diffusion factor (i.e., edge weight), gives the background weight
    
    int metamapConfidenceOn; // use metamap confidence score as P(Q) for P(D|Q) = P(Q|D)*P(D) / P(Q)
    int docRatioQueryPrior; // use the number of documents containg the query term as a P(Q)
    int arcCountQueryPrior; // use the number of arc (in & out) as query prior - P(q) = |E(q)| / |E|
    
    int mapToVisitsOn; // after retrieval is complete map individual report names to visits
    string reportToVisitMappingFile; // location of the file that contains mapping of: "reportFileName.xml VisitId"
    
    std::string rerankResultsFile; // path to a TREC results file to use for reranking
    int rerankCount; // number of top docs to rerank
    int onlyJudgedDocuments; // filter the retrieval results to include only document judged according to qrels
    
    string graphVizDir; // generate a graphViz graph for debugging and evaluation purposes. 
    
    void get() {
        // the string with quotes are the actual variable names to use for specifying the parameters
        indexPath                   = ParamGetString("index");
        queryFile                   = ParamGetString("query");
        resultCount                 = ParamGetInt("resultCount", 1000);
        resultFile                  = ParamGetString("resultsFile", LocalParameter::queryFile+".results");
        stopwordFile                = ParamGetString("stopwordFile");
        
        // calc prior: mleidf, tfidf, LM-JM, LM-Dirichlet
        weightScheme                = ParamGetString("weightScheme","LM-JM");
        LMjm_lambda                 = ParamGetDouble("LM-JM.lambda", 0.9);
        LMdirichlet_mu              = ParamGetInt("LM-Dirichlet.mu", 2000);
        indri_tfidf_k1              = ParamGetDouble("indri_tfidf_k1", 0.9);
        indri_tfidf_b               = ParamGetInt("indri_tfidf_b", 0.5);
        
        depth                       = ParamGetInt("depth", 1);
        diffusionMix                = ParamGetDouble("diffusionMix", 0.5);
        diffusionOn                 = ParamGetInt("diffusionOn", 1);
        minDiffusionFactor          = ParamGetDouble("minDiffusionFactor", -0.1);
        maxNodesVisited             = ParamGetInt("maxNodesVisited", -1);
        dfBackgroundSmoothing       = ParamGetDouble("dfBackgroundSmoothing", 0.9);
        
        metamapConfidenceOn         = ParamGetInt("metamapConfidenceOn", 0);
        docRatioQueryPrior          = ParamGetInt("docRatioQueryPrior", 0);
        arcCountQueryPrior          = ParamGetInt("arcCountQueryPrior", 0);
        
        mapToVisitsOn               = ParamGetInt("mapToVisitsOn", 0);
        reportToVisitMappingFile    = ParamGetString("reportToVisitMappingFile");
        
        rerankResultsFile           = ParamGetString("rerankResultsFile");
        rerankCount                 = ParamGetInt("rerankCount", 10);
        onlyJudgedDocuments         = ParamGetInt("onlyJudgedDocuments", 0);
        
        graphVizDir                = ParamGetString("graphVizDir");
    }
};

void GetAppParam()
{
    LocalParameter::get();
}


///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

void logger(string msg, loglevel msglvl=TRACE)
{
    if (msglvl <= level)
        cout << "Log: (" << msglvl << ") -- " << msg << endl;
}

// prints out node CUIS for a node vector (for debugging)
void print_vector(vector<Arc> &v)
{
    //cout << "\tpath vector len " << v.size() << ": ";
    for (vector<Arc>::iterator it = v.begin(); it != v.end(); it++)
    {
        cout << (*nodeToCUIMap)[g->source(*it)] << " -> " << (*nodeToCUIMap)[g->target(*it)] << ", ";
    }
    cout << endl;
}

// prints out node CUIS for a node vector (for debugging)
void print_node_vector(vector<Node> &v)
{
    //cout << "\tpath vector len " << v.size() << ": ";
    for (vector<Node>::iterator it = v.begin(); it != v.end(); it++)
    {
        cout << (*nodeToCUIMap)[*it] << ", ";
    }
    cout << endl;
}

int lookupTermId(Node node)
{
    string cui = (*nodeToCUIMap)[node];
    return idx->term(cui);
}

///////////////////////////////////////////////////////////////////////////////
// Scoring functions: prior calculation on nodes
///////////////////////////////////////////////////////////////////////////////

// calculate inverse doc freq idf(t) = log ( N / d_t )
double idf(TERMID_T termId)
{
    return log(idx->docCount() / (double)idx->docCount(termId));
}

// term freq
double tf(TERMID_T termId, DocInfo *docInfo)
{
    return docInfo->termCount();
}

// calculate maximum liklihood estimate mle(t,d) = tf(t) / |d|
double mle(TERMID_T termId, DocInfo *docInfo)
{
    return tf(termId, docInfo) / (double) idx->docLength(docInfo->docID());
}

// calculate term & collection maximum liklihood estimate mle(t,d) = [ tf(t) / |d| ] * [ C / Dt ]
// C = collection size, Dt = num documents containing term t
double cmle(TERMID_T termId, DocInfo *docInfo)
{
    return mle(termId, docInfo) / 
        (double) ( idx->docCount(termId) / idx->docCount() );
}

// bm25
double bm25(TERMID_T termId, DocInfo *docInfo)
{
    float k1 = 5.2, b = 0.4;
    
    double _tf = tf(termId, docInfo);
    double avgDocLength = idx->docLengthAvg();
    double docLength = idx->docLength(docInfo->docID());
    
    double numerator = (_tf * (k1+1));
    double denominator =  _tf + k1 * ( 1-b + b * (docLength / avgDocLength) );
    return idf(termId) * (numerator / denominator);
}

// Indri/lemur tfidf score function:
//
//                        (queryTermWeight * IDF * (K1 + 1)) * occurrences
// score = ------------------------------------------------------------------------
//          occurrences + (K1 * (1-B)) + (K1 * B * 1/avgDocLength) * documentLength

double indri_tfidf(TERMID_T termId, DocInfo *docInfo)
{
    double _tf = tf(termId, docInfo);
    double _idf = idf(termId);
    double avgDocLength = idx->docLengthAvg();
    double docLength = idx->docLength(docInfo->docID());

    double k1 = LocalParameter::indri_tfidf_k1; double b = LocalParameter::indri_tfidf_b;
    
    double queryTermWeight = (_idf * 1000 * 1) / (1 + 1000);
    
    double numerator = (queryTermWeight * _idf * (k1 + 1)) * _tf; // old
    //double numerator = (queryTermWeight * _idf * k1) * _tf;
    double denominator = _tf + (k1 * (1-b)) + (k1 * b * (1/avgDocLength)) * docLength;
    
    return numerator / denominator;
}

double p_idf(TERMID_T termId) {
    //cout << "P(" << idx->term(termId) << ") " << idx->termCount(termId) / (double) idx->termCount() << endl;
    
    return idx->docCount(termId) / (double) idx->docCount();
}

// language model with Jelinek-Mercer Smoothing
double LM_JM(TERMID_T termId, DocInfo *docInfo, double lambda)
{
    return lambda*(docInfo->termCount()/(float)idx->docLength(docInfo->docID())) + (1-lambda)*(idx->termCount(termId) / idx->termCount());
}

// mu *  tf_i_C
//       ------
//        |C|
// -------------
//  |D| + mu
double LM_Dirichlet_Background(TERMID_T termId, int docLength, double qProb) {
    double numerator = LocalParameter::LMdirichlet_mu * ( idx->termCount(termId) / (double)idx->termCount());
    double denominator = (double)LocalParameter::LMdirichlet_mu + docLength;
    
    
    //cout << "LMB: " << numerator << " / " << denominator << " / " << qProb << " = " << (numerator / denominator) / qProb;
    
    return (numerator / denominator);
    
}

// language model with Dirichlet Smoothing
//                        cf
// S(t, D) = tf + mu * ( ---- )
//                       |C|
//           -------------------
//               |D| + mu
double LM_Dirichlet(TERMID_T termId, DocInfo *docInfo, double qProb)
{
    double foreground =  (double)docInfo->termCount() / (double)( idx->docLength(docInfo->docID()) + LocalParameter::LMdirichlet_mu );
    double background = LM_Dirichlet_Background(termId, idx->docLength(docInfo->docID()), qProb);
    
    return (foreground + background ) / qProb;
        
}



// return the prior probability of this node
double prior(Node &node, DocInfo *docInfo, double qProb)
{
    int termId = idx->term((*nodeToCUIMap)[node]);
    if(termId == 0)
    {
        cout << "ERROR: invalid node with cui " << (*nodeToCUIMap)[node] << " not found in index" << endl;
    }
    
    if (LocalParameter::weightScheme == "mle")
    {
        return mle(termId, docInfo);
    }
    if (LocalParameter::weightScheme == "cmle")
    {
        return cmle(termId, docInfo);
    }
    else if (LocalParameter::weightScheme == "tfidf") {
        return tf(termId, docInfo)*idf(termId);
    }
    else if (LocalParameter::weightScheme == "bm25") {
        return bm25(termId, docInfo);
    }
    else if (LocalParameter::weightScheme == "indri_tfidf") {
        return indri_tfidf(termId, docInfo)*idf(termId);
    }
    else if (LocalParameter::weightScheme == "mleidf") {
        return mle(termId, docInfo)*idf(termId);
    }
    else if (LocalParameter::weightScheme == "LM-JM") {
        return LM_JM(termId, docInfo, LocalParameter::LMjm_lambda);
    }
    else if (LocalParameter::weightScheme == "LM-Dirichlet") {
        return LM_Dirichlet(termId, docInfo, qProb);
    }
    cout << "ERROR: unknown weighting scheme " << LocalParameter::weightScheme << endl;
}


///////////////////////////////////////////////////////////////////////////////
// Diffusion Factor and Edge Weighting
///////////////////////////////////////////////////////////////////////////////


// returns a score based on priors of nodes related to the supplied node
double edgeCount(Node &node)
{
    double numEdges = 0;
    for (ListDigraph::InArcIt arc(*g, node); arc != INVALID; ++arc)
    {
        numEdges++;
    }
    
    VLOG(2) << "edgeCount for " << (*nodeToCUIMap)[node] << " " << numEdges;
    
    return numEdges;
}

// caluculate the avg cosine between the given node and all query nodes
double query_sim(Node node, vector<string> queryConcepts)
{
    int termId = idx->term((*nodeToCUIMap)[node]);
    map<int, double> nodeDocumentVector = documentVector(termId);
    
    double simSum = 0;
    int count = 0;
    for (vector<string>::iterator it = queryConcepts.begin(); it != queryConcepts.end(); it++) {
        simSum += cosine(nodeDocumentVector, documentVector(idx->term(*it)));
        count++;
    }
    
    return simSum / count;
}

// returns the diffusionFactor to transistion trough a series of edges
double diffusionFactor(vector<Arc> pathToQuery)
{
    double total_df = 1;
    
    //cout << "calculate df " << endl;
    
    if(LocalParameter::diffusionOn) {
        
        for (vector<Arc>::iterator it = pathToQuery.begin(); it != pathToQuery.end(); it++)
        {
            Arc arc = *it;
            Node source = g->source(arc);
            Node target = g->target(arc);
            
            // the semantic sim weight
            double sim_weight = (*simArcs)[arc];
            
            // the relationship type weigth
            string reltype = (*reltypeArcs)[arc];

            
            double rel_weight;
            if(relweightsweep > 0) {
                rel_weight = snomedRelationshipTypeWeight(reltype, "../snomed_relations/reltypes.txt."+double_to_str(relweightsweep));
            } else {
                rel_weight = relationshipTypeWeight(reltype, (*nodeToCUIMap)[source].substr(0,1) == "c" || (*nodeToCUIMap)[source].substr(0,1) == "C");
            }
            
            // combined weight - the df
            double df = LocalParameter::diffusionMix * sim_weight + (1-LocalParameter::diffusionMix)*rel_weight;
            
            // apply topology connecteness
            //df /= (double)edgeCount(source);
            
            
            
            //VLOG(2) << "diffusion Factor for arc " << (*nodeToCUIMap)[source] << " -> " << (*nodeToCUIMap)[target] << ": avg(" << sim_weight << "," << reltype << "=" << rel_weight << ") =" << df;
            
            // df constanst background smoothing to deal with when df = 0.
            df = LocalParameter::dfBackgroundSmoothing*df + (1-LocalParameter::dfBackgroundSmoothing);
            total_df = total_df * df;
        }
        
        // smooth
        //total_df = LocalParameter::dfBackgroundSmoothing*total_df + (1-LocalParameter::dfBackgroundSmoothing);
    } 
    
    return total_df;
}

double calcEmptyDocLMScore(TERMID_T termId, vector<Arc> pathToQuery, double qProb) {
    return log(LM_Dirichlet_Background(termId, 0, qProb));
}

// updates the emptry document score by applying document length normalisation
double lengthNormalisedEmptyDocScore(int depth, TERMID_T termId, int docLength, double df, double qProb) {
    double score = 0;
    for(int i=0; i <= LocalParameter::depth; i++) {
        for (map<TERMID_T, double>::iterator it = emptyDocScores[i].begin(); it != emptyDocScores[i].end(); it++) {
            double prob = exp(it->second) * qProb;
            double lenNormalisedScore = log( (df * ((prob * (double)LocalParameter::LMdirichlet_mu) / ((double)docLength + LocalParameter::LMdirichlet_mu)) / qProb) );
            VLOG(2) << "\t" <<  "Length normalised background for lvl " << depth << " term " << idx->term(it->first) << ": " << lenNormalisedScore << endl;
            score += lenNormalisedScore;
        }
    }
    
    return score;
}

///////////////////////////////////////////////////////////////////////////////
// Comparison against qrels
///////////////////////////////////////////////////////////////////////////////

map<int, map<DOCID_T, int> > relevanceJudgements;
int relevanceJudgement(int queryId, DOCID_T docId) {

    if(relevanceJudgements.size() == 0) {
        ifstream file("params/medtrack-all.qrel"); // e.g. 101 0 corpus/+C3SPIppTQGH 0
        string line;
        while(getline(file, line))
        {
            stringstream linestream(line);
            string qIdStr;
            string zero;
            string docName;
            string relevance;
            
            getline(linestream, qIdStr, ' ');
            getline(linestream, zero, ' ');
            getline(linestream, docName, ' ');
            getline(linestream, relevance, ' ');
            
            int qId = string_to_int(qIdStr);
            
            
            DOCID_T docId = idx->document(docName);
            
            map<DOCID_T, int> docIds;
            if(relevanceJudgements.count(qId)) {
                docIds = relevanceJudgements[qId];
            }
            docIds[docId] = string_to_int(relevance);
            relevanceJudgements[qId] = docIds;
            
        }
        file.close();
        
    }
    
    int relevance = -1;
    if(relevanceJudgements.count(queryId) > 0 && relevanceJudgements[queryId].count(docId) > 0) {
        relevance = relevanceJudgements[queryId][docId];
    }
    
    
    return relevance;
}



int currentQueryId;
bool isRelevantDoc(int queryId, DOCID_T docId) {
    return relevanceJudgement(queryId, docId) > 0;
}

///////////////////////////////////////////////////////////////////////////////
// Graph traversal functions
///////////////////////////////////////////////////////////////////////////////

// count the number of edges connecting this node
int countNodeArcs(Node &node) {
    int edgeCount = 0;
    
    for (ListDigraph::OutArcIt i(*g, node); i!=INVALID; ++i)
        edgeCount++;
    
    for (ListDigraph::InArcIt i(*g, node); i!=INVALID; ++i)
        edgeCount++;
    
    return edgeCount;
}


// the visit node function - for each document found at this node score the document based on its prior and its diffusionFactor
void scoreNode(pair<Node, double> &nodePair, vector<Arc> &pathToQuery, map<DOCID_T, double> &scores, int depth)
{
    
    Node node = nodePair.first;
    double qProb = nodePair.second;
    
    // the term does not appear in the collection, there no document to score against it
    // note that the traverse method will still pass through this non-collection node
    TERMID_T termId = idx->term((*nodeToCUIMap)[node]);
    VLOG(1) << "hn: " << (*nodeToCUIMap)[node] << " lvl: " << pathToQuery.size() << " #docs:" << idx->docCount(termId) << endl;
    if (termId == 0)
        return;
    
    // node has previous been visited via some other path?
//    bool previously_visited = diffusionFactors.count(node) > 0;
//    if (!previously_visited) {
//        map<DOCID_T, double> e;
//        diffusionFactors[node] = e;
//    }
    

    
    double theDiffusionFactor = diffusionFactor(pathToQuery);
    //double theTopoScore = 1 /(*nodeToTopoScoreMap)[node];
    
    //VLOG(2) << " toposcore: " << (*nodeToCUIMap)[node] << " " << theTopoScore << endl;
    int relevantDocCount = 0;
    int newRelevantDocCount = 0;

    // master
    
    DocInfoList *docInfoList = idx->docInfoList(termId);
    docInfoList->startIteration();
    while (docInfoList->hasMore())
    {
        DocInfo *doc = docInfoList->nextEntry();
        
        if (LocalParameter::onlyJudgedDocuments && relevanceJudgement(currentQueryId, doc->docID()) < 0) {
            continue;
        }
                
        VLOG(2) << "Scoring document: " << idx->document(doc->docID()) << " (id=" << doc->docID() << ") against query " << (*nodeToCUIMap)[node] << " ";
        
        // MEDTRACK 2012 FILTER
        //if(idx->document(doc->docID()) == "corpus/0VppFYq+9Tjm")
        //    continue;
        
        // first if the doc has not been scored, push on the emptyDocScore
        if(LocalParameter::weightScheme == "LM-Dirichlet") {
            VLOG(2) << "\t" << "doc scored before?: " << scores.count(doc->docID()) << endl;
            if(scores.count(doc->docID()) == 0) {
                double emptyDocScore = lengthNormalisedEmptyDocScore(depth, termId, idx->docLength(doc->docID()), theDiffusionFactor, qProb);
                VLOG(2) << "\t" <<  "Writing empty document score: " << emptyDocScore << endl;
                scores[doc->docID()] += emptyDocScore;
            }
            
            VLOG(2) <<  "\t" << "2. Current accumulator score: " << scores[doc->docID()] << endl;
            
            // take off the background smoothing for this doc
            double backGroundProb = LM_Dirichlet_Background(termId, idx->docLength(doc->docID()), qProb);
            VLOG(2) << "\t" <<  "background for " << idx->term(termId) << "|" << idx->document(doc->docID()) << " " <<  backGroundProb << " * " << theDiffusionFactor << " = " << log(backGroundProb * theDiffusionFactor) << endl;
            scores[doc->docID()] -= log(backGroundProb * theDiffusionFactor);
            
            VLOG(2) << "\t" <<  "3. After background removal: " << scores[doc->docID()] << endl;
        }
        
        double thePrior = prior(node, doc, qProb);
        double theScore = log(thePrior * theDiffusionFactor);  //* theRelatedPrior; // * theTopoScore *  ((float)conceptScore / 1000) * ( 1 + log(conceptCount) ) ;
        // apply dampening for non-query nodes
        if(pathToQuery.size() > 0) {
            //theScore = related_dampener * theScore;
        }
        // apply proposition count
        //theScore = theScore * ((float)conceptScore / 1000) * ( 1 + log(conceptCount) );
        
        string relevanceIndicator = " ";
        if(isRelevantDoc(currentQueryId, doc->docID())) {
            relevanceIndicator = "*";
            if(pathToQuery.size() > 0) {
                
                if (relevantLvl0Docs.find(doc->docID()) == relevantLvl0Docs.end()) { // not already seen at lvl0
                    newRelevantDocCount++;
                }
                relevantDocCount++; 

            } else {
                relevantLvl0Docs.insert(doc->docID());
                relevantDocCount++;
                newRelevantDocCount++;
            }
        }
        
        VLOG(1) <<  relevanceIndicator << "RSV(" << (*nodeToCUIMap)[node] << "|" << idx->document(doc->docID()) << ") = prior:" << thePrior << " * df_" << pathToQuery.size() << ":" << theDiffusionFactor << "\t=> " << theScore << endl;
        
        //        VLOG(1) << "RSV(" << (*nodeToCUIMap)[node] << "|" << idx->document(doc->docID()) << ") = prior:" << thePrior << " * df_" << pathToQuery.size() << ":" << theDiffusionFactor << " * related_prior:" << theRelatedPrior << " * topo:" << theTopoScore << " * concept_score:" << ((float)conceptScore / 1000) << " * concept_count:" << ( 1 + log(conceptCount) ) << "\t=> " << theScore << endl;
        
        // we've previously visited this node via some other path:s
        // keep the path with least effor. 
//        if (previously_visited && false)
//        {
//            double prevDiffusionFactor = diffusionFactors[node][doc->docID()];
//            if (theDiffusionFactor < prevDiffusionFactor) // current path is better than prev
//            {
//                scoreAccumulator.incScore(doc->docID(), -(thePrior*prevDiffusionFactor)); // remove old path contribution
//                conceptScores[doc->docID()]
//                scoreAccumulator.incScore(doc->docID(), theScore); // add new path contribution
//                diffusionFactors[node][doc->docID()] = theDiffusionFactor;
//            }
//        } else { // first time visiting this node, score it
//        scoreAccumulator.incScore(doc->docID(), theScore);
        
        
        // WRITE THE SCORE FOR THIS DOCUMENT
        
        VLOG(2) << "\twriting the final score, before: " << scores[doc->docID()] << " + " << theScore << " = " << scores[doc->docID()]+theScore << endl;
        
        scores[doc->docID()] += theScore;
        
        VLOG(2) << "\t" <<  "4. Final recorded score " << scores[doc->docID()] << endl;
        
        
        //diffusionFactors[node][doc->docID()] = theDiffusionFactor;
        
        // DEBUG TO SEE HOW MUCH SCORE IS CHANGING
        if(pathToQuery.size() == 0) {
            lvl0Scores[doc->docID()] += theScore;
        } else {
            lvl1Scores[doc->docID()] += theScore;
        }
        
        
    }
    
    VLOG(1) << "complete hn: " << (*nodeToCUIMap)[node] << " lvl: " << pathToQuery.size() << " #docs:" << idx->docCount(termId) << " (" << relevantDocCount << " relevant)" << endl;
    
    if(LocalParameter::graphVizDir.length() > 0) {
    
        *graphViz << "\t" << (*nodeToCUIMap)[node] << " [label=\"" << concept_lookup((*nodeToCUIMap)[node].c_str()) << " (" << newRelevantDocCount << "/" << relevantDocCount << ") #" << idx->docCount(termId) << "\"";
        if(pathToQuery.size() == 0) {
            *graphViz << ", color=red";
        }
        *graphViz << "];" << endl;
    }
    
    delete docInfoList;
} 

// recursive traversal of the graph
int traverse(pair<Node, double> &node, vector<Arc> &pathToQuery, int nodeCount, map<DOCID_T, double> &scores, int &depth, bool touchNodesOnly) {

    TERMID_T termId = idx->term((*nodeToCUIMap)[node.first]);

    
    
    
    bool continueTraverse = true;
    if(LocalParameter::minDiffusionFactor > 0) { // adpative depth based on diffusion factor
        continueTraverse = (diffusionFactor(pathToQuery) >= LocalParameter::minDiffusionFactor);
        //cout << "min " << diffusionFactor(pathToQuery) << endl;
    } 
    
    if(LocalParameter::maxNodesVisited > 0) { 
        continueTraverse = (nodeCount < LocalParameter::maxNodesVisited);
        //cout << "max " << nodeCount << " / " << LocalParameter::maxNodesVisited << endl;
    } 
    
    if(LocalParameter::depth >= 0) { // static depth based on depth level
        continueTraverse = (pathToQuery.size() <= depth);
//        cout << "depth " << pathToQuery.size() << " / " << depth << endl;
    }
    
    // stopword? skip
    if(stopwordList.find((*nodeToCUIMap)[node.first]) != stopwordList.end())
        continueTraverse = false;
    
    // avoid cycles
    for (vector<Arc>::iterator it = pathToQuery.begin(); it != pathToQuery.end(); it++)
    {
        if (g->target(*it) == node.first) {
            //cout << "loop found: node:" << (*nodeToCUIMap)[node.first] << " arc " << (*nodeToCUIMap)[g->source(*it)] << " -> " << (*nodeToCUIMap)[g->target(*it)] << endl;
            continueTraverse = false;
        }
    }


        
    if(continueTraverse)
    {
        
        VLOG_IF(2, !touchNodesOnly) << "calling traverse on " << (*nodeToCUIMap)[node.first] << endl;
        
        // score all the docs attached to this node
        if(touchNodesOnly) {
            if (termId != 0) {
                double score = calcEmptyDocLMScore(termId, pathToQuery, node.second);
                emptyDocScores[pathToQuery.size()][termId] = score;
                VLOG(2) << "Empty_document_score(" <<  idx->term(termId) << ") = " << score << endl;
            }
        } else {
            scoreNode(node, pathToQuery, scores, pathToQuery.size());
        }
        nodeCount++;
    
        
        // follow children of this node
        for (ListDigraph::InArcIt arc(*g, node.first); arc != INVALID; ++arc)
        {
            if (g->source(arc) != node.first) {
                pathToQuery.push_back(arc);
                pair<Node, double> newPair (g->source(arc), node.second);
                
                VLOG_IF(1, !touchNodesOnly) << " arc: " << (*nodeToCUIMap)[g->source(arc)] << " -> " << (*nodeToCUIMap)[g->target(arc)] << " lvl:" << pathToQuery.size() << "-" << pathToQuery.size()-1 << " df:" << diffusionFactor(pathToQuery);
                
                if(!touchNodesOnly && pathToQuery.size() <= depth) {
                    string reltype = (*reltypeArcs)[arc];
                    
                    
                    if (LocalParameter::graphVizDir.length() > 0) {
                        if(isdigit(reltype[0])) {
                            reltype = concept_lookup((*reltypeArcs)[arc]);
                        }
                        *graphViz << "\t" << (*nodeToCUIMap)[g->source(arc)] << " -> " << (*nodeToCUIMap)[g->target(arc)] << " [label=\"" << reltype << " (" << diffusionFactor(pathToQuery) << ")\"];" << endl;
                    }
                }
                
                //if ( (*reltypeArcs)[arc] != "116680003" ) {
                    nodeCount = traverse(newPair, pathToQuery, nodeCount, scores, depth, touchNodesOnly);
                //}
                
                pathToQuery.pop_back();
            }
        }
        
    }
    return nodeCount;
}

///////////////////////////////////////////////////////////////////////////////
// Query priors
///////////////////////////////////////////////////////////////////////////////

// query prior - metamap score
map<string, double> calcMetaMapConfidenceQueryPrior(map<string, double> &queryConceptsWithPriors) {
    map<string, double> priors;
    int confidenceSum = 0;
    for (map<string, double>::iterator it = queryConceptsWithPriors.begin(); it != queryConceptsWithPriors.end(); it++) {
        stringstream ss(it->first);
        string query, confidenceStr;
        getline(ss, query, '_');
        getline(ss, confidenceStr, '_');
        
        int confidence = 1;
        if(confidenceStr.size() > 0 && LocalParameter::metamapConfidenceOn) {
            confidence = string_to_int(confidenceStr);
            confidenceSum += confidence;
        } else {
            confidenceSum = 1;
        }
        priors[query] = confidence; // this will index based on query with confidence score stipped out
    }
    
    queryConceptsWithPriors = priors; // gets rid of the confidence scores from the key

    for (map<string, double>::iterator it = queryConceptsWithPriors.begin(); it != queryConceptsWithPriors.end(); it++) {
        queryConceptsWithPriors[it->first] = it->second / confidenceSum;
    }
    
    return queryConceptsWithPriors;
}

// query prior - doc ratio
void calcDocRatioQueryPrior(map<string, double> &queryConceptsWithPriors) {
    int docCountSum = 0;
    for (map<string, double>::iterator it = queryConceptsWithPriors.begin(); it != queryConceptsWithPriors.end(); it++) {
        if(idx->term(it->first) > 0 ) {
            queryConceptsWithPriors[it->first] = idx->docCount(idx->term(it->first));
            docCountSum += queryConceptsWithPriors[it->first];
        }
    }
    for (map<string, double>::iterator it = queryConceptsWithPriors.begin(); it != queryConceptsWithPriors.end(); it++) {
        if(idx->term(it->first) > 0 ) {
            VLOG(2) << "calculating doc ratio query prior P(" << it->first << ") = " << it->second << "/" << docCountSum << "=" << it->second / (double)docCountSum << endl;
            //queryConceptsWithPriors[it->first] = (it->second / (double)docCountSum);
            queryConceptsWithPriors[it->first] = it->second / idx->docCount();
        }
    }
}

// query prior - arc count
void calcArcCountQueryPrior(map<string, double> &queryConceptsWithPriors) {
    for (map<string, double>::iterator it = queryConceptsWithPriors.begin(); it != queryConceptsWithPriors.end(); it++) {
        double probQ = 0;
        if (cuiToNodeId.count(it->first) > 0) {
            Node node = g->nodeFromId(cuiToNodeId[it->first]);
            probQ = countNodeArcs(node) / (double) countArcs(*g);
        } 
        queryConceptsWithPriors[it->first] = probQ;
    }
}

// apply the varies query priors
void calcQueryPriors(map<string, double> &queryConceptsWithPriors) {
    calcMetaMapConfidenceQueryPrior(queryConceptsWithPriors); // always call this to strip confidence scores (even if they are never used)
    
    if(LocalParameter::docRatioQueryPrior) {
        calcDocRatioQueryPrior(queryConceptsWithPriors);
    }
    
    if(LocalParameter::arcCountQueryPrior) {
        calcArcCountQueryPrior(queryConceptsWithPriors);
    }
    
}

///////////////////////////////////////////////////////////////////////////////
// Query processing functions
///////////////////////////////////////////////////////////////////////////////

// process a single query term/concept by locating the node in the graph and initiating traversal from that node
void processQueryTerm(string query, map<DOCID_T, double> &scores, bool touchNodesOnly, double queryProb)
{

    
    transform(query.begin(), query.end(), query.begin(), ::tolower);
    
    // query reduction by snomed connectedness (a la ADCS'12)
//    int theSnomedRelCount = (float)findRelatedConcepts(query).size();
//    if(theSnomedRelCount == 1)
//        return;
    
    if(!touchNodesOnly)
        cout << " " << query << flush;
        VLOG(2) << "P(" << query << ") = " << queryProb << endl;
    
    diffusionFactors.clear();
    
    if(cuiToNodeId.count(query) > 0) {
        Node queryNode = g->nodeFromId(cuiToNodeId[query]);
        vector<Arc> pathToQuery;
        
        // traverse the graph from this node
        pair<Node, double> nodeAndScore (queryNode, queryProb);
        int nodeCount = traverse(nodeAndScore, pathToQuery, 0, scores, LocalParameter::depth, touchNodesOnly);
        relevantLvl0Docs.clear();
        
        if(!touchNodesOnly)
            cout << "(" << nodeCount << ") ";
        
        if (idx->term(query) == 0 && !touchNodesOnly) {
            cout << " [not found in index]" << " ";
        }
    } else if(!touchNodesOnly) {
        cout << "(x)" << " ";
    }
}

// process a single query and score the results to the trec results file
void processQuery(Document *qryDoc, ResultFile &resultFile)
{
    const char *queryID = qryDoc->getID();
    currentQueryId = atoi(queryID);
    cout << "Running query: "<< queryID << "\n\t";


    string graphVizLoc = "/dev/null";
    if(LocalParameter::graphVizDir.length() > 0) {
        graphVizLoc = (LocalParameter::graphVizDir+"/"+int_to_str(currentQueryId)+".dot").c_str();    
    } 

    ofstream theGraphviz( graphVizLoc.c_str() );
    graphViz = &theGraphviz;

    *graphViz << "digraph " << currentQueryId << " {" << endl;
    *graphViz << "label=\"" << currentQueryId << "\"" << endl;


    //*graphViz << "testing" << endl;
    
    // rerank another results file
    if(LocalParameter::rerankResultsFile.length() > 0) {
        lvl0Results.findResult(qryDoc->getID(), rerankResults);
    }
    
    // used only for debug / statistics
    lvl0Scores.clear();
    lvl1Scores.clear();
    
    docNodeCount.clear();
    


    // build up the query vector
    //vector<string> queryConcepts;
    map<string, double> queryConceptsWithPriors;
    
    qryDoc->startTermIteration();
    while (qryDoc->hasMore()) {
        const Term *qryTerm = qryDoc->nextTerm();
        queryConceptsWithPriors[qryTerm->spelling()] = 1;
    }

    emptyDocScores.clear();
    
    map<DOCID_T, double> scores;
    
    calcQueryPriors(queryConceptsWithPriors);

    
    if(LocalParameter::weightScheme == "LM-Dirichlet") {
        // VISIT ALL THE NODES TO PRECOMPUTE THE EMPTY DOCUMENT SCORE
        VLOG(1) << "Calculating empty document score for query " << queryID << endl;
        for (map<string, double>::iterator it = queryConceptsWithPriors.begin(); it != queryConceptsWithPriors.end(); it++) {
            processQueryTerm(it->first, scores, true, it->second);
        }
    }
    
    // EXECUTE RETRIEVAL ON THE QUERY TERM
    VLOG(1) << "Running retrieval on query " << queryID << endl;
    for (map<string, double>::iterator it = queryConceptsWithPriors.begin(); it != queryConceptsWithPriors.end(); it++) {
        processQueryTerm(it->first, scores, false, it->second);
    }
    
    cout << endl;

    // print results    
    IndexedRealVector results(idx->docCount());
    results.clear();
    
    if(LocalParameter::rerankResultsFile.length() == 0) {
        for (map<DOCID_T, double>::iterator it=scores.begin(); it!=scores.end(); it++) {
            int queryLen = queryConceptsWithPriors.size();
            queryLen = 1;
            results.PushValue(it->first, it->second);                
        }
    } else {    // do reranking
        results = *rerankResults;
        int count = 0;
        for(IndexedRealVector::iterator it = rerankResults->begin(); it!=rerankResults->end() && count < LocalParameter::rerankCount; it++) {
            results.IncreaseValueFor((*it).ind, (scores[(*it).ind] / queryConceptsWithPriors.size())-(*it).val);
            count++;
        }
        delete rerankResults;
    }
    
    results.Sort();
    resultFile.writeResults(qryDoc->getID(), &results, LocalParameter::resultCount);
    
    //printDebugScores(queryId)

    if(LocalParameter::graphVizDir.length() > 0) {
        *graphViz << "}" << endl;
        graphViz->close();
    }
}

// print the debug query results
void printDebugScores(string queryID) {
    set<DOCID_T> docs;
    for (map<DOCID_T, double>::iterator it = lvl0Scores.begin(); it != lvl0Scores.end(); it++) {
        docs.insert(it->first);
    }
    for (map<DOCID_T, double>::iterator it = lvl1Scores.begin(); it != lvl1Scores.end(); it++) {
        docs.insert(it->first);
    }
    
    for(set<DOCID_T>::iterator it = docs.begin(); it != docs.end(); it++)
    {
        cout << queryID << "\t" << idx->document(*it) << "\t" << lvl0Scores[*it] << "\t" << lvl1Scores[*it] << "\t" << (int)((lvl1Scores[*it]/lvl0Scores[*it])*100) << "%\t" << lvl0Scores[*it]+lvl1Scores[*it] << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Setup & param reading functions
///////////////////////////////////////////////////////////////////////////////

void readStopwordFile(string &stopwordFile) {
    ifstream fh(stopwordFile.c_str());
    string line;
    while(getline(fh, line)) {
        stopwordList.insert(line);
    }
    fh.close();
}

void printParams()
{
    cout << "Running retrieval: ";

    cout << "weighting: " << LocalParameter::weightScheme << " ";

    if (LocalParameter::weightScheme == "LM-JM")
    {
        cout << "lambda=" << LocalParameter::LMjm_lambda;
    }
    else if (LocalParameter::weightScheme == "LM-Dirichlet")
    {
        cout << "mu=" << LocalParameter::LMdirichlet_mu;
    }
    if (LocalParameter::stopwordFile.length() > 0) {
        cout << " stoplist=" << LocalParameter::stopwordFile;        
    }
    
    if(LocalParameter::docRatioQueryPrior) {
        cout << " docRatioQueryPrior: on ";
    }
    
    if(relweightsweep > 0) {
        cout << " relweightsweep: " << relweightsweep;
    }
    
    if(LocalParameter::depth >= 0) {
        cout << " depth: " << LocalParameter::depth;
    }
    
    cout << endl;
    
    if(!LocalParameter::diffusionOn) {
        cout << "WARNING: diffusionFactor: OFF!! ";
    }
    
    if(LocalParameter::onlyJudgedDocuments) {
        cout << "WARNING: onlyJudgedDocuments: ON!! ";
    }
    
    if(LocalParameter::graphVizDir.length() > 0) {
        cout << "GraphViz: " << LocalParameter::graphVizDir << " ";
    }
    cout << endl;
}

///////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////

void AppMain(int argc, char *argv[])
{
    google::InitGoogleLogging(argv[0]);
    
    if (argc > 2) {
        LocalParameter::depth = string_to_int(argv[2]);
        if(LocalParameter::depth < 0 || LocalParameter::depth > 100)
        {
            cerr << "Invalid setting for depth: " << LocalParameter::depth << endl;
            exit(1);
        }
    }
    
    if (argc > 3) {
        relweightsweep = string_to_double(argv[3]);
    }
    
    if(LocalParameter::stopwordFile.length() > 0) {
        readStopwordFile(LocalParameter::stopwordFile);
    }
    
    printParams();
        
    // open index
    idx = lemur::api::IndexManager::openIndex(LocalParameter::indexPath);
    prepareDB("../snomed_relations/snomed_rel.db");
    
    // read the graph, populating required maps
    ListDigraph theGraph;
    g = &theGraph;
    ListDigraph::NodeMap<string> theNodeToCUIMap(*g);
    nodeToCUIMap = &theNodeToCUIMap;
    ListDigraph::NodeMap<float> theNodeToTopoScoreMap(*g);
    nodeToTopoScoreMap = &theNodeToTopoScoreMap;
    ListDigraph::ArcMap<double> theSimArcs(*g);
    simArcs = &theSimArcs;
    ListDigraph::ArcMap<string> theRelTypeArcs(*g);
    reltypeArcs = &theRelTypeArcs;
    
    cout << "Loading graph..." << flush;
    digraphReader(*g, LocalParameter::indexPath+"/graph-active.lgf").nodeMap("Concepts", *nodeToCUIMap).nodeMap("TopoScore", *nodeToTopoScoreMap).arcMap("SimArcWeights", *simArcs).arcMap("RelTypeArc", *reltypeArcs).run();
    for (ListDigraph::NodeIt n(*g); n != INVALID; ++n)
    {
        cuiToNodeId[theNodeToCUIMap[n]] = g->id(n);
    }
    cout << "complete, " << countNodes(*g) << " nodes; " << countArcs(*g) << " edges." << endl;
    
    
    // create a results file in trec format
    ResultFile resultFile(true); // trec format
    
    string resultFileName = LocalParameter::resultFile;
    if(relweightsweep > 0 )
        resultFileName = resultFileName + "." + int_to_str(LocalParameter::depth) + "_" + double_to_str(relweightsweep);
    
    ofstream outfile(resultFileName.c_str());
    resultFile.openForWrite(outfile, *idx);
    
    // rerank another results file
    if(LocalParameter::rerankResultsFile.length() > 0) {
        cout << "Re-ranking using: " << LocalParameter::rerankResultsFile << endl;
        ifstream infile((LocalParameter::rerankResultsFile).c_str());
        lvl0Results.load(infile, *idx);
        infile.close();
    }
    
    // cycle through queries
    DocStream *qryStream = new lemur::parse::BasicDocStream(LocalParameter::queryFile);
    qryStream->startDocIteration();
    while (qryStream->hasMore()) {
        Document *qryDoc = qryStream->nextDoc();
        processQuery(qryDoc, resultFile);
    }
    
    outfile.close();
    delete qryStream;
    closeDB();
    
    if(LocalParameter::mapToVisitsOn) {
        results_to_visit_mapper(LocalParameter::reportToVisitMappingFile, LocalParameter::resultFile);
    }
    
    cout << "Retrieval complete, results written to " << resultFileName << endl;
}



