//
// BuildNetworkQuery_Blanco
//
// This program takes a UMLS concept query and creates a graph based on the queries UMLS relationships.
// It uses the C++ Lemon Graph Library.
//
// NOTE: Currently it only supports single UMLS concept queryies.
//
// Created: 26-04-2012 by Bevan.Koopman@csiro.au

#include <iostream>

// lemon
#include <lemon/list_graph.h>

// lemur
#include "common_headers.hpp"
#include "IndexManager.hpp"

using namespace lemon;
using namespace std;
using namespace lemur::api;
using namespace indri::index;

#include "Utils.cpp"
#include "UMLS.cpp"

struct Concept {
    string cui;
    int weight;
} concept;

/////////////////////////////////////////////////////////////
// HELPERS
/////////////////////////////////////////////////////////////

void fetchNode(ListGraph::Node &sourceNode, string &cui, ListGraph* pG, map<string, int> &cuis)
{
    if (cuis.count(cui) > 0) {
        sourceNode = pG->nodeFromId(cuis[cui]);
    } else {
        sourceNode = pG->addNode();
        cuis[cui] = pG->id(sourceNode);
    }
}

void createEdgeIfNotExists(ListGraph::Node &sourceNode, ListGraph::Node &targetNode, ListGraph* pG)
{
    if (sourceNode != targetNode) {
        for (ListGraph::IncEdgeIt ed(*pG, sourceNode); ed != INVALID; ++ed)
        {

            if (pG->target(ed) == targetNode ) {
                return;
            }
        }

        pG->addEdge(sourceNode, targetNode);
    }
}

// Generate Graphviz dot code for the this graph for visualisation
void toDot(ListGraph* pG, map<string, int> &cuis)
{
    cout << "graph G {" << endl;
    for (map<string, int>::iterator it = cuis.begin(); it != cuis.end(); it++) {
        ListGraph::Node n = pG->nodeFromId(cuis[it->first]);
        for (ListGraph::IncEdgeIt e(*pG, n); e != INVALID; ++e)
        {
            ListGraph::Node uNode = pG->source(e);
            ListGraph::Node vNode = pG->target(e);

            cout << "\t" << pG->id(uNode) << " -- " << pG->id(vNode) << ";" << endl;
        }
    }
    cout << "}" << endl;
}


/////////////////////////////////////////////////////////////
// QUERY
/////////////////////////////////////////////////////////////

// void BuildQueryGraph(string queryCui)
// {
//     // setup graph, each node records the CUI associated with that node
//     ListGraph g;
//     ListGraph::NodeMap<std::string> cuis(g);
//
//     // create the root query node
//     ListGraph::Node root = g.addNode();
//     cuis[root] = queryCui;
//
//     // create nodes and edges to related concepts
//     vector<string> relatedConcepts = findRelatedConcepts(cuis[root]);
//     for (vector<string>::iterator iter = relatedConcepts.begin(); iter != relatedConcepts.end(); iter++)
//     {
//         ListGraph::Node v = g.addNode();
//         cuis[v] = *iter;
//         ListGraph::Arc  a = g.addArc(root, v);
//     }
//
//     toDot(&g, cuis, queryCui);
// }

/////////////////////////////////////////////////////////////
// DOCUMENT
/////////////////////////////////////////////////////////////

ListGraph* BuildDocumentGraph(string docName, map<string, int> &cuis, int N=3)
{

    // setup graph, each node records the CUI associated with that node
    ListGraph* g = new ListGraph;

    string doc[] = {"a", "b", "c", "b", "e", "b", "g", "b", "i", "b", "k", "x"};

    string indexLoc = "/export/data/ir.indexes/medtrack/visit_level/umls_representation/indri-5.1/visits_umls_nofields";
    //indexLoc="testdata/index";

    lemur::api::Index *idx = lemur::api::IndexManager::openIndex ( indexLoc );
    int docId = idx->document(docName);

    TermInfoList *termInfoList = idx->termInfoListSeq(docId);
    TermInfo *term;
    termInfoList->startIteration();
    list<string> window;
    while (termInfoList->hasMore())
        //for (int i = 0; i < sizeof(doc)/sizeof(doc[0]); i++)
    {
        term = termInfoList->nextEntry();
        string termStr = idx->term(term->termID());
        //string termStr = doc[i];
        //cout << "current: " << termStr << endl;

        for (list<string>::iterator iter = window.begin(); iter != window.end(); iter++)
        {
            string related = *iter;
            if (termStr != related) // avoid self edges
            {
                ListGraph::Node sourceNode;
                fetchNode(sourceNode, termStr, g, cuis);

                ListGraph::Node targetNode;
                fetchNode(targetNode, related, g, cuis);

                createEdgeIfNotExists(sourceNode, targetNode, g);
            }
        }
        window.push_back(termStr);
        //cout << endl;

        if (window.size() == N)
        {
            //cout << "remove: " << window.front() << endl;
            window.pop_front();
        }
    }

    return g;
}

////////////////////////////////////////////////////////////////
// WEIGHTING
////////////////////////////////////////////////////////////////

map<string, float>* CalculateNodeWeights(ListGraph* pG, map<string, int> &cuis, int iterations=1)
{
    // initialise all weights to 1
    ListGraph::NodeMap<float> weights(*pG, 1);

    for (int i = 0; i < iterations; i++) // iterate until (hopefully) convergence
    {
        // each CUI in the graph
        for (map<string, int>::iterator iter = cuis.begin(); iter != cuis.end(); ++iter)
        {
            ListGraph::Node node = pG->nodeFromId(iter->second);

            float weight = 0;

            //cout << iter->first << "_" << pG->id(node) << " (" << weights[node] << ") : ";

            // for each in edge
            for (ListGraph::IncEdgeIt edge(*pG, node); edge != INVALID; ++edge)
            {
                ListGraph::Node sourceNode = pG->source(edge);

                float numOutEdges = 0;
                for (ListGraph::OutEdgeIt sourceEdge(*pG, sourceNode); sourceEdge != INVALID; ++sourceEdge)
                    numOutEdges++;

                weight += weights[node] / numOutEdges;
                //cout << pG->id(sourceNode) << ":" << "(" <<  weights[sourceNode] <<"/"<< numOutEdges << "=" << weights[sourceNode] / numOutEdges << ")" << ", ";

            }
            //cout << "-> " << weight << endl;
            weights[node] = weight;
        }
    }


    map<string, float> cuiToWeight;
    for (map<string, int>::iterator it = cuis.begin(); it!=cuis.end(); it++)
    {
        cuiToWeight[it->first] = weights[pG->nodeFromId(it->second)];
    }

    return &cuiToWeight;
}

////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: retrieval_inference_using_umls <CUI>" << endl;
        exit(1);
    }
    string arg = argv[1];
    //BuildQueryGraph(arg);

    string doc = "corpus/8EmDAa+QdRvH";
    //doc = "docs/2.txt";

    cout << "Building graph for " << doc << endl;
    map<string, int> cuis;
    ListGraph* pGraph = BuildDocumentGraph(doc, cuis, string_to_int(arg));
    
    cout << "Calulating weights..." << endl;
    map<string, float> *pWeights = CalculateNodeWeights(pGraph, cuis);

    //toDot(pGraph, cuis);


    for (map<string, float>::iterator it = pWeights->begin(); it!=pWeights->end(); it++) {
        cout << it->first << ": " << it->second << endl;
    }
}
