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
ListDigraph::ArcMap<int> *reltypeArcs;
map<Node, map<DOCID_T, double> > diffusionFactors;

#include "../indexer/doccosine.cpp"

// tuning params
float related_dampener = 0.6;
float semantic_type_weight_booster = 2;
double s = 0;

// debug maps
map<DOCID_T, double> lvl0Scores;
map<DOCID_T, double> lvl1Scores;

// rerank list
ResultFile lvl0Results(true);
IndexedRealVector *rerankResults;

// stopword list
set<string> stopwordList;

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
    
    std::string rerankResultsFile; // path to a TREC results file to use for reranking
    int rerankCount; // number of top docs to rerank
    
    void get() {
        // the string with quotes are the actual variable names to use for specifying the parameters
        indexPath		= ParamGetString("index");
        queryFile		= ParamGetString("query");
        resultCount		= ParamGetInt("resultCount", 1000);
        resultFile		= ParamGetString("resultsFile", LocalParameter::queryFile+".results");
        stopwordFile	= ParamGetString("stopwordFile");
        
        // calc prior: mleidf, tfidf, LM-JM, LM-Dirichlet
        weightScheme	= ParamGetString("weightScheme","LM-JM");
        LMjm_lambda		= ParamGetDouble("LM-JM.lambda", 0.9);
        LMdirichlet_mu	= ParamGetInt("LM-Dirichlet.mu", 2000);
        indri_tfidf_k1		= ParamGetDouble("indri_tfidf_k1", 0.9);
        indri_tfidf_b	= ParamGetInt("indri_tfidf_b", 0.5);
        
        depth           = ParamGetInt("depth", 1);
        diffusionMix    = ParamGetDouble("diffusionMix", 0.5);
        
        rerankResultsFile   = ParamGetString("rerankResultsFile");
        rerankCount     = ParamGetInt("rerankCount", 10);
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

// language model with Jelinek-Mercer Smoothing
double LM_JM(TERMID_T termId, DocInfo *docInfo, double lambda)
{
    return lambda*(docInfo->termCount()/(float)idx->docLength(docInfo->docID())) + (1-lambda)*(idx->termCount(termId) / idx->termCount());
}

// language model with Dirichlet Smoothing
double LM_Dirichlet(TERMID_T termId, DocInfo *docInfo, int mu)
{
    return log(
               ( (double)docInfo->termCount() + mu*(idx->termCount(termId) / (double)idx->termCount()) )
               /
               (double)( idx->docLength(docInfo->docID()) + mu)
           );
}

// mu *  tf_i_C
//       ------
//        |C|
// -------------
//  |D| + mu
double LM_Dirichlet_Background(TERMID_T termId, int docLength) {
    return (LocalParameter::LMdirichlet_mu * ( idx->termCount(termId) / (double)idx->termCount() ))
            /
            ( (double)LocalParameter::LMdirichlet_mu + docLength );
}

double calcEmptyDocLMScore(vector<string> queryConcepts, int docLength) {
    double score = 0;
    for (vector<string>::iterator it = queryConcepts.begin(); it != queryConcepts.end(); it++) {
        score += log(LM_Dirichlet_Background(idx->term(*it), docLength));
    }
    return score;
}

// return the prior probability of this node
double prior(Node &node, DocInfo *docInfo)
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
        return LM_Dirichlet(termId, docInfo, LocalParameter::LMdirichlet_mu);
    }
    cout << "ERROR: unknown weighting scheme " << LocalParameter::weightScheme << endl;
}


///////////////////////////////////////////////////////////////////////////////
// Diffusion Factor and Edge Weighting
///////////////////////////////////////////////////////////////////////////////


// returns a score based on priors of nodes related to the supplied node
double related_prior(Node &node)
{
    double theRelatedPrior = 1;
    for (ListDigraph::InArcIt arc(*g, node); arc != INVALID; ++arc)
    {
        theRelatedPrior++;
    }
    
    VLOG(2) << "related_prior for " << (*nodeToCUIMap)[node] << " " << 1 + theRelatedPrior;
    
    return theRelatedPrior;
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
    float alpha = LocalParameter::diffusionMix;
    double total_df = 1;
    
    for (vector<Arc>::iterator it = pathToQuery.begin(); it != pathToQuery.end(); it++)
    {
        Arc arc = *it;
        Node source = g->source(arc);
        Node target = g->target(arc);
        
        // the semantic sim weight
        double sim_weight = (*simArcs)[arc];
        //double sim_weight = query_sim(source, queryConcepts);
        if(sim_weight == 0) sim_weight = 1/idx->docCount();
        
        // the relationship type weigth
        int reltype = (*reltypeArcs)[arc];
        double rel_weight = snomedRelationshipTypeWeight(int_to_str(reltype));
        
        // combined weight - the df
        double df = alpha * sim_weight + (1-alpha)*rel_weight;
        
        VLOG(2) << "diffusion Factor for arc " << (*nodeToCUIMap)[source] << " -> " << (*nodeToCUIMap)[target] << ": avg(" << sim_weight << "," << reltype << "=" << rel_weight << ") =" << df;
        
        total_df = total_df * df;
    }
    
    return total_df;
}


///////////////////////////////////////////////////////////////////////////////
// Comparison against qrels
///////////////////////////////////////////////////////////////////////////////

// isJudged checks if the doc id is contain in the qrels
// Remember, judged docs are either 2 (very relevant), 1 (relevant) and 0 (irrelevant).
// This method can be used to give and idea of the effect of unjudged documents on retrieval performance.
map<int, string> judged;
bool isJudged(int docId)
{
    if(judged.size() == 0)
    {
        ifstream file("../judged.txt");
        string line;
        while(getline(file, line))
        {
            judged[idx->document(line)] = line;
            cout << idx->document(line) << " " << judged[idx->document(line)] << endl;
        }
        file.close();
    }
    
    return (judged.count(docId) > 0);
}

map<int, vector<DOCID_T> > qrels; // queryId -> [docId, ..]
int currentQueryId;
bool isRelevantDoc(int queryId, DOCID_T docId) {
    if(qrels.size() == 0) {
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
            
//            cout << line << endl;
//            cout << qId << " " << zero << " " << docName << " " << relevance << " :" << docId << endl;
//            
            if(relevance == "0") continue;
            
            vector<DOCID_T> docIds;
            if(qrels.count(qId)) {
                docIds = qrels[qId];
//                cout << "previous size is " << docIds.size() << endl;
            }
            docIds.push_back(docId);
            qrels[qId] = docIds;
        }
        file.close();
        
//        cout << "FINISHED QRELS 154=" << qrels[154].size() << endl;
    }
    
    bool isRel = std::find ( qrels[queryId].begin(), qrels[queryId].end(), docId ) != qrels[queryId].end();
    //VLOG(2) << "query " << queryId << " contains " << qrels[queryId].size() << " relevant documents " << idx->document(docId) << " " << isRel << endl;
    return isRel;
    
}

///////////////////////////////////////////////////////////////////////////////
// Graph traversal functions
///////////////////////////////////////////////////////////////////////////////

// the visit node function - for each document found at this node score the document based on its prior and its diffusionFactor
void handleNode(Node &node, vector<Arc> &pathToQuery, map<DOCID_T, double> &scores, vector<string> queryConcepts)
{
    
    // the term does not appear in the collection, there no document to score against it
    // note that the traverse method will still pass through this non-collection node
    TERMID_T termId = idx->term((*nodeToCUIMap)[node]);
    VLOG(1) << "hn: " << (*nodeToCUIMap)[node] << " #docs:" << idx->docCount(termId) << endl;
    if (termId == 0)
        return;
    
    // node has previous been visited via some other path?
//    bool previously_visited = diffusionFactors.count(node) > 0;
//    if (!previously_visited) {
//        map<DOCID_T, double> e;
//        diffusionFactors[node] = e;
//    }
    
    double theDiffusionFactor = diffusionFactor(pathToQuery);
    double theRelatedPrior = related_prior(node);
    double theTopoScore = 1 /(*nodeToTopoScoreMap)[node];
    
    VLOG(2) << " toposcore: " << (*nodeToCUIMap)[node] << " " << theTopoScore << endl;
    int relevantDocCount = 0;
    
    DocInfoList *docInfoList = idx->docInfoList(termId);
    docInfoList->startIteration();
    while (docInfoList->hasMore())
    {
        DocInfo *doc = docInfoList->nextEntry();
        //cout << "d:" << doc->docID() << " ";
        
        // MEDTRACK 2012 FILTER
        //if(idx->document(doc->docID()) == "corpus/0VppFYq+9Tjm")
        //    continue;
        
//        if(idx->document(doc->docID()) != "corpus/jXvsweEykldF")
//            continue;
        
        //VLOG(2) << "score for corpus/jXvsweEykldF: 1) " << scores[doc->docID()] << " ";
        
        // first if the doc has not been scored, push on the emptyDocScore
        VLOG(2) << "count for doc is " << scores.count(doc->docID()) << endl;
        if(scores.count(doc->docID()) == 0) {
            scores[doc->docID()] += calcEmptyDocLMScore(queryConcepts, idx->docLength(doc->docID()));
        }
        
        VLOG(2) << "2) " << scores[doc->docID()] << endl;
        
        // take off the background smoothing for this doc
        VLOG(2) << "background for " << idx->term(termId) << " " <<  log(LM_Dirichlet_Background(termId, idx->docLength(doc->docID()))) << endl;
        scores[doc->docID()] -= log(LM_Dirichlet_Background(termId, idx->docLength(doc->docID())));
        
        VLOG(2) << "3) " << scores[doc->docID()] << endl;
        
        double thePrior = prior(node, doc);
        double theScore = thePrior * theDiffusionFactor;  //* theRelatedPrior; // * theTopoScore *  ((float)conceptScore / 1000) * ( 1 + log(conceptCount) ) ;
        // apply dampening for non-query nodes
        if(pathToQuery.size() > 0) {
            //theScore = related_dampener * theScore;
        }
        // apply proposition count
        //theScore = theScore * ((float)conceptScore / 1000) * ( 1 + log(conceptCount) );
        
        string relevanceIndicator = " ";
        if(isRelevantDoc(currentQueryId, doc->docID())) {
            relevanceIndicator = "*";
            relevantDocCount++;
        }
        
        VLOG(1) << relevanceIndicator << "RSV(" << (*nodeToCUIMap)[node] << "|" << idx->document(doc->docID()) << ") = prior:" << thePrior << " * df_" << pathToQuery.size() << ":" << theDiffusionFactor << "\t=> " << theScore << endl;
        
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
        scores[doc->docID()] += theScore;
        
        VLOG(2) << "4) " << scores[doc->docID()] << endl;
        
        
        //diffusionFactors[node][doc->docID()] = theDiffusionFactor;
        
        // DEBUG TO SEE HOW MUCH SCORE IS CHANGING
        if(pathToQuery.size() == 0) {
            lvl0Scores[doc->docID()] += theScore;
        } else {
            lvl1Scores[doc->docID()] += theScore;
        }
        
        
    }
    
    VLOG(1) << "complete hn: " << (*nodeToCUIMap)[node] << " #docs:" << idx->docCount(termId) << " (" << relevantDocCount << " relevant)" << endl;
    
    delete docInfoList;
} 

// recursive traversal of the graph
int traverse(Node node, vector<Arc> &pathToQuery, map<DOCID_T, double> &scores, int &depth, vector<string> &queryConcepts) {
    int nodeCount = 0;
    // avoid cycles
    for (vector<Arc>::iterator it = pathToQuery.begin(); it != pathToQuery.end(); it++)
    {
        if (g->target(*it) == node) {
            //cout << "loop found: node:" << (*nodeToCUIMap)[node] << " arc " << (*nodeToCUIMap)[g->source(*it)] << " -> " << (*nodeToCUIMap)[g->target(*it)] << endl;
            return 0;
        }
    }
    
    //cout << " t: " << (*nodeToCUIMap)[node] << " " ;
    //cout << "num: " << numActiveQueryInEdges.size() << endl;
    
    if (pathToQuery.size() <= depth)
    {
        //cout << "calling traverse on " << (*nodeToCUIMap)[node] << endl;
        
        // score all the docs attached to this node
        handleNode(node, pathToQuery, scores, queryConcepts);
        nodeCount++;
        
        // follow children of this node
        for (ListDigraph::InArcIt arc(*g, node); arc != INVALID; ++arc)
        {
            VLOG(1) << " arc: " << (*nodeToCUIMap)[g->source(arc)] << " -> " << (*nodeToCUIMap)[g->target(arc)];
            if (g->source(arc) != node) {
                pathToQuery.push_back(arc);
                nodeCount += traverse(g->source(arc), pathToQuery, scores, depth, queryConcepts);
                pathToQuery.pop_back();
            }
        }
    }
    return nodeCount;
}

///////////////////////////////////////////////////////////////////////////////
// Query processing functions
///////////////////////////////////////////////////////////////////////////////

// process a single query term/concept by locating the node in the graph and initiating traversal from that node
void processQueryTerm(string query, map<DOCID_T, double> &scores, vector<string> &queryConcepts)
{
    transform(query.begin(), query.end(), query.begin(), ::tolower);
    
    // stopword? skip
    if(stopwordList.find(query) != stopwordList.end())
        return;
    
    // query reduction by snomed connectedness (a la ADCS'2)
//    int theSnomedRelCount = (float)findRelatedConcepts(query).size();
//    if(theSnomedRelCount == 1)
//        return;
    
    cout << " " << query << flush;
    
    diffusionFactors.clear();
    
    if(cuiToNodeId.count(query) > 0) {
        Node queryNode = g->nodeFromId(cuiToNodeId[query]);
        vector<Arc> pathToQuery;
        
        // traverse the graph from this node
        int nodeCount = traverse(queryNode, pathToQuery, scores, LocalParameter::depth, queryConcepts);
        cout << "(" << nodeCount << ") ";
        
        if (idx->term(query) == 0) {
            cout << " [not found in index]" << " ";
        }
    } else {
        cout << " [not present in graph, ignoring]" << " ";
    }
}

// process a single query and score the results to the trec results file
void processQuery(Document *qryDoc, ResultFile &resultFile)
{
    const char *queryID = qryDoc->getID();
    currentQueryId = atoi(queryID);
    cout << "Running query: "<< queryID << "\n\t";
    
    // rerank another results file
    if(LocalParameter::rerankResultsFile.length() > 0) {
        lvl0Results.findResult(qryDoc->getID(), rerankResults);
    }
    
    // used only for debug / statistics
    lvl0Scores.clear();
    lvl1Scores.clear();
    
    
    // build up the query vector
    vector<string> queryConcepts;
    qryDoc->startTermIteration();
    while (qryDoc->hasMore()) {
        const Term *qryTerm = qryDoc->nextTerm();
        queryConcepts.push_back(qryTerm->spelling());
        
    }

    
    map<DOCID_T, double> scores;
    
    // EXECUTE RETRIEVAL ON THE QUERY TERM
    for (vector<string>::iterator it = queryConcepts.begin(); it != queryConcepts.end(); it++) {
        processQueryTerm(*it, scores, queryConcepts);
    }
    
    
    cout << endl;
    
    // print results    
    IndexedRealVector results(idx->docCount());
    results.clear();
    
    if(LocalParameter::rerankResultsFile.length() == 0) {
        for (map<DOCID_T, double>::iterator it=scores.begin(); it!=scores.end(); it++) {
            results.PushValue(it->first, it->second / queryConcepts.size());                
        }
    } else {    // do reranking
        results = *rerankResults;
        int count = 0;
        for(IndexedRealVector::iterator it = rerankResults->begin(); it!=rerankResults->end() && count < LocalParameter::rerankCount; it++) {
            results.IncreaseValueFor((*it).ind, (scores[(*it).ind] / queryConcepts.size())-(*it).val);
            count++;
        }
        delete rerankResults;
    }
    
    results.Sort();
    resultFile.writeResults(qryDoc->getID(), &results, LocalParameter::resultCount);
    
    //printDebugScores(queryId)
    
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
    cout << "depth: " << LocalParameter::depth << " ";
    cout << "weighting: " << LocalParameter::weightScheme << " ";
    cout << "dampener: " << related_dampener << " ";
    cout << "semtype boost: " << semantic_type_weight_booster << " ";
    cout << "smoothing: " << s << " ";
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
    
    cout << endl;
}

///////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////

void AppMain(int argc, char *argv[])
{
    google::InitGoogleLogging(argv[0]);
    
    if (argc > 3) {
        s= atof(argv[3]);
    }
    
    if (argc > 2) {
        LocalParameter::depth = string_to_int(argv[2]);
        if(LocalParameter::depth < 0 || LocalParameter::depth > 100)
        {
            cerr << "Invalid setting for depth: " << LocalParameter::depth << endl;
            exit(1);
        }
    }
    if (argc > 3) {
        related_dampener = string_to_float(argv[3]);
        if(related_dampener < 0 || related_dampener > 1)
        {
            cerr << "Invalid setting for related_dampener: " << related_dampener << endl;
            exit(1);
        }
    }
    if (argc > 4) {
        semantic_type_weight_booster = string_to_float(argv[4]);
        if(semantic_type_weight_booster < 0 || semantic_type_weight_booster > 100)
        {
            cerr << "Invalid setting for semantic_type_weight_booster: " << semantic_type_weight_booster << endl;
            exit(1);
        }
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
    ListDigraph::ArcMap<int> theRelTypeArcs(*g);
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
    ofstream outfile((LocalParameter::resultFile).c_str());
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
    
    cout << "Retrieval complete, results written to " << LocalParameter::resultFile << endl;
}



