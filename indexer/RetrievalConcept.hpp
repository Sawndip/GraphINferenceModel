//
//  RetrievalConcept.h
//
//
//  Created by Bevan.Koopman@csiro.au on 9/05/12.
//
//  This class represents a single UMLS concept for retrieval purposes. It captures the weight of the
//  concept and other characteristics like semantic type

#ifndef _RetrievalConcept_hpp
#define _RetrievalConcept_hpp

#include <iostream>
#include <map>
#include <vector>

using namespace std;

typedef int DOCID;
typedef string CUI;

class RetrievalConcept
{

public:
    RetrievalConcept();
    RetrievalConcept(CUI cui);
    ~RetrievalConcept();
    
    CUI cui();
    float idf();
    void setIdf(float idf);
    
    void addDocument(DOCID docId, float conceptFreqNorm);
    map<DOCID, float> conceptFreqNorm();
    
    float weight(DOCID docId);
    
private:
    CUI itsCui;
    float itsIdf;
    map<DOCID,float> itsDocConceptFreq;
};

#endif
