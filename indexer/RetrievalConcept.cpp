#include "RetrievalConcept.hpp"

RetrievalConcept::RetrievalConcept() {}

RetrievalConcept::RetrievalConcept(CUI cui)
{
    itsCui = cui;
}

RetrievalConcept::~RetrievalConcept()
{

}

// idf
void RetrievalConcept::setIdf(float idf) { itsIdf = idf; }
float RetrievalConcept::idf() { return itsIdf; }

// CUI
string RetrievalConcept::cui() {
    return itsCui;
}

void RetrievalConcept::addDocument(DOCID docId, float conceptFreqNorm)
{
    itsDocConceptFreq[docId] = conceptFreqNorm;
}

map<DOCID, float> RetrievalConcept::conceptFreqNorm() { return itsDocConceptFreq; }

// weight
float RetrievalConcept::weight(DOCID docId)
{
    return idf() * itsDocConceptFreq[docId];
}

/*
int main() 
{
    RetrievalConcept patient("C123");
    patient.setIdf(10000/20);
    patient.addDocument(1, 5);
    
    cout << patient.cui() << " " << patient.weight(1) << endl;
    
    return 0;
}
*/