//
//  ConceptLookup.cpp
//  NetworkBasedRetrievalAsInference
//
//  Created by Bevan on 4/02/13.
//  Copyright (c) 2013 CSIRO. All rights reserved.
//

#include <iostream>

map<string, string> concepts;

// SNOMED
string snomed_concept_lookup(string scui) {
    if(concepts.size() == 0) {
        ifstream file("/export/data/ir.collections/snomed/Terminology/Content/sct1_Concepts_Core_INT_20110131.txt.idx");
        string line, id, z, desc, c1, c2, c3;
        while (getline(file, line)) {
            stringstream ss(line);
            getline(ss, id, '\t');
//            getline(ss, z, '\t');
            getline(ss, desc, '\t');
//            getline(ss, c1, '\t');
//            getline(ss, c2, '\t');
//            getline(ss, c3, '\t');
            
            
            concepts.insert(pair<string, string>(id, desc));
        }
    }
    
    if(concepts.count(scui) == 0) {
        return "Unknown SNOMED " + scui;
    } else {
        return concepts.at(scui);
    }
}

// UMLS
string umls_concept_lookup(string cui) {
    if(concepts.size() == 0) {
        ifstream file("/export/data/ir.collections/umls/MRCONSO.RRF.idx");
        string line, id, desc;
        while (getline(file, line)) {
            stringstream ss(line);
            getline(ss, id, '|');
            getline(ss, desc, '|');
            
            transform(id.begin(), id.end(), id.begin(), ::tolower);
            
            concepts.insert(pair<string, string>(id, desc));
        }
    }
    
    if(concepts.count(cui) == 0) {
        return "Unknown UMLS " + cui;
    } else {
        return concepts.at(cui);
    }
}

string concept_lookup(string cui) {
    
    if(cui == "0") 
        return cui;
    
    if(cui.substr(0,1) == "c" || cui.substr(0,1) == "C") {
        return umls_concept_lookup(cui);
    } else {
        return snomed_concept_lookup(cui);
    }
}