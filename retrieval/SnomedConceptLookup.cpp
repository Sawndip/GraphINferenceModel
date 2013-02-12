//
//  SnomedConceptLookup.cpp
//  NetworkBasedRetrievalAsInference
//
//  Created by Bevan on 4/02/13.
//  Copyright (c) 2013 CSIRO. All rights reserved.
//

#include <iostream>

map<int, string> concepts;
string concept_lookup(int scui) {
    if(concepts.size() == 0) {
        ifstream file("/Users/bevan/ir.collections/snomed/Terminology/Content/sct1_Concepts_Core_INT_20110131.txt.idx");
        string line, idStr, z, desc, c1, c2, c3;
        while (getline(file, line)) {
            stringstream ss(line);
            getline(ss, idStr, '\t');
//            getline(ss, z, '\t');
            getline(ss, desc, '\t');
//            getline(ss, c1, '\t');
//            getline(ss, c2, '\t');
//            getline(ss, c3, '\t');
            
            int id = atoi(idStr.c_str());
            
            concepts.insert(pair<int, string>(id, desc));
        }
    }
    
    return concepts.at(scui);
}