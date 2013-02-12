#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

map<string, string> read_report_mappings(string mapping_file) {
    map<string, string> report_mappings;

    ifstream fh(mapping_file.c_str());
    string line;
    while(getline(fh, line)) {
        stringstream ss(line);
        string report, visit;
        getline(ss, report, ' ');
        getline(ss, visit, ' ');

        report = "corpus/"+report;

        if("NULL" != visit) {
            report_mappings[report] = "corpus/" + visit;
            //cout << report << " -> " << report_mappings[report] << endl;
        }
    }
    fh.close();

    return report_mappings;
}

void results_to_visit_mapper(string mapping_file, string results_file, bool update=true) {
     map<string, string> report_mappings = read_report_mappings(mapping_file);
     map<string, int> seenVisits; // visits processed for this qId

     ofstream outfile((results_file+".tmp").c_str());

    ifstream fh(results_file.c_str());
    string line;
    string prevQId = "";
    string qid, q0, doc, pos, score, runId;
    int count = 1;
    while(getline(fh, line)) {
        stringstream ss(line);
        getline(ss, qid, ' ');
        getline(ss, q0, ' ');
        getline(ss, doc, ' ');
        getline(ss, pos, ' ');
        getline(ss, score, ' ');
        getline(ss, runId, ' ');

        if(prevQId != qid) 
            seenVisits.clear();

        if(report_mappings.count(doc) > 0 && seenVisits.count(report_mappings[doc]) == 0) {
            outfile << qid << " " << q0 << " " << report_mappings[doc] << " " << count << " " << score << " " << runId << endl;
            seenVisits[report_mappings[doc]]++;
            count++;
        }

        prevQId = qid;
    }

    fh.close();
    outfile.close();

    if(update)
        rename((results_file+".tmp").c_str(), (results_file).c_str());
}

