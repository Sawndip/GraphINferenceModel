//
// Helper program that searches the sqlite DB that contains cui to cui mapping of related (UMLS or snomed) concepts
//
// Reads a SQLite DB of related concepts which has previously been populated by either:
//  retrieval_inference_using_umls/umls_relations/      - for UMLS relationships
//  retrieval_inference_using_umls/snomed_relations/    - for SNOMED relationships

#include <iostream>
#include <map>
#include <vector>
#include <iostream>
#include <sqlite3.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// SQLite access for concept relationship
///////////////////////////////////////////////////////////////////////////////

sqlite3 *database;
map<string, vector< vector<string> > > relatedConceptsCache;

void prepareDB(string dbFile)
{
    sqlite3_open(dbFile.c_str(), &database);
}

void closeDB()
{
    sqlite3_close(database);
}

// searches the sqlite DB that contains cui to cui mapping of related concepts
// read each from sql
vector< vector<string> > findRelatedConceptsSQL(const string &cui)
{
    
    std::vector< vector<string> > relatedConcepts;
    
   
    sqlite3_stmt *statement;
//    string query = "SELECT cui2, reltype FROM crel where cui1 = '"+cui+"';";
    string query = "SELECT cui1, reltype FROM crel where cui2 = '"+cui+"' and (relcharacteristic=0 or relcharacteristic=3);";

    if (sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
    {
        int cols = sqlite3_column_count(statement);
        int result = 0;
        while (true)
        {
            result = sqlite3_step(statement);

            if (result == SQLITE_ROW)
            {
                vector<string> relationship;
                for (int col = 0; col < cols; col++)
                {
                    char *val = (char*)sqlite3_column_text(statement, col);
                    relationship.push_back(val);
                }
                relatedConcepts.push_back(relationship);
            }
            else
            {
                break;
            }
        }

        sqlite3_finalize(statement);
    }

    string error = sqlite3_errmsg(database);
    if (error != "not an error") cout << "Error" << " " << error << endl;
    
    //sqlite3_close(database);

    return relatedConcepts;
}

// looks ALL the relationships from the DB into a map for quicker access
vector< vector<string> > findRelatedConceptsCache(const string &cui)
{
    if(relatedConceptsCache.size() == 0)
    {
        sqlite3_stmt *statement;
        string query = "SELECT * FROM crel;";
        
        if (sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
        {

            int result = 0;
            while (true)
            {
                result = sqlite3_step(statement);
                
                if (result == SQLITE_ROW)
                {
                    char *cui1 = (char*)sqlite3_column_text(statement, 0);
                    char *cui2 = (char*)sqlite3_column_text(statement, 1);
                    char *reltype = (char*)sqlite3_column_text(statement, 2);
                    
                    vector< vector<string> > mappings;
                    if (relatedConceptsCache.count(cui1) > 0)
                    {
                        mappings = relatedConceptsCache[cui1];
                    }
                    vector<string> entry;
                    entry.push_back(cui2);
                    entry.push_back(reltype);
                    mappings.push_back(entry);
                    relatedConceptsCache[cui1] = mappings;
                }
                else
                {
                    break;
                }
            }
            
            sqlite3_finalize(statement);
        }
        
        string error = sqlite3_errmsg(database);
        if (error != "not an error") cout << "Error" << " " << error << endl;
        
        sqlite3_close(database);
    }
    
    return relatedConceptsCache[cui];
}

vector<string> getAllConcepts()
{
    vector<string> concepts;
    
    sqlite3_stmt *statement;
    string query = "select distinct c from (select cui1 as c from crel union select cui2 as c from crel);";
    
    if (sqlite3_prepare_v2(database, query.c_str(), -1, &statement, 0) == SQLITE_OK)
    {
        
        int result = 0;
        while (true)
        {
            result = sqlite3_step(statement);
            
            if (result == SQLITE_ROW)
            {
                char *cui = (char*)sqlite3_column_text(statement, 0);
                concepts.push_back(cui);
            }
            else
            {
                break;
            }
        }
        
        sqlite3_finalize(statement);
    }
    
    return concepts;

}


// searches the sqlite DB that contains cui to cui mapping of related concepts
vector< vector<string> > findRelatedConcepts(const string &cui)
{
    return findRelatedConceptsSQL(cui);
}


///////////////////////////////////////////////////////////////////////////////
// Weight of different types of SNOMED relationships (ISA, PartOf, ProcedureSite, etc.)
///////////////////////////////////////////////////////////////////////////////

// lookup the weight assinged to this snomed relationship type, done by reading the
// reltypes.txt file.
map<string, double> snomedRelationshipTypeWeightMap;
double snomedRelationshipTypeWeight(string reltype, string relWeightsFile)
{
    // populate the weights on the first call
    if(snomedRelationshipTypeWeightMap.size() == 0)
    {
        ifstream file(relWeightsFile.c_str());
        
        LOG_IF(FATAL, !file.good()) << "Error: snomed relationship file " << relWeightsFile << " not found";
        
        string line;
        while(getline(file, line))
        {
            stringstream linestream(line);
            string reltype, name, weightStr;
            getline(linestream, reltype, '\t');
            getline(linestream, name, '\t');
            getline(linestream, weightStr, '\t');
            float weight = string_to_float(weightStr.c_str());
            
            snomedRelationshipTypeWeightMap[reltype] = weight;
        }
        file.close();
    }
    
    //cout << "looking up weight for " << reltype << " = " << snomedRelationshipTypeWeightMap[reltype] << endl;
    
    return snomedRelationshipTypeWeightMap[reltype];
}

double snomedRelationshipTypeWeight(string reltype) {
    return snomedRelationshipTypeWeight(reltype, "../snomed_relations/reltypes.txt");
}

///////////////////////////////////////////////////////////////////////////////
// SNOMED concept Semantic Types 
///////////////////////////////////////////////////////////////////////////////

// lookup the semantic type of a given concept, done by reading the 
// snomed_semantic_types.txt file.
map<string, string> snomedSemanticTypesMap;
string snomedSemanticType(string cui)
{
    if(snomedSemanticTypesMap.size() == 0)
    {
        ifstream file("../snomed_relations/snomed_semantic_types.txt");
        LOG_IF(FATAL, !file.good()) << "Error: snomed semantic type file snomed_semantic_types.txt not found";
        
        string line;
        while (getline(file, line)) {
            stringstream linestream(line);
            string id, type;
            getline(linestream, id, '\t');
            getline(linestream, type, '\t');
            
            snomedSemanticTypesMap[id] = type;
        }
    }
    
    if(snomedSemanticTypesMap.count(cui) > 0)
    {
        return snomedSemanticTypesMap[cui];
    } else {
        return "unknown";
    }
}

// weigting based on semantic types of two concepts
float semanticType_weight_snomed(string sourceCUI, string targetCUI, float weight_boost) {
    float weight = 1; // i.e. no change
    string sourceType = snomedSemanticType(sourceCUI);
    string targetType = snomedSemanticType(targetCUI);
    
    if (sourceType == targetType) weight = weight_boost;
    
    VLOG(2) << "semanticType_weight_snomed between " << sourceCUI << " (" << sourceType << ") and " << targetCUI << " (" << targetType << ") is: " << weight << endl; 
    
    return weight;
}


///////////////////////////////////////////////////////////////////////////////
// END
///////////////////////////////////////////////////////////////////////////////



