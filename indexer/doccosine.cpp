
// returns a map of docId -> termFreq for a given term
map<int, double> documentVector(int termId)
{
    map<int, double> docVector;
    
        DocInfoList* docList = idx->docInfoList( termId );
        DocInfo* doc;
        docList->startIteration();
        while (docList->hasMore())
        {
            doc = docList->nextEntry();
            int docId = doc->docID();
            int tf = doc->termCount();
            double idf = log(idx->docCount() / (double)idx->docCount(termId));

            // mle
            docVector[docId] = tf / (float) idx->docLength(docId);
            // tf-idf
            //docVector[docId] = tf / idf;
        }
        delete docList;
    
    return docVector;
}

// normalise a term vector to unit length norm(v) = sqrt( sum(vi^2) )
double norm(map<int, double> v)
{
    double norm = 0;
    map<int, double>::iterator it;
    for (it = v.begin(); it != v.end(); it++)
    {
        norm += pow(it->second, 2);
    }
    
    return sqrt(norm);
}

// calculate the cosine similarity between two document vectors cosine(x,y) = x . y / |x|*|y|
double cosine(map<int, double> x, map<int, double> y)
{
    double dotProd = 0;
    
    map<int,double>::iterator itx;
    for (itx=x.begin() ; itx != x.end(); itx++ )
    {
        int docId = itx->first;
        if (y.count(docId)) {
            dotProd += x[docId]*y[docId];
        }
    }
    
    return dotProd / ( norm(x)*norm(y) );
}

double cosine(int termIdOne, int termIdTwo)
{
    return cosine(documentVector(termIdOne), documentVector(termIdTwo));
}

double cosine(int termIdOne, int *termIds)
{
    return 0;
}