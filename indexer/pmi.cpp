// PMI

set<int> docSet(int termId)
{
    set<int> docSet;
    
    DocInfoList* docList = idx->docInfoList(termId);
    DocInfo* doc;
    docList->startIteration();
    while (docList->hasMore())
    {
        doc = docList->nextEntry();
        docSet.insert(doc->docID());
    }
    
    delete docList;
    
    return docSet;
}


//
// PMI(x,y) = P(x,y) / p(x) * P(y)
//
float PMI(int x, int y)
{
    set<int> xSet = docSet(x);
    set<int> ySet = docSet(y);
    
    if(xSet.size() == 0 || ySet.size() == 0)
        return 0.0;
    
    set<int> intersect;
    double totalDocs = idx->docCount();
    
    set_intersection( xSet.begin(), xSet.end(), ySet.begin(), ySet.end(),  insert_iterator< std::set<int> >( intersect, intersect.begin() ));
    
    //cout << intersect.size() << " " << totalDocs << " " << xSet.size() << " " << ySet.size() << " | ";
    
    float score = (intersect.size() / totalDocs)
            / 
        (xSet.size()*ySet.size() / totalDocs);
    
    return score;
}

float PositivePMI(int x, int y)
{
    float pmi = PMI(x, y);
    return pmi > 0 ? pmi : 0.0;
}