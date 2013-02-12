float probCheck(map<TERMID_T, float> P) 
{
    float total = 0;
    for (int i = 1; i <= idx->termCountUnique(); i++)
        total += P[i];
    
    cout << P.size() << "/" << idx->termCountUnique() << " = " << total << endl;
}

map<TERMID_T, float> background_model()
{
    map<TERMID_T, float> M;
    for (int i = 1; i <= idx->termCountUnique(); i++)
    {
        M[i] = (float)idx->termCount( i ) / idx->termCount();
    }
    return M;
}

map<TERMID_T, float> conceptLM(int termId, float lambda, map<TERMID_T, float> &M) {    
    map<TERMID_T, float> theta;
    int totalTerms = 0;
    int totalDocs = 0;
    
    DocInfoList *docInfoList = idx->docInfoList(termId);
    docInfoList->startIteration();
    while (docInfoList->hasMore())
    {
        DocInfo *docInfo = docInfoList->nextEntry();
        TermInfoList *termInfoList = idx->termInfoListSeq(docInfo->docID());
        termInfoList->startIteration();
        while (termInfoList->hasMore())
        {
            TermInfo *termInfo = termInfoList->nextEntry();
            theta[termInfo->termID()]++;
        }
        
        totalTerms += idx->docLength(docInfo->docID());
        totalDocs++;
        
        delete termInfoList;
    }
    delete docInfoList;
    
    // update counts to probability
    for (int i = 1; i <= idx->termCountUnique(); i++)
    {
        float p_term = theta.count(i) > 0 ? (float)theta[i] / totalTerms : 0.0;
        float background = M[i];
        theta[i] = lambda*(p_term) + (1-lambda)*(background);
    }
    
    return theta;
}

///////////////////////////////////////////////////////////////////////////////
// CER model from [1] D. Trieschnigg, et al. 
// 	“Measuring concept relatedness using language models,” SIGIR08
///////////////////////////////////////////////////////////////////////////////

float CER(map<DOCID_T, float> &theta1, map<DOCID_T, float> &M, map<DOCID_T, float> &theta2)
{
    float sum = 0;
    for (int i = 1; i <= idx->termCountUnique(); i++)
    {
        sum += theta2[i] * log(theta1[i]/M[i]);
    }
    return sum;
}

float dCER(map<DOCID_T, float> &theta1, map<DOCID_T, float> &M, map<DOCID_T, float> &theta2)
{
    float cer1 = CER(theta1, M, theta2);
    float cer2 = CER(theta2, M, theta1);
    return ( cer1 + cer2) /2;
}

///////////////////////////////////////////////////////////////////////////////
// standard JSD + KLD for comparison against CER
///////////////////////////////////////////////////////////////////////////////

// Kullback-Leibler divergence
// D(P||Q) = E P(i) * log (  P(i) / Q(i)  )
float kld(map<DOCID_T, float> &c_LM1, map<DOCID_T, float> &c_LM2)
{
    float sum = 0;
    for (int i = 1; i <= idx->termCountUnique(); i++)
    {
        sum += c_LM1[i] * log (c_LM1[i]/c_LM2[i]);
    }
    return sum;
}

// Jensen-Shannon divergence
float jsd(map<DOCID_T, float> &c_LM1, map<DOCID_T, float> &c_LM2)
{
    return -( kld(c_LM1,c_LM2) + kld(c_LM2,c_LM1) ) / 2;
}


///////////////////////////////////////////////////////////////////////////////
// main methods
///////////////////////////////////////////////////////////////////////////////

/*        
        map<DOCID_T, float> theta1 = conceptLM(termId1, lambda, M);
        map<DOCID_T, float> theta2 = conceptLM(termId2, lambda, M);
        
        if(method == "jsd")
        {
            cout << jsd(theta1, theta2) << endl;
        }else {
            cout << dCER(theta1, M, theta2) << endl;
		}


*/