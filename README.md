# Graph INference (GIN) Model

Implementation of:

B. Koopman, G. Zuccon, P. Bruza, L. Sitbon, and M. Lawley. Information retrieval as semantic inference: A graph inference model applied to medical search. Information Retrieval, 19(1):6–37, 2015. [https://link.springer.com/article/10.1007/s10791-015-9268-9](https://link.springer.com/article/10.1007/s10791-015-9268-9)

> This paper presents a Graph Inference retrieval model that integrates structured knowledge resources, statistical information retrieval methods and inference in a unified framework. Key components of the model are a graph-based representation of the corpus and retrieval driven by an inference mechanism achieved as a traversal over the graph. The model is proposed to tackle the semantic gap problem—the mismatch between the raw data and the way a human being interprets it. We break down the semantic gap problem into five core issues, each requiring a specific type of inference in order to be overcome. Our model and evaluation is applied to the medical domain because search within this domain is particularly challenging and, as we show, often requires inference. In addition, this domain features both structured knowledge resources as well as unstructured text. Our evaluation shows that inference can be effective, retrieving many new relevant documents that are not retrieved by state-of-the-art information retrieval models. We show that many retrieved documents were not pooled by keyword-based search methods, prompting us to perform additional relevance assessment on these new documents. A third of the newly retrieved documents judged were found to be relevant. Our analysis provides a thorough understanding of when and how to apply inference for retrieval, including a categorisation of queries according to the effect of inference. The inference mechanism promoted recall by retrieving new relevant documents not found by previous keyword-based approaches. In addition, it promoted precision by an effective reranking of documents. When inference is used, performance gains can generally be expected on hard queries. However, inference should not be applied universally: for easy, unambiguous queries and queries with few relevant documents, inference did adversely affect effectiveness. These conclusions reflect the fact that for retrieval as inference to be effective, a careful balancing act is involved. Finally, although the Graph Inference model is developed and applied to medical search, it is a general retrieval model applicable to other areas such as web search, where an emerging research trend is to utilise structured knowledge resources for more effective semantic search.


Contains the following subdirectories:

## `indexer`
 
* Creates the GIN graph and serialises it.
* **Inputs**:
	* i) Indri index ii) relation database (by default it uses `../snomed_relations/snomed_rel.db`)
* **Output**:
	* Lemin `lgf` file written to the Indri index directory.

## `query_preprocessor`
* Calls MetaMap and produces a concept-based query
* **Input**:
	* Indri stype query file
* **Output**:
	* Indri style query file, but as concepts.
	* Run `query_preprocessor.py -h` for help

## `retrieval`

Run retrieval using the GIN. First it reads the graph created by the indexer; then reads Indri query file; runs retrieval for each query.

* **Inputs**:
	* A Indri stype parameter file (see, for example: "params/cons-medical-topics101-185-active.params")
	* All the options are specified in this parameters, including weighting schemes, file locations and more.
* **Outputs**:
	* TREC style results file written to same directory as the query input file (typically "params" directory)
	* A batch script called "run_network_retrieval_lvl0-2.sh" can be called that will:
		1. Use default param file and qrels
		2. Run retrieval for depth 0,1,2
		3. Evaluate the results of each depth run
		4. Display all results to stdout

## `snomed_relations`

Takes SNOMED as input and creates an sqlite database:

```
CREATE TABLE crel (	
	cui1 varchar(20),
	cui2 varchar(20),
	reltype varchar(20),
	relcharacteristic int,
	PRIMARY KEY(cui1, cui2)
);
```

Each entry represents a SNOMED relationship (and therefore, possible edge in GIN).

The indexer reads this database to lookup SNOMED relationships.

For more info, see the separate README.txt in the directory.

