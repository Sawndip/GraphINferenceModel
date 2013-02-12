reltypes.txt:
	This file contains a list of reltypes in the form:
		<reltype_id> [tab] <name/desc> [tab] <weight>
	The file is read by network_retrieval to get a weight for each reltype it 
	encounters during retrieval

reltypes_uniform.txt:
	As above, but all reltype set to 1.0, i.e., reltypes do not affect scoring.
	To utilise this file you should "cp reltypes_uniform.txt reltypes.txt"
	
reltypes.txt.template:
	This file is read by ../retrieval/tune_on_reltype.sh, which sweeps across the 
	reltypes weights. tune_on_reltype.sh will read reltypes.txt.template, set the 
	weight for a given reltype and the write the results to reltypes.txt so that
	it can be read by network_retrieval
	
reltypes_by_number_in_snomed.txt:
	Reltypes ordered (and therefore weighted) by frequency of reltype occurrence in
	MRREL file, i.e. SNOMED/UMLS.
	
reltypes_counts.txt:
	Statistics about how many of each reltype was seen when processing the MedTrack'11
	queries.

snomed_relations.py:
	Python script of r

snomed_rel.db:
	Database of reltypes in the form: cui1, cui2, reltype. This is read by network_index
	to create the concept graph.
	
snomed_relations.py:
	Python script for creating snomed_rel.db.
	
snomed_semantic_types.txt:
	The SNOMED semantic type for each SNOMED concept. As specified in the SNOMED 
	sct1_Concepts_Core_INT_20110131.txt file, e.g. the line:
	"152309009       0       Entire submental vein (body structure)  XS066   T-48144 1"
	would be of type "body structure".
	