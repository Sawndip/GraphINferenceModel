for t in $(cut -f1 ../../snomed_relations/reltypes_counts.txt); do
	grep ^bpref ~/ir.experiments/tune_graph_reltypes/*${t}_* | grep all | sort -k3n | tail -1
	#echo $t $bpref
	#echo '--'
done