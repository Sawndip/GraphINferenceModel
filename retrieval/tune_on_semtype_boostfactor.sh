for i in `seq 1 10`; do
	./network_retrieval params/cons-medical-topics101-135.params 1 1 $i
	trec_eval -q params/medtrack_visits.qrel /Users/bevan/phd/papers/retrieval_inference_using_umls/retrieval/params/cons-medical-topics101-135.results > /Users/bevan/ir.experiments/graph-model/tune_semtype_boost/cons-medical-topics101-135.$i.eval
done