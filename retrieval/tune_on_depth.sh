for i in `seq 0 10`; do
	echo $i
	#./network_retrieval params/cons-medical-topics101-185-active.params $i
	#trec_eval -q params/medtrack-all.qrel /Users/bevan/phd/papers/network_based_retrieval_as_inference/retrieval/params/cons-medical-topics101-185-active.results > ~/phd/papers/network_based_retrieval_as_inference/retrieval/evals/depth_sweep/cons-medtrack_all-extended-$i.eval
done


for i in `seq 101 185`; do
	grep bpref ~/phd/papers/network_based_retrieval_as_inference/retrieval/evals/depth_sweep/*.eval | grep "\b$i\b" | sort -k3n | tail -1 | sed -E "s/cons-medtrack_all-extended-([0-9]+).eval:bpref/\1/g" >> ~/phd/papers/network_based_retrieval_as_inference/retrieval/evals/depth_sweep/bpref.txt
	grep "\bP_10\b" ~/phd/papers/network_based_retrieval_as_inference/retrieval/evals/depth_sweep/*.eval | grep "\b$i\b" | sort -k3n | tail -1 | sed -E "s/cons-medtrack_all-extended-([0-9]+).eval:P_10/\1/g" >> ~/phd/papers/network_based_retrieval_as_inference/retrieval/evals/depth_sweep/p10.txt
done