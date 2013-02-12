#!/bin/bash

params=$1

for k1 in `seq 0.0 0.1 10.0`; do # 0 to 2.0 with 0.1 increments
    for b in `seq 0 0.05 1.0`; do # 0 to 1 with 0.05 increments
    	echo "k1: " $k1 "b: " $b
		sed -E "s/(<indri_tfidf_k1>).*(<\/indri_tfidf_k1>)/\1"$k1"\2/g" $params > ~/ir.experiments/graph-model/sweep-indri_tfidf_params/$params-${k1}_$b.params
		sed -E -i "" "s/(<indri_tfidf_b>).*(<\/indri_tfidf_b>)/\1"$b"\2/g" ~/ir.experiments/graph-model/sweep-indri_tfidf_params/$params-${k1}_$b.params
		./network_retrieval ~/ir.experiments/graph-model/sweep-indri_tfidf_params/$params-${k1}_$b.params
		trec_eval -q params/medtrack-all.qrel params/cons-medtrack_all-extended.results > ~/ir.experiments/graph-model/sweep-indri_tfidf_params/eval/${k1}_$b.eval
	done
done