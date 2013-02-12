for s in `seq 0.0 0.1 1`; do
	./run_network_retrieval_lvl0-2.sh params/cons-medtrack_all-extended.params params/medtrack-all.qrel $s
	cp params/cons-medtrack_all-extended.eval ${s}_damp.eval
done