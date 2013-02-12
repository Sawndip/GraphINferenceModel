#!/bin/bash
#
# Runs network_retrieval with different combinations of relation_type weights. Basically, for each relationship type runs retrieval
# with weight 0.1 to 1.0. Results are written to ir.experiments. 
#

# THE NEW WAY - A FULL SWEEP
# for i in `seq 0.2 0.2 1.0`; do
# 	type="116680003"
# 	sed -e  "s/\(^${type}.*\)	[0-9].*$/\1	$i/g" ../snomed_relations/reltypes.txt.template > ../snomed_relations/reltypes.txt
# 	for j in `seq 0.2 0.2 1.0`; do
# 		type="123005000"
# 		sed -e  "s/\(^${type}.*\)	[0-9].*$/\1	$j/g" ../snomed_relations/reltypes.txt.template > ../snomed_relations/reltypes.txt
# 		for k in `seq 0.2 0.2 1.0`; do
# 			type="363698007"
# 			sed -e  "s/\(^${type}.*\)	[0-9].*$/\1	$k/g" ../snomed_relations/reltypes.txt.template > ../snomed_relations/reltypes.txt
			for l in `seq 0.2 0.2 1.0`; do
				type="116676008"
				sed -e  "s/\(^${type}.*\)	[0-9].*$/\1	$l/g" ../snomed_relations/reltypes.txt > ../snomed_relations/reltypes.txt.$l
			done

			./network_retrieval params/cons-medical-topics101-185-active.params 1 0.2 &
			./network_retrieval params/cons-medical-topics101-185-active.params 2 0.2 &

			# ./network_retrieval params/cons-medical-topics101-185-active.params 1 0.4 &
			# ./network_retrieval params/cons-medical-topics101-185-active.params 2 0.4 &
			
			# ./network_retrieval params/cons-medical-topics101-185-active.params 1 0.6 &
			# ./network_retrieval params/cons-medical-topics101-185-active.params 2 0.6 &
			
			# ./network_retrieval params/cons-medical-topics101-185-active.params 1 0.8 &
			# ./network_retrieval params/cons-medical-topics101-185-active.params 2 0.8 &
			
			./network_retrieval params/cons-medical-topics101-185-active.params 1 1.0 &
			./network_retrieval params/cons-medical-topics101-185-active.params 2 1.0

			trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.1_0.2 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.1.${i}_${j}_${k}_0.2 &
			trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.2_0.2 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.2.${i}_${j}_${k}_0.2 &

			# trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.1_0.4 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.1.${i}_${j}_${k}_0.4 &
			# trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.2_0.4 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.2.${i}_${j}_${k}_0.4 &

			# trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.1_0.6 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.1.${i}_${j}_${k}_0.6 &
			# trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.2_0.6 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.2.${i}_${j}_${k}_0.6 &

			# trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.1_0.8 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.1.${i}_${j}_${k}_0.8 &
			# trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.2_0.8 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.2.${i}_${j}_${k}_0.8 &

			trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.1_1 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.1.${i}_${j}_${k}_1.0 &
			trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results.2_1 > /export/data/ir.experiments/tune_graph_reltypes/cons-medical-topics101-185-active.results.2.${i}_${j}_${k}_1.0

# 		done
# 	done
# done

# THE OLD WAY - INDEPEDNENET WEIGHT
# for type in $(cut -f1 ../snomed_relations/reltypes.txt.template); do
# 	#type="246075003"
# 	grep -q $type ../snomed_relations/reltypes_counts.txt
# 	RETVAL=$?
# 	if [ $RETVAL -eq 0 ]; then
# 		echo $type
# 		for i in `seq 0.1 0.1 1.0`; do
# 			sed -e  "s/\(^${type}.*\)	[0-9].*$/\1	$i/g" ../snomed_relations/reltypes.txt.template > ../snomed_relations/reltypes.txt
# 			echo $type $i
# 			./network_retrieval params/cons-medical-topics101-185-active.params 1
# 			trec_eval -q params/medtrack-all.qrel params/cons-medical-topics101-185-active.results > ~/ir.experiments/graph-model/tune_graph_reltypes/cons-medical-topics101-185.${type}_$i.eval
# 		done
# 	fi
# done