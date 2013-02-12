for a in `seq 1 1 10`; do
	sed -i "" "s/<maxNodesVisited>.*<\/maxNodesVisited>/<maxNodesVisited>$a<\/maxNodesVisited>/g" params/cons-medical-topics101-185-active.params
	./network_retrieval params/cons-medical-topics101-185-active.params
	./te.sh > evals/maxNodesVisited/tuned-maxNodesVisited_$a-cons-medical-topics101-185-active.eval	
done