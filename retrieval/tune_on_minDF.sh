

for a in `seq 0.5 0.05 1.0`; do
	sed -i "" "s/<minDiffusionFactor>.*<\/minDiffusionFactor>/<minDiffusionFactor>$a<\/minDiffusionFactor>/g" params/cons-medical-topics101-185-active.params
	./network_retrieval params/cons-medical-topics101-185-active.params
	./te.sh > evals/minDF/tuned-minDF_$a-cons-medical-topics101-185-active.eval	
done