params=params/cons-reports-medical-topics101-185-active.params
base=$(echo $params | sed "s/params//g" | cut -c2-)

for a in `seq 100 100 1000`; do
	sed -i "" "s/<LM-Dirichlet.mu>.*<\/LM-Dirichlet.mu>/<LM-Dirichlet.mu>$a<\/LM-Dirichlet.mu>/g" $params
	./run_network_retrieval_lvl0-2.sh $params
	cp params/${base}eval evals/mu/reports/tuned-mu_$a-${base}eval
done