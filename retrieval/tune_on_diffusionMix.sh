for a in `seq 0.0 0.1 1.0`; do
	sed -i "" "s/<diffusionMix>.*<\/diffusionMix>/<diffusionMix>$a<\/diffusionMix>/g" params/cons-medical-topics101-185-active.params
	./run_network_retrieval_lvl0-2.sh
	cp params/cons-medical-topics101-185-active.eval evals/diffusion_mix/tuned-diffusionMix_$a-cons-medical-topics101-185-active.eval	
done
