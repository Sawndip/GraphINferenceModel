for m in `seq 0.3 0.05 1.0`; do
	sed -i "" "s/<minDiffusionFactor>.*<\/minDiffusionFactor>/<minDiffusionFactor>$m<\/minDiffusionFactor>/g" params/cons-medical-topics101-185-active.params
	

	for x in `seq 0.3 0.1 1.0`; do

		echo min=$m mix=$x

		sed -i "" "s/<diffusionMix>.*<\/diffusionMix>/<diffusionMix>$x<\/diffusionMix>/g" params/cons-medical-topics101-185-active.params
		./network_retrieval params/cons-medical-topics101-185-active.params
		./te.sh > evals/minDF_and_dfMix/tuned-minDF_and_dfMix_${m}_${x}-cons-medical-topics101-185-active.eval	
	done

done
