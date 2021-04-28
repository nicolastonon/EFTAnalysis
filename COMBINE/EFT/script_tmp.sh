for i in {1..100}
do
	python EFTFitter.py -d ./EFTWorkspace_SM.root -P cpt -n 1ktoys_v2_cpt_part$i -m scan --freeze -t 20 --batch

	#mv scan1D_${NAME}_toy${i}.png plots_$NAME
done
