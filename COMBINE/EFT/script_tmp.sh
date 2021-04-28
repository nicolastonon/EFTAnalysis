for i in {1..50}
do
	python EFTFitter.py -d ./EFTWorkspace_SM.root -P cpt -n 1ktoys_v3_cpt_part$i -m scan --exp --freeze -t 10 --batch -points 40

	#mv scan1D_${NAME}_toy${i}.png plots_$NAME
done
