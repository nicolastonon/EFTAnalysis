for j in {1..100}
do
	NAME=1ktoys_v2_cpt_part${j}.POINTS.0.49
	NAME2=1ktoys_v2_cpt_part${j}_POINTS_0_49

	dir=v2_plots_$NAME
	mkdir $dir

	for i in {1..10}
	#for i in {1..20}
	do
		filename_tmp=$(ls higgsCombine.$NAME.MultiDimFit.mH12*root)
		python EFTPlotter.py -P cpt -n $NAME -f $filename_tmp -itoy $i
		mv scan1D_${NAME2}_toy${i}.png $dir
	done
done
