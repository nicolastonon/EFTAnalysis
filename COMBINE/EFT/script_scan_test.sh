#-- To scan multiple operators, do e.g.: '-P ctz -P ctw'
#--

#--------------------------

#datacard_name_prefix=COMBINED_Datacard_TemplateFit_noShape
datacard_name_prefix=COMBINED_Datacard_TemplateFit_noShape_noStat
#datacard_name_prefix=COMBINED_Datacard_TemplateFit
operator="ctz"
datacard_name_suffix="_2017"
#datacard_name_suffix="_Run2"

#--------------------------

rm scan_EFT.log

#-- Common single workspace
#cmd="text2workspace.py ../datacards/${datacard_name_prefix}${datacard_name_suffix}.txt -P EFTModel:eftmodel --PO fits=./Parameterization_EFT.npy -o EFTWorkspace.root"
#echo "$cmd"
#eval "$cmd"

#--------------------------

#for STEP in -3 -2 -1.5 -1 -0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.5 2 3
#for STEP in -4 -3 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 3 4
for STEP in -2 -1 -1.5 0 0.5 1 2
#for STEP in 3
do

#== For each point, create EFT workspace
    cmd="text2workspace.py ../datacards/${datacard_name_prefix}_${operator}_${STEP}${datacard_name_suffix}.txt -P EFTModel:eftmodel --PO fits=./Parameterization_EFT.npy -o EFTWorkspace_${operator}_${STEP}.root"
    echo "$cmd"
    eval "$cmd"

#=== For each point, compute likelihood ratio L_EFT/L_0:
    #-- More options: --robustFit 1 -H AsymptoticLimits --cminPoiOnlyFit --keepFailures >> scan_EFT.log
    cmd="combine ./EFTWorkspace_${operator}_${STEP}.root -M MultiDimFit -P ${operator} --algo fixed --fixedPointPOIs ${operator}=${STEP} -t -1 --saveNLL --saveFitResult --setParameters ${operator}=0 -n _${operator}_${STEP} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
    echo "$cmd"
    eval "$cmd"

done
#--------------------------
