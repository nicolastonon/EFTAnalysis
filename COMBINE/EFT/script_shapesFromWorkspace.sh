#-- NB: need to read proper EFTparameterization file ! Not actually reading the coefficients, *BUT* must make sure the relevant WCs are treated as POIs when creating datacards... !

#Check first argument
# //--------------------------------------------
# if [[ $1 == "" && $1 == *"SM"* ]] ; then
#     echo "Missing argument ! Usage: ./script_shapesFromWorkspace.sh multidimfit.root [filter_region_string]"
#     exit
# fi

fitresult=""
filter=""
if [ "$#" -eq 2 ]; then
    fitresult=$1
    filter=$2
elif [ "$#" -eq 1 ]; then
    if [[ $1 == *"multidimfit"* ]]; then
        fitresult=$1
    else
        filter=$1
    fi
else
    echo "Illegal number of parameters"
    exit 2
fi

# //--------------------------------------------


#-- Options
# //--------------------------------------------
nsamples=100 #Can reduce to speed up ?

#-- Must freeze WCs to 0 to ignore their uncertainties #Will crash if trying to freeze missing variable (then should double-check)
#-- Must dump/read relevant EFTparam file, see comment above (hardcoded filepaths)
freeze="" #Default, do not freeze WCs (not needed e.g. for NN_SM prefit, etc.)
postfit=""
card_prefix="datacard_"
fit_type="prefit"

if [[ ${fitresult} == *"multidimfit"* ]]; then #If a fit result is provided, will get prefit+postfit plots (otherwise, prefit only)
    postfit="--postfit --sampling --samples ${nsamples} -f ${fitresult}:fit_mdf"
    fit_type="postfit"
fi

if [[ ${filter} != *"CR"* ]]; then #No EFT in CRs
    if [[ ${fitresult} == *"ctz"* ]]; then
        freeze="--freeze ctz=0"
        python DumpEFTParameterization.py ../templates/Templates_NN_ctz_EFT2_Run2.root
    elif [[ ${fitresult} == *"ctw"* ]]; then
       freeze="--freeze ctw=0"
        python DumpEFTParameterization.py ../templates/Templates_NN_ctw_EFT2_Run2.root
    elif [[ ${fitresult} == *"cpq3"* ]]; then
        freeze="--freeze cpq3=0"
        python DumpEFTParameterization.py ../templates/Templates_NN_cpq3_EFT2_Run2.root
    elif [[ ${fitresult} == *"cpqm"* ]]; then
        freeze="--freeze cpqm=0"
        python DumpEFTParameterization.py ../templates/Templates_NN_SM_EFT2_Run2.root
    elif [[ ${fitresult} == *"cpt"* ]]; then
        freeze="--freeze cpt=0"
        python DumpEFTParameterization.py ../templates/Templates_NN_SM_EFT2_Run2.root
    elif [[ ${fitresult} == *"5D"* ]]; then #Freeze all 5 WCs
        freeze="--freeze ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0"
        python DumpEFTParameterization.py ../templates/Templates_NN_5D_EFT2_Run2.root
    fi
fi

# regionStr=""
# if [[ $2 != "" ]] ; then
#     regionStr=$2
# fi

# if [[ ${postfit} == "" && ${filter} == *"SM"* ]]; then #Use different 'card prefix' for NN_SM prefit shapes, because not splitting into individual bins (full templates <-> SM only)
if [[ ${postfit} == "" ]]; then #Use different 'card prefix' for NN_SM prefit shapes, because not splitting into individual bins (full templates <-> SM only)
    card_prefix="COMBINED_Datacard_TemplateFit_"
fi

echo "[fitresult=${fitresult}]"
echo "[filter=${filter}]"
echo "[postfit=${postfit}]"
echo "[card_prefix=${card_prefix}]"

# //--------------------------------------------

#-- Main script
# //--------------------------------------------
mkdir dir_shapes_tmp

cp ../datacards/datacards_TemplateFit/${card_prefix}*${filter}*txt ../datacards/
echo "cp ../datacards/datacards_TemplateFit/${card_prefix}*${filter}*txt ../datacards/"

for CARD in $(ls ../datacards/${card_prefix}*${filter}*.txt)
do
    echo "INPUT CARD --> $CARD"
    BASENAME="$(basename $CARD .txt)"
    echo "OUTPUT FILE --> dir_shapes_tmp/shapes_${fit_type}_${BASENAME}.root"
    echo "... Executing: [ python EFTFitter.py -d ../datacards/$CARD --onlyworkspace -o ws_tmp.root ] ..."
    python EFTFitter.py -d ../datacards/$CARD --onlyworkspace -o ws_tmp.root

    echo "... Executing: [ PostFitShapesFromWorkspace --print --total-shapes -w ./ws_tmp.root -o dir_shapes_tmp/shapes_${fit_type}_${BASENAME}.root ${postfit} ${freeze} ] ..."
    PostFitShapesFromWorkspace --print --total-shapes -w ./ws_tmp.root -o dir_shapes_tmp/shapes_${fit_type}_${BASENAME}.root ${postfit} ${freeze}
done
# //--------------------------------------------





#-- Old
# //--------------------------------------------
#-- SRtZq
# python EFTFitter.py -d ../datacards/COMBINED_Datacard_TemplateFit_SRtZq_Run2.txt --onlyworkspace -o ws_srtzq.root
# PostFitShapesFromWorkspace --print --total-shapes --postfit --sampling --samples ${nsamples} -f $1:fit_mdf -w ./ws_srtzq.root -o shapes_postfit_SRtZq.root ${freeze}

# cp ../datacards/datacards_TemplateFit/201*/datacard_*txt ../datacards/
# sed -i 's/..\/..\/..\/templates/..\/templates/g' ../datacards/datacard_*.txt
# //--------------------------------------------
