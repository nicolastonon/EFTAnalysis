echo "Usage : ./PrintDiffNuisances.sh [path_fitDiag_file.root]"

python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -g nuisances.root $1
python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/printCorrelations.py -i $1:fit_s -p mysyst -m 20

#Plot pre/postfit pulls
expect -c 'spawn -noecho root -l -b nuisances.root
           send ".X Draw_Canvas.C\r"
           send ".q\r"
           interact'

eog nuisances.png &
eog post_fit_errs.png &
