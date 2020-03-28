# http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/#likelihood-fits-and-scans
# saturated, KS, AD
# recommended to use '--toysFreq' for 'saturated' only
#could add option '--fixedSignalStrength=1'

#Check datacard argument
if [[ $1 == *".txt" || $1 == *".root" ]] ; then
    echo "Wrong datacard argument ! You must remove the file extension ! Abort..."
    exit
fi

text2workspace.py $1.txt
combine -M GoodnessOfFit --algo=saturated --fixedSignalStrength=1 -m 125 -d $1.root -t $2 -s 1234 --saveToys --toysFreq
combine -M GoodnessOfFit --algo=saturated --fixedSignalStrength=1 -m 125 -d  $1.root
combineTool.py -M CollectGoodnessOfFit --input higgsCombineTest.GoodnessOfFit.mH125.root higgsCombineTest.GoodnessOfFit.mH125.1234.root -o GoF.json
plotGof.py --statistic saturated --mass 125.0 GoF.json -o GoF
