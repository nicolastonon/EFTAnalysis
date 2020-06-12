#Produce plot showing impact of each Systematic
#argument must be name of datacard without extension (e.g. "COMBINED_datacard_TemplateFit_tZq")

echo "--- Usage : ./Breakdown_Stat_Syst.sh + [datacard_name_without_extension]"

#Produce plot showing breakdown of statistical and systematic uncertainties
text2workspace.py -m 125 $1.txt
combine -M MultiDimFit --algo grid --points 50 --rMin -1 --rMax 4 $1.root -m 125 -n nominal
plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root
combine -M MultiDimFit --algo none --rMin -1 --rMax 4 $1.root -m 125 -n bestfit --saveWorkspace
combine -M MultiDimFit --algo grid --points 50 --rMin -1 --rMax 4 -m 125 -n stat higgsCombinebestfit.MultiDimFit.mH125.root --snapshotName MultiDimFit --freezeParameters all
plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root --others 'higgsCombinestat.MultiDimFit.mH125.root:Freeze all:2' --breakdown syst,stat
