#Produce plot showing impact of each syst on the signal strength
#First arg must be name of datacard without extension (e.g. "COMBINED_datacard")

# echo "--- Usage : ./Create_Syst_Impact_Plot.sh [datacard_name_without_extension]"
# echo "[-- NB : don't include MC stat. error (not converging) ]"

#Check datacard argument
if [[ $1 == *".txt" || $1 == *".root" ]] ; then
    echo "Wrong datacard argument ! You must remove the file extension ! Abort..."
    exit
fi

rm $1.root

bkg_only=""
if [ "$3" == "1" ]; then
    bkg_only="--setParameters r=0 --freezeParameters r "
fi

verbosity=""
if [ "$4" == "verbose" ]; then
    verbosity="--verbose 9"
fi

if [ "$2" == "obs" ]; then
    text2workspace.py -m 125 $1.txt
    combineTool.py -M Impacts -d $1.root -m 125 --doInitialFit --robustFit 1 --rMin -1 $bkg_only $verbosity
    combineTool.py -M Impacts -d $1.root -m 125 --robustFit 1 --doFits --parallel 4 --rMin -1 $bkg_only $verbosity
    combineTool.py -M Impacts -d $1.root -m 125 -o impacts.json --rMin -1 $bkg_only $verbosity
    plotImpacts.py --transparent -i impacts.json -o impacts --per-page 20 --translate rename.json --cms-label Internal
elif [ "$2" == "exp" ]; then
    text2workspace.py -m 125 $1.txt
    combineTool.py -t -1 --expectSignal 1 -M Impacts -d $1.root -m 125 --doInitialFit --robustFit 1 $verbosity
    combineTool.py -t -1 --expectSignal 1 -M Impacts -d $1.root -m 125 --robustFit 1 --doFits --parallel 4 $verbosity
    combineTool.py -t -1 --expectSignal 1 -M Impacts -d $1.root -m 125 -o impacts.json $verbosity
    plotImpacts.py --transparent -i impacts.json -o impacts --per-page 20 --translate rename.json --cms-label Internal
else
    echo "Usage : Plot_Syst_Impacts.sh [datacard name w/o extension] [exp/obs]"
    echo "\$1 => Name of the datacard, without file extension"
    echo "\$2 => 'exp' or 'obs', depending if you want to use the data or not"
fi

#Plotting
plotImpacts.py -i impacts.json -o impacts -t rename.json
