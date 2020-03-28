#Check datacard argument
if [[ $1 == *".txt" || $1 == *".root" ]] ; then
    echo "Wrong datacard argument ! You must remove the file extension ! Abort..."
    exit
fi

text2workspace.py $1.txt -m 125
combine -M MultiDimFit --algo grid --points 50 --rMin -1 --rMax 4 $1.root -m 125 -n nominal
plot1DScan.py higgsCombinenominal.MultiDimFit.mH125.root
