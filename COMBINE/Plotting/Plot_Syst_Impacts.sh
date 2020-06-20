#Produce plot showing impact of each syst on the signal strength
#First arg must be name of datacard without extension (e.g. "COMBINED_datacard")

# echo "--- Usage : ./Create_Syst_Impact_Plot.sh [datacard_name_without_extension]"
# echo "[-- NB : don't include MC stat. error (not converging) ]"

'''
See: https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools/scripts/plotImpacts.py
Available options:
parser.add_argument('--input', '-i', help='input json file')
parser.add_argument('--output', '-o', help='name of the output file to create')
parser.add_argument('--translate', '-t', help='JSON file for remapping of parameter names')
parser.add_argument('--units', default=None, help='Add units to the best-fit parameter value')
parser.add_argument('--per-page', type=int, default=30, help='Number of parameters to show per page')
parser.add_argument('--max-pages', type=int, default=None, help='Maximum number of pages to write')
parser.add_argument('--height', type=int, default=600, help='Canvas height, in pixels')
parser.add_argument('--left-margin', type=float, default=0.4, help='Left margin, expressed as a fraction')
parser.add_argument('--label-size', type=float, default=0.021, help='Parameter name label size')
parser.add_argument('--cms-label', default='Internal', help='Label next to the CMS logo')
parser.add_argument('--transparent', action='store_true', help='Draw areas as hatched lines instead of solid')
parser.add_argument('--checkboxes', action='store_true', help='Draw an extra panel with filled checkboxes')
parser.add_argument('--blind', action='store_true', help='Do not print best fit signal strength')
parser.add_argument('--color-groups', default=None, help='Comma separated list of GROUP=COLOR') #Ex.:  --color-groups xxx=2 (red)
parser.add_argument('--pullDef',  default=None, help="Choose the definition of the pull, see HiggsAnalysis/CombinedLimit/python/calculate_pulls.py for options")
parser.add_argument('--POI', default=None, help='Specify a POI to draw')

- Add colors: see https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools/scripts/plotImpacts.py#L95-L119
'''


#Check datacard argument
if [[ $1 == *".txt" || $1 == *".root" ]] ; then
    echo "Wrong datacard name argument ! You must remove the file extension ! Abort..."
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
    plotImpacts.py -i impacts.json -o impacts --per-page 20 --translate rename.json --cms-label Internal
elif [ "$2" == "exp" ]; then
    text2workspace.py -m 125 $1.txt
    combineTool.py -t -1 --expectSignal 1 -M Impacts -d $1.root -m 125 --doInitialFit --robustFit 1 $verbosity
    combineTool.py -t -1 --expectSignal 1 -M Impacts -d $1.root -m 125 --robustFit 1 --doFits --parallel 4 $verbosity
    combineTool.py -t -1 --expectSignal 1 -M Impacts -d $1.root -m 125 -o impacts.json $verbosity
    plotImpacts.py -i impacts.json -o impacts --per-page 20 --translate rename.json --cms-label Internal
else
    echo "Usage : Plot_Syst_Impacts.sh [datacard name w/o extension] [exp/obs]"
    echo "\$1 => Name of the datacard, without file extension"
    echo "\$2 => 'exp' or 'obs', depending whether you want to use data or not"
fi

#Plotting
plotImpacts.py -i impacts.json -o impacts -t rename.json
