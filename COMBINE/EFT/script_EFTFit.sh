# //--------------------------------------------

#Check first argument
if [[ $1 == "" ]] ; then
    #echo "Missing argument ! Usage: [./script_EFTFit.sh ../datacards/datacard.txt]"
    echo "Missing argument ! Usage: [./script_EFTFit.sh ctz]"
    exit
fi

#-- CREATE WORKSPACE
#python EFTFitter.py -d ../datacards/COMBINED_Datacard_TemplateFit_Run2.txt --onlyworkspace


# //--------------------------------------------

python EFTFitter.py --exp -m bestfit -d ./EFTWorkspace.root --freeze --exp -P $1

# //--------------------------------------------
