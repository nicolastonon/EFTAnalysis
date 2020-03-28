#Check compatibility between channels : outputs a chi2-like value, and a compatibility plot
#argument must be name of datacard without extension (e.g. "COMBINED_datacard_TemplateFit_tZq")

echo "--- Usage : ./Check_Compatibility_Channels.sh + [datacard_name_without_extension]"

combine -M ChannelCompatibilityCheck $1.txt --saveFitResult

expect -c 'spawn -noecho root -l -b higgsCombineTest.ChannelCompatibilityCheck.mH120.root
           send ".X cccPlot.cxx\r"
           send ".q\r"
           interact'

evince ccc.pdf &
