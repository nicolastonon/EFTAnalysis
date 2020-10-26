#Nicolas -- This script merges the ntuples from all 3 years ; these merged files will be used when processing the full run 2 dataset
#-- OBSOLETE

#New repo, store merged files
outDir="Run2"
mkdir $outDir

# dirCheck="./2016"
# sampleList2016=$(ls $dirCheck | cut -d"." -f1)
# echo $sampleList2016
# dirCheck="./2017"
# sampleList2017=$(ls $dirCheck | cut -d"." -f1)
# echo $sampleList2017

#Get list of samples from 2018 repository (contains all samples)
dirCheck="./2018"
sampleList2018=$(ls $dirCheck | cut -d"." -f1)
# echo $sampleList2018

#For each sample, merge all 3 years
for sample in $sampleList2018
do
    # echo $sample
    hadd -f $outDir/$sample.root 2016/$sample.root 2017/$sample.root 2018/$sample.root
done
