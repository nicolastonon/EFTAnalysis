#Taken from Potato code

#!/bin/bash

############ DESCRIPTION ###############
# Script to submit a condor job
# Usage: cs [options] <command>
# options could be (argument should be provided after option without space, see example below):
#   -t for max. time (in seconds)
#   -m for max. memory (in KB)
#   -p (takes no argument) to prepare submission files but do not submit the job, just print the command (pretend, or dry run)
#   -n for job name
#   -o for output file names
#   -d for output dir
#   -g for group ('MyProject')
#   -a for array of jobs in the format -amin:max (example: -a1:3 submits 3 jobs with arguments 1,2,3)
#   key=value for any other option to pass to condor_q -append
# Example to submit job with name 'xfit-MyJob' for maximum 1 hour executing command "./xfitter":
# cs -nxfit-MyJob -t3600 ./xfitter
# By default max. time and memory are 12h and 2GB
# Files used for submission as well as job output files are stored in directory 'condor' wich is automatically created in current directory
############ END OF DESCRIPTION ###############

# split arguments into two parts (2nd part starts with first argument which does not begin with '-' and does not contain '='):
# (1) options for submission script
# (2) command to run on the cluster
ArgsSub=''
ArgsExe=''
for a do
  #echo $a
  if [[ ! -z "${ArgsExe// }" ]]; then
    ArgsExe=${ArgsExe}' '$a
  elif [[ ${a:0:1} != '-' ]] && [[ *"$a"* != *"="* ]]; then
    ArgsExe=${ArgsExe}' '$a
  else
    ArgsSub=${ArgsSub}' '$a
  fi
done
#echo "ArgsSub: $ArgsSub"
#echo "ArgsExe: $ArgsExe"

# check command
# check if the command to run on the cluster was provided
if [[ -z "${ArgsExe// }" ]]; then
  echo "Error: command to run is not provided"
  echo "Usage: cs [options] <command>"
  exit 1
fi
# check whether command to run is valid
ArgsExe1=`echo $ArgsExe | awk '{print $1}'`
which ${ArgsExe1} >& /dev/null
if [ $? -ne 0 ]; then
  echo "Error: command to run is not valid, see below:"
  which ${ArgsExe1}
  exit 1
fi

# check arguments
FlagPretend=0
#JobName=`echo $ArgsExe | sed -e 's/ /-/g' | sed -e 's/\.\///g'`
#JobName=`echo $ArgsExe | sed -e 's/\.\///g' | awk '{print $1}'`
JobName=`basename "$ArgsExe" | awk '{print $1}'`
JobMem='1999'
JobTime='10799'
ExtraOptions=()
FlagOuputNameFileExplicit=0
JobArray=0
JobArrayMin=JobArrayMin=''

# directory for intermediate scripts and output
NameCondorDir='condor'

for arg in ${ArgsSub}; do
  # check for pretend
  if [ ${arg} == "--pretend" ] || [ ${arg} == "-p" ]; then
    FlagPretend=1
  # check for job name
  elif [[ "${arg}"* == "-n"* ]]; then
    JobName=${arg:2}
  # check for requested memory
  elif [[ "${arg}"* == "-m"* ]]; then
    JobMem=${arg:2}
  # check for requested run time
  elif [[ "${arg}"* == "-t"* ]]; then
    JobTime=${arg:2}
  # check for output file pattern
  elif [[ "${arg}"* == "-o"* ]]; then
    NameOutPattern=${arg:2}
    NameFileOutLog=${NameCondorDir}'/'${NameOutPattern}'.log'
    NameFileOutTxt=${NameCondorDir}'/'${NameOutPattern}'.txt'
    NameFileOutErr=${NameCondorDir}'/'${NameOutPattern}'.err'
    FlagOuputNameFileExplicit=1
  elif [[ "${arg}"* == "-d"* ]]; then
	echo 'NameOutDir='${arg:2}
    NameOutDir=${arg:2}
    NameFileOutLog=${NameCondorDir}'/'${NameOutDir}'/'${JobName}'.log'
    NameFileOutTxt=${NameCondorDir}'/'${NameOutDir}'/'${JobName}'.txt'
    NameFileOutErr=${NameCondorDir}'/'${NameOutDir}'/'${JobName}'.err'
    FlagOuputNameFileExplicit=1
  elif [[ *"${arg}"* == *"="* ]]; then
    ExtraOptions+=('-append '${arg})
  # check if this is a job array (with an array of arguments)
  elif [[ "${arg}"* == "-a"* ]]; then
    JobArray=1
    JobArrayMin=`echo ${arg:2} | awk -F: '{print $1}'`
    JobArrayMax=`echo ${arg:2} | awk -F: '{print $2}'`
  else
    echo "Warning: ignoring unknown argument ${arg}"
  fi
done

# create directory for intermediate scripts and output
mkdir -p ${NameCondorDir}

if [ $FlagOuputNameFileExplicit -eq 0 ]; then
  NameFileOutLog=${NameCondorDir}'/'${JobName}'.log'
  NameFileOutTxt=${NameCondorDir}'/'${JobName}'.txt'
  NameFileOutErr=${NameCondorDir}'/'${JobName}'.err'
fi

# for array of jobs, modified command to run and output file names by extending them with the job number
if [ ${JobArray} -eq 1 ]; then
  ArgsExe="${ArgsExe} \$1"
  NameFileOutLog=${NameFileOutLog::${#NameFileOutLog}-4}'-$(Item).log'
  NameFileOutTxt=${NameFileOutTxt::${#NameFileOutTxt}-4}'-$(Item).txt'
  NameFileOutErr=${NameFileOutErr::${#NameFileOutErr}-4}'-$(Item).err'
fi

# trick for condor: store LD_LIBRARY_PATH
export LD_LIBRARY_PATH_STORED=$LD_LIBRARY_PATH

# create script to run on  cluster side
#NameScriptFile=${NameCondorDir}'/condor-run-'${JobName}'.sh'
NameScriptFile=${NameCondorDir}'/'${NameOutDir}'/condor-run-'${JobName}'.sh'
rm -f ${NameScriptFile} #Remove any pre-existing file with same name #-f avoids warning

echo '#!/bin/bash' >> ${NameScriptFile}
echo '' >> ${NameScriptFile}
echo 'if [ ! -z $LD_LIBRARY_PATH_STORED ]; then' >> ${NameScriptFile}
echo '  echo "Using LD_LIBRARY_PATH_STORED"' >> ${NameScriptFile}
echo '  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_STORED:$LD_LIBRARY_PATH' >> ${NameScriptFile}
echo 'fi' >> ${NameScriptFile}
echo '' >> ${NameScriptFile}

#ADDED
#echo 'export SCRAM_ARCH='${SCRAM_ARCH} >> ${NameScriptFile}
#echo 'export CMSSW_VERSION='${CMSSW_VERSION} >> ${NameScriptFile}
echo 'export LD_LIBRARY_PATH='${LD_LIBRARY_PATH}':/nfs/dust/cms/user/ntonon/TopEFT/'${CMSSW_VERSION}'/lib/'${SCRAM_ARCH} >> ${NameScriptFile}
echo 'export PYTHONPATH='${PYTHONPATH}':/nfs/dust/cms/user/ntonon/TopEFT/CMSSW_11_1_2/src/EFTAnalysis/NeuralNetworks' >> ${NameScriptFile}

echo 'source ~/.profile' >> ${NameScriptFile} #Must source conda, packages, ... on HTCondor server
echo 'echo "command to run: '${ArgsExe}'"' >> ${NameScriptFile}
echo '' >> ${NameScriptFile}
echo '' >> ${NameScriptFile}
echo ${ArgsExe} >> ${NameScriptFile}
chmod +x ${NameScriptFile}

# create configureation submit file to submit job to the cluster
#NameSubmitFile='condor/condor-submit-'${JobName}
NameSubmitFile=${NameCondorDir}'/'${NameOutDir}'/condor-submit-'${JobName}
rm -f ${NameSubmitFile} #Remove any pre-existing file with same name #-f avoids warning

echo 'Executable = '`pwd`'/'${NameScriptFile} >> ${NameSubmitFile}
echo 'log = '${NameFileOutLog} >> ${NameSubmitFile}
echo 'output = '${NameFileOutTxt} >> ${NameSubmitFile}
echo 'error = '${NameFileOutErr} >> ${NameSubmitFile}
echo 'Universe = vanilla' >> ${NameSubmitFile}
echo 'notification = Error' >> ${NameSubmitFile}
echo 'Initialdir = '`pwd` >> ${NameSubmitFile}
echo '+RequestRuntime = '${JobTime} >> ${NameSubmitFile}
echo 'RequestMemory = '${JobMem} >> ${NameSubmitFile}
#echo 'getenv=True' >> ${NameSubmitFile} #CONFLICTING WITH CONDA ENV !

#echo 'Requirements = ( OpSysAndVer == "CentOS7" || OpSysAndVer == "SL6")' >> ${NameSubmitFile}
#if [[ ${SCRAM_ARCH:3:1} == 7 ]]; then
#  echo 'Requirements = (OpSysAndVer == "CentOS7")' >> ${NameSubmitFile}
#else
#  echo 'Requirements = (OpSysAndVer == "SL6")' >> ${NameSubmitFile}
#fi
echo 'Requirements = (OpSysAndVer == "CentOS7")' >> ${NameSubmitFile}

#if [ ${JobArray} -eq 0 ]; then
#  echo 'queue 1' >> ${NameSubmitFile}
#else
#  echo 'arguments = $(Item)' >> ${NameSubmitFile}
#  echo "queue from seq ${JobArrayMin} ${JobArrayMax} |" >> ${NameSubmitFile}
#fi
echo 'queue 1' >> ${NameSubmitFile} #Done describing the job, can submit it #Number defines $(Item)?

# submit job (or print corresponding command)
com="condor_submit -batch-name ${JobName} ${ExtraOptions} ${NameSubmitFile}"
if [ ${FlagPretend} -eq 1 ]; then
  com='echo '${com}
fi
${com}

exit $?
