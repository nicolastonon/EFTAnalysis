#!/bin/bash
cd $1
source /cvmfs/cms.cern.ch/cmsset_default.sh
#eval `scramv1 runtime -sh`
echo ''
echo 'Run job:'
echo ${@:2}
echo ''
${@:2}
echo ''
echo ''

