#!/bin/bash

echo "SHELL: $SHELL"
echo "SCRAM_ARCH: $SCRAM_ARCH"
echo "PATH: $PATH"
echo "PYTHONPATH: $PYTHONPATH"
#echo "source ~/.profile"
#source ~/.profile
#echo "source /cvmfs/cms.cern.ch/cmsset_default.sh"; source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "cd $1"; cd $1

echo ''
echo 'Run job:'
echo ${@:2}
echo ''
${@:2}
echo ''
echo ''

