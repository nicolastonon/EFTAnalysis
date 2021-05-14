#!/bin/sh
ulimit -s unlimited
set -e
cd /nfs/dust/cms/user/ntonon/TopEFT/CMSSW_10_2_24/src
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /nfs/dust/cms/user/ntonon/TopEFT/CMSSW_10_2_24/src/EFTAnalysis/COMBINE/EFT

if [ $1 -eq 0 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_AbsoluteMPFBias --algo impact --redefineSignalPOIs cpq3 -P AbsoluteMPFBias --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 1 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_AbsoluteScale --algo impact --redefineSignalPOIs cpq3 -P AbsoluteScale --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 2 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_AbsoluteStat2016 --algo impact --redefineSignalPOIs cpq3 -P AbsoluteStat2016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 3 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_AbsoluteStat2017 --algo impact --redefineSignalPOIs cpq3 -P AbsoluteStat2017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 4 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_AbsoluteStat2018 --algo impact --redefineSignalPOIs cpq3 -P AbsoluteStat2018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 5 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagCFerr1 --algo impact --redefineSignalPOIs cpq3 -P BtagCFerr1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 6 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagCFerr2 --algo impact --redefineSignalPOIs cpq3 -P BtagCFerr2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 7 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagHF --algo impact --redefineSignalPOIs cpq3 -P BtagHF --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 8 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagHFstats12016 --algo impact --redefineSignalPOIs cpq3 -P BtagHFstats12016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 9 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagHFstats12017 --algo impact --redefineSignalPOIs cpq3 -P BtagHFstats12017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 10 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagHFstats12018 --algo impact --redefineSignalPOIs cpq3 -P BtagHFstats12018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 11 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagHFstats22016 --algo impact --redefineSignalPOIs cpq3 -P BtagHFstats22016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 12 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagHFstats22017 --algo impact --redefineSignalPOIs cpq3 -P BtagHFstats22017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 13 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagHFstats22018 --algo impact --redefineSignalPOIs cpq3 -P BtagHFstats22018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 14 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagLF --algo impact --redefineSignalPOIs cpq3 -P BtagLF --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 15 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagLFstats12016 --algo impact --redefineSignalPOIs cpq3 -P BtagLFstats12016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 16 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagLFstats12017 --algo impact --redefineSignalPOIs cpq3 -P BtagLFstats12017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 17 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagLFstats12018 --algo impact --redefineSignalPOIs cpq3 -P BtagLFstats12018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 18 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagLFstats22016 --algo impact --redefineSignalPOIs cpq3 -P BtagLFstats22016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 19 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagLFstats22017 --algo impact --redefineSignalPOIs cpq3 -P BtagLFstats22017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 20 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_BtagLFstats22018 --algo impact --redefineSignalPOIs cpq3 -P BtagLFstats22018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 21 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_FRe_be --algo impact --redefineSignalPOIs cpq3 -P FRe_be --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 22 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_FRe_norm --algo impact --redefineSignalPOIs cpq3 -P FRe_norm --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 23 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_FRe_pt --algo impact --redefineSignalPOIs cpq3 -P FRe_pt --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 24 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_FRm_be --algo impact --redefineSignalPOIs cpq3 -P FRm_be --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 25 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_FRm_norm --algo impact --redefineSignalPOIs cpq3 -P FRm_norm --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 26 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_FRm_pt --algo impact --redefineSignalPOIs cpq3 -P FRm_pt --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 27 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_FSR --algo impact --redefineSignalPOIs cpq3 -P FSR --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 28 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_FlavorQCD --algo impact --redefineSignalPOIs cpq3 -P FlavorQCD --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 29 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_Fragmentation --algo impact --redefineSignalPOIs cpq3 -P Fragmentation --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 30 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_ISRtWZ --algo impact --redefineSignalPOIs cpq3 -P ISRtWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 31 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_ISRtZq --algo impact --redefineSignalPOIs cpq3 -P ISRtZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 32 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_ISRttZ --algo impact --redefineSignalPOIs cpq3 -P ISRttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 33 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_JER2016 --algo impact --redefineSignalPOIs cpq3 -P JER2016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 34 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_JER2017 --algo impact --redefineSignalPOIs cpq3 -P JER2017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 35 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_JER2018 --algo impact --redefineSignalPOIs cpq3 -P JER2018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 36 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_LepEff_elLoose --algo impact --redefineSignalPOIs cpq3 -P LepEff_elLoose --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 37 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_LepEff_elTight --algo impact --redefineSignalPOIs cpq3 -P LepEff_elTight --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 38 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_LepEff_muLoose --algo impact --redefineSignalPOIs cpq3 -P LepEff_muLoose --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 39 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_LepEff_muTight --algo impact --redefineSignalPOIs cpq3 -P LepEff_muTight --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 40 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_Lumi16 --algo impact --redefineSignalPOIs cpq3 -P Lumi16 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 41 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_Lumi161718 --algo impact --redefineSignalPOIs cpq3 -P Lumi161718 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 42 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_Lumi17 --algo impact --redefineSignalPOIs cpq3 -P Lumi17 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 43 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_Lumi1718 --algo impact --redefineSignalPOIs cpq3 -P Lumi1718 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 44 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_Lumi18 --algo impact --redefineSignalPOIs cpq3 -P Lumi18 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 45 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_MET2016 --algo impact --redefineSignalPOIs cpq3 -P MET2016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 46 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_MET2017 --algo impact --redefineSignalPOIs cpq3 -P MET2017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 47 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_MET2018 --algo impact --redefineSignalPOIs cpq3 -P MET2018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 48 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_MEtWZ --algo impact --redefineSignalPOIs cpq3 -P MEtWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 49 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_MEtZq --algo impact --redefineSignalPOIs cpq3 -P MEtZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 50 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_MEttZ --algo impact --redefineSignalPOIs cpq3 -P MEttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 51 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_NPL_rate --algo impact --redefineSignalPOIs cpq3 -P NPL_rate --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 52 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PDF --algo impact --redefineSignalPOIs cpq3 -P PDF --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 53 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PU --algo impact --redefineSignalPOIs cpq3 -P PU --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 54 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PileUpDataMC --algo impact --redefineSignalPOIs cpq3 -P PileUpDataMC --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 55 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PileUpPtBB --algo impact --redefineSignalPOIs cpq3 -P PileUpPtBB --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 56 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PileUpPtEC1 --algo impact --redefineSignalPOIs cpq3 -P PileUpPtEC1 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 57 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PileUpPtEC2 --algo impact --redefineSignalPOIs cpq3 -P PileUpPtEC2 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 58 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PileUpPtHF --algo impact --redefineSignalPOIs cpq3 -P PileUpPtHF --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 59 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PileUpPtRef --algo impact --redefineSignalPOIs cpq3 -P PileUpPtRef --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 60 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PrivMC_tWZ_rate --algo impact --redefineSignalPOIs cpq3 -P PrivMC_tWZ_rate --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 61 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PrivMC_tZq_rate --algo impact --redefineSignalPOIs cpq3 -P PrivMC_tZq_rate --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 62 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_PrivMC_ttZ_rate --algo impact --redefineSignalPOIs cpq3 -P PrivMC_ttZ_rate --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 63 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeBal --algo impact --redefineSignalPOIs cpq3 -P RelativeBal --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 64 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeFSR --algo impact --redefineSignalPOIs cpq3 -P RelativeFSR --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 65 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeJEREC12016 --algo impact --redefineSignalPOIs cpq3 -P RelativeJEREC12016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 66 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeJEREC12017 --algo impact --redefineSignalPOIs cpq3 -P RelativeJEREC12017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 67 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeJEREC12018 --algo impact --redefineSignalPOIs cpq3 -P RelativeJEREC12018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 68 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeJEREC22016 --algo impact --redefineSignalPOIs cpq3 -P RelativeJEREC22016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 69 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeJEREC22017 --algo impact --redefineSignalPOIs cpq3 -P RelativeJEREC22017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 70 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeJEREC22018 --algo impact --redefineSignalPOIs cpq3 -P RelativeJEREC22018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 71 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeJERHF --algo impact --redefineSignalPOIs cpq3 -P RelativeJERHF --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 72 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativePtBB --algo impact --redefineSignalPOIs cpq3 -P RelativePtBB --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 73 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativePtEC12016 --algo impact --redefineSignalPOIs cpq3 -P RelativePtEC12016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 74 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativePtEC12017 --algo impact --redefineSignalPOIs cpq3 -P RelativePtEC12017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 75 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativePtEC12018 --algo impact --redefineSignalPOIs cpq3 -P RelativePtEC12018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 76 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativePtEC22016 --algo impact --redefineSignalPOIs cpq3 -P RelativePtEC22016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 77 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativePtEC22017 --algo impact --redefineSignalPOIs cpq3 -P RelativePtEC22017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 78 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativePtEC22018 --algo impact --redefineSignalPOIs cpq3 -P RelativePtEC22018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 79 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativePtHF --algo impact --redefineSignalPOIs cpq3 -P RelativePtHF --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 80 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeSample2016 --algo impact --redefineSignalPOIs cpq3 -P RelativeSample2016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 81 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeSample2017 --algo impact --redefineSignalPOIs cpq3 -P RelativeSample2017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 82 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeSample2018 --algo impact --redefineSignalPOIs cpq3 -P RelativeSample2018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 83 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeStatEC2016 --algo impact --redefineSignalPOIs cpq3 -P RelativeStatEC2016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 84 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeStatEC2017 --algo impact --redefineSignalPOIs cpq3 -P RelativeStatEC2017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 85 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeStatEC2018 --algo impact --redefineSignalPOIs cpq3 -P RelativeStatEC2018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 86 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeStatFSR2016 --algo impact --redefineSignalPOIs cpq3 -P RelativeStatFSR2016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 87 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeStatFSR2017 --algo impact --redefineSignalPOIs cpq3 -P RelativeStatFSR2017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 88 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeStatFSR2018 --algo impact --redefineSignalPOIs cpq3 -P RelativeStatFSR2018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 89 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeStatHF2016 --algo impact --redefineSignalPOIs cpq3 -P RelativeStatHF2016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 90 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeStatHF2017 --algo impact --redefineSignalPOIs cpq3 -P RelativeStatHF2017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 91 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_RelativeStatHF2018 --algo impact --redefineSignalPOIs cpq3 -P RelativeStatHF2018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 92 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_SinglePionECAL --algo impact --redefineSignalPOIs cpq3 -P SinglePionECAL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 93 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_SinglePionHCAL --algo impact --redefineSignalPOIs cpq3 -P SinglePionHCAL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 94 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_TimePtEta2016 --algo impact --redefineSignalPOIs cpq3 -P TimePtEta2016 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 95 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_TimePtEta2017 --algo impact --redefineSignalPOIs cpq3 -P TimePtEta2017 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 96 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_TimePtEta2018 --algo impact --redefineSignalPOIs cpq3 -P TimePtEta2018 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 97 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_Trigger16 --algo impact --redefineSignalPOIs cpq3 -P Trigger16 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 98 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_Trigger17 --algo impact --redefineSignalPOIs cpq3 -P Trigger17 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 99 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_Trigger18 --algo impact --redefineSignalPOIs cpq3 -P Trigger18 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 100 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_VVV_rate --algo impact --redefineSignalPOIs cpq3 -P VVV_rate --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 101 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_WZ_HF_extrap --algo impact --redefineSignalPOIs cpq3 -P WZ_HF_extrap --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 102 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_WZ_rate --algo impact --redefineSignalPOIs cpq3 -P WZ_rate --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 103 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_XG_rate --algo impact --redefineSignalPOIs cpq3 -P XG_rate --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 104 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_alphas --algo impact --redefineSignalPOIs cpq3 -P alphas --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 105 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_cpqm --algo impact --redefineSignalPOIs cpq3 -P cpqm --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 106 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_cpt --algo impact --redefineSignalPOIs cpq3 -P cpt --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 107 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_ctw --algo impact --redefineSignalPOIs cpq3 -P ctw --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 108 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_ctz --algo impact --redefineSignalPOIs cpq3 -P ctz --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 109 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_njets_tZq --algo impact --redefineSignalPOIs cpq3 -P njets_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 110 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prefire --algo impact --redefineSignalPOIs cpq3 -P prefire --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 111 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin10_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin10_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 112 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin10_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin10_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 113 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin10_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin10_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 114 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin11_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin11_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 115 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin11_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin11_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 116 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin11_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin11_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 117 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin12_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin12_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 118 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin12_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin12_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 119 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin12_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin12_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 120 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin13_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin13_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 121 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin13_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin13_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 122 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin13_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin13_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 123 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin14_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin14_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 124 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin14_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin14_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 125 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin14_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin14_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 126 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin15_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin15_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 127 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin15_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin15_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 128 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin15_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin15_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 129 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin1_NN_cpq3_SRtZq_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin1_NN_cpq3_SRtZq_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 130 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin1_NN_cpq3_SRtZq_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin1_NN_cpq3_SRtZq_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 131 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin1_NN_cpq3_SRtZq_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin1_NN_cpq3_SRtZq_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 132 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin1_NN_cpq3_SRttZ_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin1_NN_cpq3_SRttZ_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 133 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin1_NN_cpq3_SRttZ_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin1_NN_cpq3_SRttZ_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 134 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin1_NN_cpq3_SRttZ_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin1_NN_cpq3_SRttZ_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 135 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin1_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin1_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 136 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin1_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin1_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 137 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin1_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin1_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 138 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin2_NN_cpq3_SRtZq_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin2_NN_cpq3_SRtZq_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 139 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin2_NN_cpq3_SRtZq_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin2_NN_cpq3_SRtZq_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 140 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin2_NN_cpq3_SRtZq_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin2_NN_cpq3_SRtZq_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 141 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin2_NN_cpq3_SRttZ_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin2_NN_cpq3_SRttZ_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 142 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin2_NN_cpq3_SRttZ_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin2_NN_cpq3_SRttZ_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 143 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin2_NN_cpq3_SRttZ_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin2_NN_cpq3_SRttZ_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 144 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin2_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin2_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 145 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin2_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin2_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 146 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin2_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin2_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 147 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin3_NN_cpq3_SRtZq_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin3_NN_cpq3_SRtZq_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 148 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin3_NN_cpq3_SRtZq_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin3_NN_cpq3_SRtZq_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 149 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin3_NN_cpq3_SRtZq_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin3_NN_cpq3_SRtZq_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 150 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin3_NN_cpq3_SRttZ_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin3_NN_cpq3_SRttZ_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 151 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin3_NN_cpq3_SRttZ_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin3_NN_cpq3_SRttZ_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 152 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin3_NN_cpq3_SRttZ_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin3_NN_cpq3_SRttZ_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 153 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin3_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin3_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 154 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin3_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin3_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 155 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin3_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin3_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 156 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRtZq_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRtZq_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 157 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRtZq_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRtZq_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 158 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRtZq_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRtZq_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 159 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 160 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2017_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2017_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 161 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2017_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2017_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 162 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2017_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2017_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 163 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2017_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2017_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 164 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2017_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2017_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 165 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2017_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2017_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 166 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2017_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2017_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 167 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2017_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2017_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 168 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2018_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2018_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 169 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2018_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2018_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 170 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2018_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2018_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 171 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2018_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2018_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 172 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2018_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2018_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 173 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2018_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2018_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 174 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2018_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2018_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 175 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_NN_cpq3_SRttZ_2018_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_NN_cpq3_SRttZ_2018_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 176 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 177 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 178 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin4_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin4_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 179 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2016_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2016_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 180 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2016_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2016_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 181 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2016_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2016_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 182 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2016_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2016_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 183 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2016_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2016_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 184 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2016_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2016_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 185 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2016_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2016_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 186 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2016_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2016_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 187 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 188 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRtZq_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRtZq_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 189 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2016_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2016_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 190 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2016_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2016_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 191 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2016_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2016_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 192 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2016_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2016_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 193 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2016_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2016_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 194 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2016_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2016_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 195 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2016_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2016_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 196 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2016_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2016_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 197 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2017_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2017_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 198 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2017_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2017_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 199 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2017_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2017_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 200 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2017_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2017_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 201 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2017_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2017_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 202 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2017_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2017_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 203 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2017_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2017_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 204 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2017_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2017_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 205 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_NN_cpq3_SRttZ_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_NN_cpq3_SRttZ_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 206 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 207 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 208 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin5_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin5_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 209 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2016_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2016_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 210 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2016_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2016_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 211 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2016_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2016_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 212 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2016_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2016_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 213 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2016_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2016_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 214 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2016_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2016_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 215 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2016_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2016_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 216 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2016_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2016_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 217 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2017_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2017_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 218 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2017_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2017_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 219 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2017_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2017_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 220 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2017_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2017_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 221 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2017_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2017_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 222 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2017_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2017_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 223 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2017_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2017_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 224 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2017_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2017_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 225 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2018_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2018_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 226 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2018_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2018_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 227 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2018_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2018_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 228 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2018_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2018_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 229 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2018_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2018_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 230 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2018_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2018_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 231 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2018_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2018_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 232 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRtZq_2018_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRtZq_2018_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 233 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRttZ_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRttZ_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 234 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRttZ_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRttZ_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 235 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_NN_cpq3_SRttZ_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_NN_cpq3_SRttZ_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 236 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 237 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 238 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin6_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin6_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 239 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2016_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2016_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 240 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2016_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2016_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 241 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2016_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2016_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 242 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2016_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2016_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 243 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2016_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2016_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 244 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2016_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2016_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 245 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2016_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2016_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 246 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2016_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2016_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 247 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2017_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2017_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 248 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2017_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2017_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 249 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2017_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2017_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 250 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2017_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2017_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 251 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2017_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2017_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 252 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2017_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2017_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 253 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2017_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2017_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 254 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2017_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2017_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 255 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_NN_cpq3_SRtZq_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_NN_cpq3_SRtZq_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 256 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 257 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 258 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin7_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin7_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 259 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2016_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2016_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 260 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2016_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2016_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 261 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2016_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2016_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 262 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2016_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2016_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 263 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2016_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2016_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 264 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2016_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2016_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 265 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2016_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2016_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 266 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2016_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2016_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 267 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2017_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2017_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 268 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2017_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2017_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 269 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2017_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2017_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 270 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2017_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2017_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 271 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2017_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2017_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 272 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2017_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2017_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 273 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2017_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2017_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 274 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2017_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2017_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 275 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2018_bin0_NPL --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2018_bin0_NPL --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 276 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2018_bin0_PrivMC_tWZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2018_bin0_PrivMC_tWZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 277 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2018_bin0_PrivMC_tZq --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2018_bin0_PrivMC_tZq --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 278 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2018_bin0_PrivMC_ttZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2018_bin0_PrivMC_ttZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 279 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2018_bin0_VVV --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2018_bin0_VVV --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 280 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2018_bin0_WZ --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2018_bin0_WZ --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 281 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2018_bin0_XG --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2018_bin0_XG --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 282 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_NN_cpq3_SRtZq_2018_bin0_tX --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_NN_cpq3_SRtZq_2018_bin0_tX --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 283 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 284 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 285 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin8_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin8_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 286 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin9_mTW_SRother_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin9_mTW_SRother_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 287 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin9_mTW_SRother_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin9_mTW_SRother_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 288 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_binbin9_mTW_SRother_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_binbin9_mTW_SRother_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 289 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_bincountExp_CRWZ_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_bincountExp_CRWZ_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 290 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_bincountExp_CRWZ_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_bincountExp_CRWZ_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 291 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_bincountExp_CRWZ_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_bincountExp_CRWZ_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 292 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_bincountExp_CRZZ_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_bincountExp_CRZZ_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 293 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_bincountExp_CRZZ_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_bincountExp_CRZZ_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 294 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_bincountExp_CRZZ_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_bincountExp_CRZZ_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 295 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_bincountExp_SRttZ4l_2016_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_bincountExp_SRttZ4l_2016_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 296 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_bincountExp_SRttZ4l_2017_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_bincountExp_SRttZ4l_2017_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 297 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_prop_bincountExp_SRttZ4l_2018_bin0 --algo impact --redefineSignalPOIs cpq3 -P prop_bincountExp_SRttZ4l_2018_bin0 --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi
if [ $1 -eq 298 ]; then
  combine -M MultiDimFit -n _paramFit_cpq3Obs_tX_rate --algo impact --redefineSignalPOIs cpq3 -P tX_rate --floatOtherPOIs 1 --saveInactivePOI 1 --robustFit 1 --cminPoiOnlyFit --freezeParameters ctz,ctw,cpqm,cpt --setParameterRanges ctz=-4,4:ctw=-4,4:cpq3=-8,8:cpqm=-25,25:cpt=-25,25 --autoBoundsPOIs=cpq3 --autoMaxPOIs=cpq3 -m 125 -d ./EFTWorkspace_cpq3.root --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0
fi