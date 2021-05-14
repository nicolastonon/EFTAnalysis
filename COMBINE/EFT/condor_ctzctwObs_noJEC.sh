#!/bin/sh
ulimit -s unlimited
set -e
cd /nfs/dust/cms/user/ntonon/TopEFT/CMSSW_10_2_24/src
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /nfs/dust/cms/user/ntonon/TopEFT/CMSSW_10_2_24/src/EFTAnalysis/COMBINE/EFT

if [ $1 -eq 0 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 0 --lastPoint 39 -n .ctzctwObs_noJEC.POINTS.0.39
fi
if [ $1 -eq 1 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 40 --lastPoint 79 -n .ctzctwObs_noJEC.POINTS.40.79
fi
if [ $1 -eq 2 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 80 --lastPoint 119 -n .ctzctwObs_noJEC.POINTS.80.119
fi
if [ $1 -eq 3 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 120 --lastPoint 159 -n .ctzctwObs_noJEC.POINTS.120.159
fi
if [ $1 -eq 4 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 160 --lastPoint 199 -n .ctzctwObs_noJEC.POINTS.160.199
fi
if [ $1 -eq 5 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 200 --lastPoint 239 -n .ctzctwObs_noJEC.POINTS.200.239
fi
if [ $1 -eq 6 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 240 --lastPoint 279 -n .ctzctwObs_noJEC.POINTS.240.279
fi
if [ $1 -eq 7 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 280 --lastPoint 319 -n .ctzctwObs_noJEC.POINTS.280.319
fi
if [ $1 -eq 8 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 320 --lastPoint 359 -n .ctzctwObs_noJEC.POINTS.320.359
fi
if [ $1 -eq 9 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 360 --lastPoint 399 -n .ctzctwObs_noJEC.POINTS.360.399
fi
if [ $1 -eq 10 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 400 --lastPoint 439 -n .ctzctwObs_noJEC.POINTS.400.439
fi
if [ $1 -eq 11 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 440 --lastPoint 479 -n .ctzctwObs_noJEC.POINTS.440.479
fi
if [ $1 -eq 12 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 480 --lastPoint 519 -n .ctzctwObs_noJEC.POINTS.480.519
fi
if [ $1 -eq 13 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 520 --lastPoint 559 -n .ctzctwObs_noJEC.POINTS.520.559
fi
if [ $1 -eq 14 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 560 --lastPoint 599 -n .ctzctwObs_noJEC.POINTS.560.599
fi
if [ $1 -eq 15 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 600 --lastPoint 639 -n .ctzctwObs_noJEC.POINTS.600.639
fi
if [ $1 -eq 16 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 640 --lastPoint 679 -n .ctzctwObs_noJEC.POINTS.640.679
fi
if [ $1 -eq 17 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 680 --lastPoint 719 -n .ctzctwObs_noJEC.POINTS.680.719
fi
if [ $1 -eq 18 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 720 --lastPoint 759 -n .ctzctwObs_noJEC.POINTS.720.759
fi
if [ $1 -eq 19 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 760 --lastPoint 799 -n .ctzctwObs_noJEC.POINTS.760.799
fi
if [ $1 -eq 20 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 800 --lastPoint 839 -n .ctzctwObs_noJEC.POINTS.800.839
fi
if [ $1 -eq 21 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 840 --lastPoint 879 -n .ctzctwObs_noJEC.POINTS.840.879
fi
if [ $1 -eq 22 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 880 --lastPoint 919 -n .ctzctwObs_noJEC.POINTS.880.919
fi
if [ $1 -eq 23 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 920 --lastPoint 959 -n .ctzctwObs_noJEC.POINTS.920.959
fi
if [ $1 -eq 24 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 960 --lastPoint 999 -n .ctzctwObs_noJEC.POINTS.960.999
fi
if [ $1 -eq 25 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1000 --lastPoint 1039 -n .ctzctwObs_noJEC.POINTS.1000.1039
fi
if [ $1 -eq 26 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1040 --lastPoint 1079 -n .ctzctwObs_noJEC.POINTS.1040.1079
fi
if [ $1 -eq 27 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1080 --lastPoint 1119 -n .ctzctwObs_noJEC.POINTS.1080.1119
fi
if [ $1 -eq 28 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1120 --lastPoint 1159 -n .ctzctwObs_noJEC.POINTS.1120.1159
fi
if [ $1 -eq 29 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1160 --lastPoint 1199 -n .ctzctwObs_noJEC.POINTS.1160.1199
fi
if [ $1 -eq 30 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1200 --lastPoint 1239 -n .ctzctwObs_noJEC.POINTS.1200.1239
fi
if [ $1 -eq 31 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1240 --lastPoint 1279 -n .ctzctwObs_noJEC.POINTS.1240.1279
fi
if [ $1 -eq 32 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1280 --lastPoint 1319 -n .ctzctwObs_noJEC.POINTS.1280.1319
fi
if [ $1 -eq 33 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1320 --lastPoint 1359 -n .ctzctwObs_noJEC.POINTS.1320.1359
fi
if [ $1 -eq 34 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1360 --lastPoint 1399 -n .ctzctwObs_noJEC.POINTS.1360.1399
fi
if [ $1 -eq 35 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1400 --lastPoint 1439 -n .ctzctwObs_noJEC.POINTS.1400.1439
fi
if [ $1 -eq 36 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1440 --lastPoint 1479 -n .ctzctwObs_noJEC.POINTS.1440.1479
fi
if [ $1 -eq 37 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1480 --lastPoint 1519 -n .ctzctwObs_noJEC.POINTS.1480.1519
fi
if [ $1 -eq 38 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1520 --lastPoint 1559 -n .ctzctwObs_noJEC.POINTS.1520.1559
fi
if [ $1 -eq 39 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1560 --lastPoint 1599 -n .ctzctwObs_noJEC.POINTS.1560.1599
fi
if [ $1 -eq 40 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1600 --lastPoint 1639 -n .ctzctwObs_noJEC.POINTS.1600.1639
fi
if [ $1 -eq 41 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1640 --lastPoint 1679 -n .ctzctwObs_noJEC.POINTS.1640.1679
fi
if [ $1 -eq 42 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1680 --lastPoint 1719 -n .ctzctwObs_noJEC.POINTS.1680.1719
fi
if [ $1 -eq 43 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1720 --lastPoint 1759 -n .ctzctwObs_noJEC.POINTS.1720.1759
fi
if [ $1 -eq 44 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1760 --lastPoint 1799 -n .ctzctwObs_noJEC.POINTS.1760.1799
fi
if [ $1 -eq 45 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1800 --lastPoint 1839 -n .ctzctwObs_noJEC.POINTS.1800.1839
fi
if [ $1 -eq 46 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1840 --lastPoint 1879 -n .ctzctwObs_noJEC.POINTS.1840.1879
fi
if [ $1 -eq 47 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1880 --lastPoint 1919 -n .ctzctwObs_noJEC.POINTS.1880.1919
fi
if [ $1 -eq 48 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1920 --lastPoint 1959 -n .ctzctwObs_noJEC.POINTS.1920.1959
fi
if [ $1 -eq 49 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 1960 --lastPoint 1999 -n .ctzctwObs_noJEC.POINTS.1960.1999
fi
if [ $1 -eq 50 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2000 --lastPoint 2039 -n .ctzctwObs_noJEC.POINTS.2000.2039
fi
if [ $1 -eq 51 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2040 --lastPoint 2079 -n .ctzctwObs_noJEC.POINTS.2040.2079
fi
if [ $1 -eq 52 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2080 --lastPoint 2119 -n .ctzctwObs_noJEC.POINTS.2080.2119
fi
if [ $1 -eq 53 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2120 --lastPoint 2159 -n .ctzctwObs_noJEC.POINTS.2120.2159
fi
if [ $1 -eq 54 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2160 --lastPoint 2199 -n .ctzctwObs_noJEC.POINTS.2160.2199
fi
if [ $1 -eq 55 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2200 --lastPoint 2239 -n .ctzctwObs_noJEC.POINTS.2200.2239
fi
if [ $1 -eq 56 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2240 --lastPoint 2279 -n .ctzctwObs_noJEC.POINTS.2240.2279
fi
if [ $1 -eq 57 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2280 --lastPoint 2319 -n .ctzctwObs_noJEC.POINTS.2280.2319
fi
if [ $1 -eq 58 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2320 --lastPoint 2359 -n .ctzctwObs_noJEC.POINTS.2320.2359
fi
if [ $1 -eq 59 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2360 --lastPoint 2399 -n .ctzctwObs_noJEC.POINTS.2360.2399
fi
if [ $1 -eq 60 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2400 --lastPoint 2439 -n .ctzctwObs_noJEC.POINTS.2400.2439
fi
if [ $1 -eq 61 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2440 --lastPoint 2479 -n .ctzctwObs_noJEC.POINTS.2440.2479
fi
if [ $1 -eq 62 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2480 --lastPoint 2519 -n .ctzctwObs_noJEC.POINTS.2480.2519
fi
if [ $1 -eq 63 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2520 --lastPoint 2559 -n .ctzctwObs_noJEC.POINTS.2520.2559
fi
if [ $1 -eq 64 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2560 --lastPoint 2599 -n .ctzctwObs_noJEC.POINTS.2560.2599
fi
if [ $1 -eq 65 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2600 --lastPoint 2639 -n .ctzctwObs_noJEC.POINTS.2600.2639
fi
if [ $1 -eq 66 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2640 --lastPoint 2679 -n .ctzctwObs_noJEC.POINTS.2640.2679
fi
if [ $1 -eq 67 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2680 --lastPoint 2719 -n .ctzctwObs_noJEC.POINTS.2680.2719
fi
if [ $1 -eq 68 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2720 --lastPoint 2759 -n .ctzctwObs_noJEC.POINTS.2720.2759
fi
if [ $1 -eq 69 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2760 --lastPoint 2799 -n .ctzctwObs_noJEC.POINTS.2760.2799
fi
if [ $1 -eq 70 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2800 --lastPoint 2839 -n .ctzctwObs_noJEC.POINTS.2800.2839
fi
if [ $1 -eq 71 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2840 --lastPoint 2879 -n .ctzctwObs_noJEC.POINTS.2840.2879
fi
if [ $1 -eq 72 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2880 --lastPoint 2919 -n .ctzctwObs_noJEC.POINTS.2880.2919
fi
if [ $1 -eq 73 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2920 --lastPoint 2959 -n .ctzctwObs_noJEC.POINTS.2920.2959
fi
if [ $1 -eq 74 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 2960 --lastPoint 2999 -n .ctzctwObs_noJEC.POINTS.2960.2999
fi
if [ $1 -eq 75 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3000 --lastPoint 3039 -n .ctzctwObs_noJEC.POINTS.3000.3039
fi
if [ $1 -eq 76 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3040 --lastPoint 3079 -n .ctzctwObs_noJEC.POINTS.3040.3079
fi
if [ $1 -eq 77 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3080 --lastPoint 3119 -n .ctzctwObs_noJEC.POINTS.3080.3119
fi
if [ $1 -eq 78 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3120 --lastPoint 3159 -n .ctzctwObs_noJEC.POINTS.3120.3159
fi
if [ $1 -eq 79 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3160 --lastPoint 3199 -n .ctzctwObs_noJEC.POINTS.3160.3199
fi
if [ $1 -eq 80 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3200 --lastPoint 3239 -n .ctzctwObs_noJEC.POINTS.3200.3239
fi
if [ $1 -eq 81 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3240 --lastPoint 3279 -n .ctzctwObs_noJEC.POINTS.3240.3279
fi
if [ $1 -eq 82 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3280 --lastPoint 3319 -n .ctzctwObs_noJEC.POINTS.3280.3319
fi
if [ $1 -eq 83 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3320 --lastPoint 3359 -n .ctzctwObs_noJEC.POINTS.3320.3359
fi
if [ $1 -eq 84 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3360 --lastPoint 3399 -n .ctzctwObs_noJEC.POINTS.3360.3399
fi
if [ $1 -eq 85 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3400 --lastPoint 3439 -n .ctzctwObs_noJEC.POINTS.3400.3439
fi
if [ $1 -eq 86 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3440 --lastPoint 3479 -n .ctzctwObs_noJEC.POINTS.3440.3479
fi
if [ $1 -eq 87 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3480 --lastPoint 3519 -n .ctzctwObs_noJEC.POINTS.3480.3519
fi
if [ $1 -eq 88 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3520 --lastPoint 3559 -n .ctzctwObs_noJEC.POINTS.3520.3559
fi
if [ $1 -eq 89 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3560 --lastPoint 3599 -n .ctzctwObs_noJEC.POINTS.3560.3599
fi
if [ $1 -eq 90 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3600 --lastPoint 3639 -n .ctzctwObs_noJEC.POINTS.3600.3639
fi
if [ $1 -eq 91 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3640 --lastPoint 3679 -n .ctzctwObs_noJEC.POINTS.3640.3679
fi
if [ $1 -eq 92 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3680 --lastPoint 3719 -n .ctzctwObs_noJEC.POINTS.3680.3719
fi
if [ $1 -eq 93 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3720 --lastPoint 3759 -n .ctzctwObs_noJEC.POINTS.3720.3759
fi
if [ $1 -eq 94 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3760 --lastPoint 3799 -n .ctzctwObs_noJEC.POINTS.3760.3799
fi
if [ $1 -eq 95 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3800 --lastPoint 3839 -n .ctzctwObs_noJEC.POINTS.3800.3839
fi
if [ $1 -eq 96 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3840 --lastPoint 3879 -n .ctzctwObs_noJEC.POINTS.3840.3879
fi
if [ $1 -eq 97 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3880 --lastPoint 3919 -n .ctzctwObs_noJEC.POINTS.3880.3919
fi
if [ $1 -eq 98 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3920 --lastPoint 3959 -n .ctzctwObs_noJEC.POINTS.3920.3959
fi
if [ $1 -eq 99 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 3960 --lastPoint 3999 -n .ctzctwObs_noJEC.POINTS.3960.3999
fi
if [ $1 -eq 100 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4000 --lastPoint 4039 -n .ctzctwObs_noJEC.POINTS.4000.4039
fi
if [ $1 -eq 101 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4040 --lastPoint 4079 -n .ctzctwObs_noJEC.POINTS.4040.4079
fi
if [ $1 -eq 102 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4080 --lastPoint 4119 -n .ctzctwObs_noJEC.POINTS.4080.4119
fi
if [ $1 -eq 103 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4120 --lastPoint 4159 -n .ctzctwObs_noJEC.POINTS.4120.4159
fi
if [ $1 -eq 104 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4160 --lastPoint 4199 -n .ctzctwObs_noJEC.POINTS.4160.4199
fi
if [ $1 -eq 105 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4200 --lastPoint 4239 -n .ctzctwObs_noJEC.POINTS.4200.4239
fi
if [ $1 -eq 106 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4240 --lastPoint 4279 -n .ctzctwObs_noJEC.POINTS.4240.4279
fi
if [ $1 -eq 107 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4280 --lastPoint 4319 -n .ctzctwObs_noJEC.POINTS.4280.4319
fi
if [ $1 -eq 108 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4320 --lastPoint 4359 -n .ctzctwObs_noJEC.POINTS.4320.4359
fi
if [ $1 -eq 109 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4360 --lastPoint 4399 -n .ctzctwObs_noJEC.POINTS.4360.4399
fi
if [ $1 -eq 110 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4400 --lastPoint 4439 -n .ctzctwObs_noJEC.POINTS.4400.4439
fi
if [ $1 -eq 111 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4440 --lastPoint 4479 -n .ctzctwObs_noJEC.POINTS.4440.4479
fi
if [ $1 -eq 112 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4480 --lastPoint 4519 -n .ctzctwObs_noJEC.POINTS.4480.4519
fi
if [ $1 -eq 113 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4520 --lastPoint 4559 -n .ctzctwObs_noJEC.POINTS.4520.4559
fi
if [ $1 -eq 114 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4560 --lastPoint 4599 -n .ctzctwObs_noJEC.POINTS.4560.4599
fi
if [ $1 -eq 115 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4600 --lastPoint 4639 -n .ctzctwObs_noJEC.POINTS.4600.4639
fi
if [ $1 -eq 116 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4640 --lastPoint 4679 -n .ctzctwObs_noJEC.POINTS.4640.4679
fi
if [ $1 -eq 117 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4680 --lastPoint 4719 -n .ctzctwObs_noJEC.POINTS.4680.4719
fi
if [ $1 -eq 118 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4720 --lastPoint 4759 -n .ctzctwObs_noJEC.POINTS.4720.4759
fi
if [ $1 -eq 119 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4760 --lastPoint 4799 -n .ctzctwObs_noJEC.POINTS.4760.4799
fi
if [ $1 -eq 120 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4800 --lastPoint 4839 -n .ctzctwObs_noJEC.POINTS.4800.4839
fi
if [ $1 -eq 121 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4840 --lastPoint 4879 -n .ctzctwObs_noJEC.POINTS.4840.4879
fi
if [ $1 -eq 122 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4880 --lastPoint 4919 -n .ctzctwObs_noJEC.POINTS.4880.4919
fi
if [ $1 -eq 123 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4920 --lastPoint 4959 -n .ctzctwObs_noJEC.POINTS.4920.4959
fi
if [ $1 -eq 124 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 4960 --lastPoint 4999 -n .ctzctwObs_noJEC.POINTS.4960.4999
fi
if [ $1 -eq 125 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5000 --lastPoint 5039 -n .ctzctwObs_noJEC.POINTS.5000.5039
fi
if [ $1 -eq 126 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5040 --lastPoint 5079 -n .ctzctwObs_noJEC.POINTS.5040.5079
fi
if [ $1 -eq 127 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5080 --lastPoint 5119 -n .ctzctwObs_noJEC.POINTS.5080.5119
fi
if [ $1 -eq 128 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5120 --lastPoint 5159 -n .ctzctwObs_noJEC.POINTS.5120.5159
fi
if [ $1 -eq 129 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5160 --lastPoint 5199 -n .ctzctwObs_noJEC.POINTS.5160.5199
fi
if [ $1 -eq 130 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5200 --lastPoint 5239 -n .ctzctwObs_noJEC.POINTS.5200.5239
fi
if [ $1 -eq 131 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5240 --lastPoint 5279 -n .ctzctwObs_noJEC.POINTS.5240.5279
fi
if [ $1 -eq 132 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5280 --lastPoint 5319 -n .ctzctwObs_noJEC.POINTS.5280.5319
fi
if [ $1 -eq 133 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5320 --lastPoint 5359 -n .ctzctwObs_noJEC.POINTS.5320.5359
fi
if [ $1 -eq 134 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5360 --lastPoint 5399 -n .ctzctwObs_noJEC.POINTS.5360.5399
fi
if [ $1 -eq 135 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5400 --lastPoint 5439 -n .ctzctwObs_noJEC.POINTS.5400.5439
fi
if [ $1 -eq 136 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5440 --lastPoint 5479 -n .ctzctwObs_noJEC.POINTS.5440.5479
fi
if [ $1 -eq 137 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5480 --lastPoint 5519 -n .ctzctwObs_noJEC.POINTS.5480.5519
fi
if [ $1 -eq 138 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5520 --lastPoint 5559 -n .ctzctwObs_noJEC.POINTS.5520.5559
fi
if [ $1 -eq 139 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5560 --lastPoint 5599 -n .ctzctwObs_noJEC.POINTS.5560.5599
fi
if [ $1 -eq 140 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5600 --lastPoint 5639 -n .ctzctwObs_noJEC.POINTS.5600.5639
fi
if [ $1 -eq 141 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5640 --lastPoint 5679 -n .ctzctwObs_noJEC.POINTS.5640.5679
fi
if [ $1 -eq 142 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5680 --lastPoint 5719 -n .ctzctwObs_noJEC.POINTS.5680.5719
fi
if [ $1 -eq 143 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5720 --lastPoint 5759 -n .ctzctwObs_noJEC.POINTS.5720.5759
fi
if [ $1 -eq 144 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5760 --lastPoint 5799 -n .ctzctwObs_noJEC.POINTS.5760.5799
fi
if [ $1 -eq 145 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5800 --lastPoint 5839 -n .ctzctwObs_noJEC.POINTS.5800.5839
fi
if [ $1 -eq 146 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5840 --lastPoint 5879 -n .ctzctwObs_noJEC.POINTS.5840.5879
fi
if [ $1 -eq 147 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5880 --lastPoint 5919 -n .ctzctwObs_noJEC.POINTS.5880.5919
fi
if [ $1 -eq 148 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5920 --lastPoint 5959 -n .ctzctwObs_noJEC.POINTS.5920.5959
fi
if [ $1 -eq 149 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 5960 --lastPoint 5999 -n .ctzctwObs_noJEC.POINTS.5960.5999
fi
if [ $1 -eq 150 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6000 --lastPoint 6039 -n .ctzctwObs_noJEC.POINTS.6000.6039
fi
if [ $1 -eq 151 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6040 --lastPoint 6079 -n .ctzctwObs_noJEC.POINTS.6040.6079
fi
if [ $1 -eq 152 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6080 --lastPoint 6119 -n .ctzctwObs_noJEC.POINTS.6080.6119
fi
if [ $1 -eq 153 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6120 --lastPoint 6159 -n .ctzctwObs_noJEC.POINTS.6120.6159
fi
if [ $1 -eq 154 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6160 --lastPoint 6199 -n .ctzctwObs_noJEC.POINTS.6160.6199
fi
if [ $1 -eq 155 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6200 --lastPoint 6239 -n .ctzctwObs_noJEC.POINTS.6200.6239
fi
if [ $1 -eq 156 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6240 --lastPoint 6279 -n .ctzctwObs_noJEC.POINTS.6240.6279
fi
if [ $1 -eq 157 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6280 --lastPoint 6319 -n .ctzctwObs_noJEC.POINTS.6280.6319
fi
if [ $1 -eq 158 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6320 --lastPoint 6359 -n .ctzctwObs_noJEC.POINTS.6320.6359
fi
if [ $1 -eq 159 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6360 --lastPoint 6399 -n .ctzctwObs_noJEC.POINTS.6360.6399
fi
if [ $1 -eq 160 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6400 --lastPoint 6439 -n .ctzctwObs_noJEC.POINTS.6400.6439
fi
if [ $1 -eq 161 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6440 --lastPoint 6479 -n .ctzctwObs_noJEC.POINTS.6440.6479
fi
if [ $1 -eq 162 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6480 --lastPoint 6519 -n .ctzctwObs_noJEC.POINTS.6480.6519
fi
if [ $1 -eq 163 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6520 --lastPoint 6559 -n .ctzctwObs_noJEC.POINTS.6520.6559
fi
if [ $1 -eq 164 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6560 --lastPoint 6599 -n .ctzctwObs_noJEC.POINTS.6560.6599
fi
if [ $1 -eq 165 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6600 --lastPoint 6639 -n .ctzctwObs_noJEC.POINTS.6600.6639
fi
if [ $1 -eq 166 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6640 --lastPoint 6679 -n .ctzctwObs_noJEC.POINTS.6640.6679
fi
if [ $1 -eq 167 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6680 --lastPoint 6719 -n .ctzctwObs_noJEC.POINTS.6680.6719
fi
if [ $1 -eq 168 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6720 --lastPoint 6759 -n .ctzctwObs_noJEC.POINTS.6720.6759
fi
if [ $1 -eq 169 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6760 --lastPoint 6799 -n .ctzctwObs_noJEC.POINTS.6760.6799
fi
if [ $1 -eq 170 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6800 --lastPoint 6839 -n .ctzctwObs_noJEC.POINTS.6800.6839
fi
if [ $1 -eq 171 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6840 --lastPoint 6879 -n .ctzctwObs_noJEC.POINTS.6840.6879
fi
if [ $1 -eq 172 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6880 --lastPoint 6919 -n .ctzctwObs_noJEC.POINTS.6880.6919
fi
if [ $1 -eq 173 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6920 --lastPoint 6959 -n .ctzctwObs_noJEC.POINTS.6920.6959
fi
if [ $1 -eq 174 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 6960 --lastPoint 6999 -n .ctzctwObs_noJEC.POINTS.6960.6999
fi
if [ $1 -eq 175 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7000 --lastPoint 7039 -n .ctzctwObs_noJEC.POINTS.7000.7039
fi
if [ $1 -eq 176 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7040 --lastPoint 7079 -n .ctzctwObs_noJEC.POINTS.7040.7079
fi
if [ $1 -eq 177 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7080 --lastPoint 7119 -n .ctzctwObs_noJEC.POINTS.7080.7119
fi
if [ $1 -eq 178 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7120 --lastPoint 7159 -n .ctzctwObs_noJEC.POINTS.7120.7159
fi
if [ $1 -eq 179 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7160 --lastPoint 7199 -n .ctzctwObs_noJEC.POINTS.7160.7199
fi
if [ $1 -eq 180 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7200 --lastPoint 7239 -n .ctzctwObs_noJEC.POINTS.7200.7239
fi
if [ $1 -eq 181 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7240 --lastPoint 7279 -n .ctzctwObs_noJEC.POINTS.7240.7279
fi
if [ $1 -eq 182 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7280 --lastPoint 7319 -n .ctzctwObs_noJEC.POINTS.7280.7319
fi
if [ $1 -eq 183 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7320 --lastPoint 7359 -n .ctzctwObs_noJEC.POINTS.7320.7359
fi
if [ $1 -eq 184 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7360 --lastPoint 7399 -n .ctzctwObs_noJEC.POINTS.7360.7399
fi
if [ $1 -eq 185 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7400 --lastPoint 7439 -n .ctzctwObs_noJEC.POINTS.7400.7439
fi
if [ $1 -eq 186 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7440 --lastPoint 7479 -n .ctzctwObs_noJEC.POINTS.7440.7479
fi
if [ $1 -eq 187 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7480 --lastPoint 7519 -n .ctzctwObs_noJEC.POINTS.7480.7519
fi
if [ $1 -eq 188 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7520 --lastPoint 7559 -n .ctzctwObs_noJEC.POINTS.7520.7559
fi
if [ $1 -eq 189 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7560 --lastPoint 7599 -n .ctzctwObs_noJEC.POINTS.7560.7599
fi
if [ $1 -eq 190 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7600 --lastPoint 7639 -n .ctzctwObs_noJEC.POINTS.7600.7639
fi
if [ $1 -eq 191 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7640 --lastPoint 7679 -n .ctzctwObs_noJEC.POINTS.7640.7679
fi
if [ $1 -eq 192 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7680 --lastPoint 7719 -n .ctzctwObs_noJEC.POINTS.7680.7719
fi
if [ $1 -eq 193 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7720 --lastPoint 7759 -n .ctzctwObs_noJEC.POINTS.7720.7759
fi
if [ $1 -eq 194 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7760 --lastPoint 7799 -n .ctzctwObs_noJEC.POINTS.7760.7799
fi
if [ $1 -eq 195 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7800 --lastPoint 7839 -n .ctzctwObs_noJEC.POINTS.7800.7839
fi
if [ $1 -eq 196 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7840 --lastPoint 7879 -n .ctzctwObs_noJEC.POINTS.7840.7879
fi
if [ $1 -eq 197 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7880 --lastPoint 7919 -n .ctzctwObs_noJEC.POINTS.7880.7919
fi
if [ $1 -eq 198 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7920 --lastPoint 7959 -n .ctzctwObs_noJEC.POINTS.7920.7959
fi
if [ $1 -eq 199 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 7960 --lastPoint 7999 -n .ctzctwObs_noJEC.POINTS.7960.7999
fi
if [ $1 -eq 200 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8000 --lastPoint 8039 -n .ctzctwObs_noJEC.POINTS.8000.8039
fi
if [ $1 -eq 201 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8040 --lastPoint 8079 -n .ctzctwObs_noJEC.POINTS.8040.8079
fi
if [ $1 -eq 202 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8080 --lastPoint 8119 -n .ctzctwObs_noJEC.POINTS.8080.8119
fi
if [ $1 -eq 203 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8120 --lastPoint 8159 -n .ctzctwObs_noJEC.POINTS.8120.8159
fi
if [ $1 -eq 204 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8160 --lastPoint 8199 -n .ctzctwObs_noJEC.POINTS.8160.8199
fi
if [ $1 -eq 205 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8200 --lastPoint 8239 -n .ctzctwObs_noJEC.POINTS.8200.8239
fi
if [ $1 -eq 206 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8240 --lastPoint 8279 -n .ctzctwObs_noJEC.POINTS.8240.8279
fi
if [ $1 -eq 207 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8280 --lastPoint 8319 -n .ctzctwObs_noJEC.POINTS.8280.8319
fi
if [ $1 -eq 208 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8320 --lastPoint 8359 -n .ctzctwObs_noJEC.POINTS.8320.8359
fi
if [ $1 -eq 209 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8360 --lastPoint 8399 -n .ctzctwObs_noJEC.POINTS.8360.8399
fi
if [ $1 -eq 210 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8400 --lastPoint 8439 -n .ctzctwObs_noJEC.POINTS.8400.8439
fi
if [ $1 -eq 211 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8440 --lastPoint 8479 -n .ctzctwObs_noJEC.POINTS.8440.8479
fi
if [ $1 -eq 212 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8480 --lastPoint 8519 -n .ctzctwObs_noJEC.POINTS.8480.8519
fi
if [ $1 -eq 213 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8520 --lastPoint 8559 -n .ctzctwObs_noJEC.POINTS.8520.8559
fi
if [ $1 -eq 214 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8560 --lastPoint 8599 -n .ctzctwObs_noJEC.POINTS.8560.8599
fi
if [ $1 -eq 215 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8600 --lastPoint 8639 -n .ctzctwObs_noJEC.POINTS.8600.8639
fi
if [ $1 -eq 216 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8640 --lastPoint 8679 -n .ctzctwObs_noJEC.POINTS.8640.8679
fi
if [ $1 -eq 217 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8680 --lastPoint 8719 -n .ctzctwObs_noJEC.POINTS.8680.8719
fi
if [ $1 -eq 218 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8720 --lastPoint 8759 -n .ctzctwObs_noJEC.POINTS.8720.8759
fi
if [ $1 -eq 219 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8760 --lastPoint 8799 -n .ctzctwObs_noJEC.POINTS.8760.8799
fi
if [ $1 -eq 220 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8800 --lastPoint 8839 -n .ctzctwObs_noJEC.POINTS.8800.8839
fi
if [ $1 -eq 221 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8840 --lastPoint 8879 -n .ctzctwObs_noJEC.POINTS.8840.8879
fi
if [ $1 -eq 222 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8880 --lastPoint 8919 -n .ctzctwObs_noJEC.POINTS.8880.8919
fi
if [ $1 -eq 223 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8920 --lastPoint 8959 -n .ctzctwObs_noJEC.POINTS.8920.8959
fi
if [ $1 -eq 224 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 8960 --lastPoint 8999 -n .ctzctwObs_noJEC.POINTS.8960.8999
fi
if [ $1 -eq 225 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9000 --lastPoint 9039 -n .ctzctwObs_noJEC.POINTS.9000.9039
fi
if [ $1 -eq 226 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9040 --lastPoint 9079 -n .ctzctwObs_noJEC.POINTS.9040.9079
fi
if [ $1 -eq 227 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9080 --lastPoint 9119 -n .ctzctwObs_noJEC.POINTS.9080.9119
fi
if [ $1 -eq 228 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9120 --lastPoint 9159 -n .ctzctwObs_noJEC.POINTS.9120.9159
fi
if [ $1 -eq 229 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9160 --lastPoint 9199 -n .ctzctwObs_noJEC.POINTS.9160.9199
fi
if [ $1 -eq 230 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9200 --lastPoint 9239 -n .ctzctwObs_noJEC.POINTS.9200.9239
fi
if [ $1 -eq 231 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9240 --lastPoint 9279 -n .ctzctwObs_noJEC.POINTS.9240.9279
fi
if [ $1 -eq 232 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9280 --lastPoint 9319 -n .ctzctwObs_noJEC.POINTS.9280.9319
fi
if [ $1 -eq 233 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9320 --lastPoint 9359 -n .ctzctwObs_noJEC.POINTS.9320.9359
fi
if [ $1 -eq 234 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9360 --lastPoint 9399 -n .ctzctwObs_noJEC.POINTS.9360.9399
fi
if [ $1 -eq 235 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9400 --lastPoint 9439 -n .ctzctwObs_noJEC.POINTS.9400.9439
fi
if [ $1 -eq 236 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9440 --lastPoint 9479 -n .ctzctwObs_noJEC.POINTS.9440.9479
fi
if [ $1 -eq 237 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9480 --lastPoint 9519 -n .ctzctwObs_noJEC.POINTS.9480.9519
fi
if [ $1 -eq 238 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9520 --lastPoint 9559 -n .ctzctwObs_noJEC.POINTS.9520.9559
fi
if [ $1 -eq 239 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9560 --lastPoint 9599 -n .ctzctwObs_noJEC.POINTS.9560.9599
fi
if [ $1 -eq 240 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9600 --lastPoint 9639 -n .ctzctwObs_noJEC.POINTS.9600.9639
fi
if [ $1 -eq 241 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9640 --lastPoint 9679 -n .ctzctwObs_noJEC.POINTS.9640.9679
fi
if [ $1 -eq 242 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9680 --lastPoint 9719 -n .ctzctwObs_noJEC.POINTS.9680.9719
fi
if [ $1 -eq 243 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9720 --lastPoint 9759 -n .ctzctwObs_noJEC.POINTS.9720.9759
fi
if [ $1 -eq 244 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9760 --lastPoint 9799 -n .ctzctwObs_noJEC.POINTS.9760.9799
fi
if [ $1 -eq 245 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9800 --lastPoint 9839 -n .ctzctwObs_noJEC.POINTS.9800.9839
fi
if [ $1 -eq 246 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9840 --lastPoint 9879 -n .ctzctwObs_noJEC.POINTS.9840.9879
fi
if [ $1 -eq 247 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9880 --lastPoint 9919 -n .ctzctwObs_noJEC.POINTS.9880.9919
fi
if [ $1 -eq 248 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9920 --lastPoint 9959 -n .ctzctwObs_noJEC.POINTS.9920.9959
fi
if [ $1 -eq 249 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 9960 --lastPoint 9999 -n .ctzctwObs_noJEC.POINTS.9960.9999
fi
if [ $1 -eq 250 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10000 --lastPoint 10039 -n .ctzctwObs_noJEC.POINTS.10000.10039
fi
if [ $1 -eq 251 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10040 --lastPoint 10079 -n .ctzctwObs_noJEC.POINTS.10040.10079
fi
if [ $1 -eq 252 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10080 --lastPoint 10119 -n .ctzctwObs_noJEC.POINTS.10080.10119
fi
if [ $1 -eq 253 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10120 --lastPoint 10159 -n .ctzctwObs_noJEC.POINTS.10120.10159
fi
if [ $1 -eq 254 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10160 --lastPoint 10199 -n .ctzctwObs_noJEC.POINTS.10160.10199
fi
if [ $1 -eq 255 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10200 --lastPoint 10239 -n .ctzctwObs_noJEC.POINTS.10200.10239
fi
if [ $1 -eq 256 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10240 --lastPoint 10279 -n .ctzctwObs_noJEC.POINTS.10240.10279
fi
if [ $1 -eq 257 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10280 --lastPoint 10319 -n .ctzctwObs_noJEC.POINTS.10280.10319
fi
if [ $1 -eq 258 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10320 --lastPoint 10359 -n .ctzctwObs_noJEC.POINTS.10320.10359
fi
if [ $1 -eq 259 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10360 --lastPoint 10399 -n .ctzctwObs_noJEC.POINTS.10360.10399
fi
if [ $1 -eq 260 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10400 --lastPoint 10439 -n .ctzctwObs_noJEC.POINTS.10400.10439
fi
if [ $1 -eq 261 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10440 --lastPoint 10479 -n .ctzctwObs_noJEC.POINTS.10440.10479
fi
if [ $1 -eq 262 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10480 --lastPoint 10519 -n .ctzctwObs_noJEC.POINTS.10480.10519
fi
if [ $1 -eq 263 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10520 --lastPoint 10559 -n .ctzctwObs_noJEC.POINTS.10520.10559
fi
if [ $1 -eq 264 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10560 --lastPoint 10599 -n .ctzctwObs_noJEC.POINTS.10560.10599
fi
if [ $1 -eq 265 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10600 --lastPoint 10639 -n .ctzctwObs_noJEC.POINTS.10600.10639
fi
if [ $1 -eq 266 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10640 --lastPoint 10679 -n .ctzctwObs_noJEC.POINTS.10640.10679
fi
if [ $1 -eq 267 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10680 --lastPoint 10719 -n .ctzctwObs_noJEC.POINTS.10680.10719
fi
if [ $1 -eq 268 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10720 --lastPoint 10759 -n .ctzctwObs_noJEC.POINTS.10720.10759
fi
if [ $1 -eq 269 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10760 --lastPoint 10799 -n .ctzctwObs_noJEC.POINTS.10760.10799
fi
if [ $1 -eq 270 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10800 --lastPoint 10839 -n .ctzctwObs_noJEC.POINTS.10800.10839
fi
if [ $1 -eq 271 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10840 --lastPoint 10879 -n .ctzctwObs_noJEC.POINTS.10840.10879
fi
if [ $1 -eq 272 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10880 --lastPoint 10919 -n .ctzctwObs_noJEC.POINTS.10880.10919
fi
if [ $1 -eq 273 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10920 --lastPoint 10959 -n .ctzctwObs_noJEC.POINTS.10920.10959
fi
if [ $1 -eq 274 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 10960 --lastPoint 10999 -n .ctzctwObs_noJEC.POINTS.10960.10999
fi
if [ $1 -eq 275 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11000 --lastPoint 11039 -n .ctzctwObs_noJEC.POINTS.11000.11039
fi
if [ $1 -eq 276 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11040 --lastPoint 11079 -n .ctzctwObs_noJEC.POINTS.11040.11079
fi
if [ $1 -eq 277 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11080 --lastPoint 11119 -n .ctzctwObs_noJEC.POINTS.11080.11119
fi
if [ $1 -eq 278 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11120 --lastPoint 11159 -n .ctzctwObs_noJEC.POINTS.11120.11159
fi
if [ $1 -eq 279 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11160 --lastPoint 11199 -n .ctzctwObs_noJEC.POINTS.11160.11199
fi
if [ $1 -eq 280 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11200 --lastPoint 11239 -n .ctzctwObs_noJEC.POINTS.11200.11239
fi
if [ $1 -eq 281 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11240 --lastPoint 11279 -n .ctzctwObs_noJEC.POINTS.11240.11279
fi
if [ $1 -eq 282 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11280 --lastPoint 11319 -n .ctzctwObs_noJEC.POINTS.11280.11319
fi
if [ $1 -eq 283 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11320 --lastPoint 11359 -n .ctzctwObs_noJEC.POINTS.11320.11359
fi
if [ $1 -eq 284 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11360 --lastPoint 11399 -n .ctzctwObs_noJEC.POINTS.11360.11399
fi
if [ $1 -eq 285 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11400 --lastPoint 11439 -n .ctzctwObs_noJEC.POINTS.11400.11439
fi
if [ $1 -eq 286 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11440 --lastPoint 11479 -n .ctzctwObs_noJEC.POINTS.11440.11479
fi
if [ $1 -eq 287 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11480 --lastPoint 11519 -n .ctzctwObs_noJEC.POINTS.11480.11519
fi
if [ $1 -eq 288 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11520 --lastPoint 11559 -n .ctzctwObs_noJEC.POINTS.11520.11559
fi
if [ $1 -eq 289 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11560 --lastPoint 11599 -n .ctzctwObs_noJEC.POINTS.11560.11599
fi
if [ $1 -eq 290 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11600 --lastPoint 11639 -n .ctzctwObs_noJEC.POINTS.11600.11639
fi
if [ $1 -eq 291 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11640 --lastPoint 11679 -n .ctzctwObs_noJEC.POINTS.11640.11679
fi
if [ $1 -eq 292 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11680 --lastPoint 11719 -n .ctzctwObs_noJEC.POINTS.11680.11719
fi
if [ $1 -eq 293 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11720 --lastPoint 11759 -n .ctzctwObs_noJEC.POINTS.11720.11759
fi
if [ $1 -eq 294 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11760 --lastPoint 11799 -n .ctzctwObs_noJEC.POINTS.11760.11799
fi
if [ $1 -eq 295 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11800 --lastPoint 11839 -n .ctzctwObs_noJEC.POINTS.11800.11839
fi
if [ $1 -eq 296 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11840 --lastPoint 11879 -n .ctzctwObs_noJEC.POINTS.11840.11879
fi
if [ $1 -eq 297 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11880 --lastPoint 11919 -n .ctzctwObs_noJEC.POINTS.11880.11919
fi
if [ $1 -eq 298 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11920 --lastPoint 11959 -n .ctzctwObs_noJEC.POINTS.11920.11959
fi
if [ $1 -eq 299 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 11960 --lastPoint 11999 -n .ctzctwObs_noJEC.POINTS.11960.11999
fi
if [ $1 -eq 300 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12000 --lastPoint 12039 -n .ctzctwObs_noJEC.POINTS.12000.12039
fi
if [ $1 -eq 301 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12040 --lastPoint 12079 -n .ctzctwObs_noJEC.POINTS.12040.12079
fi
if [ $1 -eq 302 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12080 --lastPoint 12119 -n .ctzctwObs_noJEC.POINTS.12080.12119
fi
if [ $1 -eq 303 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12120 --lastPoint 12159 -n .ctzctwObs_noJEC.POINTS.12120.12159
fi
if [ $1 -eq 304 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12160 --lastPoint 12199 -n .ctzctwObs_noJEC.POINTS.12160.12199
fi
if [ $1 -eq 305 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12200 --lastPoint 12239 -n .ctzctwObs_noJEC.POINTS.12200.12239
fi
if [ $1 -eq 306 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12240 --lastPoint 12279 -n .ctzctwObs_noJEC.POINTS.12240.12279
fi
if [ $1 -eq 307 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12280 --lastPoint 12319 -n .ctzctwObs_noJEC.POINTS.12280.12319
fi
if [ $1 -eq 308 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12320 --lastPoint 12359 -n .ctzctwObs_noJEC.POINTS.12320.12359
fi
if [ $1 -eq 309 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12360 --lastPoint 12399 -n .ctzctwObs_noJEC.POINTS.12360.12399
fi
if [ $1 -eq 310 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12400 --lastPoint 12439 -n .ctzctwObs_noJEC.POINTS.12400.12439
fi
if [ $1 -eq 311 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12440 --lastPoint 12479 -n .ctzctwObs_noJEC.POINTS.12440.12479
fi
if [ $1 -eq 312 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12480 --lastPoint 12519 -n .ctzctwObs_noJEC.POINTS.12480.12519
fi
if [ $1 -eq 313 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12520 --lastPoint 12559 -n .ctzctwObs_noJEC.POINTS.12520.12559
fi
if [ $1 -eq 314 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12560 --lastPoint 12599 -n .ctzctwObs_noJEC.POINTS.12560.12599
fi
if [ $1 -eq 315 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12600 --lastPoint 12639 -n .ctzctwObs_noJEC.POINTS.12600.12639
fi
if [ $1 -eq 316 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12640 --lastPoint 12679 -n .ctzctwObs_noJEC.POINTS.12640.12679
fi
if [ $1 -eq 317 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12680 --lastPoint 12719 -n .ctzctwObs_noJEC.POINTS.12680.12719
fi
if [ $1 -eq 318 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12720 --lastPoint 12759 -n .ctzctwObs_noJEC.POINTS.12720.12759
fi
if [ $1 -eq 319 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12760 --lastPoint 12799 -n .ctzctwObs_noJEC.POINTS.12760.12799
fi
if [ $1 -eq 320 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12800 --lastPoint 12839 -n .ctzctwObs_noJEC.POINTS.12800.12839
fi
if [ $1 -eq 321 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12840 --lastPoint 12879 -n .ctzctwObs_noJEC.POINTS.12840.12879
fi
if [ $1 -eq 322 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12880 --lastPoint 12919 -n .ctzctwObs_noJEC.POINTS.12880.12919
fi
if [ $1 -eq 323 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12920 --lastPoint 12959 -n .ctzctwObs_noJEC.POINTS.12920.12959
fi
if [ $1 -eq 324 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 12960 --lastPoint 12999 -n .ctzctwObs_noJEC.POINTS.12960.12999
fi
if [ $1 -eq 325 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13000 --lastPoint 13039 -n .ctzctwObs_noJEC.POINTS.13000.13039
fi
if [ $1 -eq 326 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13040 --lastPoint 13079 -n .ctzctwObs_noJEC.POINTS.13040.13079
fi
if [ $1 -eq 327 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13080 --lastPoint 13119 -n .ctzctwObs_noJEC.POINTS.13080.13119
fi
if [ $1 -eq 328 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13120 --lastPoint 13159 -n .ctzctwObs_noJEC.POINTS.13120.13159
fi
if [ $1 -eq 329 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13160 --lastPoint 13199 -n .ctzctwObs_noJEC.POINTS.13160.13199
fi
if [ $1 -eq 330 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13200 --lastPoint 13239 -n .ctzctwObs_noJEC.POINTS.13200.13239
fi
if [ $1 -eq 331 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13240 --lastPoint 13279 -n .ctzctwObs_noJEC.POINTS.13240.13279
fi
if [ $1 -eq 332 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13280 --lastPoint 13319 -n .ctzctwObs_noJEC.POINTS.13280.13319
fi
if [ $1 -eq 333 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13320 --lastPoint 13359 -n .ctzctwObs_noJEC.POINTS.13320.13359
fi
if [ $1 -eq 334 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13360 --lastPoint 13399 -n .ctzctwObs_noJEC.POINTS.13360.13399
fi
if [ $1 -eq 335 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13400 --lastPoint 13439 -n .ctzctwObs_noJEC.POINTS.13400.13439
fi
if [ $1 -eq 336 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13440 --lastPoint 13479 -n .ctzctwObs_noJEC.POINTS.13440.13479
fi
if [ $1 -eq 337 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13480 --lastPoint 13519 -n .ctzctwObs_noJEC.POINTS.13480.13519
fi
if [ $1 -eq 338 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13520 --lastPoint 13559 -n .ctzctwObs_noJEC.POINTS.13520.13559
fi
if [ $1 -eq 339 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13560 --lastPoint 13599 -n .ctzctwObs_noJEC.POINTS.13560.13599
fi
if [ $1 -eq 340 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13600 --lastPoint 13639 -n .ctzctwObs_noJEC.POINTS.13600.13639
fi
if [ $1 -eq 341 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13640 --lastPoint 13679 -n .ctzctwObs_noJEC.POINTS.13640.13679
fi
if [ $1 -eq 342 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13680 --lastPoint 13719 -n .ctzctwObs_noJEC.POINTS.13680.13719
fi
if [ $1 -eq 343 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13720 --lastPoint 13759 -n .ctzctwObs_noJEC.POINTS.13720.13759
fi
if [ $1 -eq 344 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13760 --lastPoint 13799 -n .ctzctwObs_noJEC.POINTS.13760.13799
fi
if [ $1 -eq 345 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13800 --lastPoint 13839 -n .ctzctwObs_noJEC.POINTS.13800.13839
fi
if [ $1 -eq 346 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13840 --lastPoint 13879 -n .ctzctwObs_noJEC.POINTS.13840.13879
fi
if [ $1 -eq 347 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13880 --lastPoint 13919 -n .ctzctwObs_noJEC.POINTS.13880.13919
fi
if [ $1 -eq 348 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13920 --lastPoint 13959 -n .ctzctwObs_noJEC.POINTS.13920.13959
fi
if [ $1 -eq 349 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 13960 --lastPoint 13999 -n .ctzctwObs_noJEC.POINTS.13960.13999
fi
if [ $1 -eq 350 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14000 --lastPoint 14039 -n .ctzctwObs_noJEC.POINTS.14000.14039
fi
if [ $1 -eq 351 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14040 --lastPoint 14079 -n .ctzctwObs_noJEC.POINTS.14040.14079
fi
if [ $1 -eq 352 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14080 --lastPoint 14119 -n .ctzctwObs_noJEC.POINTS.14080.14119
fi
if [ $1 -eq 353 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14120 --lastPoint 14159 -n .ctzctwObs_noJEC.POINTS.14120.14159
fi
if [ $1 -eq 354 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14160 --lastPoint 14199 -n .ctzctwObs_noJEC.POINTS.14160.14199
fi
if [ $1 -eq 355 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14200 --lastPoint 14239 -n .ctzctwObs_noJEC.POINTS.14200.14239
fi
if [ $1 -eq 356 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14240 --lastPoint 14279 -n .ctzctwObs_noJEC.POINTS.14240.14279
fi
if [ $1 -eq 357 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14280 --lastPoint 14319 -n .ctzctwObs_noJEC.POINTS.14280.14319
fi
if [ $1 -eq 358 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14320 --lastPoint 14359 -n .ctzctwObs_noJEC.POINTS.14320.14359
fi
if [ $1 -eq 359 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14360 --lastPoint 14399 -n .ctzctwObs_noJEC.POINTS.14360.14399
fi
if [ $1 -eq 360 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14400 --lastPoint 14439 -n .ctzctwObs_noJEC.POINTS.14400.14439
fi
if [ $1 -eq 361 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14440 --lastPoint 14479 -n .ctzctwObs_noJEC.POINTS.14440.14479
fi
if [ $1 -eq 362 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14480 --lastPoint 14519 -n .ctzctwObs_noJEC.POINTS.14480.14519
fi
if [ $1 -eq 363 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14520 --lastPoint 14559 -n .ctzctwObs_noJEC.POINTS.14520.14559
fi
if [ $1 -eq 364 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14560 --lastPoint 14599 -n .ctzctwObs_noJEC.POINTS.14560.14599
fi
if [ $1 -eq 365 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14600 --lastPoint 14639 -n .ctzctwObs_noJEC.POINTS.14600.14639
fi
if [ $1 -eq 366 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14640 --lastPoint 14679 -n .ctzctwObs_noJEC.POINTS.14640.14679
fi
if [ $1 -eq 367 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14680 --lastPoint 14719 -n .ctzctwObs_noJEC.POINTS.14680.14719
fi
if [ $1 -eq 368 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14720 --lastPoint 14759 -n .ctzctwObs_noJEC.POINTS.14720.14759
fi
if [ $1 -eq 369 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14760 --lastPoint 14799 -n .ctzctwObs_noJEC.POINTS.14760.14799
fi
if [ $1 -eq 370 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14800 --lastPoint 14839 -n .ctzctwObs_noJEC.POINTS.14800.14839
fi
if [ $1 -eq 371 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14840 --lastPoint 14879 -n .ctzctwObs_noJEC.POINTS.14840.14879
fi
if [ $1 -eq 372 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14880 --lastPoint 14919 -n .ctzctwObs_noJEC.POINTS.14880.14919
fi
if [ $1 -eq 373 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14920 --lastPoint 14959 -n .ctzctwObs_noJEC.POINTS.14920.14959
fi
if [ $1 -eq 374 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 14960 --lastPoint 14999 -n .ctzctwObs_noJEC.POINTS.14960.14999
fi
if [ $1 -eq 375 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15000 --lastPoint 15039 -n .ctzctwObs_noJEC.POINTS.15000.15039
fi
if [ $1 -eq 376 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15040 --lastPoint 15079 -n .ctzctwObs_noJEC.POINTS.15040.15079
fi
if [ $1 -eq 377 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15080 --lastPoint 15119 -n .ctzctwObs_noJEC.POINTS.15080.15119
fi
if [ $1 -eq 378 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15120 --lastPoint 15159 -n .ctzctwObs_noJEC.POINTS.15120.15159
fi
if [ $1 -eq 379 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15160 --lastPoint 15199 -n .ctzctwObs_noJEC.POINTS.15160.15199
fi
if [ $1 -eq 380 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15200 --lastPoint 15239 -n .ctzctwObs_noJEC.POINTS.15200.15239
fi
if [ $1 -eq 381 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15240 --lastPoint 15279 -n .ctzctwObs_noJEC.POINTS.15240.15279
fi
if [ $1 -eq 382 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15280 --lastPoint 15319 -n .ctzctwObs_noJEC.POINTS.15280.15319
fi
if [ $1 -eq 383 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15320 --lastPoint 15359 -n .ctzctwObs_noJEC.POINTS.15320.15359
fi
if [ $1 -eq 384 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15360 --lastPoint 15399 -n .ctzctwObs_noJEC.POINTS.15360.15399
fi
if [ $1 -eq 385 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15400 --lastPoint 15439 -n .ctzctwObs_noJEC.POINTS.15400.15439
fi
if [ $1 -eq 386 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15440 --lastPoint 15479 -n .ctzctwObs_noJEC.POINTS.15440.15479
fi
if [ $1 -eq 387 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15480 --lastPoint 15519 -n .ctzctwObs_noJEC.POINTS.15480.15519
fi
if [ $1 -eq 388 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15520 --lastPoint 15559 -n .ctzctwObs_noJEC.POINTS.15520.15559
fi
if [ $1 -eq 389 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15560 --lastPoint 15599 -n .ctzctwObs_noJEC.POINTS.15560.15599
fi
if [ $1 -eq 390 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15600 --lastPoint 15639 -n .ctzctwObs_noJEC.POINTS.15600.15639
fi
if [ $1 -eq 391 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15640 --lastPoint 15679 -n .ctzctwObs_noJEC.POINTS.15640.15679
fi
if [ $1 -eq 392 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15680 --lastPoint 15719 -n .ctzctwObs_noJEC.POINTS.15680.15719
fi
if [ $1 -eq 393 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15720 --lastPoint 15759 -n .ctzctwObs_noJEC.POINTS.15720.15759
fi
if [ $1 -eq 394 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15760 --lastPoint 15799 -n .ctzctwObs_noJEC.POINTS.15760.15799
fi
if [ $1 -eq 395 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15800 --lastPoint 15839 -n .ctzctwObs_noJEC.POINTS.15800.15839
fi
if [ $1 -eq 396 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15840 --lastPoint 15879 -n .ctzctwObs_noJEC.POINTS.15840.15879
fi
if [ $1 -eq 397 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15880 --lastPoint 15919 -n .ctzctwObs_noJEC.POINTS.15880.15919
fi
if [ $1 -eq 398 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15920 --lastPoint 15959 -n .ctzctwObs_noJEC.POINTS.15920.15959
fi
if [ $1 -eq 399 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 15960 --lastPoint 15999 -n .ctzctwObs_noJEC.POINTS.15960.15999
fi
if [ $1 -eq 400 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16000 --lastPoint 16039 -n .ctzctwObs_noJEC.POINTS.16000.16039
fi
if [ $1 -eq 401 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16040 --lastPoint 16079 -n .ctzctwObs_noJEC.POINTS.16040.16079
fi
if [ $1 -eq 402 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16080 --lastPoint 16119 -n .ctzctwObs_noJEC.POINTS.16080.16119
fi
if [ $1 -eq 403 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16120 --lastPoint 16159 -n .ctzctwObs_noJEC.POINTS.16120.16159
fi
if [ $1 -eq 404 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16160 --lastPoint 16199 -n .ctzctwObs_noJEC.POINTS.16160.16199
fi
if [ $1 -eq 405 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16200 --lastPoint 16239 -n .ctzctwObs_noJEC.POINTS.16200.16239
fi
if [ $1 -eq 406 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16240 --lastPoint 16279 -n .ctzctwObs_noJEC.POINTS.16240.16279
fi
if [ $1 -eq 407 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16280 --lastPoint 16319 -n .ctzctwObs_noJEC.POINTS.16280.16319
fi
if [ $1 -eq 408 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16320 --lastPoint 16359 -n .ctzctwObs_noJEC.POINTS.16320.16359
fi
if [ $1 -eq 409 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16360 --lastPoint 16399 -n .ctzctwObs_noJEC.POINTS.16360.16399
fi
if [ $1 -eq 410 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16400 --lastPoint 16439 -n .ctzctwObs_noJEC.POINTS.16400.16439
fi
if [ $1 -eq 411 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16440 --lastPoint 16479 -n .ctzctwObs_noJEC.POINTS.16440.16479
fi
if [ $1 -eq 412 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16480 --lastPoint 16519 -n .ctzctwObs_noJEC.POINTS.16480.16519
fi
if [ $1 -eq 413 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16520 --lastPoint 16559 -n .ctzctwObs_noJEC.POINTS.16520.16559
fi
if [ $1 -eq 414 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16560 --lastPoint 16599 -n .ctzctwObs_noJEC.POINTS.16560.16599
fi
if [ $1 -eq 415 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16600 --lastPoint 16639 -n .ctzctwObs_noJEC.POINTS.16600.16639
fi
if [ $1 -eq 416 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16640 --lastPoint 16679 -n .ctzctwObs_noJEC.POINTS.16640.16679
fi
if [ $1 -eq 417 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16680 --lastPoint 16719 -n .ctzctwObs_noJEC.POINTS.16680.16719
fi
if [ $1 -eq 418 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16720 --lastPoint 16759 -n .ctzctwObs_noJEC.POINTS.16720.16759
fi
if [ $1 -eq 419 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16760 --lastPoint 16799 -n .ctzctwObs_noJEC.POINTS.16760.16799
fi
if [ $1 -eq 420 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16800 --lastPoint 16839 -n .ctzctwObs_noJEC.POINTS.16800.16839
fi
if [ $1 -eq 421 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16840 --lastPoint 16879 -n .ctzctwObs_noJEC.POINTS.16840.16879
fi
if [ $1 -eq 422 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16880 --lastPoint 16919 -n .ctzctwObs_noJEC.POINTS.16880.16919
fi
if [ $1 -eq 423 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16920 --lastPoint 16959 -n .ctzctwObs_noJEC.POINTS.16920.16959
fi
if [ $1 -eq 424 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 16960 --lastPoint 16999 -n .ctzctwObs_noJEC.POINTS.16960.16999
fi
if [ $1 -eq 425 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17000 --lastPoint 17039 -n .ctzctwObs_noJEC.POINTS.17000.17039
fi
if [ $1 -eq 426 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17040 --lastPoint 17079 -n .ctzctwObs_noJEC.POINTS.17040.17079
fi
if [ $1 -eq 427 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17080 --lastPoint 17119 -n .ctzctwObs_noJEC.POINTS.17080.17119
fi
if [ $1 -eq 428 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17120 --lastPoint 17159 -n .ctzctwObs_noJEC.POINTS.17120.17159
fi
if [ $1 -eq 429 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17160 --lastPoint 17199 -n .ctzctwObs_noJEC.POINTS.17160.17199
fi
if [ $1 -eq 430 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17200 --lastPoint 17239 -n .ctzctwObs_noJEC.POINTS.17200.17239
fi
if [ $1 -eq 431 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17240 --lastPoint 17279 -n .ctzctwObs_noJEC.POINTS.17240.17279
fi
if [ $1 -eq 432 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17280 --lastPoint 17319 -n .ctzctwObs_noJEC.POINTS.17280.17319
fi
if [ $1 -eq 433 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17320 --lastPoint 17359 -n .ctzctwObs_noJEC.POINTS.17320.17359
fi
if [ $1 -eq 434 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17360 --lastPoint 17399 -n .ctzctwObs_noJEC.POINTS.17360.17399
fi
if [ $1 -eq 435 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17400 --lastPoint 17439 -n .ctzctwObs_noJEC.POINTS.17400.17439
fi
if [ $1 -eq 436 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17440 --lastPoint 17479 -n .ctzctwObs_noJEC.POINTS.17440.17479
fi
if [ $1 -eq 437 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17480 --lastPoint 17519 -n .ctzctwObs_noJEC.POINTS.17480.17519
fi
if [ $1 -eq 438 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17520 --lastPoint 17559 -n .ctzctwObs_noJEC.POINTS.17520.17559
fi
if [ $1 -eq 439 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17560 --lastPoint 17599 -n .ctzctwObs_noJEC.POINTS.17560.17599
fi
if [ $1 -eq 440 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17600 --lastPoint 17639 -n .ctzctwObs_noJEC.POINTS.17600.17639
fi
if [ $1 -eq 441 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17640 --lastPoint 17679 -n .ctzctwObs_noJEC.POINTS.17640.17679
fi
if [ $1 -eq 442 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17680 --lastPoint 17719 -n .ctzctwObs_noJEC.POINTS.17680.17719
fi
if [ $1 -eq 443 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17720 --lastPoint 17759 -n .ctzctwObs_noJEC.POINTS.17720.17759
fi
if [ $1 -eq 444 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17760 --lastPoint 17799 -n .ctzctwObs_noJEC.POINTS.17760.17799
fi
if [ $1 -eq 445 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17800 --lastPoint 17839 -n .ctzctwObs_noJEC.POINTS.17800.17839
fi
if [ $1 -eq 446 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17840 --lastPoint 17879 -n .ctzctwObs_noJEC.POINTS.17840.17879
fi
if [ $1 -eq 447 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17880 --lastPoint 17919 -n .ctzctwObs_noJEC.POINTS.17880.17919
fi
if [ $1 -eq 448 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17920 --lastPoint 17959 -n .ctzctwObs_noJEC.POINTS.17920.17959
fi
if [ $1 -eq 449 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 17960 --lastPoint 17999 -n .ctzctwObs_noJEC.POINTS.17960.17999
fi
if [ $1 -eq 450 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18000 --lastPoint 18039 -n .ctzctwObs_noJEC.POINTS.18000.18039
fi
if [ $1 -eq 451 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18040 --lastPoint 18079 -n .ctzctwObs_noJEC.POINTS.18040.18079
fi
if [ $1 -eq 452 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18080 --lastPoint 18119 -n .ctzctwObs_noJEC.POINTS.18080.18119
fi
if [ $1 -eq 453 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18120 --lastPoint 18159 -n .ctzctwObs_noJEC.POINTS.18120.18159
fi
if [ $1 -eq 454 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18160 --lastPoint 18199 -n .ctzctwObs_noJEC.POINTS.18160.18199
fi
if [ $1 -eq 455 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18200 --lastPoint 18239 -n .ctzctwObs_noJEC.POINTS.18200.18239
fi
if [ $1 -eq 456 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18240 --lastPoint 18279 -n .ctzctwObs_noJEC.POINTS.18240.18279
fi
if [ $1 -eq 457 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18280 --lastPoint 18319 -n .ctzctwObs_noJEC.POINTS.18280.18319
fi
if [ $1 -eq 458 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18320 --lastPoint 18359 -n .ctzctwObs_noJEC.POINTS.18320.18359
fi
if [ $1 -eq 459 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18360 --lastPoint 18399 -n .ctzctwObs_noJEC.POINTS.18360.18399
fi
if [ $1 -eq 460 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18400 --lastPoint 18439 -n .ctzctwObs_noJEC.POINTS.18400.18439
fi
if [ $1 -eq 461 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18440 --lastPoint 18479 -n .ctzctwObs_noJEC.POINTS.18440.18479
fi
if [ $1 -eq 462 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18480 --lastPoint 18519 -n .ctzctwObs_noJEC.POINTS.18480.18519
fi
if [ $1 -eq 463 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18520 --lastPoint 18559 -n .ctzctwObs_noJEC.POINTS.18520.18559
fi
if [ $1 -eq 464 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18560 --lastPoint 18599 -n .ctzctwObs_noJEC.POINTS.18560.18599
fi
if [ $1 -eq 465 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18600 --lastPoint 18639 -n .ctzctwObs_noJEC.POINTS.18600.18639
fi
if [ $1 -eq 466 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18640 --lastPoint 18679 -n .ctzctwObs_noJEC.POINTS.18640.18679
fi
if [ $1 -eq 467 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18680 --lastPoint 18719 -n .ctzctwObs_noJEC.POINTS.18680.18719
fi
if [ $1 -eq 468 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18720 --lastPoint 18759 -n .ctzctwObs_noJEC.POINTS.18720.18759
fi
if [ $1 -eq 469 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18760 --lastPoint 18799 -n .ctzctwObs_noJEC.POINTS.18760.18799
fi
if [ $1 -eq 470 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18800 --lastPoint 18839 -n .ctzctwObs_noJEC.POINTS.18800.18839
fi
if [ $1 -eq 471 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18840 --lastPoint 18879 -n .ctzctwObs_noJEC.POINTS.18840.18879
fi
if [ $1 -eq 472 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18880 --lastPoint 18919 -n .ctzctwObs_noJEC.POINTS.18880.18919
fi
if [ $1 -eq 473 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18920 --lastPoint 18959 -n .ctzctwObs_noJEC.POINTS.18920.18959
fi
if [ $1 -eq 474 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 18960 --lastPoint 18999 -n .ctzctwObs_noJEC.POINTS.18960.18999
fi
if [ $1 -eq 475 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19000 --lastPoint 19039 -n .ctzctwObs_noJEC.POINTS.19000.19039
fi
if [ $1 -eq 476 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19040 --lastPoint 19079 -n .ctzctwObs_noJEC.POINTS.19040.19079
fi
if [ $1 -eq 477 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19080 --lastPoint 19119 -n .ctzctwObs_noJEC.POINTS.19080.19119
fi
if [ $1 -eq 478 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19120 --lastPoint 19159 -n .ctzctwObs_noJEC.POINTS.19120.19159
fi
if [ $1 -eq 479 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19160 --lastPoint 19199 -n .ctzctwObs_noJEC.POINTS.19160.19199
fi
if [ $1 -eq 480 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19200 --lastPoint 19239 -n .ctzctwObs_noJEC.POINTS.19200.19239
fi
if [ $1 -eq 481 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19240 --lastPoint 19279 -n .ctzctwObs_noJEC.POINTS.19240.19279
fi
if [ $1 -eq 482 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19280 --lastPoint 19319 -n .ctzctwObs_noJEC.POINTS.19280.19319
fi
if [ $1 -eq 483 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19320 --lastPoint 19359 -n .ctzctwObs_noJEC.POINTS.19320.19359
fi
if [ $1 -eq 484 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19360 --lastPoint 19399 -n .ctzctwObs_noJEC.POINTS.19360.19399
fi
if [ $1 -eq 485 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19400 --lastPoint 19439 -n .ctzctwObs_noJEC.POINTS.19400.19439
fi
if [ $1 -eq 486 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19440 --lastPoint 19479 -n .ctzctwObs_noJEC.POINTS.19440.19479
fi
if [ $1 -eq 487 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19480 --lastPoint 19519 -n .ctzctwObs_noJEC.POINTS.19480.19519
fi
if [ $1 -eq 488 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19520 --lastPoint 19559 -n .ctzctwObs_noJEC.POINTS.19520.19559
fi
if [ $1 -eq 489 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19560 --lastPoint 19599 -n .ctzctwObs_noJEC.POINTS.19560.19599
fi
if [ $1 -eq 490 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19600 --lastPoint 19639 -n .ctzctwObs_noJEC.POINTS.19600.19639
fi
if [ $1 -eq 491 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19640 --lastPoint 19679 -n .ctzctwObs_noJEC.POINTS.19640.19679
fi
if [ $1 -eq 492 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19680 --lastPoint 19719 -n .ctzctwObs_noJEC.POINTS.19680.19719
fi
if [ $1 -eq 493 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19720 --lastPoint 19759 -n .ctzctwObs_noJEC.POINTS.19720.19759
fi
if [ $1 -eq 494 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19760 --lastPoint 19799 -n .ctzctwObs_noJEC.POINTS.19760.19799
fi
if [ $1 -eq 495 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19800 --lastPoint 19839 -n .ctzctwObs_noJEC.POINTS.19800.19839
fi
if [ $1 -eq 496 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19840 --lastPoint 19879 -n .ctzctwObs_noJEC.POINTS.19840.19879
fi
if [ $1 -eq 497 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19880 --lastPoint 19919 -n .ctzctwObs_noJEC.POINTS.19880.19919
fi
if [ $1 -eq 498 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19920 --lastPoint 19959 -n .ctzctwObs_noJEC.POINTS.19920.19959
fi
if [ $1 -eq 499 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 19960 --lastPoint 19999 -n .ctzctwObs_noJEC.POINTS.19960.19999
fi
if [ $1 -eq 500 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20000 --lastPoint 20039 -n .ctzctwObs_noJEC.POINTS.20000.20039
fi
if [ $1 -eq 501 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20040 --lastPoint 20079 -n .ctzctwObs_noJEC.POINTS.20040.20079
fi
if [ $1 -eq 502 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20080 --lastPoint 20119 -n .ctzctwObs_noJEC.POINTS.20080.20119
fi
if [ $1 -eq 503 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20120 --lastPoint 20159 -n .ctzctwObs_noJEC.POINTS.20120.20159
fi
if [ $1 -eq 504 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20160 --lastPoint 20199 -n .ctzctwObs_noJEC.POINTS.20160.20199
fi
if [ $1 -eq 505 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20200 --lastPoint 20239 -n .ctzctwObs_noJEC.POINTS.20200.20239
fi
if [ $1 -eq 506 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20240 --lastPoint 20279 -n .ctzctwObs_noJEC.POINTS.20240.20279
fi
if [ $1 -eq 507 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20280 --lastPoint 20319 -n .ctzctwObs_noJEC.POINTS.20280.20319
fi
if [ $1 -eq 508 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20320 --lastPoint 20359 -n .ctzctwObs_noJEC.POINTS.20320.20359
fi
if [ $1 -eq 509 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20360 --lastPoint 20399 -n .ctzctwObs_noJEC.POINTS.20360.20399
fi
if [ $1 -eq 510 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20400 --lastPoint 20439 -n .ctzctwObs_noJEC.POINTS.20400.20439
fi
if [ $1 -eq 511 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20440 --lastPoint 20479 -n .ctzctwObs_noJEC.POINTS.20440.20479
fi
if [ $1 -eq 512 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20480 --lastPoint 20519 -n .ctzctwObs_noJEC.POINTS.20480.20519
fi
if [ $1 -eq 513 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20520 --lastPoint 20559 -n .ctzctwObs_noJEC.POINTS.20520.20559
fi
if [ $1 -eq 514 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20560 --lastPoint 20599 -n .ctzctwObs_noJEC.POINTS.20560.20599
fi
if [ $1 -eq 515 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20600 --lastPoint 20639 -n .ctzctwObs_noJEC.POINTS.20600.20639
fi
if [ $1 -eq 516 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20640 --lastPoint 20679 -n .ctzctwObs_noJEC.POINTS.20640.20679
fi
if [ $1 -eq 517 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20680 --lastPoint 20719 -n .ctzctwObs_noJEC.POINTS.20680.20719
fi
if [ $1 -eq 518 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20720 --lastPoint 20759 -n .ctzctwObs_noJEC.POINTS.20720.20759
fi
if [ $1 -eq 519 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20760 --lastPoint 20799 -n .ctzctwObs_noJEC.POINTS.20760.20799
fi
if [ $1 -eq 520 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20800 --lastPoint 20839 -n .ctzctwObs_noJEC.POINTS.20800.20839
fi
if [ $1 -eq 521 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20840 --lastPoint 20879 -n .ctzctwObs_noJEC.POINTS.20840.20879
fi
if [ $1 -eq 522 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20880 --lastPoint 20919 -n .ctzctwObs_noJEC.POINTS.20880.20919
fi
if [ $1 -eq 523 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20920 --lastPoint 20959 -n .ctzctwObs_noJEC.POINTS.20920.20959
fi
if [ $1 -eq 524 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 20960 --lastPoint 20999 -n .ctzctwObs_noJEC.POINTS.20960.20999
fi
if [ $1 -eq 525 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21000 --lastPoint 21039 -n .ctzctwObs_noJEC.POINTS.21000.21039
fi
if [ $1 -eq 526 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21040 --lastPoint 21079 -n .ctzctwObs_noJEC.POINTS.21040.21079
fi
if [ $1 -eq 527 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21080 --lastPoint 21119 -n .ctzctwObs_noJEC.POINTS.21080.21119
fi
if [ $1 -eq 528 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21120 --lastPoint 21159 -n .ctzctwObs_noJEC.POINTS.21120.21159
fi
if [ $1 -eq 529 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21160 --lastPoint 21199 -n .ctzctwObs_noJEC.POINTS.21160.21199
fi
if [ $1 -eq 530 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21200 --lastPoint 21239 -n .ctzctwObs_noJEC.POINTS.21200.21239
fi
if [ $1 -eq 531 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21240 --lastPoint 21279 -n .ctzctwObs_noJEC.POINTS.21240.21279
fi
if [ $1 -eq 532 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21280 --lastPoint 21319 -n .ctzctwObs_noJEC.POINTS.21280.21319
fi
if [ $1 -eq 533 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21320 --lastPoint 21359 -n .ctzctwObs_noJEC.POINTS.21320.21359
fi
if [ $1 -eq 534 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21360 --lastPoint 21399 -n .ctzctwObs_noJEC.POINTS.21360.21399
fi
if [ $1 -eq 535 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21400 --lastPoint 21439 -n .ctzctwObs_noJEC.POINTS.21400.21439
fi
if [ $1 -eq 536 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21440 --lastPoint 21479 -n .ctzctwObs_noJEC.POINTS.21440.21479
fi
if [ $1 -eq 537 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21480 --lastPoint 21519 -n .ctzctwObs_noJEC.POINTS.21480.21519
fi
if [ $1 -eq 538 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21520 --lastPoint 21559 -n .ctzctwObs_noJEC.POINTS.21520.21559
fi
if [ $1 -eq 539 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21560 --lastPoint 21599 -n .ctzctwObs_noJEC.POINTS.21560.21599
fi
if [ $1 -eq 540 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21600 --lastPoint 21639 -n .ctzctwObs_noJEC.POINTS.21600.21639
fi
if [ $1 -eq 541 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21640 --lastPoint 21679 -n .ctzctwObs_noJEC.POINTS.21640.21679
fi
if [ $1 -eq 542 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21680 --lastPoint 21719 -n .ctzctwObs_noJEC.POINTS.21680.21719
fi
if [ $1 -eq 543 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21720 --lastPoint 21759 -n .ctzctwObs_noJEC.POINTS.21720.21759
fi
if [ $1 -eq 544 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21760 --lastPoint 21799 -n .ctzctwObs_noJEC.POINTS.21760.21799
fi
if [ $1 -eq 545 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21800 --lastPoint 21839 -n .ctzctwObs_noJEC.POINTS.21800.21839
fi
if [ $1 -eq 546 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21840 --lastPoint 21879 -n .ctzctwObs_noJEC.POINTS.21840.21879
fi
if [ $1 -eq 547 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21880 --lastPoint 21919 -n .ctzctwObs_noJEC.POINTS.21880.21919
fi
if [ $1 -eq 548 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21920 --lastPoint 21959 -n .ctzctwObs_noJEC.POINTS.21920.21959
fi
if [ $1 -eq 549 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 21960 --lastPoint 21999 -n .ctzctwObs_noJEC.POINTS.21960.21999
fi
if [ $1 -eq 550 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22000 --lastPoint 22039 -n .ctzctwObs_noJEC.POINTS.22000.22039
fi
if [ $1 -eq 551 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22040 --lastPoint 22079 -n .ctzctwObs_noJEC.POINTS.22040.22079
fi
if [ $1 -eq 552 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22080 --lastPoint 22119 -n .ctzctwObs_noJEC.POINTS.22080.22119
fi
if [ $1 -eq 553 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22120 --lastPoint 22159 -n .ctzctwObs_noJEC.POINTS.22120.22159
fi
if [ $1 -eq 554 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22160 --lastPoint 22199 -n .ctzctwObs_noJEC.POINTS.22160.22199
fi
if [ $1 -eq 555 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22200 --lastPoint 22239 -n .ctzctwObs_noJEC.POINTS.22200.22239
fi
if [ $1 -eq 556 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22240 --lastPoint 22279 -n .ctzctwObs_noJEC.POINTS.22240.22279
fi
if [ $1 -eq 557 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22280 --lastPoint 22319 -n .ctzctwObs_noJEC.POINTS.22280.22319
fi
if [ $1 -eq 558 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22320 --lastPoint 22359 -n .ctzctwObs_noJEC.POINTS.22320.22359
fi
if [ $1 -eq 559 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22360 --lastPoint 22399 -n .ctzctwObs_noJEC.POINTS.22360.22399
fi
if [ $1 -eq 560 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22400 --lastPoint 22439 -n .ctzctwObs_noJEC.POINTS.22400.22439
fi
if [ $1 -eq 561 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22440 --lastPoint 22479 -n .ctzctwObs_noJEC.POINTS.22440.22479
fi
if [ $1 -eq 562 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22480 --lastPoint 22519 -n .ctzctwObs_noJEC.POINTS.22480.22519
fi
if [ $1 -eq 563 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22520 --lastPoint 22559 -n .ctzctwObs_noJEC.POINTS.22520.22559
fi
if [ $1 -eq 564 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22560 --lastPoint 22599 -n .ctzctwObs_noJEC.POINTS.22560.22599
fi
if [ $1 -eq 565 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22600 --lastPoint 22639 -n .ctzctwObs_noJEC.POINTS.22600.22639
fi
if [ $1 -eq 566 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22640 --lastPoint 22679 -n .ctzctwObs_noJEC.POINTS.22640.22679
fi
if [ $1 -eq 567 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22680 --lastPoint 22719 -n .ctzctwObs_noJEC.POINTS.22680.22719
fi
if [ $1 -eq 568 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22720 --lastPoint 22759 -n .ctzctwObs_noJEC.POINTS.22720.22759
fi
if [ $1 -eq 569 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22760 --lastPoint 22799 -n .ctzctwObs_noJEC.POINTS.22760.22799
fi
if [ $1 -eq 570 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22800 --lastPoint 22839 -n .ctzctwObs_noJEC.POINTS.22800.22839
fi
if [ $1 -eq 571 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22840 --lastPoint 22879 -n .ctzctwObs_noJEC.POINTS.22840.22879
fi
if [ $1 -eq 572 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22880 --lastPoint 22919 -n .ctzctwObs_noJEC.POINTS.22880.22919
fi
if [ $1 -eq 573 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22920 --lastPoint 22959 -n .ctzctwObs_noJEC.POINTS.22920.22959
fi
if [ $1 -eq 574 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 22960 --lastPoint 22999 -n .ctzctwObs_noJEC.POINTS.22960.22999
fi
if [ $1 -eq 575 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23000 --lastPoint 23039 -n .ctzctwObs_noJEC.POINTS.23000.23039
fi
if [ $1 -eq 576 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23040 --lastPoint 23079 -n .ctzctwObs_noJEC.POINTS.23040.23079
fi
if [ $1 -eq 577 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23080 --lastPoint 23119 -n .ctzctwObs_noJEC.POINTS.23080.23119
fi
if [ $1 -eq 578 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23120 --lastPoint 23159 -n .ctzctwObs_noJEC.POINTS.23120.23159
fi
if [ $1 -eq 579 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23160 --lastPoint 23199 -n .ctzctwObs_noJEC.POINTS.23160.23199
fi
if [ $1 -eq 580 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23200 --lastPoint 23239 -n .ctzctwObs_noJEC.POINTS.23200.23239
fi
if [ $1 -eq 581 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23240 --lastPoint 23279 -n .ctzctwObs_noJEC.POINTS.23240.23279
fi
if [ $1 -eq 582 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23280 --lastPoint 23319 -n .ctzctwObs_noJEC.POINTS.23280.23319
fi
if [ $1 -eq 583 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23320 --lastPoint 23359 -n .ctzctwObs_noJEC.POINTS.23320.23359
fi
if [ $1 -eq 584 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23360 --lastPoint 23399 -n .ctzctwObs_noJEC.POINTS.23360.23399
fi
if [ $1 -eq 585 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23400 --lastPoint 23439 -n .ctzctwObs_noJEC.POINTS.23400.23439
fi
if [ $1 -eq 586 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23440 --lastPoint 23479 -n .ctzctwObs_noJEC.POINTS.23440.23479
fi
if [ $1 -eq 587 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23480 --lastPoint 23519 -n .ctzctwObs_noJEC.POINTS.23480.23519
fi
if [ $1 -eq 588 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23520 --lastPoint 23559 -n .ctzctwObs_noJEC.POINTS.23520.23559
fi
if [ $1 -eq 589 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23560 --lastPoint 23599 -n .ctzctwObs_noJEC.POINTS.23560.23599
fi
if [ $1 -eq 590 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23600 --lastPoint 23639 -n .ctzctwObs_noJEC.POINTS.23600.23639
fi
if [ $1 -eq 591 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23640 --lastPoint 23679 -n .ctzctwObs_noJEC.POINTS.23640.23679
fi
if [ $1 -eq 592 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23680 --lastPoint 23719 -n .ctzctwObs_noJEC.POINTS.23680.23719
fi
if [ $1 -eq 593 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23720 --lastPoint 23759 -n .ctzctwObs_noJEC.POINTS.23720.23759
fi
if [ $1 -eq 594 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23760 --lastPoint 23799 -n .ctzctwObs_noJEC.POINTS.23760.23799
fi
if [ $1 -eq 595 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23800 --lastPoint 23839 -n .ctzctwObs_noJEC.POINTS.23800.23839
fi
if [ $1 -eq 596 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23840 --lastPoint 23879 -n .ctzctwObs_noJEC.POINTS.23840.23879
fi
if [ $1 -eq 597 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23880 --lastPoint 23919 -n .ctzctwObs_noJEC.POINTS.23880.23919
fi
if [ $1 -eq 598 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23920 --lastPoint 23959 -n .ctzctwObs_noJEC.POINTS.23920.23959
fi
if [ $1 -eq 599 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 23960 --lastPoint 23999 -n .ctzctwObs_noJEC.POINTS.23960.23999
fi
if [ $1 -eq 600 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24000 --lastPoint 24039 -n .ctzctwObs_noJEC.POINTS.24000.24039
fi
if [ $1 -eq 601 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24040 --lastPoint 24079 -n .ctzctwObs_noJEC.POINTS.24040.24079
fi
if [ $1 -eq 602 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24080 --lastPoint 24119 -n .ctzctwObs_noJEC.POINTS.24080.24119
fi
if [ $1 -eq 603 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24120 --lastPoint 24159 -n .ctzctwObs_noJEC.POINTS.24120.24159
fi
if [ $1 -eq 604 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24160 --lastPoint 24199 -n .ctzctwObs_noJEC.POINTS.24160.24199
fi
if [ $1 -eq 605 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24200 --lastPoint 24239 -n .ctzctwObs_noJEC.POINTS.24200.24239
fi
if [ $1 -eq 606 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24240 --lastPoint 24279 -n .ctzctwObs_noJEC.POINTS.24240.24279
fi
if [ $1 -eq 607 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24280 --lastPoint 24319 -n .ctzctwObs_noJEC.POINTS.24280.24319
fi
if [ $1 -eq 608 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24320 --lastPoint 24359 -n .ctzctwObs_noJEC.POINTS.24320.24359
fi
if [ $1 -eq 609 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24360 --lastPoint 24399 -n .ctzctwObs_noJEC.POINTS.24360.24399
fi
if [ $1 -eq 610 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24400 --lastPoint 24439 -n .ctzctwObs_noJEC.POINTS.24400.24439
fi
if [ $1 -eq 611 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24440 --lastPoint 24479 -n .ctzctwObs_noJEC.POINTS.24440.24479
fi
if [ $1 -eq 612 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24480 --lastPoint 24519 -n .ctzctwObs_noJEC.POINTS.24480.24519
fi
if [ $1 -eq 613 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24520 --lastPoint 24559 -n .ctzctwObs_noJEC.POINTS.24520.24559
fi
if [ $1 -eq 614 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24560 --lastPoint 24599 -n .ctzctwObs_noJEC.POINTS.24560.24599
fi
if [ $1 -eq 615 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24600 --lastPoint 24639 -n .ctzctwObs_noJEC.POINTS.24600.24639
fi
if [ $1 -eq 616 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24640 --lastPoint 24679 -n .ctzctwObs_noJEC.POINTS.24640.24679
fi
if [ $1 -eq 617 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24680 --lastPoint 24719 -n .ctzctwObs_noJEC.POINTS.24680.24719
fi
if [ $1 -eq 618 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24720 --lastPoint 24759 -n .ctzctwObs_noJEC.POINTS.24720.24759
fi
if [ $1 -eq 619 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24760 --lastPoint 24799 -n .ctzctwObs_noJEC.POINTS.24760.24799
fi
if [ $1 -eq 620 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24800 --lastPoint 24839 -n .ctzctwObs_noJEC.POINTS.24800.24839
fi
if [ $1 -eq 621 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24840 --lastPoint 24879 -n .ctzctwObs_noJEC.POINTS.24840.24879
fi
if [ $1 -eq 622 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24880 --lastPoint 24919 -n .ctzctwObs_noJEC.POINTS.24880.24919
fi
if [ $1 -eq 623 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24920 --lastPoint 24959 -n .ctzctwObs_noJEC.POINTS.24920.24959
fi
if [ $1 -eq 624 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 24960 --lastPoint 24999 -n .ctzctwObs_noJEC.POINTS.24960.24999
fi
if [ $1 -eq 625 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25000 --lastPoint 25039 -n .ctzctwObs_noJEC.POINTS.25000.25039
fi
if [ $1 -eq 626 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25040 --lastPoint 25079 -n .ctzctwObs_noJEC.POINTS.25040.25079
fi
if [ $1 -eq 627 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25080 --lastPoint 25119 -n .ctzctwObs_noJEC.POINTS.25080.25119
fi
if [ $1 -eq 628 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25120 --lastPoint 25159 -n .ctzctwObs_noJEC.POINTS.25120.25159
fi
if [ $1 -eq 629 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25160 --lastPoint 25199 -n .ctzctwObs_noJEC.POINTS.25160.25199
fi
if [ $1 -eq 630 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25200 --lastPoint 25239 -n .ctzctwObs_noJEC.POINTS.25200.25239
fi
if [ $1 -eq 631 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25240 --lastPoint 25279 -n .ctzctwObs_noJEC.POINTS.25240.25279
fi
if [ $1 -eq 632 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25280 --lastPoint 25319 -n .ctzctwObs_noJEC.POINTS.25280.25319
fi
if [ $1 -eq 633 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25320 --lastPoint 25359 -n .ctzctwObs_noJEC.POINTS.25320.25359
fi
if [ $1 -eq 634 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25360 --lastPoint 25399 -n .ctzctwObs_noJEC.POINTS.25360.25399
fi
if [ $1 -eq 635 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25400 --lastPoint 25439 -n .ctzctwObs_noJEC.POINTS.25400.25439
fi
if [ $1 -eq 636 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25440 --lastPoint 25479 -n .ctzctwObs_noJEC.POINTS.25440.25479
fi
if [ $1 -eq 637 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25480 --lastPoint 25519 -n .ctzctwObs_noJEC.POINTS.25480.25519
fi
if [ $1 -eq 638 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25520 --lastPoint 25559 -n .ctzctwObs_noJEC.POINTS.25520.25559
fi
if [ $1 -eq 639 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25560 --lastPoint 25599 -n .ctzctwObs_noJEC.POINTS.25560.25599
fi
if [ $1 -eq 640 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25600 --lastPoint 25639 -n .ctzctwObs_noJEC.POINTS.25600.25639
fi
if [ $1 -eq 641 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25640 --lastPoint 25679 -n .ctzctwObs_noJEC.POINTS.25640.25679
fi
if [ $1 -eq 642 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25680 --lastPoint 25719 -n .ctzctwObs_noJEC.POINTS.25680.25719
fi
if [ $1 -eq 643 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25720 --lastPoint 25759 -n .ctzctwObs_noJEC.POINTS.25720.25759
fi
if [ $1 -eq 644 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25760 --lastPoint 25799 -n .ctzctwObs_noJEC.POINTS.25760.25799
fi
if [ $1 -eq 645 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25800 --lastPoint 25839 -n .ctzctwObs_noJEC.POINTS.25800.25839
fi
if [ $1 -eq 646 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25840 --lastPoint 25879 -n .ctzctwObs_noJEC.POINTS.25840.25879
fi
if [ $1 -eq 647 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25880 --lastPoint 25919 -n .ctzctwObs_noJEC.POINTS.25880.25919
fi
if [ $1 -eq 648 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25920 --lastPoint 25959 -n .ctzctwObs_noJEC.POINTS.25920.25959
fi
if [ $1 -eq 649 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 25960 --lastPoint 25999 -n .ctzctwObs_noJEC.POINTS.25960.25999
fi
if [ $1 -eq 650 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26000 --lastPoint 26039 -n .ctzctwObs_noJEC.POINTS.26000.26039
fi
if [ $1 -eq 651 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26040 --lastPoint 26079 -n .ctzctwObs_noJEC.POINTS.26040.26079
fi
if [ $1 -eq 652 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26080 --lastPoint 26119 -n .ctzctwObs_noJEC.POINTS.26080.26119
fi
if [ $1 -eq 653 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26120 --lastPoint 26159 -n .ctzctwObs_noJEC.POINTS.26120.26159
fi
if [ $1 -eq 654 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26160 --lastPoint 26199 -n .ctzctwObs_noJEC.POINTS.26160.26199
fi
if [ $1 -eq 655 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26200 --lastPoint 26239 -n .ctzctwObs_noJEC.POINTS.26200.26239
fi
if [ $1 -eq 656 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26240 --lastPoint 26279 -n .ctzctwObs_noJEC.POINTS.26240.26279
fi
if [ $1 -eq 657 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26280 --lastPoint 26319 -n .ctzctwObs_noJEC.POINTS.26280.26319
fi
if [ $1 -eq 658 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26320 --lastPoint 26359 -n .ctzctwObs_noJEC.POINTS.26320.26359
fi
if [ $1 -eq 659 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26360 --lastPoint 26399 -n .ctzctwObs_noJEC.POINTS.26360.26399
fi
if [ $1 -eq 660 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26400 --lastPoint 26439 -n .ctzctwObs_noJEC.POINTS.26400.26439
fi
if [ $1 -eq 661 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26440 --lastPoint 26479 -n .ctzctwObs_noJEC.POINTS.26440.26479
fi
if [ $1 -eq 662 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26480 --lastPoint 26519 -n .ctzctwObs_noJEC.POINTS.26480.26519
fi
if [ $1 -eq 663 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26520 --lastPoint 26559 -n .ctzctwObs_noJEC.POINTS.26520.26559
fi
if [ $1 -eq 664 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26560 --lastPoint 26599 -n .ctzctwObs_noJEC.POINTS.26560.26599
fi
if [ $1 -eq 665 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26600 --lastPoint 26639 -n .ctzctwObs_noJEC.POINTS.26600.26639
fi
if [ $1 -eq 666 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26640 --lastPoint 26679 -n .ctzctwObs_noJEC.POINTS.26640.26679
fi
if [ $1 -eq 667 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26680 --lastPoint 26719 -n .ctzctwObs_noJEC.POINTS.26680.26719
fi
if [ $1 -eq 668 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26720 --lastPoint 26759 -n .ctzctwObs_noJEC.POINTS.26720.26759
fi
if [ $1 -eq 669 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26760 --lastPoint 26799 -n .ctzctwObs_noJEC.POINTS.26760.26799
fi
if [ $1 -eq 670 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26800 --lastPoint 26839 -n .ctzctwObs_noJEC.POINTS.26800.26839
fi
if [ $1 -eq 671 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26840 --lastPoint 26879 -n .ctzctwObs_noJEC.POINTS.26840.26879
fi
if [ $1 -eq 672 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26880 --lastPoint 26919 -n .ctzctwObs_noJEC.POINTS.26880.26919
fi
if [ $1 -eq 673 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26920 --lastPoint 26959 -n .ctzctwObs_noJEC.POINTS.26920.26959
fi
if [ $1 -eq 674 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 26960 --lastPoint 26999 -n .ctzctwObs_noJEC.POINTS.26960.26999
fi
if [ $1 -eq 675 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27000 --lastPoint 27039 -n .ctzctwObs_noJEC.POINTS.27000.27039
fi
if [ $1 -eq 676 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27040 --lastPoint 27079 -n .ctzctwObs_noJEC.POINTS.27040.27079
fi
if [ $1 -eq 677 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27080 --lastPoint 27119 -n .ctzctwObs_noJEC.POINTS.27080.27119
fi
if [ $1 -eq 678 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27120 --lastPoint 27159 -n .ctzctwObs_noJEC.POINTS.27120.27159
fi
if [ $1 -eq 679 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27160 --lastPoint 27199 -n .ctzctwObs_noJEC.POINTS.27160.27199
fi
if [ $1 -eq 680 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27200 --lastPoint 27239 -n .ctzctwObs_noJEC.POINTS.27200.27239
fi
if [ $1 -eq 681 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27240 --lastPoint 27279 -n .ctzctwObs_noJEC.POINTS.27240.27279
fi
if [ $1 -eq 682 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27280 --lastPoint 27319 -n .ctzctwObs_noJEC.POINTS.27280.27319
fi
if [ $1 -eq 683 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27320 --lastPoint 27359 -n .ctzctwObs_noJEC.POINTS.27320.27359
fi
if [ $1 -eq 684 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27360 --lastPoint 27399 -n .ctzctwObs_noJEC.POINTS.27360.27399
fi
if [ $1 -eq 685 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27400 --lastPoint 27439 -n .ctzctwObs_noJEC.POINTS.27400.27439
fi
if [ $1 -eq 686 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27440 --lastPoint 27479 -n .ctzctwObs_noJEC.POINTS.27440.27479
fi
if [ $1 -eq 687 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27480 --lastPoint 27519 -n .ctzctwObs_noJEC.POINTS.27480.27519
fi
if [ $1 -eq 688 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27520 --lastPoint 27559 -n .ctzctwObs_noJEC.POINTS.27520.27559
fi
if [ $1 -eq 689 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27560 --lastPoint 27599 -n .ctzctwObs_noJEC.POINTS.27560.27599
fi
if [ $1 -eq 690 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27600 --lastPoint 27639 -n .ctzctwObs_noJEC.POINTS.27600.27639
fi
if [ $1 -eq 691 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27640 --lastPoint 27679 -n .ctzctwObs_noJEC.POINTS.27640.27679
fi
if [ $1 -eq 692 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27680 --lastPoint 27719 -n .ctzctwObs_noJEC.POINTS.27680.27719
fi
if [ $1 -eq 693 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27720 --lastPoint 27759 -n .ctzctwObs_noJEC.POINTS.27720.27759
fi
if [ $1 -eq 694 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27760 --lastPoint 27799 -n .ctzctwObs_noJEC.POINTS.27760.27799
fi
if [ $1 -eq 695 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27800 --lastPoint 27839 -n .ctzctwObs_noJEC.POINTS.27800.27839
fi
if [ $1 -eq 696 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27840 --lastPoint 27879 -n .ctzctwObs_noJEC.POINTS.27840.27879
fi
if [ $1 -eq 697 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27880 --lastPoint 27919 -n .ctzctwObs_noJEC.POINTS.27880.27919
fi
if [ $1 -eq 698 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27920 --lastPoint 27959 -n .ctzctwObs_noJEC.POINTS.27920.27959
fi
if [ $1 -eq 699 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 27960 --lastPoint 27999 -n .ctzctwObs_noJEC.POINTS.27960.27999
fi
if [ $1 -eq 700 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28000 --lastPoint 28039 -n .ctzctwObs_noJEC.POINTS.28000.28039
fi
if [ $1 -eq 701 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28040 --lastPoint 28079 -n .ctzctwObs_noJEC.POINTS.28040.28079
fi
if [ $1 -eq 702 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28080 --lastPoint 28119 -n .ctzctwObs_noJEC.POINTS.28080.28119
fi
if [ $1 -eq 703 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28120 --lastPoint 28159 -n .ctzctwObs_noJEC.POINTS.28120.28159
fi
if [ $1 -eq 704 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28160 --lastPoint 28199 -n .ctzctwObs_noJEC.POINTS.28160.28199
fi
if [ $1 -eq 705 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28200 --lastPoint 28239 -n .ctzctwObs_noJEC.POINTS.28200.28239
fi
if [ $1 -eq 706 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28240 --lastPoint 28279 -n .ctzctwObs_noJEC.POINTS.28240.28279
fi
if [ $1 -eq 707 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28280 --lastPoint 28319 -n .ctzctwObs_noJEC.POINTS.28280.28319
fi
if [ $1 -eq 708 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28320 --lastPoint 28359 -n .ctzctwObs_noJEC.POINTS.28320.28359
fi
if [ $1 -eq 709 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28360 --lastPoint 28399 -n .ctzctwObs_noJEC.POINTS.28360.28399
fi
if [ $1 -eq 710 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28400 --lastPoint 28439 -n .ctzctwObs_noJEC.POINTS.28400.28439
fi
if [ $1 -eq 711 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28440 --lastPoint 28479 -n .ctzctwObs_noJEC.POINTS.28440.28479
fi
if [ $1 -eq 712 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28480 --lastPoint 28519 -n .ctzctwObs_noJEC.POINTS.28480.28519
fi
if [ $1 -eq 713 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28520 --lastPoint 28559 -n .ctzctwObs_noJEC.POINTS.28520.28559
fi
if [ $1 -eq 714 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28560 --lastPoint 28599 -n .ctzctwObs_noJEC.POINTS.28560.28599
fi
if [ $1 -eq 715 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28600 --lastPoint 28639 -n .ctzctwObs_noJEC.POINTS.28600.28639
fi
if [ $1 -eq 716 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28640 --lastPoint 28679 -n .ctzctwObs_noJEC.POINTS.28640.28679
fi
if [ $1 -eq 717 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28680 --lastPoint 28719 -n .ctzctwObs_noJEC.POINTS.28680.28719
fi
if [ $1 -eq 718 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28720 --lastPoint 28759 -n .ctzctwObs_noJEC.POINTS.28720.28759
fi
if [ $1 -eq 719 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28760 --lastPoint 28799 -n .ctzctwObs_noJEC.POINTS.28760.28799
fi
if [ $1 -eq 720 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28800 --lastPoint 28839 -n .ctzctwObs_noJEC.POINTS.28800.28839
fi
if [ $1 -eq 721 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28840 --lastPoint 28879 -n .ctzctwObs_noJEC.POINTS.28840.28879
fi
if [ $1 -eq 722 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28880 --lastPoint 28919 -n .ctzctwObs_noJEC.POINTS.28880.28919
fi
if [ $1 -eq 723 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28920 --lastPoint 28959 -n .ctzctwObs_noJEC.POINTS.28920.28959
fi
if [ $1 -eq 724 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 28960 --lastPoint 28999 -n .ctzctwObs_noJEC.POINTS.28960.28999
fi
if [ $1 -eq 725 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29000 --lastPoint 29039 -n .ctzctwObs_noJEC.POINTS.29000.29039
fi
if [ $1 -eq 726 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29040 --lastPoint 29079 -n .ctzctwObs_noJEC.POINTS.29040.29079
fi
if [ $1 -eq 727 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29080 --lastPoint 29119 -n .ctzctwObs_noJEC.POINTS.29080.29119
fi
if [ $1 -eq 728 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29120 --lastPoint 29159 -n .ctzctwObs_noJEC.POINTS.29120.29159
fi
if [ $1 -eq 729 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29160 --lastPoint 29199 -n .ctzctwObs_noJEC.POINTS.29160.29199
fi
if [ $1 -eq 730 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29200 --lastPoint 29239 -n .ctzctwObs_noJEC.POINTS.29200.29239
fi
if [ $1 -eq 731 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29240 --lastPoint 29279 -n .ctzctwObs_noJEC.POINTS.29240.29279
fi
if [ $1 -eq 732 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29280 --lastPoint 29319 -n .ctzctwObs_noJEC.POINTS.29280.29319
fi
if [ $1 -eq 733 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29320 --lastPoint 29359 -n .ctzctwObs_noJEC.POINTS.29320.29359
fi
if [ $1 -eq 734 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29360 --lastPoint 29399 -n .ctzctwObs_noJEC.POINTS.29360.29399
fi
if [ $1 -eq 735 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29400 --lastPoint 29439 -n .ctzctwObs_noJEC.POINTS.29400.29439
fi
if [ $1 -eq 736 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29440 --lastPoint 29479 -n .ctzctwObs_noJEC.POINTS.29440.29479
fi
if [ $1 -eq 737 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29480 --lastPoint 29519 -n .ctzctwObs_noJEC.POINTS.29480.29519
fi
if [ $1 -eq 738 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29520 --lastPoint 29559 -n .ctzctwObs_noJEC.POINTS.29520.29559
fi
if [ $1 -eq 739 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29560 --lastPoint 29599 -n .ctzctwObs_noJEC.POINTS.29560.29599
fi
if [ $1 -eq 740 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29600 --lastPoint 29639 -n .ctzctwObs_noJEC.POINTS.29600.29639
fi
if [ $1 -eq 741 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29640 --lastPoint 29679 -n .ctzctwObs_noJEC.POINTS.29640.29679
fi
if [ $1 -eq 742 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29680 --lastPoint 29719 -n .ctzctwObs_noJEC.POINTS.29680.29719
fi
if [ $1 -eq 743 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29720 --lastPoint 29759 -n .ctzctwObs_noJEC.POINTS.29720.29759
fi
if [ $1 -eq 744 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29760 --lastPoint 29799 -n .ctzctwObs_noJEC.POINTS.29760.29799
fi
if [ $1 -eq 745 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29800 --lastPoint 29839 -n .ctzctwObs_noJEC.POINTS.29800.29839
fi
if [ $1 -eq 746 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29840 --lastPoint 29879 -n .ctzctwObs_noJEC.POINTS.29840.29879
fi
if [ $1 -eq 747 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29880 --lastPoint 29919 -n .ctzctwObs_noJEC.POINTS.29880.29919
fi
if [ $1 -eq 748 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29920 --lastPoint 29959 -n .ctzctwObs_noJEC.POINTS.29920.29959
fi
if [ $1 -eq 749 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 29960 --lastPoint 29999 -n .ctzctwObs_noJEC.POINTS.29960.29999
fi
if [ $1 -eq 750 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30000 --lastPoint 30039 -n .ctzctwObs_noJEC.POINTS.30000.30039
fi
if [ $1 -eq 751 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30040 --lastPoint 30079 -n .ctzctwObs_noJEC.POINTS.30040.30079
fi
if [ $1 -eq 752 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30080 --lastPoint 30119 -n .ctzctwObs_noJEC.POINTS.30080.30119
fi
if [ $1 -eq 753 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30120 --lastPoint 30159 -n .ctzctwObs_noJEC.POINTS.30120.30159
fi
if [ $1 -eq 754 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30160 --lastPoint 30199 -n .ctzctwObs_noJEC.POINTS.30160.30199
fi
if [ $1 -eq 755 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30200 --lastPoint 30239 -n .ctzctwObs_noJEC.POINTS.30200.30239
fi
if [ $1 -eq 756 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30240 --lastPoint 30279 -n .ctzctwObs_noJEC.POINTS.30240.30279
fi
if [ $1 -eq 757 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30280 --lastPoint 30319 -n .ctzctwObs_noJEC.POINTS.30280.30319
fi
if [ $1 -eq 758 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30320 --lastPoint 30359 -n .ctzctwObs_noJEC.POINTS.30320.30359
fi
if [ $1 -eq 759 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30360 --lastPoint 30399 -n .ctzctwObs_noJEC.POINTS.30360.30399
fi
if [ $1 -eq 760 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30400 --lastPoint 30439 -n .ctzctwObs_noJEC.POINTS.30400.30439
fi
if [ $1 -eq 761 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30440 --lastPoint 30479 -n .ctzctwObs_noJEC.POINTS.30440.30479
fi
if [ $1 -eq 762 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30480 --lastPoint 30519 -n .ctzctwObs_noJEC.POINTS.30480.30519
fi
if [ $1 -eq 763 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30520 --lastPoint 30559 -n .ctzctwObs_noJEC.POINTS.30520.30559
fi
if [ $1 -eq 764 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30560 --lastPoint 30599 -n .ctzctwObs_noJEC.POINTS.30560.30599
fi
if [ $1 -eq 765 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30600 --lastPoint 30639 -n .ctzctwObs_noJEC.POINTS.30600.30639
fi
if [ $1 -eq 766 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30640 --lastPoint 30679 -n .ctzctwObs_noJEC.POINTS.30640.30679
fi
if [ $1 -eq 767 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30680 --lastPoint 30719 -n .ctzctwObs_noJEC.POINTS.30680.30719
fi
if [ $1 -eq 768 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30720 --lastPoint 30759 -n .ctzctwObs_noJEC.POINTS.30720.30759
fi
if [ $1 -eq 769 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30760 --lastPoint 30799 -n .ctzctwObs_noJEC.POINTS.30760.30799
fi
if [ $1 -eq 770 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30800 --lastPoint 30839 -n .ctzctwObs_noJEC.POINTS.30800.30839
fi
if [ $1 -eq 771 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30840 --lastPoint 30879 -n .ctzctwObs_noJEC.POINTS.30840.30879
fi
if [ $1 -eq 772 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30880 --lastPoint 30919 -n .ctzctwObs_noJEC.POINTS.30880.30919
fi
if [ $1 -eq 773 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30920 --lastPoint 30959 -n .ctzctwObs_noJEC.POINTS.30920.30959
fi
if [ $1 -eq 774 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 30960 --lastPoint 30999 -n .ctzctwObs_noJEC.POINTS.30960.30999
fi
if [ $1 -eq 775 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31000 --lastPoint 31039 -n .ctzctwObs_noJEC.POINTS.31000.31039
fi
if [ $1 -eq 776 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31040 --lastPoint 31079 -n .ctzctwObs_noJEC.POINTS.31040.31079
fi
if [ $1 -eq 777 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31080 --lastPoint 31119 -n .ctzctwObs_noJEC.POINTS.31080.31119
fi
if [ $1 -eq 778 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31120 --lastPoint 31159 -n .ctzctwObs_noJEC.POINTS.31120.31159
fi
if [ $1 -eq 779 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31160 --lastPoint 31199 -n .ctzctwObs_noJEC.POINTS.31160.31199
fi
if [ $1 -eq 780 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31200 --lastPoint 31239 -n .ctzctwObs_noJEC.POINTS.31200.31239
fi
if [ $1 -eq 781 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31240 --lastPoint 31279 -n .ctzctwObs_noJEC.POINTS.31240.31279
fi
if [ $1 -eq 782 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31280 --lastPoint 31319 -n .ctzctwObs_noJEC.POINTS.31280.31319
fi
if [ $1 -eq 783 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31320 --lastPoint 31359 -n .ctzctwObs_noJEC.POINTS.31320.31359
fi
if [ $1 -eq 784 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31360 --lastPoint 31399 -n .ctzctwObs_noJEC.POINTS.31360.31399
fi
if [ $1 -eq 785 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31400 --lastPoint 31439 -n .ctzctwObs_noJEC.POINTS.31400.31439
fi
if [ $1 -eq 786 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31440 --lastPoint 31479 -n .ctzctwObs_noJEC.POINTS.31440.31479
fi
if [ $1 -eq 787 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31480 --lastPoint 31519 -n .ctzctwObs_noJEC.POINTS.31480.31519
fi
if [ $1 -eq 788 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31520 --lastPoint 31559 -n .ctzctwObs_noJEC.POINTS.31520.31559
fi
if [ $1 -eq 789 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31560 --lastPoint 31599 -n .ctzctwObs_noJEC.POINTS.31560.31599
fi
if [ $1 -eq 790 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31600 --lastPoint 31639 -n .ctzctwObs_noJEC.POINTS.31600.31639
fi
if [ $1 -eq 791 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31640 --lastPoint 31679 -n .ctzctwObs_noJEC.POINTS.31640.31679
fi
if [ $1 -eq 792 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31680 --lastPoint 31719 -n .ctzctwObs_noJEC.POINTS.31680.31719
fi
if [ $1 -eq 793 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31720 --lastPoint 31759 -n .ctzctwObs_noJEC.POINTS.31720.31759
fi
if [ $1 -eq 794 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31760 --lastPoint 31799 -n .ctzctwObs_noJEC.POINTS.31760.31799
fi
if [ $1 -eq 795 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31800 --lastPoint 31839 -n .ctzctwObs_noJEC.POINTS.31800.31839
fi
if [ $1 -eq 796 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31840 --lastPoint 31879 -n .ctzctwObs_noJEC.POINTS.31840.31879
fi
if [ $1 -eq 797 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31880 --lastPoint 31919 -n .ctzctwObs_noJEC.POINTS.31880.31919
fi
if [ $1 -eq 798 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31920 --lastPoint 31959 -n .ctzctwObs_noJEC.POINTS.31920.31959
fi
if [ $1 -eq 799 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 31960 --lastPoint 31999 -n .ctzctwObs_noJEC.POINTS.31960.31999
fi
if [ $1 -eq 800 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32000 --lastPoint 32039 -n .ctzctwObs_noJEC.POINTS.32000.32039
fi
if [ $1 -eq 801 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32040 --lastPoint 32079 -n .ctzctwObs_noJEC.POINTS.32040.32079
fi
if [ $1 -eq 802 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32080 --lastPoint 32119 -n .ctzctwObs_noJEC.POINTS.32080.32119
fi
if [ $1 -eq 803 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32120 --lastPoint 32159 -n .ctzctwObs_noJEC.POINTS.32120.32159
fi
if [ $1 -eq 804 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32160 --lastPoint 32199 -n .ctzctwObs_noJEC.POINTS.32160.32199
fi
if [ $1 -eq 805 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32200 --lastPoint 32239 -n .ctzctwObs_noJEC.POINTS.32200.32239
fi
if [ $1 -eq 806 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32240 --lastPoint 32279 -n .ctzctwObs_noJEC.POINTS.32240.32279
fi
if [ $1 -eq 807 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32280 --lastPoint 32319 -n .ctzctwObs_noJEC.POINTS.32280.32319
fi
if [ $1 -eq 808 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32320 --lastPoint 32359 -n .ctzctwObs_noJEC.POINTS.32320.32359
fi
if [ $1 -eq 809 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32360 --lastPoint 32399 -n .ctzctwObs_noJEC.POINTS.32360.32399
fi
if [ $1 -eq 810 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32400 --lastPoint 32439 -n .ctzctwObs_noJEC.POINTS.32400.32439
fi
if [ $1 -eq 811 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32440 --lastPoint 32479 -n .ctzctwObs_noJEC.POINTS.32440.32479
fi
if [ $1 -eq 812 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32480 --lastPoint 32519 -n .ctzctwObs_noJEC.POINTS.32480.32519
fi
if [ $1 -eq 813 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32520 --lastPoint 32559 -n .ctzctwObs_noJEC.POINTS.32520.32559
fi
if [ $1 -eq 814 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32560 --lastPoint 32599 -n .ctzctwObs_noJEC.POINTS.32560.32599
fi
if [ $1 -eq 815 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32600 --lastPoint 32639 -n .ctzctwObs_noJEC.POINTS.32600.32639
fi
if [ $1 -eq 816 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32640 --lastPoint 32679 -n .ctzctwObs_noJEC.POINTS.32640.32679
fi
if [ $1 -eq 817 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32680 --lastPoint 32719 -n .ctzctwObs_noJEC.POINTS.32680.32719
fi
if [ $1 -eq 818 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32720 --lastPoint 32759 -n .ctzctwObs_noJEC.POINTS.32720.32759
fi
if [ $1 -eq 819 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32760 --lastPoint 32799 -n .ctzctwObs_noJEC.POINTS.32760.32799
fi
if [ $1 -eq 820 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32800 --lastPoint 32839 -n .ctzctwObs_noJEC.POINTS.32800.32839
fi
if [ $1 -eq 821 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32840 --lastPoint 32879 -n .ctzctwObs_noJEC.POINTS.32840.32879
fi
if [ $1 -eq 822 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32880 --lastPoint 32919 -n .ctzctwObs_noJEC.POINTS.32880.32919
fi
if [ $1 -eq 823 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32920 --lastPoint 32959 -n .ctzctwObs_noJEC.POINTS.32920.32959
fi
if [ $1 -eq 824 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 32960 --lastPoint 32999 -n .ctzctwObs_noJEC.POINTS.32960.32999
fi
if [ $1 -eq 825 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33000 --lastPoint 33039 -n .ctzctwObs_noJEC.POINTS.33000.33039
fi
if [ $1 -eq 826 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33040 --lastPoint 33079 -n .ctzctwObs_noJEC.POINTS.33040.33079
fi
if [ $1 -eq 827 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33080 --lastPoint 33119 -n .ctzctwObs_noJEC.POINTS.33080.33119
fi
if [ $1 -eq 828 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33120 --lastPoint 33159 -n .ctzctwObs_noJEC.POINTS.33120.33159
fi
if [ $1 -eq 829 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33160 --lastPoint 33199 -n .ctzctwObs_noJEC.POINTS.33160.33199
fi
if [ $1 -eq 830 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33200 --lastPoint 33239 -n .ctzctwObs_noJEC.POINTS.33200.33239
fi
if [ $1 -eq 831 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33240 --lastPoint 33279 -n .ctzctwObs_noJEC.POINTS.33240.33279
fi
if [ $1 -eq 832 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33280 --lastPoint 33319 -n .ctzctwObs_noJEC.POINTS.33280.33319
fi
if [ $1 -eq 833 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33320 --lastPoint 33359 -n .ctzctwObs_noJEC.POINTS.33320.33359
fi
if [ $1 -eq 834 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33360 --lastPoint 33399 -n .ctzctwObs_noJEC.POINTS.33360.33399
fi
if [ $1 -eq 835 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33400 --lastPoint 33439 -n .ctzctwObs_noJEC.POINTS.33400.33439
fi
if [ $1 -eq 836 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33440 --lastPoint 33479 -n .ctzctwObs_noJEC.POINTS.33440.33479
fi
if [ $1 -eq 837 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33480 --lastPoint 33519 -n .ctzctwObs_noJEC.POINTS.33480.33519
fi
if [ $1 -eq 838 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33520 --lastPoint 33559 -n .ctzctwObs_noJEC.POINTS.33520.33559
fi
if [ $1 -eq 839 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33560 --lastPoint 33599 -n .ctzctwObs_noJEC.POINTS.33560.33599
fi
if [ $1 -eq 840 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33600 --lastPoint 33639 -n .ctzctwObs_noJEC.POINTS.33600.33639
fi
if [ $1 -eq 841 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33640 --lastPoint 33679 -n .ctzctwObs_noJEC.POINTS.33640.33679
fi
if [ $1 -eq 842 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33680 --lastPoint 33719 -n .ctzctwObs_noJEC.POINTS.33680.33719
fi
if [ $1 -eq 843 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33720 --lastPoint 33759 -n .ctzctwObs_noJEC.POINTS.33720.33759
fi
if [ $1 -eq 844 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33760 --lastPoint 33799 -n .ctzctwObs_noJEC.POINTS.33760.33799
fi
if [ $1 -eq 845 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33800 --lastPoint 33839 -n .ctzctwObs_noJEC.POINTS.33800.33839
fi
if [ $1 -eq 846 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33840 --lastPoint 33879 -n .ctzctwObs_noJEC.POINTS.33840.33879
fi
if [ $1 -eq 847 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33880 --lastPoint 33919 -n .ctzctwObs_noJEC.POINTS.33880.33919
fi
if [ $1 -eq 848 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33920 --lastPoint 33959 -n .ctzctwObs_noJEC.POINTS.33920.33959
fi
if [ $1 -eq 849 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 33960 --lastPoint 33999 -n .ctzctwObs_noJEC.POINTS.33960.33999
fi
if [ $1 -eq 850 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34000 --lastPoint 34039 -n .ctzctwObs_noJEC.POINTS.34000.34039
fi
if [ $1 -eq 851 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34040 --lastPoint 34079 -n .ctzctwObs_noJEC.POINTS.34040.34079
fi
if [ $1 -eq 852 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34080 --lastPoint 34119 -n .ctzctwObs_noJEC.POINTS.34080.34119
fi
if [ $1 -eq 853 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34120 --lastPoint 34159 -n .ctzctwObs_noJEC.POINTS.34120.34159
fi
if [ $1 -eq 854 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34160 --lastPoint 34199 -n .ctzctwObs_noJEC.POINTS.34160.34199
fi
if [ $1 -eq 855 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34200 --lastPoint 34239 -n .ctzctwObs_noJEC.POINTS.34200.34239
fi
if [ $1 -eq 856 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34240 --lastPoint 34279 -n .ctzctwObs_noJEC.POINTS.34240.34279
fi
if [ $1 -eq 857 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34280 --lastPoint 34319 -n .ctzctwObs_noJEC.POINTS.34280.34319
fi
if [ $1 -eq 858 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34320 --lastPoint 34359 -n .ctzctwObs_noJEC.POINTS.34320.34359
fi
if [ $1 -eq 859 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34360 --lastPoint 34399 -n .ctzctwObs_noJEC.POINTS.34360.34399
fi
if [ $1 -eq 860 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34400 --lastPoint 34439 -n .ctzctwObs_noJEC.POINTS.34400.34439
fi
if [ $1 -eq 861 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34440 --lastPoint 34479 -n .ctzctwObs_noJEC.POINTS.34440.34479
fi
if [ $1 -eq 862 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34480 --lastPoint 34519 -n .ctzctwObs_noJEC.POINTS.34480.34519
fi
if [ $1 -eq 863 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34520 --lastPoint 34559 -n .ctzctwObs_noJEC.POINTS.34520.34559
fi
if [ $1 -eq 864 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34560 --lastPoint 34599 -n .ctzctwObs_noJEC.POINTS.34560.34599
fi
if [ $1 -eq 865 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34600 --lastPoint 34639 -n .ctzctwObs_noJEC.POINTS.34600.34639
fi
if [ $1 -eq 866 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34640 --lastPoint 34679 -n .ctzctwObs_noJEC.POINTS.34640.34679
fi
if [ $1 -eq 867 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34680 --lastPoint 34719 -n .ctzctwObs_noJEC.POINTS.34680.34719
fi
if [ $1 -eq 868 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34720 --lastPoint 34759 -n .ctzctwObs_noJEC.POINTS.34720.34759
fi
if [ $1 -eq 869 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34760 --lastPoint 34799 -n .ctzctwObs_noJEC.POINTS.34760.34799
fi
if [ $1 -eq 870 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34800 --lastPoint 34839 -n .ctzctwObs_noJEC.POINTS.34800.34839
fi
if [ $1 -eq 871 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34840 --lastPoint 34879 -n .ctzctwObs_noJEC.POINTS.34840.34879
fi
if [ $1 -eq 872 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34880 --lastPoint 34919 -n .ctzctwObs_noJEC.POINTS.34880.34919
fi
if [ $1 -eq 873 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34920 --lastPoint 34959 -n .ctzctwObs_noJEC.POINTS.34920.34959
fi
if [ $1 -eq 874 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 34960 --lastPoint 34999 -n .ctzctwObs_noJEC.POINTS.34960.34999
fi
if [ $1 -eq 875 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35000 --lastPoint 35039 -n .ctzctwObs_noJEC.POINTS.35000.35039
fi
if [ $1 -eq 876 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35040 --lastPoint 35079 -n .ctzctwObs_noJEC.POINTS.35040.35079
fi
if [ $1 -eq 877 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35080 --lastPoint 35119 -n .ctzctwObs_noJEC.POINTS.35080.35119
fi
if [ $1 -eq 878 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35120 --lastPoint 35159 -n .ctzctwObs_noJEC.POINTS.35120.35159
fi
if [ $1 -eq 879 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35160 --lastPoint 35199 -n .ctzctwObs_noJEC.POINTS.35160.35199
fi
if [ $1 -eq 880 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35200 --lastPoint 35239 -n .ctzctwObs_noJEC.POINTS.35200.35239
fi
if [ $1 -eq 881 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35240 --lastPoint 35279 -n .ctzctwObs_noJEC.POINTS.35240.35279
fi
if [ $1 -eq 882 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35280 --lastPoint 35319 -n .ctzctwObs_noJEC.POINTS.35280.35319
fi
if [ $1 -eq 883 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35320 --lastPoint 35359 -n .ctzctwObs_noJEC.POINTS.35320.35359
fi
if [ $1 -eq 884 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35360 --lastPoint 35399 -n .ctzctwObs_noJEC.POINTS.35360.35399
fi
if [ $1 -eq 885 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35400 --lastPoint 35439 -n .ctzctwObs_noJEC.POINTS.35400.35439
fi
if [ $1 -eq 886 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35440 --lastPoint 35479 -n .ctzctwObs_noJEC.POINTS.35440.35479
fi
if [ $1 -eq 887 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35480 --lastPoint 35519 -n .ctzctwObs_noJEC.POINTS.35480.35519
fi
if [ $1 -eq 888 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35520 --lastPoint 35559 -n .ctzctwObs_noJEC.POINTS.35520.35559
fi
if [ $1 -eq 889 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35560 --lastPoint 35599 -n .ctzctwObs_noJEC.POINTS.35560.35599
fi
if [ $1 -eq 890 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35600 --lastPoint 35639 -n .ctzctwObs_noJEC.POINTS.35600.35639
fi
if [ $1 -eq 891 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35640 --lastPoint 35679 -n .ctzctwObs_noJEC.POINTS.35640.35679
fi
if [ $1 -eq 892 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35680 --lastPoint 35719 -n .ctzctwObs_noJEC.POINTS.35680.35719
fi
if [ $1 -eq 893 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35720 --lastPoint 35759 -n .ctzctwObs_noJEC.POINTS.35720.35759
fi
if [ $1 -eq 894 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35760 --lastPoint 35799 -n .ctzctwObs_noJEC.POINTS.35760.35799
fi
if [ $1 -eq 895 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35800 --lastPoint 35839 -n .ctzctwObs_noJEC.POINTS.35800.35839
fi
if [ $1 -eq 896 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35840 --lastPoint 35879 -n .ctzctwObs_noJEC.POINTS.35840.35879
fi
if [ $1 -eq 897 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35880 --lastPoint 35919 -n .ctzctwObs_noJEC.POINTS.35880.35919
fi
if [ $1 -eq 898 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35920 --lastPoint 35959 -n .ctzctwObs_noJEC.POINTS.35920.35959
fi
if [ $1 -eq 899 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 35960 --lastPoint 35999 -n .ctzctwObs_noJEC.POINTS.35960.35999
fi
if [ $1 -eq 900 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36000 --lastPoint 36039 -n .ctzctwObs_noJEC.POINTS.36000.36039
fi
if [ $1 -eq 901 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36040 --lastPoint 36079 -n .ctzctwObs_noJEC.POINTS.36040.36079
fi
if [ $1 -eq 902 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36080 --lastPoint 36119 -n .ctzctwObs_noJEC.POINTS.36080.36119
fi
if [ $1 -eq 903 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36120 --lastPoint 36159 -n .ctzctwObs_noJEC.POINTS.36120.36159
fi
if [ $1 -eq 904 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36160 --lastPoint 36199 -n .ctzctwObs_noJEC.POINTS.36160.36199
fi
if [ $1 -eq 905 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36200 --lastPoint 36239 -n .ctzctwObs_noJEC.POINTS.36200.36239
fi
if [ $1 -eq 906 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36240 --lastPoint 36279 -n .ctzctwObs_noJEC.POINTS.36240.36279
fi
if [ $1 -eq 907 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36280 --lastPoint 36319 -n .ctzctwObs_noJEC.POINTS.36280.36319
fi
if [ $1 -eq 908 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36320 --lastPoint 36359 -n .ctzctwObs_noJEC.POINTS.36320.36359
fi
if [ $1 -eq 909 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36360 --lastPoint 36399 -n .ctzctwObs_noJEC.POINTS.36360.36399
fi
if [ $1 -eq 910 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36400 --lastPoint 36439 -n .ctzctwObs_noJEC.POINTS.36400.36439
fi
if [ $1 -eq 911 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36440 --lastPoint 36479 -n .ctzctwObs_noJEC.POINTS.36440.36479
fi
if [ $1 -eq 912 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36480 --lastPoint 36519 -n .ctzctwObs_noJEC.POINTS.36480.36519
fi
if [ $1 -eq 913 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36520 --lastPoint 36559 -n .ctzctwObs_noJEC.POINTS.36520.36559
fi
if [ $1 -eq 914 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36560 --lastPoint 36599 -n .ctzctwObs_noJEC.POINTS.36560.36599
fi
if [ $1 -eq 915 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36600 --lastPoint 36639 -n .ctzctwObs_noJEC.POINTS.36600.36639
fi
if [ $1 -eq 916 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36640 --lastPoint 36679 -n .ctzctwObs_noJEC.POINTS.36640.36679
fi
if [ $1 -eq 917 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36680 --lastPoint 36719 -n .ctzctwObs_noJEC.POINTS.36680.36719
fi
if [ $1 -eq 918 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36720 --lastPoint 36759 -n .ctzctwObs_noJEC.POINTS.36720.36759
fi
if [ $1 -eq 919 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36760 --lastPoint 36799 -n .ctzctwObs_noJEC.POINTS.36760.36799
fi
if [ $1 -eq 920 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36800 --lastPoint 36839 -n .ctzctwObs_noJEC.POINTS.36800.36839
fi
if [ $1 -eq 921 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36840 --lastPoint 36879 -n .ctzctwObs_noJEC.POINTS.36840.36879
fi
if [ $1 -eq 922 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36880 --lastPoint 36919 -n .ctzctwObs_noJEC.POINTS.36880.36919
fi
if [ $1 -eq 923 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36920 --lastPoint 36959 -n .ctzctwObs_noJEC.POINTS.36920.36959
fi
if [ $1 -eq 924 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 36960 --lastPoint 36999 -n .ctzctwObs_noJEC.POINTS.36960.36999
fi
if [ $1 -eq 925 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37000 --lastPoint 37039 -n .ctzctwObs_noJEC.POINTS.37000.37039
fi
if [ $1 -eq 926 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37040 --lastPoint 37079 -n .ctzctwObs_noJEC.POINTS.37040.37079
fi
if [ $1 -eq 927 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37080 --lastPoint 37119 -n .ctzctwObs_noJEC.POINTS.37080.37119
fi
if [ $1 -eq 928 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37120 --lastPoint 37159 -n .ctzctwObs_noJEC.POINTS.37120.37159
fi
if [ $1 -eq 929 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37160 --lastPoint 37199 -n .ctzctwObs_noJEC.POINTS.37160.37199
fi
if [ $1 -eq 930 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37200 --lastPoint 37239 -n .ctzctwObs_noJEC.POINTS.37200.37239
fi
if [ $1 -eq 931 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37240 --lastPoint 37279 -n .ctzctwObs_noJEC.POINTS.37240.37279
fi
if [ $1 -eq 932 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37280 --lastPoint 37319 -n .ctzctwObs_noJEC.POINTS.37280.37319
fi
if [ $1 -eq 933 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37320 --lastPoint 37359 -n .ctzctwObs_noJEC.POINTS.37320.37359
fi
if [ $1 -eq 934 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37360 --lastPoint 37399 -n .ctzctwObs_noJEC.POINTS.37360.37399
fi
if [ $1 -eq 935 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37400 --lastPoint 37439 -n .ctzctwObs_noJEC.POINTS.37400.37439
fi
if [ $1 -eq 936 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37440 --lastPoint 37479 -n .ctzctwObs_noJEC.POINTS.37440.37479
fi
if [ $1 -eq 937 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37480 --lastPoint 37519 -n .ctzctwObs_noJEC.POINTS.37480.37519
fi
if [ $1 -eq 938 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37520 --lastPoint 37559 -n .ctzctwObs_noJEC.POINTS.37520.37559
fi
if [ $1 -eq 939 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37560 --lastPoint 37599 -n .ctzctwObs_noJEC.POINTS.37560.37599
fi
if [ $1 -eq 940 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37600 --lastPoint 37639 -n .ctzctwObs_noJEC.POINTS.37600.37639
fi
if [ $1 -eq 941 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37640 --lastPoint 37679 -n .ctzctwObs_noJEC.POINTS.37640.37679
fi
if [ $1 -eq 942 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37680 --lastPoint 37719 -n .ctzctwObs_noJEC.POINTS.37680.37719
fi
if [ $1 -eq 943 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37720 --lastPoint 37759 -n .ctzctwObs_noJEC.POINTS.37720.37759
fi
if [ $1 -eq 944 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37760 --lastPoint 37799 -n .ctzctwObs_noJEC.POINTS.37760.37799
fi
if [ $1 -eq 945 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37800 --lastPoint 37839 -n .ctzctwObs_noJEC.POINTS.37800.37839
fi
if [ $1 -eq 946 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37840 --lastPoint 37879 -n .ctzctwObs_noJEC.POINTS.37840.37879
fi
if [ $1 -eq 947 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37880 --lastPoint 37919 -n .ctzctwObs_noJEC.POINTS.37880.37919
fi
if [ $1 -eq 948 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37920 --lastPoint 37959 -n .ctzctwObs_noJEC.POINTS.37920.37959
fi
if [ $1 -eq 949 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 37960 --lastPoint 37999 -n .ctzctwObs_noJEC.POINTS.37960.37999
fi
if [ $1 -eq 950 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38000 --lastPoint 38039 -n .ctzctwObs_noJEC.POINTS.38000.38039
fi
if [ $1 -eq 951 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38040 --lastPoint 38079 -n .ctzctwObs_noJEC.POINTS.38040.38079
fi
if [ $1 -eq 952 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38080 --lastPoint 38119 -n .ctzctwObs_noJEC.POINTS.38080.38119
fi
if [ $1 -eq 953 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38120 --lastPoint 38159 -n .ctzctwObs_noJEC.POINTS.38120.38159
fi
if [ $1 -eq 954 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38160 --lastPoint 38199 -n .ctzctwObs_noJEC.POINTS.38160.38199
fi
if [ $1 -eq 955 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38200 --lastPoint 38239 -n .ctzctwObs_noJEC.POINTS.38200.38239
fi
if [ $1 -eq 956 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38240 --lastPoint 38279 -n .ctzctwObs_noJEC.POINTS.38240.38279
fi
if [ $1 -eq 957 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38280 --lastPoint 38319 -n .ctzctwObs_noJEC.POINTS.38280.38319
fi
if [ $1 -eq 958 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38320 --lastPoint 38359 -n .ctzctwObs_noJEC.POINTS.38320.38359
fi
if [ $1 -eq 959 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38360 --lastPoint 38399 -n .ctzctwObs_noJEC.POINTS.38360.38399
fi
if [ $1 -eq 960 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38400 --lastPoint 38439 -n .ctzctwObs_noJEC.POINTS.38400.38439
fi
if [ $1 -eq 961 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38440 --lastPoint 38479 -n .ctzctwObs_noJEC.POINTS.38440.38479
fi
if [ $1 -eq 962 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38480 --lastPoint 38519 -n .ctzctwObs_noJEC.POINTS.38480.38519
fi
if [ $1 -eq 963 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38520 --lastPoint 38559 -n .ctzctwObs_noJEC.POINTS.38520.38559
fi
if [ $1 -eq 964 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38560 --lastPoint 38599 -n .ctzctwObs_noJEC.POINTS.38560.38599
fi
if [ $1 -eq 965 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38600 --lastPoint 38639 -n .ctzctwObs_noJEC.POINTS.38600.38639
fi
if [ $1 -eq 966 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38640 --lastPoint 38679 -n .ctzctwObs_noJEC.POINTS.38640.38679
fi
if [ $1 -eq 967 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38680 --lastPoint 38719 -n .ctzctwObs_noJEC.POINTS.38680.38719
fi
if [ $1 -eq 968 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38720 --lastPoint 38759 -n .ctzctwObs_noJEC.POINTS.38720.38759
fi
if [ $1 -eq 969 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38760 --lastPoint 38799 -n .ctzctwObs_noJEC.POINTS.38760.38799
fi
if [ $1 -eq 970 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38800 --lastPoint 38839 -n .ctzctwObs_noJEC.POINTS.38800.38839
fi
if [ $1 -eq 971 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38840 --lastPoint 38879 -n .ctzctwObs_noJEC.POINTS.38840.38879
fi
if [ $1 -eq 972 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38880 --lastPoint 38919 -n .ctzctwObs_noJEC.POINTS.38880.38919
fi
if [ $1 -eq 973 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38920 --lastPoint 38959 -n .ctzctwObs_noJEC.POINTS.38920.38959
fi
if [ $1 -eq 974 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 38960 --lastPoint 38999 -n .ctzctwObs_noJEC.POINTS.38960.38999
fi
if [ $1 -eq 975 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39000 --lastPoint 39039 -n .ctzctwObs_noJEC.POINTS.39000.39039
fi
if [ $1 -eq 976 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39040 --lastPoint 39079 -n .ctzctwObs_noJEC.POINTS.39040.39079
fi
if [ $1 -eq 977 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39080 --lastPoint 39119 -n .ctzctwObs_noJEC.POINTS.39080.39119
fi
if [ $1 -eq 978 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39120 --lastPoint 39159 -n .ctzctwObs_noJEC.POINTS.39120.39159
fi
if [ $1 -eq 979 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39160 --lastPoint 39199 -n .ctzctwObs_noJEC.POINTS.39160.39199
fi
if [ $1 -eq 980 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39200 --lastPoint 39239 -n .ctzctwObs_noJEC.POINTS.39200.39239
fi
if [ $1 -eq 981 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39240 --lastPoint 39279 -n .ctzctwObs_noJEC.POINTS.39240.39279
fi
if [ $1 -eq 982 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39280 --lastPoint 39319 -n .ctzctwObs_noJEC.POINTS.39280.39319
fi
if [ $1 -eq 983 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39320 --lastPoint 39359 -n .ctzctwObs_noJEC.POINTS.39320.39359
fi
if [ $1 -eq 984 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39360 --lastPoint 39399 -n .ctzctwObs_noJEC.POINTS.39360.39399
fi
if [ $1 -eq 985 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39400 --lastPoint 39439 -n .ctzctwObs_noJEC.POINTS.39400.39439
fi
if [ $1 -eq 986 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39440 --lastPoint 39479 -n .ctzctwObs_noJEC.POINTS.39440.39479
fi
if [ $1 -eq 987 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39480 --lastPoint 39519 -n .ctzctwObs_noJEC.POINTS.39480.39519
fi
if [ $1 -eq 988 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39520 --lastPoint 39559 -n .ctzctwObs_noJEC.POINTS.39520.39559
fi
if [ $1 -eq 989 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39560 --lastPoint 39599 -n .ctzctwObs_noJEC.POINTS.39560.39599
fi
if [ $1 -eq 990 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39600 --lastPoint 39639 -n .ctzctwObs_noJEC.POINTS.39600.39639
fi
if [ $1 -eq 991 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39640 --lastPoint 39679 -n .ctzctwObs_noJEC.POINTS.39640.39679
fi
if [ $1 -eq 992 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39680 --lastPoint 39719 -n .ctzctwObs_noJEC.POINTS.39680.39719
fi
if [ $1 -eq 993 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39720 --lastPoint 39759 -n .ctzctwObs_noJEC.POINTS.39720.39759
fi
if [ $1 -eq 994 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39760 --lastPoint 39799 -n .ctzctwObs_noJEC.POINTS.39760.39799
fi
if [ $1 -eq 995 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39800 --lastPoint 39839 -n .ctzctwObs_noJEC.POINTS.39800.39839
fi
if [ $1 -eq 996 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39840 --lastPoint 39879 -n .ctzctwObs_noJEC.POINTS.39840.39879
fi
if [ $1 -eq 997 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39880 --lastPoint 39919 -n .ctzctwObs_noJEC.POINTS.39880.39919
fi
if [ $1 -eq 998 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39920 --lastPoint 39959 -n .ctzctwObs_noJEC.POINTS.39920.39959
fi
if [ $1 -eq 999 ]; then
  combine --algo grid --cminPreScan --cminDefaultMinimizerStrategy=0 -P ctz -P ctw --setParameters ctz=0,ctw=0,cpq3=0,cpqm=0,cpt=0 --freezeParameters cpq3,cpqm,cpt --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --saveNLL --freezeNuisanceGroups jes -M MultiDimFit -d ./EFTWorkspace_5D.root --setParameterRanges ctz=-2,2:ctw=-1.5,1.5:cpq3=-2,2:cpqm=-15,25:cpt=-28,20 --points 40000 --firstPoint 39960 --lastPoint 39999 -n .ctzctwObs_noJEC.POINTS.39960.39999
fi