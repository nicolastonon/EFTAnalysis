#NB: not important to read correct EFTparam, since freezing WCs (?)

#-----------------------
#-- Not needed for paper

#./script_shapesFromWorkspace.sh multidimfit.ctzExp.root NN_ctz_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.ctzExp.root NN_ctz_SRttZ
#./script_shapesFromWorkspace.sh multidimfit.ctwExp.root NN_ctw_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.ctwExp.root NN_ctw_SRttZ
#./script_shapesFromWorkspace.sh multidimfit.cpq3Exp.root NN_cpq3_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.cpq3Exp.root NN_cpq3_SRttZ
#./script_shapesFromWorkspace.sh multidimfit.cpqmExp.root NN_cpqm_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.cpqmExp.root NN_cpqm_SRttZ
#./script_shapesFromWorkspace.sh multidimfit.cptExp.root NN_cpt_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.cptExp.root NN_cpt_SRttZ
#./script_shapesFromWorkspace.sh multidimfit.5DExp.root NN_5D_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.5DExp.root NN_5D_SRttZ

#./script_shapesFromWorkspace.sh multidimfit.cpqmObs.root NN_cpqm_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.cpqmObs.root NN_cpqm_SRttZ
#./script_shapesFromWorkspace.sh multidimfit.cptObs.root NN_cpt_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.cptObs.root NN_cpt_SRttZ

#./script_shapesFromWorkspace.sh multidimfit.5DExp.root CRWZ
#./script_shapesFromWorkspace.sh multidimfit.5DExp.root CRZZ
#./script_shapesFromWorkspace.sh multidimfit.5DExp.root SRttZ4l
#./script_shapesFromWorkspace.sh multidimfit.5DExp.root SRother

#-----------------------
#-- POSTFIT 
# /!\ /!\ Should run twice (and move output files in corresponding subdirs): once without freezing any parameter (to get correct central values), and once freezing WCs and split JEC (to get correct uncertainties) #NB: only for postfit /!\ /!\

#SR OBS

#./script_shapesFromWorkspace.sh multidimfit.ctzObs.root NN_ctz_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.ctzObs.root NN_ctz_SRttZ
#./script_shapesFromWorkspace.sh multidimfit.ctwObs.root NN_ctw_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.ctwObs.root NN_ctw_SRttZ
#./script_shapesFromWorkspace.sh multidimfit.cpq3Obs.root NN_cpq3_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.cpq3Obs.root NN_cpq3_SRttZ
#./script_shapesFromWorkspace.sh multidimfit.5DObs.root NN_5D_SRtZq
#./script_shapesFromWorkspace.sh multidimfit.5DObs.root NN_5D_SRttZ

#COMMON REGIONS OBS

#./script_shapesFromWorkspace.sh multidimfit.5DObs.root CRWZ
#./script_shapesFromWorkspace.sh multidimfit.5DObs.root CRZZ
#./script_shapesFromWorkspace.sh multidimfit.5DObs.root SRttZ4l
#./script_shapesFromWorkspace.sh multidimfit.5DObs.root SRother

#-----------------------
#-- PREFIT

#-- NN_SM PREFIT <-> hardcoded prefit shapes for all 3 NN-SM output nodes #NB: fit result not needed, but must use keyword 'SM' as first arg #NB: for these, must read template file 'Templates_NN_SM_SM_Run2.root' (because the file Templates_NN_SM_EFT2_Run2.root used for cpqm/cpt 1D limits uses mTW distribution in SRother !)
#./script_shapesFromWorkspace.sh NN_SM_SRtZq
#./script_shapesFromWorkspace.sh NN_SM_SRttZ
#./script_shapesFromWorkspace.sh NN_SM_SRother

#CONTROL PLOTS (prefit) #No keyword to filter #Prefit #Make sure relevant datacards were generated

#./script_shapesFromWorkspace.sh

#-----------------------
