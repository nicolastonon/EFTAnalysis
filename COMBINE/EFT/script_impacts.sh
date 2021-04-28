#-- CREATE WORKSPACE
#python EFTFitter.py -d ../datacards/COMBINED_Datacard_TemplateFit_Run2.txt --onlyworkspace

#-- PRODUCE IMPACTS

#python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctz -n ctzExp -d ./EFTWorkspace_ctz.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctw -n ctwExp -d ./EFTWorkspace_ctw.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpq3 -n cpq3Exp -d ./EFTWorkspace_cpq3.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpqm -n cpqmExp -d ./EFTWorkspace_SM.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpt -n cptExp -d ./EFTWorkspace_SM.root -m impacts --batch condor

#python Plot_Syst_Impacts_EFT.py --freeze -P ctz -n ctzObs -d ./EFTWorkspace_ctz.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --freeze -P ctw -n ctwObs -d ./EFTWorkspace_ctw.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --freeze -P cpq3 -n cpq3Obs -d ./EFTWorkspace_cpq3.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --freeze -P cpqm -n cpqmObs -d ./EFTWorkspace_SM.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --freeze -P cpt -n cptObs -d ./EFTWorkspace_SM.root -m impacts --batch condor

#-- COLLECT OUTPUTS

#python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctz -n ctzExp -d ./EFTWorkspace_ctz.root -m impacts --collect
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctw -n ctwExp -d ./EFTWorkspace_ctw.root -m impacts --collect
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpq3 -n cpq3Exp -d ./EFTWorkspace_cpq3.root -m impacts --collect
python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpqm -n cpqmExp -d ./EFTWorkspace_SM.root -m impacts --collect
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpt -n cptExp -d ./EFTWorkspace_SM.root -m impacts --collect

#python Plot_Syst_Impacts_EFT.py --freeze -P ctz -n ctzObs -d ./EFTWorkspace_ctz.root -m impacts --collect
#python Plot_Syst_Impacts_EFT.py --freeze -P ctw -n ctwObs -d ./EFTWorkspace_ctw.root -m impacts --collect
#python Plot_Syst_Impacts_EFT.py --freeze -P cpq3 -n cpq3Obs -d ./EFTWorkspace_cpq3.root -m impacts --collect
#python Plot_Syst_Impacts_EFT.py --freeze -P cpqm -n cpqmObs -d ./EFTWorkspace_SM.root -m impacts --collect
#python Plot_Syst_Impacts_EFT.py --freeze -P cpt -n cptObs -d ./EFTWorkspace_SM.root -m impacts --collect

#-- 5D FIT
#NB: no sense to make impact plot for 5D fit... ? (else would need to adapt impact script similarly to David, for multiple POIs) -- but could check impacts on single POI

##python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctz ctw cpq3 cpqm cpt -n 5DExp -d ./EFTWorkspace_5D.root -m impacts --batch condor
##python Plot_Syst_Impacts_EFT.py --freeze -P ctz ctw cpq3 cpqm cpt -n 5DObs -d ./EFTWorkspace_5D.root -m impacts --batch condor

##python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctz -n 5DExp -d ./EFTWorkspace_5D.root -m impacts --collect
