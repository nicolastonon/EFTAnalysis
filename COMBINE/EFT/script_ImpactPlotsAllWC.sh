#-- CREATE WORKSPACE
#python EFTFitter.py -d ../datacards/COMBINED_Datacard_TemplateFit_Run2.txt --onlyworkspace

#-- PRODUCE IMPACTS

#python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctz -n ctzexp -d ./EFTWorkspace.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctw -n ctwexp -d ./EFTWorkspace.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpqm -n cpqmexp -d ./EFTWorkspace.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpq3 -n cpq3exp -d ./EFTWorkspace.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpt -n cptexp -d ./EFTWorkspace.root -m impacts --batch condor

#python Plot_Syst_Impacts_EFT.py --freeze -P ctz -n ctzobs -d ./EFTWorkspace.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --freeze -P ctw -n ctwobs -d ./EFTWorkspace.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --freeze -P cpqm -n cpqmobs -d ./EFTWorkspace.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --freeze -P cpq3 -n cpq3obs -d ./EFTWorkspace.root -m impacts --batch condor
#python Plot_Syst_Impacts_EFT.py --freeze -P cpt -n cptobs -d ./EFTWorkspace.root -m impacts --batch condor

#-- COLLECT OUTPUTS

#python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctz -n ctzexp -d ./EFTWorkspace.root -m impacts --batch condor --collect
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P ctw -n ctwexp -d ./EFTWorkspace.root -m impacts --batch condor --collect
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpqm -n cpqmexp -d ./EFTWorkspace.root -m impacts --batch condor --collect
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpq3 -n cpq3exp -d ./EFTWorkspace.root -m impacts --batch condor --collect
#python Plot_Syst_Impacts_EFT.py --exp --freeze -P cpt -n cptexp -d ./EFTWorkspace.root -m impacts --batch condor --collect

#python Plot_Syst_Impacts_EFT.py --freeze -P ctz -n ctzobs -d ./EFTWorkspace.root -m impacts --batch condor --collect
#python Plot_Syst_Impacts_EFT.py --freeze -P ctw -n ctwobs -d ./EFTWorkspace.root -m impacts --batch condor --collect
#python Plot_Syst_Impacts_EFT.py --freeze -P cpqm -n cpqmobs -d ./EFTWorkspace.root -m impacts --batch condor --collect
#python Plot_Syst_Impacts_EFT.py --freeze -P cpq3 -n cpq3obs -d ./EFTWorkspace.root -m impacts --batch condor --collect
#python Plot_Syst_Impacts_EFT.py --freeze -P cpt -n cptobs -d ./EFTWorkspace.root -m impacts --batch condor --collect
