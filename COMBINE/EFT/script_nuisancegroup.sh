#NB: takes ~45min per operator !

# //--------------------------------------------
# ctz

if [[ $1 == "" || $1 == "ctz" ]] ; then
    rm log_ctz.out

    FILE=./EFTWorkspace_ctz_group.root
    if [ -f "$FILE" ]; then
        echo "$FILE already exists."
    else
        echo "$FILE does not exist yet, will create it."
        python DumpEFTParameterization.py ../templates/Templates_NN_ctz_EFT2_Run2.root
        python EFTFitter.py -d ../datacards/COMBINED_Datacard_TemplateFit_NN_ctz_Run2_splitBin.txt --onlyworkspace -o EFTWorkspace_ctz_group.root
    fi

    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups tzq_rate 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups ttz_rate 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups twz_rate 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups bkg_rate 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups jes 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups jer 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups met 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups btag 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups other_exp 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups fr 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups lep_eff 2>&1 | tee -a log_ctz.out
    python EFTFitter.py -d ./EFTWorkspace_ctz_group.root -P ctz --freeze -m bestfit --freezenuisancegroups theory 2>&1 | tee -a log_ctz.out
fi

# //--------------------------------------------
# ctw

if [[ $1 == "" || $1 == "ctw" ]] ; then
    rm log_ctw.out

    FILE=./EFTWorkspace_ctw_group.root
    if [ -f "$FILE" ]; then
        echo "$FILE already exists."
    else
        echo "$FILE does not exist yet, will create it."
        python DumpEFTParameterization.py ../templates/Templates_NN_ctw_EFT2_Run2.root
        python EFTFitter.py -d ../datacards/COMBINED_Datacard_TemplateFit_NN_ctw_Run2_splitBin.txt --onlyworkspace -o EFTWorkspace_ctw_group.root
    fi
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups tzq_rate 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups ttz_rate 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups twz_rate 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups bkg_rate 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups jes 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups jer 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups met 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups btag 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups other_exp 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups fr 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups lep_eff 2>&1 | tee -a log_ctw.out
    python EFTFitter.py -d ./EFTWorkspace_ctw_group.root -P ctw --freeze -m bestfit --freezenuisancegroups theory 2>&1 | tee -a log_ctw.out
fi

# //--------------------------------------------
# cpq3

if [[ $1 == "" || $1 == "cpq3" ]] ; then
    rm log_cpq3.out

    FILE=./EFTWorkspace_cpq3_group.root
    if [ -f "$FILE" ]; then
        echo "$FILE already exists."
    else
        echo "$FILE does not exist yet, will create it."
        python DumpEFTParameterization.py ../templates/Templates_NN_cpq3_EFT2_Run2.root
        python EFTFitter.py -d ../datacards/COMBINED_Datacard_TemplateFit_NN_cpq3_Run2_splitBin.txt --onlyworkspace -o EFTWorkspace_cpq3_group.root
    fi

    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups tzq_rate 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups ttz_rate 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups twz_rate 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups bkg_rate 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups jes 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups jer 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups met 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups btag 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups other_exp 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups fr 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups lep_eff 2>&1 | tee -a log_cpq3.out
    python EFTFitter.py -d ./EFTWorkspace_cpq3_group.root -P cpq3 --freeze -m bestfit --freezenuisancegroups theory 2>&1 | tee -a log_cpq3.out
fi

# //--------------------------------------------
# cpqm

if [[ $1 == "" || $1 == "cpqm" ]] ; then
    rm log_cpqm.out

    FILE=./EFTWorkspace_SM_group.root
    if [ -f "$FILE" ]; then
        echo "$FILE already exists."
    else
        echo "$FILE does not exist yet, will create it."
        python DumpEFTParameterization.py ../templates/Templates_NN_SM_EFT2_Run2.root
        python EFTFitter.py -d ../datacards/COMBINED_Datacard_TemplateFit_NN_SM_Run2_splitBin.txt --onlyworkspace -o EFTWorkspace_SM_group.root
    fi

    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups tzq_rate 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups ttz_rate 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups twz_rate 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups bkg_rate 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups jes 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups jer 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups met 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups btag 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups other_exp 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups fr 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups lep_eff 2>&1 | tee -a log_cpqm.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpqm --freeze -m bestfit --freezenuisancegroups theory 2>&1 | tee -a log_cpqm.out
fi

# //--------------------------------------------
# cpt

if [[ $1 == "" || $1 == "cpt" ]] ; then
    rm log_cpt.out

    FILE=./EFTWorkspace_SM_group.root
    if [ -f "$FILE" ]; then
        echo "$FILE already exists."
    else
        echo "$FILE does not exist yet, will create it."
        python DumpEFTParameterization.py ../templates/Templates_NN_SM_EFT2_Run2.root
        python EFTFitter.py -d ../datacards/COMBINED_Datacard_TemplateFit_NN_SM_Run2_splitBin.txt --onlyworkspace -o EFTWorkspace_SM_group.root
    fi

    #-- Bestfit
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups tzq_rate 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups ttz_rate 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups twz_rate 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups bkg_rate 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups jes 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups jer 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups met 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups btag 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups other_exp 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups fr 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups lep_eff 2>&1 | tee -a log_cpt.out
    # python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m bestfit --freezenuisancegroups theory 2>&1 | tee -a log_cpt.out

    #-- Scans (more reliable because of double-minimum)
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_tzq_rate --freezenuisancegroups tzq_rate 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_ttz_rate --freezenuisancegroups ttz_rate 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_twz_rate --freezenuisancegroups twz_rate 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_bkg_rate --freezenuisancegroups bkg_rate 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_jes --freezenuisancegroups jes 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_jer --freezenuisancegroups jer 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_met --freezenuisancegroups met 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_btag --freezenuisancegroups btag 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_other_exp --freezenuisancegroups other_exp 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_fr --freezenuisancegroups fr 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_lep_eff --freezenuisancegroups lep_eff 2>&1 | tee -a log_cpt.out
    python EFTFitter.py -d ./EFTWorkspace_SM_group.root -P cpt --freeze -m scan -n cptgroup_theory --freezenuisancegroups theory 2>&1 | tee -a log_cpt.out
    python EFTPlotter.py -n cptgroup_tzq_rate -P cpt
    python EFTPlotter.py -n cptgroup_ttz_rate -P cpt
    python EFTPlotter.py -n cptgroup_twz_rate -P cpt
    python EFTPlotter.py -n cptgroup_bkg_rate -P cpt
    python EFTPlotter.py -n cptgroup_jes -P cpt
    python EFTPlotter.py -n cptgroup_jer -P cpt
    python EFTPlotter.py -n cptgroup_met -P cpt
    python EFTPlotter.py -n cptgroup_btag -P cpt
    python EFTPlotter.py -n cptgroup_other_exp -P cpt
    python EFTPlotter.py -n cptgroup_fr -P cpt
    python EFTPlotter.py -n cptgroup_lep_eff -P cpt
    python EFTPlotter.py -n cptgroup_theory -P cpt
fi

# //--------------------------------------------
