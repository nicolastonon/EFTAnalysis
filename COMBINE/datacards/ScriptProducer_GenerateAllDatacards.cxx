/* BASH CUSTOM */
#define RST   "\e[0m"
#define KRED  "\e[31m"
#define KGRN  "\e[32m"
#define KYEL  "\e[33m"
#define KBLU  "\e[34m"
#define KMAG  "\e[35m"
#define KCYN  "\e[36m"
#define KWHT  "\e[37m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST
#define BOLD(x) "\e[1m" x RST
#define ITAL(x) "\e[3m" x RST
#define UNDL(x) "\e[4m" x RST
#define STRIKE(x) "\e[9m" x RST
#define DIM(x) "\e[2m" x RST
#define DOUBLEUNDERLINE(x) "\e[21m" x RST
#define CURLYUNDERLINE(x) "\e[4:3m" x RST
#define BLINK(x) "\e[5m" x RST
#define REVERSE(x) "\e[7m" x RST
#define INVISIBLE(x) "\e[8m" x RST
#define OVERLINE(x) "\e[53m" x RST
#define TURQUOISE  "\e[38;5;42m"
#define SALMON  "\e[38;2;240;143;104m"
#define FTURQUOISE(x) TURQUOISE x RST
#define FSALMON(x) SALMON x RST
#define YELBKG(x) "\e[43m" x RST
#define CYANBKG(x) "\e[46m" x RST

#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include <sys/stat.h> // to be able to check file existence
#include <sstream>
#include <iomanip>

using namespace std;


//--------------------------------------------
// ##     ## ######## ##       ########  ######## ########
// ##     ## ##       ##       ##     ## ##       ##     ##
// ##     ## ##       ##       ##     ## ##       ##     ##
// ######### ######   ##       ########  ######   ########
// ##     ## ##       ##       ##        ##       ##   ##
// ##     ## ##       ##       ##        ##       ##    ##
// ##     ## ######## ######## ##        ######## ##     ##
//--------------------------------------------

//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Existence(const TString& name)
{
    struct stat buffer;
    bool found = (stat (name.Data(), &buffer) == 0); //true if file exists
    return found;
}


//Convert a double into a TString
// precision --> can choose if TString how many digits the TString should display
TString Convert_Number_To_TString(double number, int precision=3)
{
	stringstream ss;
	ss << std::setprecision(precision) << number;
	TString ts = ss.str();
	return ts;
}


//--------------------------------------------
//  ######  ########  ########    ###    ######## ########
// ##    ## ##     ## ##         ## ##      ##    ##
// ##       ##     ## ##        ##   ##     ##    ##
// ##       ########  ######   ##     ##    ##    ######
// ##       ##   ##   ##       #########    ##    ##
// ##    ## ##    ##  ##       ##     ##    ##    ##
//  ######  ##     ## ######## ##     ##    ##    ########

//  ######   ######  ########  #### ########  ########    ########  #######  ########
// ##    ## ##    ## ##     ##  ##  ##     ##    ##       ##       ##     ## ##     ##
// ##       ##       ##     ##  ##  ##     ##    ##       ##       ##     ## ##     ##
//  ######  ##       ########   ##  ########     ##       ######   ##     ## ########
//       ## ##       ##   ##    ##  ##           ##       ##       ##     ## ##   ##
// ##    ## ##    ## ##    ##   ##  ##           ##       ##       ##     ## ##    ##
//  ######   ######  ##     ## #### ##           ##       ##        #######  ##     ##

// ######## ######## ##     ## ########  ##          ###    ######## ########    ######## #### ########
//    ##    ##       ###   ### ##     ## ##         ## ##      ##    ##          ##        ##     ##
//    ##    ##       #### #### ##     ## ##        ##   ##     ##    ##          ##        ##     ##
//    ##    ######   ## ### ## ########  ##       ##     ##    ##    ######      ######    ##     ##
//    ##    ##       ##     ## ##        ##       #########    ##    ##          ##        ##     ##
//    ##    ##       ##     ## ##        ##       ##     ##    ##    ##          ##        ##     ##
//    ##    ######## ##     ## ##        ######## ##     ##    ##    ########    ##       ####    ##
//--------------------------------------------

/**
 * Produce script containing the commands to produce the datacards (single and combination) automatically
 */
void Script_Datacards_TemplateFit(char include_systematics, char include_statistical, TString template_name, TString region, vector<TString> v_templates, vector<TString> v_channel, vector<TString> v_regions, TString lumiName, int mode_histoBins, bool scan_operator_hardcoded, char use_SM_setup, TString selection, TString filename_template_suffix, char include_otherRegions, bool use_rph, TString path_tmp_workspace, TString scriptname, bool hardcoded_spectator_variables=false)
{
//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

    //-- Check if use shape syst or not
	TString systChoice;
	if(include_systematics == 'y') {systChoice = "withShape";}
	else if(include_systematics == 'n') {systChoice = "noShape";}
    else {cout<<"Wrong arguments ! Abort !"<<endl; return;}

    //-- Check if use stat. uncert. or not
	TString statChoice;
	if(include_statistical == 'y') {statChoice = "withStat";}
	else if(include_statistical == 'n') {statChoice = "noStat";}
    else {cout<<"Wrong arguments ! Abort !"<<endl; return;}

    //-- If specific arguments were chosen at command line, modify the vectors defined in the main
    if(template_name != "0") {v_templates.resize(0); v_templates.push_back(template_name);}
    else if(!v_templates.size() && include_otherRegions != '1') {cout<<"Template name not set ! Abort !"<<endl; return;}

    //-- Trick: in case there is ==1 template and >=1 regions defined, we want to use this template for all the regions; but to allow to then add additional pairs of regions/templates, expand the current lists so that they match 1 to 1
    if(v_templates.size() == 1 && v_regions.size() >= 1)
    {
        for(int iregion=1; iregion<v_regions.size(); iregion++)
        {
            v_templates.push_back(v_templates[0]); //Duplicate first element for each additional region
        }
    }

    //-- Special case: for cpq3, we use NN_cpq3 in SRtZq *but* different distribution in SRttZ --> Hardcode it here //Obsolete
    /*
    for(int itemplate=0; itemplate<v_templates.size(); itemplate++)
    {
        for(int iregion=0; iregion<v_regions.size(); iregion++)
        {
            if(v_templates[itemplate] == "NN_cpq3" && v_regions[iregion] == "SRttZ") {v_templates[iregion] = "NN_5D";}
        }
    }*/

    //-- HARDCODED
    vector<TString> v_WCs_operator_scan1 = {"-999"}; //HARDCODED
    TString operator_scan1 = "";
    if(scan_operator_hardcoded)
    {
        mode_histoBins = 1; //Scan on parametrized NN --> must treat each histogram bin separately
        filename_template_suffix+= "param"; //Default identifier

        operator_scan1 = "ctw";
        //v_WCs_operator_scan1 = {"-4","-2","-1","0","1","2","4"};
        // v_WCs_operator_scan1 = {"-3", "-2", "-1.5", "-1", "-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6", "0.8", "1", "1.5", "2", "3"};
        v_WCs_operator_scan1 = {"-1.5", "-1", "-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6", "0.8", "1", "1.5"};

        v_regions.resize(0); v_templates.resize(0);
        v_regions.push_back("SRtZq"); v_templates.push_back(v_templates[0]);
        v_regions.push_back("SRttZ"); v_templates.push_back(v_templates[0]);
        v_regions.push_back("SRother"); v_templates.push_back("mTW");
    }

    //-- Add predefined additional regions (NB: corresponding templates must have been produced and merged into main template file !)
    else if(include_otherRegions == 'y' || include_otherRegions == '1')
    {
        v_regions.push_back("SRttZ4l"); v_templates.push_back("countExp");
        v_regions.push_back("CRWZ"); v_templates.push_back("countExp");
        v_regions.push_back("CRZZ"); v_templates.push_back("countExp");
        // v_regions.push_back("CRWZ"); v_templates.push_back("mTW"); //Obsolete

        //v_regions.push_back("CRDY"); v_templates.push_back("countExp"); //-- remove ?
    }

    //-- Overrides some option to perform a fit following the setup of the main SM tZq differential analysis
    //NB: just add all region and template names together here; then call dedicated hard-coded function to sort out relevant combinations
    else if(use_SM_setup == 'y')
    {
        v_regions.resize(0); v_templates.resize(0);
        v_regions.push_back("tZq"); v_templates.push_back("NN");
        v_regions.push_back("ttZ"); v_templates.push_back("NN");
        v_regions.push_back("wz"); v_templates.push_back("mTW");
        v_regions.push_back("zz"); v_templates.push_back("countExp");
        v_regions.push_back("Vg"); v_templates.push_back("channel");
    }

    if(hardcoded_spectator_variables)
    {
        v_regions.clear(); v_templates.clear();
        v_templates.push_back("njets");
        v_templates.push_back("nbjets");
        v_templates.push_back("jPrimeAbsEta");
        v_templates.push_back("jprime_Pt");
        v_templates.push_back("metEt");
        v_templates.push_back("mTW");
        v_templates.push_back("lAsymmetry");
        v_templates.push_back("recoZ_Pt");
        v_templates.push_back("recoZ_dPhill");
        v_templates.push_back("dR_tZ");
        v_templates.push_back("maxDeepJet");

        for(int ivar=0; ivar<v_templates.size(); ivar++)
        {
            v_regions.push_back("");
        }
    }
    // else if(template_name == "1")
    // {
    //     v_regions.resize(0); v_templates.resize(0);
    //     v_regions.push_back("SR"); v_templates.push_back("njets");
    //     v_regions.push_back("SR"); v_templates.push_back("nbjets");
    //     v_regions.push_back("SR"); v_templates.push_back("jPrimeAbsEta");
    //     v_regions.push_back("SR"); v_templates.push_back("jprime_Pt");
    //     v_regions.push_back("SR"); v_templates.push_back("metEt");
    //     v_regions.push_back("SR"); v_templates.push_back("mTW");
    //     v_regions.push_back("SR"); v_templates.push_back("lAsymmetry");
    //     v_regions.push_back("SR"); v_templates.push_back("recoZ_Pt");
    //     v_regions.push_back("SR"); v_templates.push_back("recoZ_dPhill");
    //     v_regions.push_back("SR"); v_templates.push_back("dR_tZ");
    // }

    vector<TString> v_lumiYears;
    if(lumiName == "2016") {v_lumiYears.push_back("2016");}
    else if(lumiName == "2017") {v_lumiYears.push_back("2017");}
    else if(lumiName == "2018") {v_lumiYears.push_back("2018");}
    else if(lumiName == "201617") {v_lumiYears.push_back("2016"); v_lumiYears.push_back("2017");}
    else if(lumiName == "201618") {v_lumiYears.push_back("2016"); v_lumiYears.push_back("2018");}
    else if(lumiName == "201718") {v_lumiYears.push_back("2017"); v_lumiYears.push_back("2018");}
    else if(lumiName == "Run2") {v_lumiYears.push_back("2016"); v_lumiYears.push_back("2017"); v_lumiYears.push_back("2018");}

	ofstream file_out; //output script

    if(use_rph) {file_out.open(scriptname.Data(), std::ios_base::app);} //Append
    else {file_out.open(scriptname.Data(), std::ofstream::out);} //Overwrite

	TString dir = "./datacards_TemplateFit/";

    int nbins = -1; //If 'mode_histoBins=1', will read histogram binning in rootfile and create 1 datacard per histo bin
    vector<vector<vector<vector<int>>>> v_nbins_EFT_year_region_template(v_WCs_operator_scan1.size()); //Binnings can differ for each template --> keep track of all different binnins per year/region/template/...
    for(int i=0; i<v_WCs_operator_scan1.size(); i++) //Reserve memory
    {
        v_nbins_EFT_year_region_template[i].resize(v_lumiYears.size());
        for(int j=0; j<v_lumiYears.size(); j++)
        {
            v_nbins_EFT_year_region_template[i][j].resize(v_regions.size());
            for(int k=0; k<v_regions.size(); k++)
            {
                v_nbins_EFT_year_region_template[i][j][k].resize(v_templates.size());
                for(int l=0; l<v_templates.size(); l++)
                {
                    v_nbins_EFT_year_region_template[i][j][k][l] = 1;
                }
            }
        }
    }

    file_out <<"mkdir "<<dir<<endl<<endl;
    for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
    {
        file_out <<"mkdir "<<dir+v_lumiYears[iyear]<<endl<<endl;
    }


//  ######  ########  ########    ###    ######## ########
// ##    ## ##     ## ##         ## ##      ##    ##
// ##       ##     ## ##        ##   ##     ##    ##
// ##       ########  ######   ##     ##    ##    ######
// ##       ##   ##   ##       #########    ##    ##
// ##    ## ##    ##  ##       ##     ##    ##    ##
//  ######  ##     ## ######## ##     ##    ##    ########

// # #    # #####  # #    # # #####  #    #   ##   #
// # ##   # #    # # #    # # #    # #    #  #  #  #
// # # #  # #    # # #    # # #    # #    # #    # #
// # #  # # #    # # #    # # #    # #    # ###### #
// # #   ## #    # #  #  #  # #    # #    # #    # #
// # #    # #####  #   ##   # #####   ####  #    # ######

// #####    ##   #####   ##    ####    ##   #####  #####   ####
// #    #  #  #    #    #  #  #    #  #  #  #    # #    # #
// #    # #    #   #   #    # #      #    # #    # #    #  ####
// #    # ######   #   ###### #      ###### #####  #    #      #
// #    # #    #   #   #    # #    # #    # #   #  #    # #    #
// #####  #    #   #   #    #  ####  #    # #    # #####   ####

//--------------------------------------------
//--- First loop over years/regions/templates
//===> Commands to produce individual datacards

    for(int ipt_EFT=0; ipt_EFT<v_WCs_operator_scan1.size(); ipt_EFT++) //1 datacard per EFT scan point
    {
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
            for(int iregion=0; iregion<v_regions.size(); iregion++) //Loop over regions
            {
                // cout<<"v_regions[iregion] "<<v_regions[iregion]<<endl;

                for(int itemplate=0; itemplate<v_templates.size(); itemplate++) //Loop over templates
                {
                    if(itemplate != iregion) {continue;} //Associate templates/regions 1 to 1 for more flexibility

                    // cout<<"v_templates[itemplate] "<<v_templates[itemplate]<<endl;
                    bool isOtherRegion = false;
                    if((include_otherRegions=='y' || include_otherRegions=='1') && (v_regions[iregion].Contains("SRttZ4l") || v_regions[iregion].Contains("CRWZ") || v_regions[iregion].Contains("CRZZ") || v_regions[iregion].Contains("CRDY")) ) {isOtherRegion = true;}

                    TString var = v_templates[itemplate] + (v_regions[iregion]!=""? "_"+v_regions[iregion]:"");
                    if(scan_operator_hardcoded && v_templates[itemplate].Contains("NN")) {var+= "_" + operator_scan1 + "_" + v_WCs_operator_scan1[ipt_EFT];} //NN scan: adapt template name
                    else if(v_regions[iregion] == "SRother" && filename_template_suffix.Contains("EFT")) {var = "mTW_SRother";} //Use mTW for now in CR //Hard-coded

                    //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                    var.ReplaceAll('-', 'm');

                    TString file_histos_pathFromHere = "./../templates/Templates_"+ v_templates[0] + (filename_template_suffix != ""? "_"+filename_template_suffix:"") + (selection != ""? "_"+selection:"") + "_"; //For use within this code
                    if(isOtherRegion) {file_histos_pathFromHere = "./../templates/Templates_otherRegions"+(selection != ""? "_"+selection:"")+"_";} //Read a different file for templates with 'fixed' observables (only change SR templates)
                    else if(v_templates[itemplate] == "NN_SM") {file_histos_pathFromHere = "./../templates/Templates_NN_SM" + (filename_template_suffix != ""? "_"+filename_template_suffix:"") + (selection != ""? "_"+selection:"") + "_";} //Trick: when reading NN_cpq3 file (for SRtZq), we need to read the NN_SM file for SRttZ (<-> NN_SM template) and SRother (<-> mTW template) !
			        else if(v_templates[itemplate] == "NN_5D") {file_histos_pathFromHere = "./../templates/Templates_NN_5D" + (filename_template_suffix != ""? "_"+filename_template_suffix:"") + (selection != ""? "_"+selection:"") + "_";} //Idem for NN_5D
                    TString file_histos_pathFromHere_Run2 = file_histos_pathFromHere + "Run2.root"; //In case year-dependent file is not found, will look for Run2 file by default
                    file_histos_pathFromHere+= lumiName+".root";
                    if(hardcoded_spectator_variables) {file_histos_pathFromHere = "./../templates/ControlPlotsPaper.root"; file_histos_pathFromHere_Run2 = file_histos_pathFromHere;} //Hardcoded filename (for control plots)
                    // if(hardcoded_spectator_variables) {file_histos_pathFromHere = "./../templates/ControlHistograms_signal_";} //Hardcoded path for now
                    // if(scan_operator_hardcoded) {file_histos_pathFromHere = "./../templates/Templates_NN_EFT2param_Run2.root";} //HARD-CODED

                    cout<<DIM("Trying to open input file "<<file_histos_pathFromHere<<" ... ");
                    if(Check_File_Existence(file_histos_pathFromHere)) {cout<<DIM("FOUND !")<<endl;}
            		else
                    {
                        cout<<endl<<DIM("Trying to open input file "<<file_histos_pathFromHere_Run2<<" ... "); //Try Run2 file
                        if(Check_File_Existence(file_histos_pathFromHere_Run2)) {cout<<DIM("FOUND !")<<endl;}
                        else {cout<<BOLD(FRED("ERROR: input template file not found ! Abort !"))<<endl;}
                        file_histos_pathFromHere = file_histos_pathFromHere_Run2; //Update main filepath
                    }

                    TString file_histos = "../." + file_histos_pathFromHere; //Path to write into datacard
        			cout<<endl<<FMAG("---> Will use filepath : ")<<file_histos<<endl<<endl;

                    if((mode_histoBins==1 && !isOtherRegion) || use_rph) //Need to infer the number of bins of the considered histograms, in order to create 1 card per bin
                    // if(true)
                    {
                        TFile* f_tmp = TFile::Open(file_histos_pathFromHere);
                        TString hname_tmp = var + "_" + v_lumiYears[iyear ] + "__data_obs"; //Hard-coded: look for data histo to infer binning
						// cout<<"Reading histogram: "<<hname_tmp<<endl;

                        if(!f_tmp->GetListOfKeys()->Contains(hname_tmp) ) {cout<<BOLD(FRED("ERROR: histogram "<<hname_tmp<<" not found in file "<<file_histos_pathFromHere<<" ! Can not infer histogram binning !"))<<endl; return;}
                        // cout<<"hname_tmp "<<hname_tmp<<endl;
                        TH1F* h_tmp = (TH1F*) f_tmp->Get(hname_tmp);
                        nbins = h_tmp->GetNbinsX();
						// cout<<"--> nbins = "<<nbins<<endl;
                        v_nbins_EFT_year_region_template[ipt_EFT][iyear][iregion][itemplate] = nbins;
                        delete h_tmp; h_tmp = NULL; f_tmp->Close();
                        if(nbins == -1) {cout<<BOLD(FRED("ERROR: histogram "<<hname_tmp<<" not found ! Can not infer histogram binning !"))<<endl; return;}
                    }
                    else {nbins = 1;} //Default: run next loop only once

                    if(nbins<0) {nbins=1;} //Need to loop at least once by default
        			for(int ilepchan=0; ilepchan<v_channel.size(); ilepchan++)
        			{
                        for(int ibin=1; ibin<nbins+1; ibin++)
                        {
                            if(mode_histoBins != 1 && ibin>1) {break;}

                            TString var_tmp = var;
                            if(mode_histoBins==1 && nbins > 1 && !isOtherRegion) {var_tmp = (TString) "bin" + Form("%d",ibin) + "_" + var;} //Also include bin number in naming scheme (--> will read single bin histos instead of full histos)
                            else if(mode_histoBins==3 && !isOtherRegion) {var_tmp = "countExp_" + var;} //Counting experiments

                            TString template_path = "Template_Datacard.txt";
                            if(use_rph) {template_path = "Template_Datacard_RPH.txt";}

            				file_out<<"python Parser_Datacard_Template.py "
                            + var_tmp + " "
                            + v_channel[ilepchan] + " "
                            + v_lumiYears[iyear] + " "
            				+ file_histos + " "
            				+ systChoice + " "
            				+ statChoice + " "
                            + dir + v_lumiYears[iyear] + " "
                            + path_tmp_workspace + " "
                            + template_path + " "
                            <<nbins<<" "
            				<<endl;
                        } //bin loop
        			} //channel loop

        			file_out<<endl<<endl;
            	} //template loop
            } //region loop
        } //year loop
    } //EFT loop

    // for(int i=0; i<v_nbins.size(); i++) {cout<<"v_nbins[i] "<<v_nbins[i]<<endl;}


//  ####   ####  #    # #####  # #    # ###### #####
// #    # #    # ##  ## #    # # ##   # #      #    #
// #      #    # # ## # #####  # # #  # #####  #    #
// #      #    # #    # #    # # #  # # #      #    #
// #    # #    # #    # #    # # #   ## #      #    #
//  ####   ####  #    # #####  # #    # ###### #####

// #####    ##   #####   ##    ####    ##   #####  #####
// #    #  #  #    #    #  #  #    #  #  #  #    # #    #
// #    # #    #   #   #    # #      #    # #    # #    #
// #    # ######   #   ###### #      ###### #####  #    #
// #    # #    #   #   #    # #    # #    # #   #  #    #
// #####  #    #   #   #    #  ####  #    # #    # #####

//--------------------------------------------
//--- Second loop over years/regions/templates
//===> Give all the single datacards as arguments to the script, for combination

    if(!hardcoded_spectator_variables) //Not needed for control plots
    {
        // int idx_v_nbins = 0; //Vector storing template binning info needs to be read in same order (same loops) as before; increment some index to keep track of the current element to read
        for(int ipt_EFT=0; ipt_EFT<v_WCs_operator_scan1.size(); ipt_EFT++)
        {
            //Command to execute the script which combines the datacards
            file_out<<"python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py ";

            for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
            {
                for(int iregion=0; iregion<v_regions.size(); iregion++) //Loop over regions
                {
                    for(int itemplate=0; itemplate<v_templates.size(); itemplate++) //Loop over templates
                    {
                        if(itemplate != iregion) {continue;} //Associate templates/regions 1 to 1 for more flexibility

                        bool isOtherRegion = false;
                        if((include_otherRegions=='y' || include_otherRegions=='1') && (v_regions[iregion].Contains("SRttZ4l") || v_regions[iregion].Contains("CRWZ") || v_regions[iregion].Contains("CRZZ") || v_regions[iregion].Contains("CRDY")) ) {isOtherRegion = true;}

                        TString var = v_templates[itemplate] + (v_regions[iregion]!=""? "_"+v_regions[iregion]:"");
                        if(scan_operator_hardcoded && v_templates[itemplate].Contains("NN")) {var+= "_" + operator_scan1 + "_" + v_WCs_operator_scan1[ipt_EFT];}
                        else if(v_regions[iregion] == "SRother" && filename_template_suffix.Contains("EFT")) {var = "mTW_SRother";} //Use mTW for now in CR //Hard-coded

                        //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                        var.ReplaceAll('-', 'm');

                        for(int ilepchan=0; ilepchan<v_channel.size(); ilepchan++) //Loop over channels
                        {
                            int nbins_tmp = 1;
                            if(mode_histoBins==1 && !isOtherRegion)
                            {
                                // nbins_tmp = v_nbins[idx_v_nbins]; idx_v_nbins++; //Read current binning; increment index to stay in sync
                                nbins_tmp = v_nbins_EFT_year_region_template[ipt_EFT][iyear][iregion][itemplate];
                                // cout<<"nbins_tmp "<<nbins_tmp<<endl;
                            }
                            for(int ibin=1; ibin<nbins_tmp+1; ibin++)
                            {
                                // if(v_templates[itemplate]=="categ" and ibin==6) {continue;} //HARDCODED TMP FIX (empty bin)

                                TString var_tmp = var;
                                if(mode_histoBins==1 && nbins_tmp > 1 && !isOtherRegion) {var_tmp = (TString) "bin" + Form("%d",ibin) + "_" + var;} //Also include bin number in naming scheme (--> will read single bin histos instead of full histos)
                                else if(mode_histoBins==3 && !isOtherRegion) {var_tmp = "countExp_" + var;}

                                file_out<<var_tmp;
                                if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
                                file_out<<"_"+v_lumiYears[iyear];
                                file_out<<"=" + dir + v_lumiYears[iyear] + "/"
                                + "datacard_"+var_tmp;
                                if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
                                if(use_rph) {file_out<<"_rph";}
                                file_out<<".txt ";
                            } //bin loop
                        } //channel loop

                    } //template loop
                } //region loop
            } //years loop

        	TString output_name = "COMBINED_Datacard_TemplateFit";
            if(systChoice == "noShape") output_name+= "_noShape";
            if(statChoice == "noStat") output_name+= "_noStat";
            if(mode_histoBins==3) output_name+= "_countExp";
            output_name+= "_" + v_templates[0]; //-- Make name template-specific //CHANGED
            if(scan_operator_hardcoded) {output_name+= "_" + operator_scan1 + "_" + v_WCs_operator_scan1[ipt_EFT];}
            output_name+= "_" + lumiName;
            if(use_rph) {output_name+= "_rph";}
            if(mode_histoBins==1) {output_name+= "_splitBin";}
            output_name+= ".txt";

        	file_out<<"> "<<output_name<<endl<<endl;

            file_out<<"echo \"\""<<endl;
            file_out<<"echo \"\e[33m===> Created combined datacard: \e[0m"<<output_name<<"\""<<endl<<endl<<endl;
        } //EFT loop
    } //If condition


//  ######   #######  ##     ## ########  #### ##    ## ######## ########
// ##    ## ##     ## ###   ### ##     ##  ##  ###   ## ##       ##     ##
// ##       ##     ## #### #### ##     ##  ##  ####  ## ##       ##     ##
// ##       ##     ## ## ### ## ########   ##  ## ## ## ######   ##     ##
// ##       ##     ## ##     ## ##     ##  ##  ##  #### ##       ##     ##
// ##    ## ##     ## ##     ## ##     ##  ##  ##   ### ##       ##     ##
//  ######   #######  ##     ## ########  #### ##    ## ######## ########

// ########     ###    ########    ###     ######     ###    ########  ########   ######
// ##     ##   ## ##      ##      ## ##   ##    ##   ## ##   ##     ## ##     ## ##    ##
// ##     ##  ##   ##     ##     ##   ##  ##        ##   ##  ##     ## ##     ## ##
// ##     ## ##     ##    ##    ##     ## ##       ##     ## ########  ##     ##  ######
// ##     ## #########    ##    ######### ##       ######### ##   ##   ##     ##       ##
// ##     ## ##     ##    ##    ##     ## ##    ## ##     ## ##    ##  ##     ## ##    ##
// ########  ##     ##    ##    ##     ##  ######  ##     ## ##     ## ########   ######

// ########  ######## ########     ########  ########  ######   ####  #######  ##    ##
// ##     ## ##       ##     ##    ##     ## ##       ##    ##   ##  ##     ## ###   ##
// ##     ## ##       ##     ##    ##     ## ##       ##         ##  ##     ## ####  ##
// ########  ######   ########     ########  ######   ##   ####  ##  ##     ## ## ## ##
// ##        ##       ##   ##      ##   ##   ##       ##    ##   ##  ##     ## ##  ####
// ##        ##       ##    ##     ##    ##  ##       ##    ##   ##  ##     ## ##   ###
// ##        ######## ##     ##    ##     ## ########  ######   ####  #######  ##    ##

//-- Also create datacards for individual regions (cross-checks)
//-- NB: datacards split by years may already be created by selecting a specific year via command-line
//-- NB: these cards are stored in 'datacards_TemplateFit', but must be run from 'datacards' from relative paths to work (--> need to copy card first)
    bool create_datacards_indivRegions = true;
    if(create_datacards_indivRegions && !scan_operator_hardcoded)
    {
        for(int iregion=0; iregion<v_regions.size(); iregion++) //Loop over regions
        {
            //Command to execute the script which combines the datacards
            file_out<<"python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py ";

            for(int itemplate=0; itemplate<v_templates.size(); itemplate++) //Loop over templates
            {
                if(itemplate != iregion) {continue;} //-- Associate templates/regions 1 to 1 for more flexibility

                for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
                {
                    TString var = v_templates[itemplate] + (v_regions[iregion]!=""? "_"+v_regions[iregion]:"");
                    if(v_regions[iregion] == "SRother" && filename_template_suffix.Contains("EFT")) {var = "mTW_SRother";} //Use mTW for now in CR //Hard-coded

                    bool isOtherRegion = false;
                    if((include_otherRegions=='y' || include_otherRegions=='1') && (v_regions[iregion].Contains("SRttZ4l") || v_regions[iregion].Contains("CRWZ") || v_regions[iregion].Contains("CRZZ") || v_regions[iregion].Contains("CRDY")) ) {isOtherRegion = true;}

                    //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                    var.ReplaceAll('-', 'm');

                    for(int ilepchan=0; ilepchan<v_channel.size(); ilepchan++) //Loop over channels
                    {
                        int nbins_tmp = 1;
                        if(mode_histoBins==1 && !isOtherRegion)
                        {
                            // nbins_tmp = v_nbins[idx_v_nbins]; idx_v_nbins++; //Read current binning; increment index to stay in sync
                            nbins_tmp = v_nbins_EFT_year_region_template[0][iyear][iregion][itemplate];
                        }
                        for(int ibin=1; ibin<nbins_tmp+1; ibin++)
                        {
                            // if(v_templates[itemplate]=="categ" and ibin==6) {continue;} //HARDCODED TMP FIX (empty bin)

                            TString var_tmp = var;
                            if(mode_histoBins==1 && nbins_tmp > 1 && !isOtherRegion) {var_tmp = (TString) "bin" + Form("%d",ibin) + "_" + var;} //Also include bin number in naming scheme (--> will read single bin histos instead of full histos)
                            else if(mode_histoBins==3 && !isOtherRegion) {var_tmp = "countExp_" + var;}

                            file_out<<var_tmp;
                            if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
                            file_out<<"_"+v_lumiYears[iyear];
                            file_out<<"=" + dir + v_lumiYears[iyear] + "/"
                            + "datacard_"+var_tmp;
                            if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
                            file_out<<".txt ";
                        } //bin loop
                    } //channel loop
                } //years loop
            } //template loop

            TString output_name = "COMBINED_Datacard_TemplateFit";
            if(systChoice == "noShape") output_name+= "_noShape";
            if(statChoice == "noStat") output_name+= "_noStat";
            if(mode_histoBins==3) output_name+= "_countExp";
            output_name+= "_" + v_templates[iregion]; //-- Make region card name template-specific
            // if(hardcoded_spectator_variables) {}
            // output_name+= "_" + v_templates[0]; //-- Make region card name template-specific
            output_name+= "_" + v_regions[iregion] + "_" + lumiName + ".txt";
            file_out<<"> "<<output_name<<endl<<endl;
            file_out<<"mv "<<output_name<<" datacards_TemplateFit/"<<endl<<endl;
            cout<<"--> Generated per-region card: "<<output_name<<endl;
        } //region loop
    } //create_datacards_indivRegions


//  ######   #######  ##     ## ########  #### ##    ## ######## ########
// ##    ## ##     ## ###   ### ##     ##  ##  ###   ## ##       ##     ##
// ##       ##     ## #### #### ##     ##  ##  ####  ## ##       ##     ##
// ##       ##     ## ## ### ## ########   ##  ## ## ## ######   ##     ##
// ##       ##     ## ##     ## ##     ##  ##  ##  #### ##       ##     ##
// ##    ## ##     ## ##     ## ##     ##  ##  ##   ### ##       ##     ##
//  ######   #######  ##     ## ########  #### ##    ## ######## ########

// ########     ###    ########    ###     ######     ###    ########  ########   ######
// ##     ##   ## ##      ##      ## ##   ##    ##   ## ##   ##     ## ##     ## ##    ##
// ##     ##  ##   ##     ##     ##   ##  ##        ##   ##  ##     ## ##     ## ##
// ##     ## ##     ##    ##    ##     ## ##       ##     ## ########  ##     ##  ######
// ##     ## #########    ##    ######### ##       ######### ##   ##   ##     ##       ##
// ##     ## ##     ##    ##    ##     ## ##    ## ##     ## ##    ##  ##     ## ##    ##
// ########  ##     ##    ##    ##     ##  ######  ##     ## ##     ## ########   ######

// ########  ######## ########     ########  #### ##    ##
// ##     ## ##       ##     ##    ##     ##  ##  ###   ##
// ##     ## ##       ##     ##    ##     ##  ##  ####  ##
// ########  ######   ########     ########   ##  ## ## ##
// ##        ##       ##   ##      ##     ##  ##  ##  ####
// ##        ##       ##    ##     ##     ##  ##  ##   ###
// ##        ######## ##     ##    ########  #### ##    ##

    //-- To get correctr prefit/postfit uncertainties from combine, must split bins just like for nominal fits (so that we can access total errors easily) --> only combine all years together, for each channel (e.g. histo bin) included in the fit
    bool create_datacards_indivBins = true;
    if(create_datacards_indivBins && !scan_operator_hardcoded && mode_histoBins == 1)
    {
        for(int iregion=0; iregion<v_regions.size(); iregion++) //Loop over regions
        {
            for(int itemplate=0; itemplate<v_templates.size(); itemplate++) //Loop over templates
            {
                if(itemplate != iregion) {continue;} //NEW -- associate templates/regions 1 to 1 for more flexibility

                for(int ilepchan=0; ilepchan<v_channel.size(); ilepchan++) //Loop over channels
                {
                    TString var = v_templates[itemplate] + (v_regions[iregion]!=""? "_"+v_regions[iregion]:"");
                    if(v_regions[iregion] == "SRother" && filename_template_suffix.Contains("EFT")) {var = "mTW_SRother";} //Use mTW for now in CR //Hard-coded

                    bool isOtherRegion = false;
                    if((include_otherRegions=='y' || include_otherRegions=='1') && (v_regions[iregion].Contains("SRttZ4l") || v_regions[iregion].Contains("CRWZ") || v_regions[iregion].Contains("CRZZ") || v_regions[iregion].Contains("CRDY")) ) {isOtherRegion = true;}

                    //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                    var.ReplaceAll('-', 'm');

                    int nbins_tmp = 1;
                    if(!isOtherRegion)
                    {
                        // nbins_tmp = v_nbins[idx_v_nbins]; idx_v_nbins++; //Read current binning; increment index to stay in sync
                        nbins_tmp = v_nbins_EFT_year_region_template[0][0][iregion][itemplate];
                    }
                    for(int ibin=1; ibin<nbins_tmp+1; ibin++)
                    {
                        //Command to execute the script which combines the datacards
                        file_out<<"python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py ";

                        TString var_tmp = var;
                        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
                        {
                            if(mode_histoBins==1 && nbins_tmp > 1 && !isOtherRegion) {var_tmp = (TString) "bin" + Form("%d",ibin) + "_" + var;} //Also include bin number in naming scheme (--> will read single bin histos instead of full histos)
                            else if(mode_histoBins==3 && !isOtherRegion) {var_tmp = "countExp_" + var;}

                            file_out<<var_tmp;
                            if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
                            file_out<<"_"+v_lumiYears[iyear];
                            file_out<<"=" + dir + v_lumiYears[iyear] + "/"
                            + "datacard_"+var_tmp;
                            if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
                            file_out<<".txt ";
                        } //years loop
                        TString output_name = "datacard_" + var_tmp + ".txt";
                        file_out<<"> "<<output_name<<endl<<endl;
                        file_out<<"mv "<<output_name<<" datacards_TemplateFit/"<<endl<<endl;
                    } //bin loop
                } //channel loop
            } //template loop
        } //region loop
    } //create_datacards_indivRegions


	file_out.close();
}


void Run_Commands(TString scriptname)
{
    system(("chmod 755 "+scriptname).Data());

    cout<<FGRN("... Created script ./"<<scriptname<<" !")<<endl;

    cout<<endl<<endl<<FGRN("-- Executing script ./"<<scriptname<<" !")<<endl;

    system(("./"+scriptname).Data());
}


//--------------------------------------------
//  ######   #######  ##     ## ##     ##    ###    ##    ## ########
// ##    ## ##     ## ###   ### ###   ###   ## ##   ###   ## ##     ##
// ##       ##     ## #### #### #### ####  ##   ##  ####  ## ##     ##
// ##       ##     ## ## ### ## ## ### ## ##     ## ## ## ## ##     ##
// ##       ##     ## ##     ## ##     ## ######### ##  #### ##     ##
// ##    ## ##     ## ##     ## ##     ## ##     ## ##   ### ##     ##
//  ######   #######  ##     ## ##     ## ##     ## ##    ## ########

// ##       #### ##    ## ########       ###    ########   ######    ######
// ##        ##  ###   ## ##            ## ##   ##     ## ##    ##  ##    ##
// ##        ##  ####  ## ##           ##   ##  ##     ## ##        ##
// ##        ##  ## ## ## ######      ##     ## ########  ##   ####  ######
// ##        ##  ##  #### ##          ######### ##   ##   ##    ##        ##
// ##        ##  ##   ### ##          ##     ## ##    ##  ##    ##  ##    ##
// ######## #### ##    ## ########    ##     ## ##     ##  ######    ######
//--------------------------------------------

//Ask user to choose options at command line for script generation
void Choose_Arguments_From_CommandLine(char& include_systematics, char& include_statistical, TString& template_name, TString& region, TString& lumiName, int& mode_histoBins, char& use_SM_setup, TString& selection, TString& filename_template_suffix, char& include_otherRegions)
{
    char chartmp = 'a'; int inttmp = -1; TString tstringtmp = ""; //Buffers

    cout<<endl<<BOLD(ITAL(DIM(<<"['0' <-> use value set in main()] ... ")))<<endl<<endl;

	//-- Choose whether to include shape syst or not
	cout<<endl<<FYEL("=== Include "<<UNDL(<<"shape"))<<FYEL(" systematics in the datacards ? ===")<<endl;
    cout<<ITAL(DIM(<<"['y'/'n'] ... "));
	cin>>chartmp;
	while(chartmp != 'y' && chartmp != 'n' && chartmp != '0')
	{
		cin.clear();
		cin.ignore(1000, '\n');

		cout<<" Wrong answer ! Need to type 'y' or 'n' ! Retry :"<<endl;
		cin>>chartmp;
	}
    if(chartmp == 'y' || chartmp == 'n') {include_systematics = chartmp;}

    //-- Choose whether to include statistic uncert. or not (e.g. for pulls of NPs, need to remove them, else too long)
    cout<<endl<<FYEL("=== Include "<<UNDL(<<"statistical"))<<FYEL(" uncertainties in the datacards ? ===")<<endl;
    cout<<ITAL(DIM(<<"['y'/'n'] ... "));
    cin>>chartmp;
    while(chartmp != 'y' && chartmp != 'n' && chartmp != '0')
    {
        cin.clear();
        cin.ignore(1000, '\n');

        cout<<" Wrong answer ! Need to type 'y' or 'n' ! Retry :"<<endl;
        cin>>chartmp;
    }
    if(chartmp == 'y' || chartmp == 'n') {include_statistical = chartmp;}

    //-- Choose whether to perform a fit following the main SM differential tZq analysis (overrides some options)
    // cout<<endl<<FYEL("=== Do you wish to perform a fit in the SM differential tZq setup ? ===")<<endl;
    // cout<<ITAL(DIM(<<"['y' for yes (else no)]"));
    // cout<<ITAL(DIM(<<"... "));
    // cin>>use_SM_setup;
    // if(use_SM_setup == 'y') {return;} //Overrides following options

    //-- Choose whether to include templates from 'Other regions' (ttZ4l, CRs -- not impacted by NNs)
    cout<<endl<<FYEL("=== Do you want to include 'other regions' in the fit ? (<-> ttZ 4l SR + WZ CR + ZZ CR) ===")<<endl;
    cout<<ITAL(DIM(<<"['y' = yes / 'n' = no / 1 = only include other regions / 0 = keep default]"));
    cout<<ITAL(DIM(<<"... "));
    cin>>chartmp;
    while(chartmp != 'y' && chartmp != 'n' && chartmp != '0' && chartmp != '1')
    {
        cin.clear();
        cin.ignore(1000, '\n');

        cout<<" Wrong answer ! Need to type 'y' or 'n' ! Retry :"<<endl;
        cin>>chartmp;
    }
    if(chartmp == 'y' || chartmp == 'n' || chartmp == '1') {include_otherRegions = chartmp;}

    //Choose whether to create a single datacard for entire histo, or separate datacards for each histo bin (allows to parametrize each bin independently)
    cout<<endl<<FYEL("=== Create separate datacards for each histogram bin ? ===")<<endl;
    cout<<ITAL(DIM("1 <-> treat each histogram bin individually (separate datacards & histos) for individual parametrizations [EFT default]"))<<endl;
    cout<<ITAL(DIM("2 <-> use full MVA distribution [SM default]"))<<endl;
    cout<<ITAL(DIM("3 <-> entire histogram treated as single bin (counting experiment)"))<<endl;
    cout<<ITAL(DIM(<<"... "));
    cin>>inttmp;
    while(inttmp != 0 && inttmp != 1 && inttmp != 2 && inttmp != 3)
    {
        cin.clear();
        cin.ignore(1000, '\n');

        cout<<" Wrong answer ! Need to type '1' / '2' / '3' or '0' ! Retry :"<<endl;
        cin>>inttmp;
    }
    if(inttmp == 1 || inttmp == 2 || inttmp == 3) {mode_histoBins = inttmp;}

    //Set the template name (e.g.'NN') to be looked for in the rootfiles //If ignored, use the value set in the main
    cout<<endl<<FYEL("=== Set the template name (to read corresponding histograms in SRs) ===")<<endl;
    cout<<ITAL(DIM("'Zpt' / 'NN' / 'BDT' / ..."))<<endl;
    cout<<ITAL(DIM("(0 <-> use the template name defined in the main)"))<<endl;
    // cout<<ITAL(DIM("(1 <-> create datacards for control plot in SR-3l (hardcoded))"))<<endl;
    cout<<ITAL(DIM("(NB: if you set a template name, it will be used for all regions defined in the main)"))<<endl;
    cout<<ITAL(DIM(<<"..."));
    cin>>tstringtmp;
    if(tstringtmp != "0") {template_name = tstringtmp;}

	//Choose the luminosity
    cout<<endl<<FYEL("=== Choose the luminosity ===")<<endl;
    cout<<ITAL(DIM(<<"['Run2'/'2016'/'2017'/'2018'/'201617'/'201618'/'201718'] ... "));
    cin>>tstringtmp;
    while(tstringtmp != "0" && tstringtmp != "Run2" && tstringtmp != "2016" && tstringtmp != "2017" && tstringtmp != "2018" && tstringtmp != "201617" && tstringtmp != "201618" && tstringtmp != "201718")
    {
        cin.clear();
        cin.ignore(1000, '\n');

        cout<<" Wrong answer ! Retry :"<<endl;
        cin>>tstringtmp;
    }
    if(tstringtmp == "Run2" || tstringtmp == "2016" || tstringtmp == "2017" || tstringtmp == "2018" || tstringtmp == "201617" || tstringtmp == "201618" || tstringtmp == "201718") {lumiName = tstringtmp;}

    //Set a 'filename suffix' if needed (if present in the filename)
    cout<<endl<<FYEL("=== Set the template name suffix for the filename ===")<<endl;
    cout<<ITAL(DIM("'EFT1' / 'EFT2' / 'SM' / ..."))<<endl;
    cout<<ITAL(DIM(<<"(NB: 1 <-> '')"));
    cout<<ITAL(DIM(<<"..."));
    cin>>tstringtmp;
    if(tstringtmp != "0") {filename_template_suffix = tstringtmp;}
    if(tstringtmp == "1") {filename_template_suffix = "";} //No suffix

    //Set a 'selection flag' if necessary (if present in the histo/file names)
    cout<<endl<<FYEL("=== Set the event selection flag ===")<<endl;
    cout<<ITAL(DIM("'signal' / 'tZq' / 'ttZ' / ..."))<<endl;
    cout<<ITAL(DIM(<<"..."));
    cin>>tstringtmp;
    if(tstringtmp != "0") {selection = tstringtmp;}

    return;
}










//--------------------------------------------
// ##     ##    ###    #### ##    ##
// ###   ###   ## ##    ##  ###   ##
// #### ####  ##   ##   ##  ####  ##
// ## ### ## ##     ##  ##  ## ## ##
// ##     ## #########  ##  ##  ####
// ##     ## ##     ##  ##  ##   ###
// ##     ## ##     ## #### ##    ##
//--------------------------------------------

int main()
{
// Can set options here
//--------------------------------------------
    vector<TString> v_channel; //'all', 'uuu', 'eeu', 'uue', 'eee' //Used only for leptonic channels (use v_regions for 'regions')
    v_channel.push_back("all");
    // v_channel.push_back("uuu");
    // v_channel.push_back("eeu");
    // v_channel.push_back("uue");
    // v_channel.push_back("eee");

    //-- Two ways to combine regions (e.g. 'SRtZq') and templates (e.g. 'Zpt') together
    // a) Define here the same number of regions & templates --> will associate them 1 to 1
    // b) Define ==1 template (possibly via command-line) and >= 1 region --> will associate the template with all considered regions
    // NB: ORDER MATTERS /!\

    vector<TString> v_templates; //'NN', 'BDT', ...
    v_templates.push_back("Zpt");
    // v_templates.push_back("NN");
    // v_templates.push_back("categ");

    vector<TString> v_regions; //'SR', 'CR_xx', ... (must reflect bin names)
    v_regions.push_back("SRtZq");
    v_regions.push_back("SRttZ");
    v_regions.push_back("SRother");

    //v_regions.push_back("signal");

    TString selection = ""; //Main event selection, before sub-categorization
    TString filename_template_suffix = "EFT2"; //Specify extension in histo filename
    bool scan_operator_hardcoded = false; //true <-> will generate datacards for several different bin names (scan steps) to be used in a script

    bool use_rph = false; //true <-> will also produce datacards for using RooParametric
    TString path_tmp_workspace = "../../../EFT/WS.root"; //Hard-coded path to the temporary workspace potentially containing RooParametricHists, etc. //Must be provided to card parser so that the WS contents get read

    bool hardcoded_spectator_variables = false; //true <-> some params are hardcoded to create cards for spectator variables (to make prefit/postfit plots with total uncertainties)

    TString scriptname = "makeDatacardsForTemplateFit.sh"; //Name of output script to create and run

//-- Modified at command-line
//--------------------------------------------
    TString lumiName = "Run2"; //'2016','2017','2018','201617','201618','201718','Run2'

    char include_systematics = 'y';//'y' <-> datacards will include syst. uncertainties (as specified in template datacard)
    char include_statistical = 'y';//'y' <-> datacards will include stat. uncertainty (as specified in template datacard)

//--------------------------------------------

//Automated
//--------------------------------------------
    int mode_histoBins = 1; //1 <-> treat each histogram bin individually (separate datacards & histos) for individual parametrizations (default); 2 <-> use full MVA distribution; 3 <-> entire histogram treated as single bin (counting experiment)
    char use_SM_setup = 'n'; //'y' <-> overrides some option to perform a fit following the setup of the main SM tZq differential analysis
    TString template_name = "0", region = "0";
    char include_otherRegions = 'y'; //'y' <-> will include 'other regions' in the fit (hardcoded: ttZ 4l SR / WZ CR / ZZ CR)
    if(hardcoded_spectator_variables) {use_rph=false; mode_histoBins=2;}
    else{Choose_Arguments_From_CommandLine(include_systematics, include_statistical, template_name, region, lumiName, mode_histoBins, use_SM_setup, selection, filename_template_suffix, include_otherRegions);}
    if(mode_histoBins==1) {use_rph = false;} //Don't split per bin if using RPHs
    if(include_otherRegions == '1') {v_regions.clear(); v_templates.clear();} //Only include other regions in the fit, not SR (e.g. to get postfit place only in CRs, etc.)

    Script_Datacards_TemplateFit(include_systematics, include_statistical, template_name, region, v_templates, v_channel, v_regions, lumiName, mode_histoBins, scan_operator_hardcoded, use_SM_setup, selection, filename_template_suffix, include_otherRegions, false, path_tmp_workspace, scriptname, hardcoded_spectator_variables);

    if(use_rph && mode_histoBins==2) {Script_Datacards_TemplateFit(include_systematics, include_statistical, template_name, region, v_templates, v_channel, v_regions, lumiName, mode_histoBins, scan_operator_hardcoded, use_SM_setup, selection, filename_template_suffix, include_otherRegions, true, path_tmp_workspace, scriptname);}

    Run_Commands(scriptname);

	return 0;
}
