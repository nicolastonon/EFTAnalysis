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

//HARD-CODED: define here whether the combination of region/template is valid (according to SM tZq differential analysis)
//FIXME remove ?
/*
bool Is_Template_Matching_Region(TString templatename, TString region, char use_SM_setup, char include_otherRegions)
{
    if(use_SM_setup == 'y')
    {
        if(templatename == "NN" && (region == "tZq" || region == "ttZ")) {return true;}
        if(templatename == "channel" && region == "Vg") {return true;}
        if(templatename == "mTW" && region == "wz") {return true;}
        if(templatename == "countExp" && region == "zz") {return true;}
    }
    else if(include_otherRegions == 'y')
    {
        if(templatename == "NN" && (region == "tZq" || region == "ttZ")) {return true;}
        if(templatename == "channel" && region == "Vg") {return true;}
        if(templatename == "mTW" && (region == "wz" || region == "CRWZ")) {return true;}
        if(templatename == "countExp" && (region == "zz" || region == "SRttZ4l" || region == "CRZZ")) {return true;}
    }
    else {return true;} //Incorrect case


    return false;
}
*/

//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Existence(const TString& name)
{
    struct stat buffer;
    bool found = (stat (name.Data(), &buffer) == 0); //true if file exists
    return found;
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
void Script_Datacards_TemplateFit(char include_systematics, char include_statistical, TString template_name, TString region, vector<TString> v_templates, vector<TString> v_channel, vector<TString> v_regions, TString lumiName, int mode_histoBins, bool scan_operator_hardcoded, char use_SM_setup, TString selection, TString filename_template_suffix, char include_otherRegions)
{
//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

    //Check if use shape syst or not
	TString systChoice;
	if(include_systematics == 'y') {systChoice = "withShape";}
	else if(include_systematics == 'n') {systChoice = "noShape";}
    else {cout<<"Wrong arguments ! Abort !"<<endl; return;}

    //Check if use stat. uncert. or not
	TString statChoice;
	if(include_statistical == 'y') {statChoice = "withStat";}
	else if(include_statistical == 'n') {statChoice = "noStat";}
    else {cout<<"Wrong arguments ! Abort !"<<endl; return;}

    // If specific arguments were chosen at command line, modify the vectors defined in the main
    if(template_name != "0") {v_templates.resize(0); v_templates.push_back(template_name);}
    else if(!v_templates.size()) {cout<<"Template name not set ! Abort !"<<endl; return;}

    // Trick: in case there is ==1 template and >=1 regions defined, we want to use this template for all the regions; but to allow to then add additional pairs of regions/templates, expand the current lists so that they match 1 to 1
    if(v_templates.size() == 1 && v_regions.size() >= 1)
    {
        for(int iregion=1; iregion<v_regions.size(); iregion++)
        {
            v_templates.push_back(v_templates[0]); //Duplicate first element for each additional region
        }
    }

    //FIXME -- HARDCODED
    vector<TString> v_WCs_operator_scan1 = {"-999"}; //HARDCODED
    TString operator_scan1 = "";
    if(scan_operator_hardcoded)
    {
        mode_histoBins = 1; //Scan on parametrized NN --> must treat each histogram bin separately

        operator_scan1 = "ctz";
        //v_WCs_operator_scan1 = {"-4","-2","-1","0","1","2","4"};
        v_WCs_operator_scan1 = {"-3", "-2", "-1.5", "-1", "-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6", "0.8", "1", "1.5", "2", "3"};

        v_regions.resize(0); v_templates.resize(0);
        v_regions.push_back("SRtZq"); v_templates.push_back("NN");
        v_regions.push_back("SRttZ"); v_templates.push_back("NN");
    }
    //Add predefined additional regions (NB: corresponding templates must have been produced and merged into main template file !)
    else if(include_otherRegions == 'y')
    {
        v_regions.push_back("SRttZ4l"); v_templates.push_back("countExp");
        v_regions.push_back("CRWZ"); v_templates.push_back("mTW");
        v_regions.push_back("CRZZ"); v_templates.push_back("countExp");
    }
    //Overrides some option to perform a fit following the setup of the main SM tZq differential analysis
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

    vector<TString> v_lumiYears;
    if(lumiName == "2016") {v_lumiYears.push_back("2016");}
    else if(lumiName == "2017") {v_lumiYears.push_back("2017");}
    else if(lumiName == "2018") {v_lumiYears.push_back("2018");}
    else if(lumiName == "201617") {v_lumiYears.push_back("2016"); v_lumiYears.push_back("2017");}
    else if(lumiName == "201618") {v_lumiYears.push_back("2016"); v_lumiYears.push_back("2018");}
    else if(lumiName == "201718") {v_lumiYears.push_back("2017"); v_lumiYears.push_back("2018");}
    else if(lumiName == "Run2") {v_lumiYears.push_back("2016"); v_lumiYears.push_back("2017"); v_lumiYears.push_back("2018");}


	ofstream file_out("makeDatacardsForTemplateFit.sh"); //output script

	TString dir = "./datacards_TemplateFit/";

    int nbins = -1; //If 'mode_histoBins=1', will read histogram binning in rootfile and create 1 datacard per histo bin
    vector<int> v_nbins; //Actually binnings can differ for each template, so keep track of all of them (need them twice: first to create individual cards, then to combine them)


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


    for(int ipt_EFT=0; ipt_EFT<v_WCs_operator_scan1.size(); ipt_EFT++)
    {
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
            for(int iregion=0; iregion<v_regions.size(); iregion++) //Loop over regions
            {
                // cout<<"v_regions[iregion] "<<v_regions[iregion]<<endl;

                for(int itemplate=0; itemplate<v_templates.size(); itemplate++) //Loop over templates
                {
                    if(itemplate != iregion) {continue;} //NEW -- associate templates/regions 1 to 1 for more flexibility

                    // if((use_SM_setup == 'y' || include_otherRegions == 'y') && !Is_Template_Matching_Region(v_templates[itemplate], v_regions[iregion], use_SM_setup, include_otherRegions)) {continue;}

                    // cout<<"v_templates[itemplate] "<<v_templates[itemplate]<<endl;

                    TString var = v_templates[itemplate] + "_" + v_regions[iregion];
                    if(scan_operator_hardcoded) {var+= "_" + operator_scan1 + "_" + v_WCs_operator_scan1[ipt_EFT];}
                    if(v_regions[iregion] == "CR" && filename_template_suffix.Contains("EFT")) {var = "mTW_CR";} //Use mTW for now in CR //Hard-coded

                    //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                    var.ReplaceAll('-', 'm');

        			//Make subdir for single datacards
        			file_out <<"mkdir "<<dir<<endl<<endl;
                    file_out <<"mkdir "<<dir+v_lumiYears[iyear]<<endl<<endl;

        		    // TString file_histos = "../templates/Combine_Input.root";
                    TString file_histos_pathFromHere = "./../templates/Templates_"+ (include_otherRegions? v_templates[0]:v_templates[itemplate]) + (filename_template_suffix? "_"+filename_template_suffix:"")+(selection != ""? "_"+selection:"")+"_"+lumiName+".root"; //For use within this code
                    if(scan_operator_hardcoded) {file_histos_pathFromHere = "./../templates/Templates_NN_EFT2param_Run2.root";} //HARD-CODED

                    cout<<DIM("Trying to open input file "<<file_histos_pathFromHere<<" ... ");
                    if(Check_File_Existence(file_histos_pathFromHere)) {cout<<DIM("FOUND !")<<endl;}
            		else
                    {
                        file_histos_pathFromHere = "./../templates/Templates_"+v_templates[itemplate]+(filename_template_suffix? "_"+filename_template_suffix:"")+(selection != ""? "_"+selection:"")+"_Run2.root"; //Try Run2 file

                        cout<<endl<<DIM("Trying to open input file "<<file_histos_pathFromHere<<" ... ");
                        if(Check_File_Existence(file_histos_pathFromHere)) {cout<<DIM("FOUND !")<<endl;}
                        else {cout<<BOLD(FRED("ERROR: input template file not found ! Abort !"))<<endl; return;}
                    }

                    TString file_histos = "../." + file_histos_pathFromHere; //Path to write into datacard
        			cout<<endl<<FMAG("---> Will use filepath : ")<<file_histos<<endl<<endl;

                    if(mode_histoBins==1) //Need to infer the number of bins of the considered histograms, in order to create 1 card per bin
                    {
                        TFile* f_tmp = TFile::Open(file_histos_pathFromHere);
                        TString hname_tmp = var + "_" + v_lumiYears[iyear ] + "__data_obs"; //Hard-coded: look for data histo to infer binning
						// cout<<"Reading histogram: "<<hname_tmp<<endl;

                        if(!f_tmp->GetListOfKeys()->Contains(hname_tmp) ) {cout<<BOLD(FRED("ERROR: histogram "<<hname_tmp<<" not found ! Can not infer histogram binning !"))<<endl; return;}
                        // cout<<"hname_tmp "<<hname_tmp<<endl;
                        TH1F* h_tmp = (TH1F*) f_tmp->Get(hname_tmp);
                        nbins = h_tmp->GetNbinsX();
						// cout<<"--> nbins = "<<nbins<<endl;
                        v_nbins.push_back(nbins); //Store the binning for this template, to be used again for combined card (--> will need to read vector in same order !)
                        delete h_tmp; h_tmp = NULL; f_tmp->Close();
                        if(nbins == -1) {cout<<BOLD(FRED("ERROR: histogram "<<hname_tmp<<" not found ! Can not infer histogram binning !"))<<endl; return;}
                    }

                    if(nbins<0) {nbins=1;} //Need to loop at least once by default
        			for(int ilepchan=0; ilepchan<v_channel.size(); ilepchan++)
        			{
                        for(int ibin=1; ibin<nbins+1; ibin++)
                        {
                            // if(v_templates[itemplate]=="categ" and ibin==6) {continue;} //HARDCODED TMP FIX (empty bin)

                            TString var_tmp = var;
                            if(mode_histoBins==1 && nbins > 1) {var_tmp = (TString) "bin" + Form("%d",ibin) + "_" + var;} //Also include bin number in naming scheme (--> will read single bin histos instead of full histos)
                            else if(mode_histoBins==2) {var_tmp = "countExp_" + var;}

            				file_out<<"python Parser_Datacard_Template.py "
                            + var_tmp + " "
                            + v_channel[ilepchan] + " "
                            + v_lumiYears[iyear] + " "
            				+ file_histos + " "
            				+ systChoice + " "
            				+ statChoice + " "
                            + dir + v_lumiYears[iyear] + " "
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

    int idx_v_nbins = 0; //Vector storing template binning info needs to be read in same order (same loops) as before; increment some index to keep track of the current element to read
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
                    if(itemplate != iregion) {continue;} //NEW -- associate templates/regions 1 to 1 for more flexibility

                    // if((use_SM_setup == 'y' || include_otherRegions == 'y') && !Is_Template_Matching_Region(v_templates[itemplate], v_regions[iregion], use_SM_setup, include_otherRegions)) {continue;}

                    TString var = v_templates[itemplate] + "_" + v_regions[iregion];
                    if(v_regions[iregion] == "CR" && filename_template_suffix.Contains("EFT")) {var = "mTW_CR";} //Use mTW for now in CR //Hard-coded

                    if(scan_operator_hardcoded) {var+= "_" + operator_scan1 + "_" + v_WCs_operator_scan1[ipt_EFT];}

                    //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                    var.ReplaceAll('-', 'm');

        			for(int ilepchan=0; ilepchan<v_channel.size(); ilepchan++) //Loop over channels
        			{
                        int nbins_tmp = 1;
                        if(mode_histoBins==1) {nbins_tmp = v_nbins[idx_v_nbins]; idx_v_nbins++;} //Read current binning; increment index to stay in sync
                        for(int ibin=1; ibin<nbins_tmp+1; ibin++)
                        {
                            // if(v_templates[itemplate]=="categ" and ibin==6) {continue;} //HARDCODED TMP FIX (empty bin)

                            TString var_tmp = var;
                            if(mode_histoBins && nbins_tmp > 1) {var_tmp = (TString) "bin" + Form("%d",ibin) + "_" + var;} //Also include bin number in naming scheme (--> will read single bin histos instead of full histos)
                            else if(mode_histoBins==2) {var_tmp = "countExp_" + var;}

            				file_out<<var_tmp;
            				if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
                            file_out<<"_"+v_lumiYears[iyear];
                            file_out<<"=" + dir + v_lumiYears[iyear] + "/"
            				+ "datacard_"+var_tmp;
            				if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
            				file_out<<".txt ";
                        } //bin loop
        			} //channel loop

        		} //template loop
        	} //region loop
        } //years loop

    	TString output_name = "COMBINED_Datacard_TemplateFit";
        if(systChoice == "noShape") output_name+= "_noShape";
        if(statChoice == "noStat") output_name+= "_noStat";
        if(scan_operator_hardcoded) {output_name+= "_" + operator_scan1 + "_" + v_WCs_operator_scan1[ipt_EFT];}
    	output_name+= "_" + lumiName + ".txt";

    	file_out<<"> "<<output_name<<endl<<endl;
    } //EFT loop

//datacard for single channels
/*
for(int ichan=0; ichan<v_channel.size(); ichan++)
{
	file_out<<"python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py ";

	for(int ivar=0; ivar<var_list.size(); ivar++)
	{
		file_out<<var_list[ivar] + "_" + v_channel[ichan] + "=datacard_"+v_channel[ichan]+"_"+var_list[ivar]+".txt ";
	}

	output_name = "COMBINED_Datacard_TemplateFit_";
	if(systChoice == "noShape") output_name+= "noShape_";
	output_name+= v_channel[ichan]+".txt";

	file_out<<"> "<<output_name<<endl<<endl;

	file_out<<"mv "<<output_name<<" datacards_TemplateFit/"<<endl<<endl;
}
*/

	file_out.close();

	system("chmod 755 makeDatacardsForTemplateFit.sh");

	cout<<FGRN("... Created script ./makeDatacardsForTemplateFit.sh !")<<endl;

	cout<<endl<<endl<<FGRN("-- Executing script ./makeDatacardsForTemplateFit.sh !")<<endl;

	system("./makeDatacardsForTemplateFit.sh");

	return;
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
	//Choose whether to include shape syst or not
	cout<<endl<<FYEL("=== Include "<<UNDL(<<"shape"))<<FYEL(" systematics in the datacards ? ===")<<endl;
    cout<<ITAL(DIM(<<"['y'/'n'] ... "));
	cin>>include_systematics;
	while(include_systematics != 'y' && include_systematics != 'n')
	{
		cin.clear();
		cin.ignore(1000, '\n');

		cout<<" Wrong answer ! Need to type 'y' or 'n' ! Retry :"<<endl;
		cin>>include_systematics;
	}

    //Choose whether to include statistic uncert. or not (e.g. for pulls of NPs, need to remove them, else too long)
    cout<<endl<<FYEL("=== Include "<<UNDL(<<"statistical"))<<FYEL(" uncertainties in the datacards ? ===")<<endl;
    cout<<ITAL(DIM(<<"['y'/'n'] ... "));
    cin>>include_statistical;
    while(include_statistical != 'y' && include_statistical != 'n')
    {
        cin.clear();
        cin.ignore(1000, '\n');

        cout<<" Wrong answer ! Need to type 'y' or 'n' ! Retry :"<<endl;
        cin>>include_statistical;
    }

    //Choose whether to perform a fit following the main SM differential tZq analysis (overrides some options)
    // cout<<endl<<FYEL("=== Do you wish to perform a fit in the SM differential tZq setup ? ===")<<endl;
    // cout<<ITAL(DIM(<<"['y' for yes (else no)]"));
    // cout<<ITAL(DIM(<<"... "));
    // cin>>use_SM_setup;
    // if(use_SM_setup == 'y') {return;} //Overrides following options

    //Choose whether to perform a fit following the main SM differential tZq analysis (overrides some options)
    cout<<endl<<FYEL("=== Do you want to include 'other regions' in the fit ? (ttZ 4l SR / WZ CR / ZZ CR) ===")<<endl;
    cout<<ITAL(DIM("(NB: requires you to have produced the relevant templates and merged them with SR templates)"))<<endl;
    cout<<ITAL(DIM(<<"['y' for yes (else no)]"));
    cout<<ITAL(DIM(<<"... "));
    cin>>include_otherRegions;
    if(include_otherRegions != 'y') {include_otherRegions = 'n';} //Default

    //Choose whether to create a single datacard for entire histo, or separate datacards for each histo bin (allows to parametrize each bin independently)
    cout<<endl<<FYEL("=== Create separate datacards for each histogram bin ? ===")<<endl;
    cout<<ITAL(DIM("0 <-> use full MVA distribution (default)"))<<endl;
    cout<<ITAL(DIM("1 <-> treat each histogram bin individually (separate datacards & histos) for individual parametrizations"))<<endl;
    cout<<ITAL(DIM("2 <-> entire histogram treated as single bin (counting experiment)"))<<endl;
    cout<<ITAL(DIM(<<"... "));
    cin>>mode_histoBins;
    while(mode_histoBins != 0 && mode_histoBins != 1 && mode_histoBins != 2)
    {
        cin.clear();
        cin.ignore(1000, '\n');

        cout<<" Wrong answer ! Need to type '0' / '1' / '2' ! Retry :"<<endl;
        cin>>mode_histoBins;
    }

	//Choose the luminosity
    cout<<endl<<FYEL("=== Choose the luminosity ===")<<endl;
    cout<<ITAL(DIM(<<"['Run2'/'2016'/'2017'/'2018'/'201617'/'201618'/'201718'] ... "));
    cin>>lumiName;
    while(lumiName != "0" && lumiName != "Run2" && lumiName != "2016" && lumiName != "2017" && lumiName != "2018" && lumiName != "201617" && lumiName != "201618" && lumiName != "201718")
    {
        cin.clear();
        cin.ignore(1000, '\n');

        cout<<" Wrong answer ! Retry :"<<endl;
        cin>>lumiName;
    }
    if(lumiName == "0") {lumiName = "Run2";}

    //Set the template name (e.g.'NN') to be looked for in the rootfiles //If ignored, use the value set in the main
    TString template_name_tmp;
    cout<<endl<<FYEL("=== Set the template name (to read corresponding histograms in SRs) ===")<<endl;
    cout<<ITAL(DIM("'0' <-> use value set in main() / 'Zpt' / 'NN' / 'BDT' / ..."))<<endl;
    cout<<ITAL(DIM("(NB: if you set a template name, it will be used for all regions defined in the main')"))<<endl;
    cout<<ITAL(DIM(<<"..."));
    cin>>template_name_tmp;
    if(template_name_tmp != "0") {template_name = template_name_tmp;}

    //Set a 'selection flag' if necessary (if present in the histo/file names)
    TString selection_tmp;
    cout<<endl<<FYEL("=== Set the event selection flag ===")<<endl;
    cout<<ITAL(DIM("'0' <-> use value set in main() / 'signal' / 'tZq' / 'ttZ' / ..."))<<endl;
    cout<<ITAL(DIM(<<"..."));
    cin>>selection_tmp;
    if(selection_tmp != "0") {selection = selection_tmp;}

    //Set a 'filename suffix' if needed (if present in the filename)
    TString filename_template_suffix_tmp;
    cout<<endl<<FYEL("=== Set the template name suffix for the filename ===")<<endl;
    cout<<ITAL(DIM("'0' <-> use value set in main() / 'EFT1' / 'EFT2' / 'SM' / ..."))<<endl;
    cout<<ITAL(DIM(<<"..."));
    cin>>filename_template_suffix_tmp;
    if(filename_template_suffix_tmp != "0") {filename_template_suffix = filename_template_suffix_tmp;}

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
    vector<TString> v_channel; //'all', 'uuu', 'eeu', 'uue', 'eee'
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
    v_regions.push_back("CR");

    TString selection = ""; //Main event selection, before sub-categorization
    TString filename_template_suffix = ""; //Specify extension in histo filename
    bool scan_operator_hardcoded = false; //true <-> will generate datacards for several different bin names (scan steps) to be used in a script

// Modified at command-line
//--------------------------------------------
    TString lumiName = "Run2"; //'2016','2017','2018','201617','201618','201718','Run2'

    char include_systematics = 'n';//'y' <-> datacards will include syst. uncertainties (as specified in template datacard)
    char include_statistical = 'n';//'y' <-> datacards will include stat. uncertainty (as specified in template datacard)

//--------------------------------------------

//Automated
//--------------------------------------------
    int mode_histoBins = 0; //0 <-> use full MVA distribution (default); 1 <-> treat each histogram bin individually (separate datacards & histos) for individual parametrizations; 2 <-> entire histogram treated as single bin (counting experiment)
    char use_SM_setup = 'n'; //'y' <-> overrides some option to perform a fit following the setup of the main SM tZq differential analysis
    TString template_name = "0", region = "0";
    char include_otherRegions = 'n'; //'y' <-> expect to find merged template files containing both SR templates and templates from other regions (hardcoded: ttZ 4l SR / WZ CR / ZZ CR)
    Choose_Arguments_From_CommandLine(include_systematics, include_statistical, template_name, region, lumiName, mode_histoBins, use_SM_setup, selection, filename_template_suffix, include_otherRegions);

	Script_Datacards_TemplateFit(include_systematics, include_statistical, template_name, region, v_templates, v_channel, v_regions, lumiName, mode_histoBins, scan_operator_hardcoded, use_SM_setup, selection, filename_template_suffix, include_otherRegions);

	return 0;
}
