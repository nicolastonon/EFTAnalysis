/* BASH COLORS */
#define RST   "[0m"
#define KRED  "[31m"
#define KGRN  "[32m"
#define KYEL  "[33m"
#define KBLU  "[34m"
#define KMAG  "[35m"
#define KCYN  "[36m"
#define KWHT  "[37m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST
#define BOLD(x) "[1m" x RST
#define UNDL(x) "[4m" x RST

#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TTree.h"
#include "TString.h"
#include "TColor.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLine.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>

#include <cassert> //Can be used to terminate program if argument is not true //Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> // to be able to use mkdir

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

//Convert a double into a TString
// precision --> can choose if TString how many digits the TString should display
TString Convert_Number_To_TString(double number, int precision/*=3*/)
{
	stringstream ss;
	ss << std::setprecision(precision) << number;
	TString ts = ss.str();
	return ts;
}

//Convert a TString into a float
double Convert_TString_To_Number(TString ts)
{
	double number = 0;
	string s = ts.Data();
	stringstream ss(s);
	ss >> number;
	return number;
}

//Can set here protections : return false if a given syst does not apply to a given sample
bool Is_Syst_Match_Sample(TString syst, TString sample, bool use_rph)
{
	// cout<<"syst "<<syst<<endl;
	// cout<<"sample "<<sample<<endl;

    //FIXME -- testing
    // if(sample.Contains("PrivMC")) {return false;} //Combine does not yet support interpolation for RooParametricHists --> nuisances incorporated directly in bin parameterizations
    if(use_rph && sample.Contains("PrivMC")) {return false;} //Combine does not yet support interpolation for RooParametricHists --> nuisances incorporated directly in bin parameterizations

    if( (syst.Contains("Fake", TString::kIgnoreCase) || syst.BeginsWith("FR") || syst.Contains("NPL")) && !sample.Contains("NPL")) {return false;}
    else if(sample.Contains("NPL") && !syst.Contains("Fake", TString::kIgnoreCase) && !syst.BeginsWith("FR") && !syst.Contains("NPL") && !syst.Contains("CRDY") ) {return false;}

    else if(syst.Contains("njets_tZq", TString::kIgnoreCase) && sample != "PrivMC_tZq") {return false;}

    else if(syst.Contains("tZq", TString::kIgnoreCase) && !sample.Contains("tZq")) {return false;}
    else if(syst.Contains("ttZ", TString::kIgnoreCase) && !sample.Contains("ttZ")) {return false;}
    else if(syst.Contains("tWZ", TString::kIgnoreCase) && !sample.Contains("tWZ")) {return false;}

    else if(syst.Contains("CRWZ", TString::kIgnoreCase) && !sample.Contains("WZ")) {return false;}
    else if(syst.Contains("CRZZ", TString::kIgnoreCase) && !sample.Contains("VVV") && !sample.Contains("ZZ")) {return false;}
    else if(syst.Contains("CRDY", TString::kIgnoreCase) && !sample.Contains("WZ") && !sample.Contains("VVV") && !sample.Contains("ZZ") && !sample.Contains("XG") && !sample.Contains("NPL") && !sample.Contains("DY") && !sample.Contains("TTbar")) {return false;}

    //else if((syst == "PDF" || syst == "alphas" || syst.BeginsWith("ME")) && !sample.Contains("PrivMC") && sample != "tZq" && sample != "ttZ" && sample != "tWZ" ) {return false;} //Hardcoded: PDF uncertainty considered/correlated for signals only
    else if((syst == "PDF" || syst == "alphas" || syst.BeginsWith("ME")) && !sample.Contains("PrivMC") && sample != "tZq" && sample != "ttZ") {return false;} //Some missing for tWZ

    //else if((syst=="JES" || syst.Contains("JER") || syst.Contains("MET")) && !sample.Contains("PrivMC")) {return false;}

	return true;
}

//Ask user to choose options at command line for script generation
void Choose_Arguments_From_CommandLine(TString& signal, bool& use_rph)
{
    //-- Choose whether to include shape syst or not
	cout<<endl<<FYEL("--- What is your SIGNAL ?")<<endl;
    cout<<"* 'eft'   \t<-> Signals are SMEFT tZq+ttZ+tWZ"<<endl;
    cout<<"* 'efttzq'   <-> Signal is SMEFT tZq only"<<endl;
    cout<<"* 'eftttz'   <-> Signal is SMEFT ttZ only"<<endl;
    cout<<"* 'efttwz'   <-> Signal is SMEFT tWZ only"<<endl;
    cout<<"* '0'   <-> Signals are tZq+ttZ+tWZ"<<endl;
    // cout<<"* 'thq' <-> Signals are tHq + tHW"<<endl;
    cout<<"* 'tzq' <-> Signal is tZq"<<endl;
    cout<<"* 'ttz' <-> Signal is ttZ"<<endl;
    cout<<"* 'twz' <-> Signal is tWZ"<<endl;
	cin>>signal;
	while(signal != "tzq" && signal != "ttz" && signal != "twz" && signal != "0" && signal != "eft" && signal != "efttzq" && signal != "eftttz" && signal != "efttwz")
	{
		cin.clear();
		cin.ignore(1000, '\n');

		cout<<" Wrong answer ! Retry :"<<endl;
		cin>>signal;
	}

    //-- Choose whether to treat EFT signals as RooParametricHists
    /*
    char chartmp = 'a';
	cout<<endl<<FYEL("--- Use RooParametricHists instead of TH1s for EFT signals ? (<-> else need to split each hist bin)")<<endl;
    cout<<"* y/n, 0 <-> keep default value in code"<<endl;
	cin>>chartmp;
	while(chartmp != 'y' && chartmp != 'n' && chartmp != '0')
	{
		cin.clear();
		cin.ignore(1000, '\n');

		cout<<" Wrong answer ! Retry :"<<endl;
		cin>>chartmp;
	}
    if(chartmp == 'y') {use_rph = true;}
    else if(chartmp == 'n') {use_rph = false;}
    */

	return;
}



//--------------------------------------------
// ##     ##    ###    ##    ## ########    ######## ######## ##     ## ########  ##          ###    ######## ########
// ###   ###   ## ##   ##   ##  ##             ##    ##       ###   ### ##     ## ##         ## ##      ##    ##
// #### ####  ##   ##  ##  ##   ##             ##    ##       #### #### ##     ## ##        ##   ##     ##    ##
// ## ### ## ##     ## #####    ######         ##    ######   ## ### ## ########  ##       ##     ##    ##    ######
// ##     ## ######### ##  ##   ##             ##    ##       ##     ## ##        ##       #########    ##    ##
// ##     ## ##     ## ##   ##  ##             ##    ##       ##     ## ##        ##       ##     ##    ##    ##
// ##     ## ##     ## ##    ## ########       ##    ######## ##     ## ##        ######## ##     ##    ##    ########
//--------------------------------------------

/**
 * Create output text file containing skeleton/template of COMBINE datacard
 * NB : must make sure that arguments given to function (e.g. sample names, syst names, etc.) are in sync with what is in the template !

 * @param outfile_name   output name
 * @param v_samples      list of samples in datacard
 * @param v_isSignal     "0" <-> sample is bkg, "1" <-> sample is signal
 * @param v_sampleUncert uncertainty associated to sample, in % (only for bkg)
 * @param v_normSyst         list of systematics applying to all samples in %(e.g. lumi)
 * @param v_normSystValue    value of these systematics
 * @param v_shapeSyst    list of shape systematics
 */
void Generate_Datacard(vector<TString> v_samples, vector<int> v_isSignal, vector<float> v_sampleUncert, vector<TString> v_normSyst, vector<TString> v_normSystValue, vector<TString> v_shapeSyst, TString signal, vector<bool> v_shapeSyst_isCorrelYears, bool use_rph, TString outfile_name="Template_Datacard.txt")
{
    //TString outfile_name = "Template_Datacard.txt";
    ofstream outfile(outfile_name.Data());
    // bool include_njets_syst = false; //Obsolete

	//-- Protections
	if(v_shapeSyst.size() != v_shapeSyst_isCorrelYears.size()) {cout<<"ERROR: incorrect size for vector 'v_shapeSyst_isCorrelYears' !"<<endl; return;}


 // #    # ######   ##   #####  ###### #####
 // #    # #       #  #  #    # #      #    #
 // ###### #####  #    # #    # #####  #    #
 // #    # #      ###### #    # #      #####
 // #    # #      #    # #    # #      #   #
 // #    # ###### #    # #####  ###### #    #
//--------------------------------------------

    //-- Add keyword at b.o.f. to indicate whether this template corresponds to the use of RooParametricHists, or only regular TH1s
    if (use_rph && signal.Contains("eft"))
    {
        outfile<<"#RPH"<<endl;
    }

    //--- Nof observables, bkgs, nuisance params
    outfile<<"imax"<<"\t"<<1<<"\t"<<"number of categories"<<endl;
    outfile<<"jmax"<<"\t"<<"*"<<"\t"<<"number of backgrounds"<<endl;
    outfile<<"kmax"<<"\t"<<"*"<<"\t"<<"number of nuisance parameters (sources of systematic uncertainties)"<<endl;

//--------------------------------------------
    //--- Filepath, naming convention
    outfile<<"---------------------------------------------------"<<endl;

    if(use_rph && signal.Contains("eft"))
    {
        outfile<<"[SR] shapes * [VAR]_[CHAN]_[YEAR] workspacetoread w:rdh_$CHANNEL__$PROCESS w:rdh_$CHANNEL__$PROCESS__$SYSTEMATIC"<<endl; //RDH

        outfile<<"[SR]shapes PrivMC_tZq [VAR]_[CHAN]_[YEAR] workspacetoread w:rph_$CHANNEL__$PROCESS"<<endl;
        outfile<<"[SR]shapes PrivMC_ttZ [VAR]_[CHAN]_[YEAR] workspacetoread w:rph_$CHANNEL__$PROCESS"<<endl;
        outfile<<"[SR]shapes PrivMC_tWZ [VAR]_[CHAN]_[YEAR] workspacetoread w:rph_$CHANNEL__$PROCESS"<<endl;

        outfile<<"[CR]shapes * [VAR]_[CHAN]_[YEAR] filetoread $CHANNEL__$PROCESS $CHANNEL__$PROCESS__$SYSTEMATIC"<<endl;
    }
    else
    {
        outfile<<"shapes * [VAR]_[CHAN]_[YEAR] filetoread $CHANNEL__$PROCESS $CHANNEL__$PROCESS__$SYSTEMATIC"<<endl;
    }

//--------------------------------------------
    //--- Var name, get yields from templates
    outfile<<"---------------------------------------------------"<<endl;
    outfile<<"bin        "<<"\t"<<"[VAR]_[CHAN]_[YEAR]"<<endl;
    outfile<<"observation"<<"\t"<<-1<<endl;

//--------------------------------------------
    //--- Processes names & indices
    outfile<<"---------------------------------------------------"<<endl;
    outfile<<"bin    ";
    for(int isample=0; isample<v_samples.size(); isample++)
    {
        outfile<<"\t";

        outfile<<"[VAR]_[CHAN]_[YEAR]";
    }
    outfile<<endl;

    outfile<<"process";
    for(int isample=0; isample<v_samples.size(); isample++)
    {
        outfile<<"\t";
        outfile<<v_samples[isample];
    }
    outfile<<endl;

    outfile<<"process";
    for(int isample=0; isample<v_samples.size(); isample++)
    {
        int index_tmp = isample;
        if(v_isSignal[isample] == 1 && isample > 0) {index_tmp = -isample;}
        else if(v_isSignal[isample] != 0 && v_isSignal[isample] != 1) {cout<<"Error ! Wrong value of v_isSignal !"<<endl; return;}

        outfile<<"\t";
        outfile<<index_tmp;
    }
    outfile<<endl;

    outfile<<"rate ";
    for(int isample=0; isample<v_samples.size(); isample++)
    {
        int index_tmp = isample;
        if(v_isSignal[isample] == 1 && isample > 0) {index_tmp = -isample;}
        else if(v_isSignal[isample] != 0 && v_isSignal[isample] != 1) {cout<<"Error ! Wrong value of v_isSignal !"<<endl; return;}

        outfile<<"\t";
        if(use_rph && v_isSignal[isample] == 1 && v_samples[isample].Contains("PrivMC")) {outfile<<"[RATE_SIG]";} //NEW -- if using RooParametricHist, must use 'rate 1' (multiplies norm)
        else {outfile<<-1;}
    }
    outfile<<endl;


 //                      #     #
 // #       ####   ####  ##    #     ####  #   #  ####  #####
 // #      #    # #    # # #   #    #       # #  #        #
 // #      #    # #      #  #  #     ####    #    ####    #
 // #      #    # #  ### #   # #         #   #        #   #
 // #      #    # #    # #    ##    #    #   #   #    #   #
 // ######  ####   ####  #     #     ####    #    ####    #

  //   ##   #      #
  //  #  #  #      #
  // #    # #      #
  // ###### #      #
  // #    # #      #
  // #    # ###### ######
//--------------------------------------------

//--- Systematics applied to all processes (e.g. lumi)
    outfile<<"---------------------------------------------------"<<endl;
    for(int isyst=0; isyst<v_normSyst.size(); isyst++)
    {
        //-- Year-specific markers
        if(v_normSyst[isyst].EndsWith("1617")) {outfile<<"[201617]";}
        else if(v_normSyst[isyst].EndsWith("1718")) {outfile<<"[201718]";}
        else if(v_normSyst[isyst].EndsWith("1618")) {outfile<<"[201618]";}
        else if(v_normSyst[isyst].EndsWith("16")) {outfile<<"[2016]";}
        else if(v_normSyst[isyst].EndsWith("17")) {outfile<<"[2017]";}
        else if(v_normSyst[isyst].EndsWith("18")) {outfile<<"[2018]";}

        //-- Region-specific markers
        if(v_normSyst[isyst].Contains("CRWZ")) {outfile<<"[CRWZ]";}
        else if(v_normSyst[isyst].Contains("CRZZ")) {outfile<<"[CRZZ]";}
        else if(v_normSyst[isyst].Contains("CRDY")) {outfile<<"[CRDY]";}

        outfile<<v_normSyst[isyst];
        // if(!v_normSyst_isCorrelYears[isyst]) {outfile<<"[YEAR]";} //Uncorrelated for each year

        outfile<<"\t"<<"lnN";

        for(int isample=0; isample<v_samples.size(); isample++)
        {
			// cout<<"sample "<<v_samples[isample]<<endl;
			// cout<<"syst "<<v_normSyst[isyst]<<endl;

            outfile<<"\t";

            if(Is_Syst_Match_Sample(v_normSyst[isyst], v_samples[isample], false)) //Other syst : check if applies to current samples (from name) //Set 'use_rph' arg to false because we want to apply lnN systs to EFT signals even if using RPHs (supported)
			{
                //-- Hard-coded special cases: e.g. if a lN syst. is correlated between years with different values, use a marker replaced with year-specific values by parsing code
                if(v_normSyst[isyst] == "Lumi1617") {outfile<<"[Lumi1617]";}
                else if(v_normSyst[isyst] == "Lumi1718") {outfile<<"[Lumi1718]";}
                else if(v_normSyst[isyst] == "LumiXY") {outfile<<"[LumiXY]";}

				else {outfile<<v_normSystValue[isyst];} //Normal cases
			}
			else {outfile<<"-";}
        }
        outfile<<endl;
    }

//                      #     #
// #       ####   ####  ##    #     ####  #   #  ####  #####
// #      #    # #    # # #   #    #       # #  #        #
// #      #    # #      #  #  #     ####    #    ####    #
// #      #    # #  ### #   # #         #   #        #   #
// #      #    # #    # #    ##    #    #   #   #    #   #
// ######  ####   ####  #     #     ####    #    ####    #

 //  ####  # #    #  ####  #      ######
 // #      # ##   # #    # #      #
 //  ####  # # #  # #      #      #####
 //      # # #  # # #  ### #      #
 // #    # # #   ## #    # #      #
  // ####  # #    #  ####  ###### ######
//--------------------------------------------

    vector<TString> v_rateParam; //Store names of samples for which we want to use a rate parameter (specified at end of datacard)

    //--- lnN Systematics applied to only 1 process (e.g. background uncert.)
    outfile<<"---------------------------------------------------"<<endl;
    for(int isample=0; isample<v_samples.size(); isample++)
    {
        if(v_isSignal[isample] == 1) {continue;} //No norm. syst for signals
		else if(v_sampleUncert[isample] == -1) //Rate param
        {
            v_rateParam.push_back(v_samples[isample]);
            continue;
        } //Don't apply lnN rate syst for some samples

		// cout<<"Sample "<<v_samples[isample]<<" / Uncert = "<<v_sampleUncert[isample]<<endl;

        outfile<<v_samples[isample] + "_rate"<<"\t"<<"lnN";

        for(int jsample=0; jsample<v_samples.size(); jsample++)
        {
            outfile<<"\t";

			if(isample == jsample) {outfile<<1.+v_sampleUncert[jsample]/100.;} //in %
            else {outfile<<"-";}
        }
        outfile<<endl;
    }


 //  ####  #    #   ##   #####  ######
 // #      #    #  #  #  #    # #
 //  ####  ###### #    # #    # #####
 //      # #    # ###### #####  #
 // #    # #    # #    # #      #
  // ####  #    # #    # #      ######

//--------------------------------------------

    //--- Shape systematics
    outfile<<"---------------------------------------------------"<<endl;
    for(int isyst=0; isyst<v_shapeSyst.size(); isyst++)
    {
        //FIXME
        // if(v_shapeSyst[isyst]=="MET" || v_shapeSyst[isyst]=="JES" || v_shapeSyst[isyst]=="JER") outfile<<"[SR]";

		//-- Markers at beginning of line :
        //[SHAPE] can be used later to easily disactivate all shape systs, at parsing
        //[201617] can be used to disactivate the prefiring syst for 2018, etc.
        outfile<<"[SHAPE]";
        if(v_shapeSyst[isyst].EndsWith("1617") || v_shapeSyst[isyst].BeginsWith("prefir")) {outfile<<"[201617]";} //Hardcoded: prefire for 16/17 only
        else if(v_shapeSyst[isyst].EndsWith("1718")) {outfile<<"[201718]";}
        else if(v_shapeSyst[isyst].EndsWith("1618")) {outfile<<"[201618]";}
        else if(v_shapeSyst[isyst].EndsWith("16")) {outfile<<"[2016]";}
        else if(v_shapeSyst[isyst].EndsWith("17")) {outfile<<"[2017]";}
        else if(v_shapeSyst[isyst].EndsWith("18")) {outfile<<"[2018]";}

        outfile<<v_shapeSyst[isyst];
        if(!v_shapeSyst_isCorrelYears[isyst]) {outfile<<"[YEAR]";} //Uncorrelated for different year --> Modify systematic name itself
        outfile<<"\t"<<"shape";

        for(int isample=0; isample<v_samples.size(); isample++)
        {
            outfile<<"\t";

            if (Is_Syst_Match_Sample(v_shapeSyst[isyst], v_samples[isample], use_rph))
            {
                if(v_samples[isample].Contains("PrivMC")) {outfile<<"[PrivMC_CR]";} //Special case: we want to disactivate systs for PrivMC samples in CRs (since nominal is dummy anyway) //in %
                else {outfile<<"1";} //in %
            }
			else {outfile<<"-";}
        }
        outfile<<endl;

        // if(v_shapeSyst[isyst] == "njets_tZq") {include_njets_syst = true;} //Obsolete
    }

    //-- If using RooDataHists (and RooParametricHists), can't use autoMCstats --> Declare 1 shape systematic for each histobin and process
    if(use_rph && signal.Contains("eft"))
    {
        for(int isample=0; isample<v_samples.size(); isample++)
        {
            if(v_samples[isample].Contains("PrivMC")) {continue;} //MC stat directly embeded in RPH parmeterizations

            TString nuis_mcstat_name = "#MCstat_[BIN]_[VAR]_[YEAR]_" + v_samples[isample]; //Ex: 'MCstat_countExp_SRttZ4l_2018_data_obs'

            outfile<<"[STAT]"<<nuis_mcstat_name<<"\t"<<"shape";
            for(int jsample=0; jsample<v_samples.size(); jsample++)
            {
                outfile<<"\t";

                //-- MCstat for given process only applies to this process
                if(isample==jsample) {outfile<<"1";}
                else {outfile<<"-";}
            }
            outfile<<endl;
        }
    }


 // #####    ##   ##### ######    #####    ##   #####    ##   #    #
 // #    #  #  #    #   #         #    #  #  #  #    #  #  #  ##  ##
 // #    # #    #   #   #####     #    # #    # #    # #    # # ## #
 // #####  ######   #   #         #####  ###### #####  ###### #    #
 // #   #  #    #   #   #         #      #    # #   #  #    # #    #
 // #    # #    #   #   ######    #      #    # #    # #    # #    #

//--------------------------------------------
//Modify normalization of any process from the datacard (e.g. FCNC)
//See : https://cms-hcomb.gitbook.io/combine/setting-up-the-analysis/preparing-the-datacard#rate-parameters
    outfile<<"---------------------------------------------------"<<endl;
	// outfile<<"[ratePar]rate_modif"<<"\t"<<"rateParam"<<"\t"<<"[VAR]_[CHAN]_[YEAR]"<<"\t"<<"sigPar"<<"\t"<<"rateVal";

    //-- May use rate parameters for some samples (specified at end of datacard)
    for(int isample=0; isample<v_rateParam.size(); isample++)
    {
        outfile<<"rate_"+v_rateParam[isample]<<"\t"<<"rateParam"<<"\t"<<"*"<<"\t"<<v_rateParam[isample]<<"\t"<<"1.0 [0.0,3.0]"<<endl;
    }


 //  ####  #####   ##   #####
 // #        #    #  #    #
 //  ####    #   #    #   #
 //      #   #   ######   #   ###
 // #    #   #   #    #   #   ###
 //  ####    #   #    #   #   ###

//--------------------------------------------
//--- MC statistical uncert.
//See : https://cms-hcomb.gitbooks.io/combine/content/part2/bin-wise-stats.html#usage-instructions
// Usage : [channel] autoMCStats [threshold] [include-signal = 0] [hist-mode = 1]
// [threshold] : A positive value sets the threshold on the effective number of unweighted events above which the uncertainty will be modeled with the Barlow-Beeston-lite approach described above. Below the threshold an individual uncertainty per-process will be created.

    outfile<<"---------------------------------------------------"<<endl;

    if(!use_rph || !signal.Contains("eft")) //Default : use autoMCstats
    {
        outfile<<"[STAT]"<<"\t"<<"*"<<"\t"<<"autoMCStats"<<"\t"<<"10"<<endl;
    }
    else if(use_rph) //Can't use autoMCstats in regions wirh rph/rdh
    {
        outfile<<"[CR][STAT]"<<"\t"<<"*"<<"\t"<<"autoMCStats"<<"\t"<<"10"<<endl;
    }


 //  ####  ##### #    # ###### #####
 // #    #   #   #    # #      #    #
 // #    #   #   ###### #####  #    #
 // #    #   #   #    # #      #####
 // #    #   #   #    # #      #   #
 //  ####    #   #    # ###### #    #

//--------------------------------------------
//-- Additional parameters

    outfile<<"---------------------------------------------------"<<endl;

    //-- If using RPHs (whose bins are directly parameterized on WCs), better to declare the WCs explicitely in the card
    if(use_rph && signal.Contains("eft"))
    {
        // if(include_njets_syst) {outfile<<"njets_tZq"<<"\t"<<"param"<<"\t"<<0<<"\t"<<1<<endl;} //Only impacts RPHs

        //-- Add explicit gaussian constraints with mean=0 and sigma=1 (this is the default constraint used for shape systematics)
        for(int isyst=0; isyst<v_shapeSyst.size(); isyst++)
        {
            outfile<<"[SHAPE]";
            if(v_shapeSyst[isyst].EndsWith("1617") || v_shapeSyst[isyst].BeginsWith("prefir")) {outfile<<"[201617]";} //Hardcoded: prefire for 16/17 only
            else if(v_shapeSyst[isyst].EndsWith("1718")) {outfile<<"[201718]";}
            else if(v_shapeSyst[isyst].EndsWith("1618")) {outfile<<"[201618]";}
            else if(v_shapeSyst[isyst].EndsWith("16")) {outfile<<"[2016]";}
            else if(v_shapeSyst[isyst].EndsWith("17")) {outfile<<"[2017]";}
            else if(v_shapeSyst[isyst].EndsWith("18")) {outfile<<"[2018]";}
            outfile<<v_shapeSyst[isyst];
            if(!v_shapeSyst_isCorrelYears[isyst]) {outfile<<"[YEAR]";} //Uncorrelated for different year --> Modify systematic name itself
            outfile<<"\t"<<"param"<<"\t"<<0<<"\t"<<1<<endl;
        }

        for(int isample=0; isample<v_samples.size(); isample++)
        {
            if(!v_samples[isample].Contains("PrivMC")) {continue;}

            TString nuis_mcstat_name = "#MCstat_[BIN]_[VAR]_[YEAR]_" + v_samples[isample]; //Ex: 'MCstat_countExp_SRttZ4l_2018_PrivMC_tZq'
            outfile<<"[STAT]"<<nuis_mcstat_name<<"\t"<<"param"<<"\t"<<0<<"\t"<<1<<endl;
        }

        outfile<<"---------------------------------------------------"<<endl;

        outfile<<"ctz"<<"\t"<<"flatParam"<<endl;
        outfile<<"ctw"<<"\t"<<"flatParam"<<endl;
        outfile<<"cpqm"<<"\t"<<"flatParam"<<endl;
        outfile<<"cpq3"<<"\t"<<"flatParam"<<endl;
        outfile<<"cpt"<<"\t"<<"flatParam"<<endl;
    }

    outfile<<"---------------------------------------------------"<<endl;

//--------------------------------------------

    cout<<endl<<endl<<"---> File ./"<<outfile_name<<" created..."<<endl<<endl<<endl;

    return;
}








//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------
//--------------------------------------------










//--------------------------------------------
// ##     ##    ###    #### ##    ##
// ###   ###   ## ##    ##  ###   ##
// #### ####  ##   ##   ##  ####  ##
// ## ### ## ##     ##  ##  ## ## ##
// ##     ## #########  ##  ##  ####
// ##     ## ##     ##  ##  ##   ###
// ##     ## ##     ## #### ##    ##
//--------------------------------------------

//Define all arguments needed by generator function (see function description for details about args)
int main()
{
    bool use_rph = true; //true <-> new, use RooParametricHist //Still debugging

//-- Read command line arguments
//--------------------------------------------
    TString signal = "";
    Choose_Arguments_From_CommandLine(signal, use_rph);


//  ####    ##   #    # #####  #      ######  ####
// #       #  #  ##  ## #    # #      #      #
//  ####  #    # # ## # #    # #      #####   ####
//      # ###### #    # #####  #      #           #
// #    # #    # #    # #      #      #      #    #
//  ####  #    # #    # #      ###### ######  ####

//-- Sample / 0 = sig, 1 = bkg / -1 = rateParam, else = lnN uncertainty
//--------------------------------------------
    vector<TString> v_samples; vector<int> v_isSignal; vector<float> v_sampleUncert;
    if(signal == "efttzq") //Signal : SMEFT tZq
    {
        v_samples.push_back("PrivMC_tZq"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("ttZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
        v_samples.push_back("tWZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    else if(signal == "eftttz") //Signal : SMEFT ttZ
    {
        v_samples.push_back("PrivMC_ttZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("tZq"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
        v_samples.push_back("tWZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    else if(signal == "efttwz") //Signal : SMEFT tWZ
    {
        v_samples.push_back("PrivMC_tWZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("tZq"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
        v_samples.push_back("ttZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    else if(signal == "eft") //Signals : SMEFT tZq+ttZ+tWZ
    {
        v_samples.push_back("PrivMC_tZq"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("PrivMC_ttZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("PrivMC_tWZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
    }
    else if(signal == "0") //Signals : SM tZq+ttZ+tWZ
    {
        v_samples.push_back("tZq"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("ttZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("tWZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
    }
    else if(signal == "tzq") //Signal : SM tZq
    {
        v_samples.push_back("tZq"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("ttZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
        v_samples.push_back("tWZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    else if(signal == "ttz") //Signal : SM ttZ
    {
        v_samples.push_back("ttZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("tZq"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
        v_samples.push_back("tWZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    else if(signal == "twz") //Signal : SM tWZ
    {
        v_samples.push_back("tWZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("ttZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
        v_samples.push_back("tZq"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    //else {cout<<FRED("Wrong arg ! Abort !")<<endl; return 0;}

    v_samples.push_back("tX"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);
    v_samples.push_back("WZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    v_samples.push_back("VVV"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    v_samples.push_back("XG"); v_isSignal.push_back(0); v_sampleUncert.push_back(10);
    v_samples.push_back("NPL"); v_isSignal.push_back(0); v_sampleUncert.push_back(30);


 //               #     #
 // #      #    # ##    #     ####  #   #  ####  #####
 // #      ##   # # #   #    #       # #  #        #
 // #      # #  # #  #  #     ####    #    ####    #
 // #      #  # # #   # #         #   #        #   #
 // #      #   ## #    ##    #    #   #   #    #   #
 // ###### #    # #     #     ####    #    ####    #

//lnN systematics
//Write "(1+X)%". E.g for lnN symmetric of 10% => "1.10"
//For a 5%/10% lnN asymmetric syst, write : "1.05/1.10"
//-1 <-> values must be hardcoded (to allow for correlations with different values per year)
//--------------------------------------------
    vector<TString> v_normSyst; vector<TString> v_normSystValue;
    v_normSyst.push_back("Lumi16"); v_normSystValue.push_back("1.022");
    v_normSyst.push_back("Lumi17"); v_normSystValue.push_back("1.020");
    v_normSyst.push_back("Lumi18"); v_normSystValue.push_back("1.015");
    v_normSyst.push_back("Lumi1617"); v_normSystValue.push_back("-1");
    v_normSyst.push_back("Lumi1718"); v_normSystValue.push_back("-1");
    v_normSyst.push_back("LumiXY"); v_normSystValue.push_back("-1");

    v_normSyst.push_back("Trigger16"); v_normSystValue.push_back("1.02");
    v_normSyst.push_back("Trigger17"); v_normSystValue.push_back("1.02");
    v_normSyst.push_back("Trigger18"); v_normSystValue.push_back("1.02");

    v_normSyst.push_back("N_CRWZ_extrap"); v_normSystValue.push_back("1.05");
    v_normSyst.push_back("N_CRZZ_extrap"); v_normSystValue.push_back("1.05");
    v_normSyst.push_back("N_CRDY_extrap"); v_normSystValue.push_back("1.05");


//  ####  #    #   ##   #####  ######     ####  #   #  ####  #####
// #      #    #  #  #  #    # #         #       # #  #        #
//  ####  ###### #    # #    # #####      ####    #    ####    #
//      # #    # ###### #####  #              #   #        #   #
// #    # #    # #    # #      #         #    #   #   #    #   #
//  ####  #    # #    # #      ######     ####    #    ####    #

//-- Name of shape systematic / Whether syst is correlated between years (single nuisance) or not (1 per year)
//--------------------------------------------
    vector<TString> v_shapeSyst; vector<bool> v_shapeSyst_isCorrelYears;
    v_shapeSyst.push_back("PU"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("prefire"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("BtagHF"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("BtagLF"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("BtagHFstats1"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("BtagHFstats2"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("BtagLFstats1"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("BtagLFstats2"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("BtagCFerr1"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("BtagCFerr2"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("jetPUIDEff"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("jetPUIDMT"); v_shapeSyst_isCorrelYears.push_back(true);

    v_shapeSyst.push_back("FRm_norm"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("FRm_pt"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("FRm_be"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("FRe_norm"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("FRe_pt"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("FRe_be"); v_shapeSyst_isCorrelYears.push_back(true);

    v_shapeSyst.push_back("LepEff_muLoose"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("LepEff_muTight"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("LepEff_elLoose"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("LepEff_elTight"); v_shapeSyst_isCorrelYears.push_back(true);

    v_shapeSyst.push_back("njets_tZq"); v_shapeSyst_isCorrelYears.push_back(true);

    v_shapeSyst.push_back("PDF"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("alphas"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("MEtZq"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("MEttZ"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("MEtWZ"); v_shapeSyst_isCorrelYears.push_back(true);

    //FIXME -- TESTING
    // v_shapeSyst.push_back("JES"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("JER"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("MET"); v_shapeSyst_isCorrelYears.push_back(false);

    //-- JEC split sources //FIXME -- TESTING
    v_shapeSyst.push_back("AbsoluteStat"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("AbsoluteScale"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("AbsoluteMPFBias"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("Fragmentation"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("SinglePionECAL"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("SinglePionHCAL"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("FlavorQCD"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("TimePtEta"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("RelativeJEREC1"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("RelativeJEREC2"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("RelativeJERHF"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("RelativePtBB"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("RelativePtEC1"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("RelativePtEC2"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("RelativePtHF"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("RelativeBal"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("RelativeSample"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("RelativeFSR"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("RelativeStatFSR"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("RelativeStatEC"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("RelativeStatHF"); v_shapeSyst_isCorrelYears.push_back(true); //false
    v_shapeSyst.push_back("PileUpDataMC"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("PileUpPtRef"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("PileUpPtBB"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("PileUpPtEC1"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("PileUpPtEC2"); v_shapeSyst_isCorrelYears.push_back(true); //true
    v_shapeSyst.push_back("PileUpPtHF"); v_shapeSyst_isCorrelYears.push_back(true); //true

    //-- Missing
    // v_shapeSyst.push_back("ISRtZq"); v_shapeSyst_isCorrelYears.push_back(true);
    // v_shapeSyst.push_back("ISRttZ"); v_shapeSyst_isCorrelYears.push_back(true);
    // v_shapeSyst.push_back("ISRtWZ"); v_shapeSyst_isCorrelYears.push_back(true);
    // v_shapeSyst.push_back("FSR"); v_shapeSyst_isCorrelYears.push_back(true);

    //-- Obsolete
    // v_shapeSyst.push_back("FakeFactor"); v_shapeSyst_isCorrelYears.push_back(true);


//  ####    ##   #      #       ####
// #    #  #  #  #      #      #
// #      #    # #      #       ####
// #      ###### #      #           #
// #    # #    # #      #      #    #
//  ####  #    # ###### ######  ####

//Function calls
//--------------------------------------------
    //Generate the template datacard
    Generate_Datacard(v_samples, v_isSignal, v_sampleUncert, v_normSyst, v_normSystValue, v_shapeSyst, signal, v_shapeSyst_isCorrelYears, false);

    if(use_rph) {Generate_Datacard(v_samples, v_isSignal, v_sampleUncert, v_normSyst, v_normSystValue, v_shapeSyst, signal, v_shapeSyst_isCorrelYears, true, "Template_Datacard_RPH.txt");}

    return 0;
}
