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
bool Is_Syst_Match_Sample(TString syst, TString sample)
{
	// cout<<"syst "<<syst<<endl;
	// cout<<"sample "<<sample<<endl;

	// if(syst.Contains("ttH", TString::kIgnoreCase) && !sample.Contains("ttH")) {return false;}

	return true;
}

//Ask user to choose options at command line for script generation
void Choose_Arguments_From_CommandLine(TString& signal)
{
    //Choose whether to include shape syst or not
	cout<<endl<<FYEL("--- What is your SIGNAL ?")<<endl;
    cout<<"* 'eft'   \t<-> Signals are SMEFT tZq+ttZ"<<endl;
    cout<<"* 'efttzq'   <-> Signal is SMEFT tZq only"<<endl;
    cout<<"* 'eftttz'   <-> Signal is SMEFT ttZ only"<<endl;
    cout<<"* '0'   <-> Signals are tZq + ttZ"<<endl;
    // cout<<"* 'thq' <-> Signals are tHq + tHW"<<endl;
    cout<<"* 'tzq' <-> Signal is tZq"<<endl;
    cout<<"* 'ttz' <-> Signal is ttZ"<<endl;
	cin>>signal;
	while(signal != "tzq" && signal != "ttz" && signal != "0" && signal != "eft" && signal != "efttzq" && signal != "eftttz")
	{
		cin.clear();
		cin.ignore(1000, '\n');

		cout<<" Wrong answer ! Retry :"<<endl;
		cin>>signal;
	}

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
void Generate_Datacard(vector<TString> v_samples, vector<int> v_isSignal, vector<float> v_sampleUncert, vector<TString> v_normSyst, vector<TString> v_normSystValue, vector<TString> v_shapeSyst, TString signal, vector<bool> v_normSyst_isCorrelYears, vector<bool> v_shapeSyst_isCorrelYears)
{
    TString outfile_name = "Template_Datacard.txt";
    ofstream outfile(outfile_name.Data());

 // #    # ######   ##   #####  ###### #####
 // #    # #       #  #  #    # #      #    #
 // ###### #####  #    # #    # #####  #    #
 // #    # #      ###### #    # #      #####
 // #    # #      #    # #    # #      #   #
 // #    # ###### #    # #####  ###### #    #
//--------------------------------------------

    //--- Nof observables, bkgs, nuisance params
    outfile<<"imax"<<"\t"<<1<<"\t"<<"number of categories"<<endl;
    outfile<<"jmax"<<"\t"<<"*"<<"\t"<<"number of backgrounds"<<endl;
    outfile<<"kmax"<<"\t"<<"*"<<"\t"<<"number of nuisance parameters (sources of systematic uncertainties)"<<endl;

//--------------------------------------------
    //--- Filepath, naming convention
    outfile<<"---------------------------------------------------"<<endl;
    outfile<<"shapes * [VAR]_[CHAN]_[YEAR] filetoread $CHANNEL__$PROCESS $CHANNEL__$PROCESS__$SYSTEMATIC"<<endl;

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
        outfile<<-1;
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
        outfile<<v_normSyst[isyst];
        if(!v_normSyst_isCorrelYears[isyst]) {outfile<<"[YEAR]";} //Uncorrelated for each year
        outfile<<"\t"<<"lnN";

        for(int isample=0; isample<v_samples.size(); isample++)
        {
			// cout<<"sample "<<v_samples[isample]<<endl;
			// cout<<"syst "<<v_normSyst[isyst]<<endl;

            outfile<<"\t";

            if(Is_Syst_Match_Sample(v_normSyst[isyst], v_samples[isample])) //Other syst : check if applies to current samples (from name)
			{
				outfile<<v_normSystValue[isyst];
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

    //--- lnN Systematics applied to only 1 process (e.g. background uncert.)
    outfile<<"---------------------------------------------------"<<endl;
    for(int isample=0; isample<v_samples.size(); isample++)
    {
        // if(v_isSignal[isample] == 1) {continue;} //No norm. syst for signals
		if(v_sampleUncert[isample] == -1) {continue;} //Don't apply lnN rate syst for some samples (ttZ, Fakes, etc.)

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
		//Markers at beginning of line :
        //-- the [SHAPE] symbol can be used later to easily disactivate all shape systs, at parsing
        //-- idem, [201617] can be used to disactivate the prefiring syst for 2018 !
        outfile<<"[SHAPE]";
        if(v_shapeSyst[isyst].BeginsWith("prefir")) {outfile<<"[201617]";}
        outfile<<v_shapeSyst[isyst]; //the [SHAPE] symbol can be used later to easily disactivate all shape systs, at parsing
        if(!v_shapeSyst_isCorrelYears[isyst]) {outfile<<"[YEAR]";} //Uncorrelated for each year
        outfile<<"\t"<<"shape";

        for(int isample=0; isample<v_samples.size(); isample++)
        {
            outfile<<"\t";

            if (Is_Syst_Match_Sample(v_shapeSyst[isyst], v_samples[isample])) {outfile<<"1";} //in %
			else {outfile<<"-";}
        }
        outfile<<endl;
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

	//The [STAT] symbol can be used later to easily disactivate all shape systs, at parsing
	if(signal == "0") {outfile<<"[STAT]"<<"\t"<<"*"<<"\t"<<"autoMCStats"<<"\t"<<"0 0 1"<<endl;} //do as THQ for now...
	else {outfile<<"[STAT]"<<"\t"<<"*"<<"\t"<<"autoMCStats"<<"\t"<<"10"<<endl;}

//--------------------------------------------

    cout<<endl<<endl<<"---> File ./Template_datacard.txt created..."<<endl<<endl<<endl;

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

//Read command line arguments
//--------------------------------------------
    TString signal = "";
    Choose_Arguments_From_CommandLine(signal);


//  ####    ##   #    # #####  #      ######  ####
// #       #  #  ##  ## #    # #      #      #
//  ####  #    # # ## # #    # #      #####   ####
//      # ###### #    # #####  #      #           #
// #    # #    # #    # #      #      #      #    #
//  ####  #    # #    # #      ###### ######  ####
//--------------------------------------------

    vector<TString> v_samples; vector<int> v_isSignal; vector<float> v_sampleUncert;
    if(signal == "efttzq") //Signal : tZq
    {
        v_samples.push_back("PrivMC_tZq"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("ttZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    else if(signal == "eftttz") //Signal : tZq
    {
        v_samples.push_back("PrivMC_ttZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("tZq"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    else if(signal == "eft") //Signals : tZq+ttZ
    {
        v_samples.push_back("PrivMC_tZq"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("PrivMC_ttZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
    }
    else if(signal == "0") //Signals : tZq+ttZ
    {
        v_samples.push_back("tZq"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("ttZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
    }
    else if(signal == "tzq") //Signals : tZq
    {
        v_samples.push_back("tZq"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("ttZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    else if(signal == "ttz") //Signals : ttZ
    {
        v_samples.push_back("ttZ"); v_isSignal.push_back(1); v_sampleUncert.push_back(-1);
        v_samples.push_back("tZq"); v_isSignal.push_back(0); v_sampleUncert.push_back(15);
    }
    else {cout<<FRED("Wrong arg ! Abort !")<<endl; return 0;}

    v_samples.push_back("tX"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);
    v_samples.push_back("VVV"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);
    v_samples.push_back("NPL"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);
    v_samples.push_back("WZ"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);
    v_samples.push_back("Xg"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);

    // v_samples.push_back("ttX"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);
    // v_samples.push_back("tX"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);
    // v_samples.push_back("VV"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);
    // v_samples.push_back("DY"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);
    // v_samples.push_back("TTbar_DiLep"); v_isSignal.push_back(0); v_sampleUncert.push_back(20);



 //               #     #
 // #      #    # ##    #     ####  #   #  ####  #####
 // #      ##   # # #   #    #       # #  #        #
 // #      # #  # #  #  #     ####    #    ####    #
 // #      #  # # #   # #         #   #        #   #
 // #      #   ## #    ##    #    #   #   #    #   #
 // ###### #    # #     #     ####    #    ####    #

//lnN systematics. Write "1+X%". E.g for lnN symmetric of 10% => "1.10"
//For a 5%/10% lnN asymmetric syst, write : "1.05/1.10"
//NB : if syst contains e.g. "ttH" in the name, it will only apply to ttH sample ! hardcoded for thz "pdf" systs
//--------------------------------------------
    vector<TString> v_normSyst; vector<TString> v_normSystValue; vector<bool> v_normSyst_isCorrelYears;
    v_normSyst.push_back("Lumi"); v_normSystValue.push_back("1.023"); v_normSyst_isCorrelYears.push_back(false);

	//-- Combine convention : "x/y" for asymmetric norm errors
	// v_normSyst.push_back("QCDscale_ttH"); v_normSystValue.push_back("0.907/1.058");
    // v_normSyst.push_back("Clos_e_bt_norm"); v_normSystValue.push_back("Closebtnorm");


//  ####  #    #   ##   #####  ######     ####  #   #  ####  #####
// #      #    #  #  #  #    # #         #       # #  #        #
//  ####  ###### #    # #    # #####      ####    #    ####    #
//      # #    # ###### #####  #              #   #        #   #
// #    # #    # #    # #      #         #    #   #   #    #   #
//  ####  #    # #    # #      ######     ####    #    ####    #

//--------------------------------------------
    vector<TString> v_shapeSyst; vector<bool> v_shapeSyst_isCorrelYears;
    // v_shapeSyst.push_back("PU"); v_shapeSyst_isCorrelYears.push_back(true);
    v_shapeSyst.push_back("prefiringWeight"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("BtagH"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("BtagL"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("LepEff_muLoose"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("LepEff_muTight"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("LepEff_elLoose"); v_shapeSyst_isCorrelYears.push_back(false);
    v_shapeSyst.push_back("LepEff_elTight"); v_shapeSyst_isCorrelYears.push_back(false);


//  ####    ##   #      #       ####
// #    #  #  #  #      #      #
// #      #    # #      #       ####
// #      ###### #      #           #
// #    # #    # #      #      #    #
//  ####  #    # ###### ######  ####

//Function calls
//--------------------------------------------
    //Generate the template datacards for ele and mu channels (some differences)
    Generate_Datacard(v_samples, v_isSignal, v_sampleUncert, v_normSyst, v_normSystValue, v_shapeSyst, signal, v_normSyst_isCorrelYears, v_shapeSyst_isCorrelYears);

    return 0;
}
