/* BASH CUSTOM */
#define RST   "\e[0m"
#define KBLK  "\e[30m"
#define KRED  "\e[31m"
#define KGRN  "\e[32m"
#define KYEL  "\e[33m"
#define KBLU  "\e[34m"
#define KMAG  "\e[35m"
#define KCYN  "\e[36m"
#define KGRAY  "\e[37m"
#define KWHT  "\e[97m"
#define KLRED  "\e[91m"
#define KLBLU  "\e[94m"
#define KBRED  "\e[41m"
#define KBGRN  "\e[42m"
#define KBYEL  "\e[43m"
#define KBBLU  "\e[44m"
#define KBCYN  "\e[46m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FGRAY(x) KGRAY x RST
#define FWHT(x) KWHT x RST
#define FLRED(x) KLRED x RST
#define FLBLU(x) KLBLU x RST
#define BRED(x) KBRED x RST
#define BGRN(x) KBGRN x RST
#define BYEL(x) KBYEL x RST
#define BBLU(x) KBBLU x RST
#define BCYN(x) KBCYN x RST
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

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TTree.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TString.h"
#include "TColor.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLine.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TObject.h"
#include "TRandom3.h"
#include "TRandom1.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TLegendEntry.h"
#include "TGaxis.h"
#include "TLeaf.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include <iostream>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <cassert>     //Can be used to terminate program if argument is not true.
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

//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Existence(const TString& name)
{
    struct stat buffer;
    return (stat (name.Data(), &buffer) == 0); //true if file exists
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

/**
 * Return the name of the directory where to store the corresponding subsample
 * NB : training events (looser selection) are also stored under the "SR" directory
 */
TString Get_Directory(TString cat, TString sample, TString analysis_type)
{
	TString ts_lep = "", ts_cat = "";

    if(cat.Contains("tZq") ) {ts_cat = "SR_tZq";}
    else if(cat.Contains("ttZ") ) {ts_cat = "SR_ttZ";}
    else if(cat.Contains("signal") ) {ts_cat = "SR";}
    else {return "";}

	return ts_cat;
}







//--------------------------------------------
//  ######   #######  ########  ##    ##     #######  ##    ## ########
// ##    ## ##     ## ##     ##  ##  ##     ##     ## ###   ## ##
// ##       ##     ## ##     ##   ####      ##     ## ####  ## ##
// ##       ##     ## ########     ##       ##     ## ## ## ## ######
// ##       ##     ## ##           ##       ##     ## ##  #### ##
// ##    ## ##     ## ##           ##       ##     ## ##   ### ##
//  ######   #######  ##           ##        #######  ##    ## ########

//  ######  ##     ## ########           ######     ###    ##     ## ########  ##       ########
// ##    ## ##     ## ##     ##         ##    ##   ## ##   ###   ### ##     ## ##       ##
// ##       ##     ## ##     ##         ##        ##   ##  #### #### ##     ## ##       ##
//  ######  ##     ## ########  #######  ######  ##     ## ## ### ## ########  ##       ######
//       ## ##     ## ##     ##               ## ######### ##     ## ##        ##       ##
// ##    ## ##     ## ##     ##         ##    ## ##     ## ##     ## ##        ##       ##
//  ######   #######  ########           ######  ##     ## ##     ## ##        ######## ########
//--------------------------------------------

/*
Create new subsample based on full sample, copying only events satisfying one selection at a time
==> Faster processing of ntuples
 */
void Create_Subsample_fromSample(TString fullsample_name, TString newsample_name, TString selection, TString sample, TString open_mode, vector<TString> v_TTrees)
{
    cout<<endl<<BYEL("                                  ")<<endl;
    cout<<FYEL("--- WILL EXTRACT SUBSET OF EVENTS SATISFYING '"<<selection<<"'")<<endl;
    cout<<FYEL("FROM ["<<fullsample_name<<"]")<<endl;
    cout<<FYEL("--> TO ["<<newsample_name<<"] ---")<<endl;

    //Input file
    if(!Check_File_Existence(fullsample_name) ) {cout<<FRED("Sample "<<fullsample_name<<" not found !")<<endl; return;}
    TFile* f_data = new TFile(fullsample_name);

    //Output file
    TFile *f_new = 0;
    f_new = new TFile(newsample_name, open_mode);

    for(int itree=0; itree<v_TTrees.size(); itree++)
    {
        cout<<"-- Tree : '"<<v_TTrees[itree]<<"'..."<<endl;
        if(itree != 0 && (!selection.Contains("SR") || newsample_name.Contains("DATA")) ) {continue;} //Don't need JES/JER for Fakes/QFlip/Train/... //-- also keep JES for training ! (for MEM, ...)

		TTree* t_data = 0; //input tree
        t_data = (TTree*) f_data->Get(v_TTrees[itree]);
        if(!t_data) {cout<<FRED("Null tree "<<v_TTrees[itree]<<" !")<<endl; return;}
        Int_t nentries = t_data->GetEntries();

    	f_new->cd();
        TTree *t_new = 0; //output tree
        t_new = (TTree*) t_data->CopyTree(selection);
        if(!t_new) {cout<<FRED("t_new is NULL !")<<endl; return;}

		//Not sure why, but seems necessary for tHq sample from MEM... else all called "tree" (metadata remaining ?)
		// if(v_TTrees[itree] == "JESUp") {t_new->SetName("JESUp");}
		// else if(v_TTrees[itree] == "JESDown") {t_new->SetName("JESDown");}

        t_new->Write();

        cout<<"("<<t_new->GetEntries()<<" entries)"<<endl;

        delete t_data; t_data = 0;
        delete t_new; t_new = 0;
    }

    f_data->Close();
    f_new->Close();

    cout<<BYEL("                                  ")<<endl<<endl;

    return;
}








//--------------------------------------------
//  ######   #######  ########  ##    ##     ######  ##     ## ##     ## ##      ## ######## ####  ######   ##     ## ########
// ##    ## ##     ## ##     ##  ##  ##     ##    ## ##     ## ###   ### ##  ##  ## ##        ##  ##    ##  ##     ##    ##
// ##       ##     ## ##     ##   ####      ##       ##     ## #### #### ##  ##  ## ##        ##  ##        ##     ##    ##
// ##       ##     ## ########     ##        ######  ##     ## ## ### ## ##  ##  ## ######    ##  ##   #### #########    ##
// ##       ##     ## ##           ##             ## ##     ## ##     ## ##  ##  ## ##        ##  ##    ##  ##     ##    ##
// ##    ## ##     ## ##           ##       ##    ## ##     ## ##     ## ##  ##  ## ##        ##  ##    ##  ##     ##    ##
//  ######   #######  ##           ##        ######   #######  ##     ##  ###  ###  ######## ####  ######   ##     ##    ##
//--------------------------------------------

/**
 * Save the histograms containing the sum of weights (for event rescaling, scale variation studies, etc.) into the split sample too !
 */
void Copy_SumWeight_Histogram_Into_SplitSample(TString filepath, TString outfile_path, TString sample, TString analysis_type)
{
    //Input file
    if(!Check_File_Existence(filepath) ) {cout<<FRED("Sample "<<filepath<<" not found !")<<endl; return;}
    TFile* f_input = TFile::Open(filepath);

    //Output file
    TFile *f_output = 0;
    if(!Check_File_Existence(outfile_path) ) {cout<<FRED("Sample "<<outfile_path<<" not found !")<<endl; return;}
    f_output = TFile::Open(outfile_path, "UPDATE");

    TH1D* EFT_SumWeights = (TH1D*) f_input->Get("EFT_SumWeights");

    f_output->cd();
    if(!f_input->GetListOfKeys()->Contains("EFT_SumWeights") ) {cout<<"Histogram EFT_SumWeights not found in file "<<filepath<<" !"<<endl; return;}
    else {EFT_SumWeights->Write("EFT_SumWeights");}

    f_input->Close();
    f_output->Close();

    return;
}










//--------------------------------------------
//  ######   ######## ##    ## ######## ########     ###    ##       #### ########    ###    ######## ####  #######  ##    ##
// ##    ##  ##       ###   ## ##       ##     ##   ## ##   ##        ##       ##    ## ##      ##     ##  ##     ## ###   ##
// ##        ##       ####  ## ##       ##     ##  ##   ##  ##        ##      ##    ##   ##     ##     ##  ##     ## ####  ##
// ##   #### ######   ## ## ## ######   ########  ##     ## ##        ##     ##    ##     ##    ##     ##  ##     ## ## ## ##
// ##    ##  ##       ##  #### ##       ##   ##   ######### ##        ##    ##     #########    ##     ##  ##     ## ##  ####
// ##    ##  ##       ##   ### ##       ##    ##  ##     ## ##        ##   ##      ##     ##    ##     ##  ##     ## ##   ###
//  ######   ######## ##    ## ######## ##     ## ##     ## ######## #### ######## ##     ##    ##    ####  #######  ##    ##
//--------------------------------------------
/*
Automatically split all samples into different subsamples (based on categories) *
* --> Call to Create_Subsample_fromSample() for all samples/selections
 */
void Split_AllNtuples_ByCategory(vector<TString> v_samples, vector<TString> v_sel, vector<TString> v_years, TString analysis_type, bool make_nominal_samples, vector<TString> v_TTrees)
{
	if(analysis_type.Contains("njet") || analysis_type.Contains("fwd")) {return;}

    cout<<endl<<endl<<FBLU("== START OF NTUPLES SPLITTING ==")<<endl;
    cout<<"This can be quite long. Make sure you have correctly selected in the code :"<<endl;
    cout<<"- The samples & categories to split"<<endl;
    cout<<"- The TTrees you want to split : Nominal/JES/JER/... TTrees"<<endl;
    cout<<"---------------------------"<<endl<<endl<<endl;
    usleep(4000000); //Pause for 4s

    TString prefix = "./input_ntuples/";

    for(int iyear=0; iyear<v_years.size(); iyear++)
    {
        TString full_dirname = prefix + "SR/";
        mkdir(full_dirname.Data(), 0777);
        full_dirname+= v_years[iyear];
        mkdir(full_dirname.Data(), 0777);

        for(int isel=0; isel<v_sel.size(); isel++)
        {
            for(int isample=0; isample<v_samples.size(); isample++)
            {
    			// cout<<"prefix = "<<prefix<<endl;
                TString filepath = prefix + v_years[iyear] + "/" + v_samples[isample] + ".root";
    			if( !Check_File_Existence(filepath) ) {cout<<BOLD(FRED("Sample "<<filepath<<" not found !") )<<endl; continue;}

    			TString dir = Get_Directory(v_sel[isel], v_samples[isample], analysis_type);
    			if(dir == "") {continue;}

    			TString outfile_path = prefix + dir + "/" + v_years[iyear] + "/" + v_samples[isample] + ".root";

                // cout<<"outfile_path = "<<outfile_path<<endl;

    			TString selection = v_sel[isel] + "==1";

                TString open_mode = "RECREATE";
            	Create_Subsample_fromSample(filepath, outfile_path, selection, v_samples[isample], open_mode, v_TTrees);
            	if(v_samples[isample].Contains("PrivMC")) {Copy_SumWeight_Histogram_Into_SplitSample(filepath, outfile_path, v_samples[isample], analysis_type);}

            } //sample loop
    	} //selections loop
    } //year loop

    // Merge_Many_TTrees_Into_One(v_sel, v_samples, v_TTrees, prefix, analysis_type);

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

int main(int argc, char **argv)
{
//--- Options---------------------------------
//--------------------------------------------
    //GammaConv and Fakes_MC require running over all MC samples, but only for the relevant categories
    bool make_nominal_samples = true; //Keep =true for now

    bool tree_name = "result";

    TString analysis_type = ""; //obsolete
//--------------------------------------------


 //  ####  #   #  ####  #####    ##### ##### #####  ###### ######  ####
 // #       # #  #        #        #     #   #    # #      #      #
 //  ####    #    ####    #        #     #   #    # #####  #####   ####
 //      #   #        #   #        #     #   #####  #      #           #
 // #    #   #   #    #   #        #     #   #   #  #      #      #    #
//  ####    #    ####    #        #     #   #    # ###### ######  ####

    //-- Copy nominal + JES/JER TTrees
    //NB: first must be nominal tree name
    vector<TString> v_TTrees;
    v_TTrees.push_back("result");
    // v_TTrees.push_back("JESDown"); v_TTrees.push_back("JESUp");
    // v_TTrees.push_back("JERDown"); v_TTrees.push_back("JERUp"); //not used by ttH

 //  ####    ##   #    # #####  #      ######  ####
 // #       #  #  ##  ## #    # #      #      #
 //  ####  #    # # ## # #    # #      #####   ####
 //      # ###### #    # #####  #      #           #
 // #    # #    # #    # #      #      #      #    #
 //  ####  #    # #    # #      ###### ######  ####

    //--- Sample list
    vector<TString> v_samples;
	v_samples.push_back("DATA");
    v_samples.push_back("PrivMC_tZq");
    v_samples.push_back("PrivMC_ttZ");
    v_samples.push_back("tZq");
    v_samples.push_back("ttZ");
    v_samples.push_back("tWZ");
	v_samples.push_back("tHq");
    v_samples.push_back("tHW");
    v_samples.push_back("ttH");
    v_samples.push_back("ttW");
    v_samples.push_back("ttZZ");
    v_samples.push_back("ttWW");
    v_samples.push_back("ttWZ");
    v_samples.push_back("ttZH");
    v_samples.push_back("ttWH");
    v_samples.push_back("tttt");
    v_samples.push_back("ttHH");
    v_samples.push_back("ZZ4l");
    v_samples.push_back("ggToZZTo4l");
    v_samples.push_back("ZZZ");
    v_samples.push_back("WZZ");
    v_samples.push_back("WWW");
    v_samples.push_back("WWZ");
    v_samples.push_back("WZ");
	v_samples.push_back("TTGamma_Dilep");
	v_samples.push_back("tGJets");
	v_samples.push_back("WGToLNuG");
	v_samples.push_back("ZGToLLG_01J");
	v_samples.push_back("DY");
    v_samples.push_back("TTbar_DiLep");


 //  ####  ###### #      ######  ####  ##### #  ####  #    #  ####
 // #      #      #      #      #    #   #   # #    # ##   # #
 //  ####  #####  #      #####  #        #   # #    # # #  #  ####
 //      # #      #      #      #        #   # #    # #  # #      #
 // #    # #      #      #      #    #   #   # #    # #   ## #    #
 //  ####  ###### ###### ######  ####    #   #  ####  #    #  ####

    //--- Will divide samples based on categories
    vector<TString> v_sel;
    v_sel.push_back("is_signal_SR");


 // #   # ######   ##   #####   ####
 //  # #  #       #  #  #    # #
 //   #   #####  #    # #    #  ####
 //   #   #      ###### #####       #
 //   #   #      #    # #   #  #    #
 //   #   ###### #    # #    #  ####

    //--- Define the data-taking years
    vector<TString> v_years;
    v_years.push_back("2016");
    v_years.push_back("2017");
    v_years.push_back("2018");


 // ###### #    # #    #  ####      ####    ##   #      #
 // #      #    # ##   # #    #    #    #  #  #  #      #
 // #####  #    # # #  # #         #      #    # #      #
 // #      #    # #  # # #         #      ###### #      #
 // #      #    # #   ## #    #    #    # #    # #      #
 // #       ####  #    #  ####      ####  #    # ###### ######

    //-- Make split ntuples per sub-category
    Split_AllNtuples_ByCategory(v_samples, v_sel, v_years, analysis_type, make_nominal_samples, v_TTrees);

    return 0;
}
