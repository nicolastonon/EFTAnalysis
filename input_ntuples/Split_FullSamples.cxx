//Could merge full ntuples per groups, then delete individual ntuples locally ?
//Could run this code on NAF ? And only DL final ntuples ? must activate corresponding bool in AN code (run on groups, ...)

//-- by Nicolas Tonon (DESY) --//

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

#include "../Utils/Helper.h"

//Custom classes for EFT (see https://github.com/Andrew42/EFTGenReader/blob/maste)
// #include "../Utils/TH1EFT.h"
#include "../Utils/WCPoint.h"
#include "../Utils/WCFit.h"

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

/**
 * Return the name of the directory where to store the corresponding subsample
 * NB : training events (looser selection) are also stored under the "SR" directory
 */
TString Get_Directory(TString cat)
{
    cat.ToLower(); //Case-insensitive comparisons
	TString ts_cat = ""; //Return value

    if(cat.Contains("tzq") ) {ts_cat = "SR_tZq";}
    else if(cat.Contains("ttz") ) {ts_cat = "SR_ttZ";}
    else if(cat.Contains("signal") ) {ts_cat = "SR";}
    else if(cat.Contains("_vg") ) {ts_cat = "CR_vg";}
    else if(cat.Contains("_zz") ) {ts_cat = "CR_zz";}
    else if(cat.Contains("_tx") ) {ts_cat = "CR_tz";}
    else if(cat.Contains("_tt") ) {ts_cat = "CR_tt";}
    else if(cat.Contains("_wz") ) {ts_cat = "CR_wz";}
    else if(cat.Contains("_dy") ) {ts_cat = "CR_dy";}

    else {cout<<BOLD(FRED("ERROR: category "<<cat<<" not recognized !"))<<endl;; return "";}

	return ts_cat + "/";
}


//--------------------------------------------
// ######## ######## ########    ########     ###    ########     ###    ##     ## ######## ######## ######## ########  #### ########    ###    ######## ####  #######  ##    ##
// ##       ##          ##       ##     ##   ## ##   ##     ##   ## ##   ###   ### ##          ##    ##       ##     ##  ##       ##    ## ##      ##     ##  ##     ## ###   ##
// ##       ##          ##       ##     ##  ##   ##  ##     ##  ##   ##  #### #### ##          ##    ##       ##     ##  ##      ##    ##   ##     ##     ##  ##     ## ####  ##
// ######   ######      ##       ########  ##     ## ########  ##     ## ## ### ## ######      ##    ######   ########   ##     ##    ##     ##    ##     ##  ##     ## ## ## ##
// ##       ##          ##       ##        ######### ##   ##   ######### ##     ## ##          ##    ##       ##   ##    ##    ##     #########    ##     ##  ##     ## ##  ####
// ##       ##          ##       ##        ##     ## ##    ##  ##     ## ##     ## ##          ##    ##       ##    ##   ##   ##      ##     ##    ##     ##  ##     ## ##   ###
// ######## ##          ##       ##        ##     ## ##     ## ##     ## ##     ## ########    ##    ######## ##     ## #### ######## ##     ##    ##    ####  #######  ##    ##
//--------------------------------------------

//For private SMEFT samples, can compute and store the per-event EFT parameterizations directly (then will only need to read it when processing the samples)
void Store_EFTparameterization(TString filepath, vector<TString> v_TTrees, TString nominal_tree_name)
{
    if(!filepath.Contains("PrivMC") || filepath.Contains("_c")) {return;} //For private SMEFT samples only

    TFile* f = new TFile(filepath, "UPDATE");

    for(int itree=0; itree<v_TTrees.size(); itree++)
    {
        if(v_TTrees[itree] != nominal_tree_name) {cout<<DIM("Skipping syst. tree "<<v_TTrees[itree]<<" !")<<endl; continue;} //Don't need to store parameterization (very slow) for systematics

        TTree* t = (TTree*) f->Get(v_TTrees[itree]);
        if(!t) {cout<<BOLD(FRED("ERROR: tree "<<v_TTrees[itree]<<" not found (filepath: "<<filepath<<") ! Skip !"))<<endl; continue;}
        int nentries = t->GetEntries();

        //-- Enforce naming convention //Obsolete in next ntuple prod
        TString systTree_name = v_TTrees[itree];
        if(v_TTrees[itree] == "TotalDown") {systTree_name = "JESDown";}
        else if(v_TTrees[itree] == "TotalUp") {systTree_name = "JESUp";}
        else if(v_TTrees[itree] == "UnclEnDown") {systTree_name = "METDown";}
        else if(v_TTrees[itree] == "UnclEnUp") {systTree_name = "METUp";}

        TString newTree_name = "EFTparameterization";
        if(v_TTrees[itree] != nominal_tree_name) {newTree_name+= "_" +systTree_name;}

        cout<<FBLU("--- STORING PER-EVENT EFT PARAMETERIZATIONS IN TTREE '"<<newTree_name<<"', IN FILE: '"<<filepath<<"' ...")<<endl;
        cout<<"("<<nentries<<" entries) "<<endl;

        // TH1EFT* th1eft = new TH1EFT("EFTparameterization", "EFTparameterization", 1, 0, 1);

        //--- Branch address
        double eventWeight; //Corresponds to 'nominalMEWeight_' in TopAnalysis code
        float eventMCFactor; //Sample-specific SF (xsec*lumi/SWE)
        Float_t weightMENominal;
        t->SetBranchStatus("eventWeight", 1);
        t->SetBranchAddress("eventWeight", &eventWeight);
        t->SetBranchStatus("eventMCFactor", 1);
        t->SetBranchAddress("eventMCFactor", &eventMCFactor);
        t->SetBranchStatus("weightMENominal", 1);
        t->SetBranchAddress("weightMENominal", &weightMENominal);

        vector<float>* v_wgts = new vector<float>;
        vector<string>* v_ids = new vector<string>;
        t->SetBranchStatus("mc_EFTweights", 1);
        t->SetBranchAddress("mc_EFTweights", &v_wgts);
        t->SetBranchStatus("mc_EFTweightIDs", 1);
        t->SetBranchAddress("mc_EFTweightIDs", &v_ids);

        vector<float> v_SWE; //Store Sums of Weights (SWE) for all reweight points -- for private MC samples only
        TString hSWE_name = "EFT_SumWeights";
        if(!f->GetListOfKeys()->Contains(hSWE_name)) {cout<<"ERROR ! Histogram "<<hSWE_name<<" containing the sums of weights not found... Abort !"<<endl; return;}

        //Read and store sums of weights (SWE)
        TH1F* h_SWE = (TH1F*) f->Get(hSWE_name);
        for(int ibin=0; ibin<h_SWE->GetNbinsX(); ibin++)
        {
            v_SWE.push_back(h_SWE->GetBinContent(ibin+1)); //1 SWE stored for each stored weight
        }
        delete h_SWE;

        //For private EFT samples, get and store index of SM reweight
        int idx_sm = -1;
        t->GetEntry(0); //Read 1 entry
        for(int iwgt=0; iwgt<v_ids->size(); iwgt++)
        {
            if(ToLower(v_ids->at(iwgt)).Contains("_sm") ) {idx_sm = iwgt; break;} //SM weight found
            if(((TString) v_ids->at(iwgt)).Contains("EFTrwgt183_") ) {idx_sm = iwgt; break;} //SM weight found //TOP19001 sample
        }
        if(idx_sm == -1) {cout<<BOLD(FRED("Error: SM reweight not found in private sample ! Abort ! "))<<endl; return;}

        //-- Create new additional branch
        WCFit* eft_fit = new WCFit;
        f->cd();
        TTree* new_tree = new TTree(newTree_name, "");
        new_tree->Branch("eft_fit", &eft_fit);

        for(int ientry=0; ientry<nentries; ientry++)
        {
            t->GetEntry(ientry);

            if(ientry%5000==0) {cout<<DIM(" --- "<<ientry<<" / "<<nentries<<"")<<endl;}

            float w_SMpoint = eventWeight * eventMCFactor * v_wgts->at(idx_sm) / (weightMENominal * v_SWE[idx_sm]);
            Get_WCFit(eft_fit, v_ids, v_wgts, v_SWE, (eventWeight*eventMCFactor), weightMENominal, w_SMpoint, idx_sm);
            new_tree->Fill();
        }

        new_tree->Write(newTree_name, TObject::kOverwrite);

        delete eft_fit; eft_fit = NULL;
        delete new_tree; new_tree = NULL;
    } //TTrees
    f->Close();

    cout<<endl<<FYEL("... Done !")<<endl<<endl;

    return;
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
void Create_Subsample_fromSample(TString fullsample_name, TString newsample_name, TString selection, TString sample, vector<TString> v_TTrees, TString nominal_tree_name, TString open_mode = "UPDATE")
{
    cout<<endl<<BYEL("                                  ")<<endl;
    cout<<FYEL("--- WILL EXTRACT SUBSET OF EVENTS SATISFYING '"<<selection<<"'")<<endl;
    cout<<FYEL("FROM ["<<fullsample_name<<"]")<<endl;
    cout<<FYEL("--> TO ["<<newsample_name<<"] ---")<<endl;

    //Input file
    if(!Check_File_Existence(fullsample_name) ) {cout<<FRED("Sample "<<fullsample_name<<" not found !")<<endl; return;}
    TFile* f_data = new TFile(fullsample_name);

    //Output file
    TFile *f_new = new TFile(newsample_name, open_mode);

    for(int itree=0; itree<v_TTrees.size(); itree++)
    {
        cout<<"-- Tree : '"<<v_TTrees[itree]<<"'..."<<endl;
        if(v_TTrees[itree] != nominal_tree_name && (newsample_name.Contains("DATA") || newsample_name.Contains("NPL"))) {continue;}

		TTree* t_input = NULL; //input tree
        t_input = (TTree*) f_data->Get(v_TTrees[itree]);
        if(!t_input) {cout<<FRED("Null tree "<<v_TTrees[itree]<<" !")<<endl; return;}
        Int_t nentries = t_input->GetEntries();

    	f_new->cd();
        TTree *t_new = NULL; //output tree
        t_new = (TTree*) t_input->CopyTree(selection);
        if(!t_new) {cout<<FRED("t_new is NULL !")<<endl; return;}

        //-- Enforce naming convention
        TString t_new_name = v_TTrees[itree];
        if(v_TTrees[itree] == "TotalDown") {t_new_name = "JESDown";}
        else if(v_TTrees[itree] == "TotalUp") {t_new_name = "JESUp";}

        //MC Prompt fake contributions : merge all contributions in same file
        f_new->cd();
        if(newsample_name.Contains("NPL_MC") ) {t_new->Write(t_new_name + "_" + sample, TObject::kOverwrite);}
        else {t_new->Write();}

        cout<<"("<<t_new->GetEntries()<<" entries)"<<endl;

        delete t_input; t_input = NULL;
        delete t_new; t_new = NULL;
    } //TTrees

    f_data->Close();
    f_new->Close();
    cout<<BYEL("                                  ")<<endl<<endl;

    return;
}


//--------------------------------------------
//  ######  ########  ##       #### ########    ########  ########  ######     ###    ##    ##
// ##    ## ##     ## ##        ##     ##       ##     ## ##       ##    ##   ## ##    ##  ##
// ##       ##     ## ##        ##     ##       ##     ## ##       ##        ##   ##    ####
//  ######  ########  ##        ##     ##       ##     ## ######   ##       ##     ##    ##
//       ## ##        ##        ##     ##       ##     ## ##       ##       #########    ##
// ##    ## ##        ##        ##     ##       ##     ## ##       ##    ## ##     ##    ##
//  ######  ##        ######## ####    ##       ########  ########  ######  ##     ##    ##
//--------------------------------------------

/**
 * Special case: may want to split the full WZ sample according to the flavours of the additional jets
 */
void Split_WZ_sample_byJetFlavour(TString prefix, TString dir, TString filepath, TString selection, TString sample, vector<TString> v_TTrees, TString nominal_tree_name)
{
    if(sample != "WZ") {return;} //Only for WZ sample

    TString outfile_path = "";
    TString selection_tmp = "";

    outfile_path = prefix + dir + "/WZ_b.root";
    selection_tmp = "fragmentationInfo==2";
	if(selection != "") {selection_tmp+= " && " + selection;}
	Create_Subsample_fromSample(filepath, outfile_path, selection_tmp, sample, v_TTrees, nominal_tree_name);

    outfile_path = prefix + dir + "/WZ_c.root";
    selection_tmp = "fragmentationInfo==1";
	if(selection != "") {selection_tmp+= " && " + selection;}
	Create_Subsample_fromSample(filepath, outfile_path, selection_tmp, sample, v_TTrees, nominal_tree_name);

    outfile_path = prefix + dir + "/WZ_l.root";
    selection_tmp = "fragmentationInfo==0";
	if(selection != "") {selection_tmp+= " && " + selection;}
	Create_Subsample_fromSample(filepath, outfile_path, selection_tmp, sample, v_TTrees, nominal_tree_name);

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
void Copy_SumWeight_Histogram_Into_SplitSample(TString filepath, TString outfile_path, TString sample)
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
// ##     ## ######## ########   ######   ########
// ###   ### ##       ##     ## ##    ##  ##
// #### #### ##       ##     ## ##        ##
// ## ### ## ######   ########  ##   #### ######
// ##     ## ##       ##   ##   ##    ##  ##
// ##     ## ##       ##    ##  ##    ##  ##
// ##     ## ######## ##     ##  ######   ########

// ########    ###    ##    ## ########  ######
// ##         ## ##   ##   ##  ##       ##    ##
// ##        ##   ##  ##  ##   ##       ##
// ######   ##     ## #####    ######    ######
// ##       ######### ##  ##   ##             ##
// ##       ##     ## ##   ##  ##       ##    ##
// ##       ##     ## ##    ## ########  ######
//--------------------------------------------


 //  ####  #    # #####   ####    ##   #    # #####  #      ######  ####
 // #      #    # #    # #       #  #  ##  ## #    # #      #      #
 //  ####  #    # #####   ####  #    # # ## # #    # #      #####   ####
 //      # #    # #    #      # ###### #    # #####  #      #           #
 // #    # #    # #    # #    # #    # #    # #      #      #      #    #
 //  ####   ####  #####   ####  #    # #    # #      ###### ######  ####

//-- MC Prompt fake contributions: for each subcategory, merge all contributions (prompt MC fakes from all prompt processes) into 1 single file
void Merge_Many_TTrees_Into_One(vector<TString> v_years, vector<TString> v_sel, vector<TString> v_samples, vector<TString> v_TTrees, TString prefix, TString nominal_tree_name)
{
    for(int iyear=0; iyear<v_years.size(); iyear++)
    {
        for(int itree=0; itree<v_TTrees.size(); itree++)
        {
            cout<<"-- Merging TTree : '"<<v_TTrees[itree]<<"' from several samples (for 'NPL_MC' samples)"<<endl;

            for(int isel=0; isel<v_sel.size(); isel++)
            {
                if(v_TTrees[itree] != nominal_tree_name) {continue;} //Don't consider JES/JER variations for Fakes(MC) sample
                if(!v_sel[isel].Contains("Fake")) {continue;} //Only care about fake events

                TString dir = Get_Directory(v_sel[isel]);
                if(dir == "") {continue;} //Sub-directories only

                //-- Fakes/Flip/Conv samples
                TString f_MCPrompt_path = prefix + dir + v_years[iyear] + "/NPL_MC.root";
                if(!Check_File_Existence(f_MCPrompt_path) ) {cout<<BOLD(FRED("Sample "<<f_MCPrompt_path<<" not found ! Can not merge multiple TTrees inside this file into a single one !") )<<endl; continue;}
                TFile* f_MCPrompt = new TFile(f_MCPrompt_path, "UPDATE");
                vector<TTree*> v_new_trees(0);

                //Get list of all TTrees
                for(int isample=0; isample<v_samples.size(); isample++)
                {
                    if(v_samples[isample].Contains("PrivMC") || v_samples[isample].Contains("TTbar") || v_samples[isample].Contains("DY")) {continue;} //Don't consider prompt fake contributions from: private samples / ttbar / DY / ...

                    // cout<<"Sample : "<<v_samples[isample]<<endl;

                    TTree* t_tmp = (TTree*) f_MCPrompt->Get(v_TTrees[itree] + "_" + v_samples[isample]);
                    // cout<<"t_tmp->GetEntries() "<<t_tmp->GetEntries()<<endl;
                    // if(t_tmp != 0 && t_tmp->GetEntries() > 0) {v_new_trees.push_back(t_tmp);}
                    if(t_tmp != 0) {v_new_trees.push_back(t_tmp);}
                    else {cout<<DIM("Tree "<<v_TTrees[itree] + "_" + v_samples[isample]<<" not found in file "<<f_MCPrompt_path<<" !")<<endl;}

                    //DEBUG inf/isnan values
                    // float test;
                    // t_tmp->SetBranchAddress("weight", &test);
                    // for(int ientry=0; ientry<t_tmp->GetEntries(); ientry++)
                    // {
                    //     t_tmp->GetEntry(ientry);
                    //     if(isinf(test) || isnan(test))
                    //     {
                    //         cout<<"Sample "<<v_samples[isample]<<endl;
                    //         cout<<endl<<endl<<"ERROR : ifnan or isinf value : "<<test<<" ! Did you properly rescale the sample ? (<-> did you produce the 'Sum of Weights' file using the code 'Get_Merged_Histograms_From_FlatTrees.cxx' ? ) ABORT !"<<endl<<endl;
                    //         return;
                    //     }
                    // }
                }

                if(!v_new_trees.size() ) {continue;}

                TList *list = new TList;
                for(int jtree=0; jtree<v_new_trees.size(); jtree++)
                {
                    // cout<<"v_new_trees[jtree] "<<v_new_trees[jtree]<<endl;
                    list->Add(v_new_trees[jtree]);
                }

                //Merge TTrees in one //WARNING : all empty TTrees are skipped !
                TTree *newtree = 0;
    			newtree = TTree::MergeTrees(list);
                // cout<<"newtree "<<newtree<<endl;

                if(list) {delete list;}

                //-- Enforce naming convention
                TString t_new_name = v_TTrees[itree];
                if(v_TTrees[itree] == "TotalDown") {t_new_name = "JESDown";}
                else if(v_TTrees[itree] == "TotalUp") {t_new_name = "JESUp";}

                //-- Write merged TTree
                if(newtree) {newtree->Write(t_new_name, TObject::kOverwrite);}
    			else {continue;}
                // cout<<"Write "<<v_TTrees[itree]<<endl;

                //Delete previous TTrees
                for(int isample=0; isample<v_samples.size(); isample++)
                {
                    if(v_samples[isample].Contains("PrivMC") || v_samples[isample].Contains("TTbar") || v_samples[isample].Contains("DY")) {continue;} //Don't consider prompt fake contributions from: private samples / ttbar / DY / ...

                    TString delete_name = v_TTrees[itree] + "_" + v_samples[isample]+";1";
                    f_MCPrompt->Delete(delete_name); //Delete (first cycle of) TTree for each sample
                    // cout<<"Deleted : "<<delete_name<<" in "<<f_MCPrompt_path<<endl;
                }

                f_MCPrompt->Close(); //Also closes associated TList/TTree !

    			cout<<FYEL("-> Updated file "<<f_MCPrompt_path<<"")<<endl;
    		} //sel loop
        } //tree loop
    } //year loop

    return;
}


// ###### #    # #      #
// #      #    # #      #
// #####  #    # #      #
// #      #    # #      #
// #      #    # #      #
// #       ####  ###### ######

/**
 * Make merged ntuples (NPL_DATA, NPL_MC, ...), but this time *without* any subcategorization
 * Arg 'datadriven' determines whether we are creating the 'NPL_DATA' or 'NPL_MC' sample; in both cases, store all events satisfying the 'NPL_flag' option (i.e. we don't care about any sub-selection)
 */
void Make_Full_Merged_Ntuples(vector<TString> v_years, vector<TString> v_TTrees, vector<TString> v_samples, bool make_FakesMC_samples, TString prefix, TString NPL_flag, bool datadriven)
{
    for(int iyear=0; iyear<v_years.size(); iyear++)
    {
        TString output_path = prefix + v_years[iyear] + "/NPL_MC.root";
        if(datadriven)
        {
            bool dataSampleFound = false;
            for(int isample=0; isample<v_samples.size(); isample++)
            {
                if(v_samples[isample] == "DATA") {dataSampleFound = true;}
            }
            if(!dataSampleFound) {continue;}
            output_path = prefix + v_years[iyear] + "/NPL_DATA.root";
        }

        cout<<endl<<endl<<FYEL("-- Creating merged ntuple : "<<output_path<<"")<<endl<<endl;

        //-- Output file
        TFile* f_new = new TFile(output_path, "RECREATE"); //NB: will erase previous file even if it is not reproduced...

        for(int itree=0; itree<v_TTrees.size(); itree++)
        {
            if(itree>0) {break;} //Don't consider JEV systematics for NPL

    		cout<<endl<<FBLU("* Tree : "<<v_TTrees[itree]<<"")<<endl;

            f_new->cd();
            vector<TTree*> v_new_trees; //Vector containing new TTree(s) (--> merged if >=1)

            for(int isample=0; isample<v_samples.size(); isample++)
            {
                if(datadriven && v_samples[isample] != "DATA") {continue;}
                else if(!datadriven && v_samples[isample] == "DATA") {continue;}
                else if(v_samples[isample].Contains("PrivMC") || v_samples[isample].Contains("TTbar") || v_samples[isample].Contains("DY")) {continue;} //Don't consider prompt fake contributions from: private samples / ttbar / DY / ...

                cout<<endl<<"Sample "<<v_samples[isample]<<" ..."<<endl;

                TString input_filepath = prefix + v_years[iyear] + "/" + v_samples[isample] + ".root";
                if(!Check_File_Existence(input_filepath) ) {cout<<BOLD(FRED("Sample "<<input_filepath<<" not found !") )<<endl; continue;}
                TFile* f_input = TFile::Open(input_filepath, "READ");

                TTree* t_input = 0; //input tree
                t_input = (TTree*) f_input->Get(v_TTrees[itree]);
                if(!t_input) {cout<<FRED("Null tree "<<v_TTrees[itree]<<" !")<<endl; return;}

                f_new->cd();
                TTree *t_new = 0; //output tree
                t_new = (TTree*) t_input->CopyTree(NPL_flag); //Copy subset of events satisfying Fakes flag
                if(!t_new) {cout<<FRED("t_new is NULL !")<<endl; return;}
    			else {cout<<DIM("("<<t_new->GetEntries()<<" entries)")<<endl;}

                f_new->cd();
                v_new_trees.push_back((TTree*) t_new->Clone());

    			delete t_input; delete t_new;

                f_input->Close();
            } //samples

            if(!v_new_trees.size()) {continue;}

            TList *list = new TList;
            for(int jtree=0; jtree<v_new_trees.size(); jtree++)
            {
                // cout<<"v_new_trees[jtree] "<<v_new_trees[jtree]<<endl;
                list->Add(v_new_trees[jtree]);
            }

            //-- Merge TTrees in one //WARNING : all empty TTrees are skipped !
            TTree *newtree = TTree::MergeTrees(list);
            // cout<<"newtree "<<newtree<<endl;

            //-- Enforce naming convention
            TString t_new_name = v_TTrees[itree];
            if(v_TTrees[itree] == "TotalDown") {t_new_name = "JESDown";}
            else if(v_TTrees[itree] == "TotalUp") {t_new_name = "JESUp";}

            //-- Write merged TTree
            f_new->cd();
            newtree->Write(t_new_name, TObject::kOverwrite);
            // cout<<"Write "<<v_TTrees[itree]<<endl;

            delete newtree;
            delete list;
        } //trees

        f_new->Close();

    	cout<<FYEL("-> Wrote file "<<output_path<<"")<<endl;
    } //years

    return;
}


//--------------------------------------------
// ##     ## ######## ########   ######   ########     ######   ########   #######  ##     ## ########   ######
// ###   ### ##       ##     ## ##    ##  ##          ##    ##  ##     ## ##     ## ##     ## ##     ## ##    ##
// #### #### ##       ##     ## ##        ##          ##        ##     ## ##     ## ##     ## ##     ## ##
// ## ### ## ######   ########  ##   #### ######      ##   #### ########  ##     ## ##     ## ########   ######
// ##     ## ##       ##   ##   ##    ##  ##          ##    ##  ##   ##   ##     ## ##     ## ##              ##
// ##     ## ##       ##    ##  ##    ##  ##          ##    ##  ##    ##  ##     ## ##     ## ##        ##    ##
// ##     ## ######## ##     ##  ######   ########     ######   ##     ##  #######   #######  ##         ######
//--------------------------------------------

void Merge_Samples_byGroups(vector<TString> v_samples, vector<TString> v_sampleGroups, vector<TString> v_sel, vector<TString> v_years, TString prefix, bool merge_full_samples)
{
    cout<<endl<<endl<<"---------------------------"<<endl;
    cout<<FBLU("== START MERGING SAMPLES BY GROUPS ==")<<endl;
    cout<<"Make sure the sample groups defined in [Split_FullSamples] are the same as in main AN code [TopEFT_analysis] !"<<endl;
    cout<<"---------------------------"<<endl<<endl<<endl;
    usleep(5000000); //Pause in ms

    //-- Trick: NPL_DATA & NPL_MC samples are not in the list (they are produced by this code from other samples). But once they have been produced, add them to the sample lists so that we merge NPL=(NPL_DATA+NPL_MC) //NB: NPL_MC events have negative weights, so can simply hadd to substract !
    v_samples.push_back("NPL_DATA"); v_sampleGroups.push_back("NPL");
    v_samples.push_back("NPL_MC"); v_sampleGroups.push_back("NPL");

    for(int iyear=0; iyear<v_years.size(); iyear++)
    {
        for(int isel=0; isel<v_sel.size(); isel++)
        {
            if(!v_sel[isel].Contains("Fake")) {continue;} //Only consider the 'real' SR/CR regions, not the convenience 'Fake' flags

            for(int igroup=0; igroup<v_samples.size(); igroup++) //Loop on sample groups
            {
                if(igroup > 0 && v_sampleGroups[igroup] == v_sampleGroups[igroup-1]) {continue;} //If this sample group was already processed, skip

                TString dir = Get_Directory(v_sel[isel]);
                if(merge_full_samples) {dir = "";} //Merging 'full' samples by groups <-> perform merging in main per-year ntuple repos, not in sub-category repos !
                TString full_dirname = prefix + dir + v_years[iyear] + "/";
                TString outfile_path = full_dirname + v_sampleGroups[igroup] + ".root";
                // cout<<"outfile_path = "<<outfile_path<<endl;

                TString command_merge = "hadd -f " + outfile_path;
                int nsamples_toMerge = 0; //Check how many samples belong to this sample group
                for(int isample=0; isample<v_samples.size(); isample++)
                {
                    TString filepath = full_dirname + v_samples[isample] + ".root";

                    if(v_samples[isample].Contains("PrivMC")) {continue;} //Private SMEFT samples should not be merged by group (for now)
                    if(v_sampleGroups[isample] != v_sampleGroups[igroup]) {continue;} //Only consider samples which belong to the current sample group
                    if(!Check_File_Existence(filepath)) {cout<<DIM("NB: can't merge file: "<<filepath<<" (not found) !")<<endl; continue;} //Can only hadd these TFile if they were already produced !

                    nsamples_toMerge++;
                    command_merge+= " " + filepath; //Merge this sample
                }

                if(nsamples_toMerge >= 2) // There are at least 2 samples to hadd
                {
                    cout<<endl<<FYEL("--- Merging sub-samples by group as such:")<<endl;
                    cout<<"--> "<<DIM(<<command_merge<<"")<<endl<<endl;

                    system(command_merge.Data());
                }
            } //groups
        } //sels
    } //years

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
void Split_AllNtuples_ByCategory(vector<TString> v_samples, vector<TString> v_sampleGroups, vector<TString> v_sel, vector<TString> v_years, bool make_nominal_samples, bool make_FakesMC_samples, bool make_FakesDATA_fullSample, vector<TString> v_TTrees, TString NPL_flag, TString nominal_tree_name, bool store_WCFit_inSMEFTsubsamples, bool split_WZ_byJetFlavour, bool hadd_subsamples_byGroup, bool hadd_fullSamples_byGroup, bool update_fullSMEFTSamples_withWCFit)
{
    cout<<endl<<endl<<FBLU("== START OF NTUPLES SPLITTING ==")<<endl;
    cout<<"This can be quite long. Make sure you have correctly selected in the code :"<<endl;
    cout<<"- The samples & categories to split"<<endl;
    cout<<"- The TTrees you want to split : Nominal/JES/JER/... TTrees"<<endl;
    cout<<"---------------------------"<<endl<<endl<<endl;
    usleep(4000000); //Pause for 4s

    TString prefix = NTUPLEDIR; //Defined in Utils/Helper.h

    //-- Split the full ntuples by sub-categories
    for(int iyear=0; iyear<v_years.size(); iyear++)
    {
        for(int isel=0; isel<v_sel.size(); isel++)
        {
            // cout<<"v_sel[isel] "<<v_sel[isel]<<endl;

            TString opening_mode_FakesMC = "RECREATE"; //Each subcategory --> 1 NPL_MC sample. Must RECREATE it at beginning, then UPDATE it with each sample

            for(int isample=0; isample<v_samples.size(); isample++)
            {
                // cout<<"v_samples[isample] "<<v_samples[isample]<<endl;

                //-- Dir. path
                TString dir = Get_Directory(v_sel[isel]);
                if(dir == "") {continue;} //Sub-directories only
                TString full_dirname = prefix + dir;
                mkdir(full_dirname.Data(), 0777);
                full_dirname+= v_years[iyear] + "/";
                mkdir(full_dirname.Data(), 0777);
                TString filepath = prefix + v_years[iyear] + "/" + v_samples[isample] + ".root";
                // cout<<"prefix = "<<prefix<<endl;

                //-- Store EFT param. in 'full' ntuples
                if((update_fullSMEFTSamples_withWCFit && isel==0) && v_samples[isample].Contains("PrivMC") && !v_samples[isample].Contains("_c")) {Store_EFTparameterization(filepath, v_TTrees, nominal_tree_name);}

                //-- Skip unwanted selection/sample combinations
                // if(!make_nominal_samples && !v_sel[isel].Contains("Fake")) {continue;} //Only fake categories
                if(!make_nominal_samples) {continue;}
                if(!make_FakesMC_samples && v_sel[isel].Contains("Fake") && v_samples[isample] != "DATA") {continue;} //No fake MC
                else if(v_sel[isel].Contains("Fake") && (v_samples[isample].Contains("PrivMC") || v_samples[isample].Contains("TTbar") || v_samples[isample].Contains("DY"))) {continue;} //Don't consider prompt fake contributions from: private samples / ttbar / DY / ...
                else if(!make_nominal_samples && v_sel[isel].Contains("Fake") && v_samples[isample] == "DATA") {continue;} //If 'make_nominal_samples=False', dont make NPL_DATA sub-ntuples neither !

    			if(!Check_File_Existence(filepath)) {cout<<BOLD(FRED("Sample "<<filepath<<" not found !") )<<endl; continue;}
    			TString outfile_path = full_dirname + v_samples[isample] + ".root";
                if(v_sel[isel].Contains("Fake"))
                {
                    if(v_samples[isample] == "DATA") {outfile_path = full_dirname + "NPL_DATA.root";}
                    else {outfile_path = full_dirname + "NPL_MC.root";} //NPL MC prompt summed contribution
                }
                // cout<<"outfile_path = "<<outfile_path<<endl;

                TString open_mode = "RECREATE"; //1 file per promtp sample
                if(v_sel[isel].Contains("Fake")) {open_mode = opening_mode_FakesMC;} //1 common file for all fake MC samples

            	Create_Subsample_fromSample(filepath, outfile_path, v_sel[isel], v_samples[isample], v_TTrees, nominal_tree_name, open_mode);

                if(v_samples[isample].Contains("PrivMC") && !v_samples[isample].Contains("_c")) {Copy_SumWeight_Histogram_Into_SplitSample(filepath, outfile_path, v_samples[isample]);} //'_c' <-> identifier for pure-EFT samples (no parameterization)
                if(store_WCFit_inSMEFTsubsamples && v_samples[isample].Contains("PrivMC") && !v_samples[isample].Contains("_c")) {Store_EFTparameterization(outfile_path, v_TTrees, nominal_tree_name);}
                if(v_sel[isel].Contains("Fake")) {opening_mode_FakesMC = "UPDATE";} //Will update the TFile with next samples
                if(split_WZ_byJetFlavour) {Split_WZ_sample_byJetFlavour(prefix, dir, filepath, v_sel[isel], v_samples[isample], v_TTrees, nominal_tree_name);}
            } //sample loop
    	} //selections loop
    } //year loop

    //-- Also create a 'full' NPL_DATA sample (no sub-cat.)
    if(make_FakesDATA_fullSample) {Make_Full_Merged_Ntuples(v_years, v_TTrees, v_samples, make_FakesMC_samples, prefix, NPL_flag, true);} //Data-driven NPL

    //-- Make the NPL_MC samples (full & split by sub-categories), by merging the 'fake' contributions from multiple prompt MC samples
    if(make_FakesMC_samples)
    {
        if(v_samples.size() < 20) //Protection
        {
            cout<<FRED("WARNING: v_samples.size()="<<v_samples.size()<<"... Are you sure you are considering the complete sample list ? (necessary to produce correct NPL_MC ntuples !)")<<endl;
            cout<<endl<<"=== ENTER 'y' TO CONTINUE ==="<<endl<<endl;
            TString answer = "n"; cin>>answer;
            while(answer != 'y') {cout<<"Try again... "<<endl; cin>>answer;}
        }

        if(make_nominal_samples) {Merge_Many_TTrees_Into_One(v_years, v_sel, v_samples, v_TTrees, prefix, nominal_tree_name);} //Make subregion/NPL_MC sample only if make_nominal_samples==make_FakesMC_samples==true
        Make_Full_Merged_Ntuples(v_years, v_TTrees, v_samples, make_FakesMC_samples, prefix, NPL_flag, false); //NPL MC
    }

    //-- Merge sub-ntuples (split by sub-categories) into ntuples grouped by 'sample groups'
    if(hadd_subsamples_byGroup) {Merge_Samples_byGroups(v_samples, v_sampleGroups, v_sel, v_years, prefix, false);} //Merge by sample group

    //-- Merge 'full' ntuples (*not* split by sub-categories) into ntuples grouped by 'sample groups'
    if(hadd_fullSamples_byGroup) {Merge_Samples_byGroups(v_samples, v_sampleGroups, v_sel, v_years, prefix, true);} //Merge by sample group
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
    bool make_nominal_samples = false; //true <-> create sub-samples satisfying given category flags
    bool make_FakesMC_samples = false; //true <-> merge the MC prompt+fake contribution into a single "NPL_MC" sample (for full ntuples, and also for sub-ntuples in sub-categories if 'make_nominal_samples=true')
    bool make_FakesDATA_fullSample = false; //true <-> make 'full' sample (no subcat.) for data-driven NPL contribution
    bool hadd_subsamples_byGroup = false; //true <-> hadd the ntuples (split by sub-categories) into 'sample group' ntuples (e.g. tX, ...)
    bool hadd_fullSamples_byGroup = false; //true <-> hadd the 'full' ntuples (*not* split by sub-categories) into 'sample group' ntuples (e.g. tX, ...)
    bool update_fullSMEFTSamples_withWCFit = false; //true <-> update the 'full' private SMEFT samples, compute+store the WCFit objects for all events in the files (faster to read the EFT parameterization later in the analysis) //Extremely slow when considering many TTrees (few hours!) -- but makes it all the more necessary
    TString NPL_flag = "isFake"; //Flag defining fake events
    bool store_WCFit_inSMEFTsubsamples = true; //true <-> also store per-event EFT parameterization for SMEFT samples (so that it can be then read directly when processing the sample)

    TString nominal_tree_name = "result"; //Hard-coded nominal tree name (special case)
    bool split_WZ_byJetFlavour = false; //NOT TESTED ! //true <-> also split WZ sample depending on flavour of additional jet
//--------------------------------------------


 //  ####  #   #  ####  #####    ##### ##### #####  ###### ######  ####
 // #       # #  #        #        #     #   #    # #      #      #
 //  ####    #    ####    #        #     #   #    # #####  #####   ####
 //      #   #        #   #        #     #   #####  #      #           #
 // #    #   #   #    #   #        #     #   #   #  #      #      #    #
//  ####    #    ####    #        #     #   #    # ###### ######  ####

    //-- Copy nominal + JES/JER TTrees
    //-- NB: first must be nominal tree name
    vector<TString> v_TTrees;
    v_TTrees.push_back("result");
    v_TTrees.push_back("JESDown"); v_TTrees.push_back("JESUp");
    v_TTrees.push_back("JERDown"); v_TTrees.push_back("JERUp");
    v_TTrees.push_back("METDown"); v_TTrees.push_back("METUp");


 //  ####    ##   #    # #####  #      ######  ####
 // #       #  #  ##  ## #    # #      #      #
 //  ####  #    # # ## # #    # #      #####   ####
 //      # ###### #    # #####  #      #           #
 // #    # #    # #    # #      #      #      #    #
 //  ####  #    # #    # #      ###### ######  ####

    //--- Sample list
    vector<TString> v_samples; vector<TString> v_sample_groups;

    //-- MAIN ANALYSIS SAMPLES
	v_samples.push_back("DATA"); v_sample_groups.push_back("DATA");
    v_samples.push_back("PrivMC_tZq"); v_sample_groups.push_back("PrivMC_tZq");
    v_samples.push_back("PrivMC_ttZ"); v_sample_groups.push_back("PrivMC_ttZ");
    v_samples.push_back("PrivMC_tWZ"); v_sample_groups.push_back("PrivMC_tWZ");
    v_samples.push_back("tZq"); v_sample_groups.push_back("tZq");
    v_samples.push_back("ttZ"); v_sample_groups.push_back("ttZ");
    v_samples.push_back("tWZ"); v_sample_groups.push_back("tX");
    v_samples.push_back("ttZ_M1to10"); v_sample_groups.push_back("tX");
	v_samples.push_back("tHq"); v_sample_groups.push_back("tX");
    v_samples.push_back("tHW"); v_sample_groups.push_back("tX");
    v_samples.push_back("ttH"); v_sample_groups.push_back("tX");
    v_samples.push_back("ttW"); v_sample_groups.push_back("tX");
    v_samples.push_back("ttZZ"); v_sample_groups.push_back("tX");
    v_samples.push_back("ttWW"); v_sample_groups.push_back("tX");
    v_samples.push_back("ttWZ"); v_sample_groups.push_back("tX");
    v_samples.push_back("ttZH"); v_sample_groups.push_back("tX");
    v_samples.push_back("ttWH"); v_sample_groups.push_back("tX");
    v_samples.push_back("tttt"); v_sample_groups.push_back("tX");
    v_samples.push_back("ttHH"); v_sample_groups.push_back("tX");
    v_samples.push_back("ZZ4l"); v_sample_groups.push_back("VVV");
    v_samples.push_back("ZZZ"); v_sample_groups.push_back("VVV");
    v_samples.push_back("WZZ"); v_sample_groups.push_back("VVV");
    v_samples.push_back("WWW"); v_sample_groups.push_back("VVV");
    v_samples.push_back("WWZ"); v_sample_groups.push_back("VVV");
    v_samples.push_back("WZ"); v_sample_groups.push_back("WZ");
	v_samples.push_back("TTGamma_Dilep"); v_sample_groups.push_back("XG");
	v_samples.push_back("ZGToLLG_01J"); v_sample_groups.push_back("XG");

    //-- OBSOLETE
    // v_samples.push_back("tGJets"); v_sample_groups.push_back("XG");
    // v_samples.push_back("WGToLNuG"); v_sample_groups.push_back("XG");
    //v_samples.push_back("ggToZZTo4l"); v_sample_groups.push_back("VVV");
	// v_samples.push_back("DY"); v_sample_groups.push_back("DY");
    // v_samples.push_back("TTbar_DiLep"); v_sample_groups.push_back("TTbar_DiLep");

    //-- OTHER PRIVATE SAMPLES
    // v_samples.push_back("PrivMC_tZq_TOP19001"); v_sample_groups.push_back("xxx");
    // v_samples.push_back("PrivMC_ttZ_TOP19001"); v_sample_groups.push_back("xxx");
    // v_samples.push_back("PrivMC_tZq_ctz"); v_sample_groups.push_back("xxx");
    // v_samples.push_back("PrivMC_ttZ_ctz"); v_sample_groups.push_back("xxx");
    // v_samples.push_back("PrivMC_tZq_v3"); v_sample_groups.push_back("xxx");


 //  ####  ###### #      ######  ####  ##### #  ####  #    #  ####
 // #      #      #      #      #    #   #   # #    # ##   # #
 //  ####  #####  #      #####  #        #   # #    # # #  #  ####
 //      # #      #      #      #        #   # #    # #  # #      #
 // #    # #      #      #      #    #   #   # #    # #   ## #    #
 //  ####  ###### ###### ######  ####    #   #  ####  #    #  ####

    //--- Will divide samples based on these subcategories
    //-- NB: also include corresponding 'Fake' flags to produce NPL samples (NPL_DATA if ["DATA" is in v_samples], NPL_MC if [make_FakesMC_samples==true], and NPL both NPL_DATA & NPL_MC are produced)
    vector<TString> v_sel;
    v_sel.push_back("is_signal_SR");
    v_sel.push_back("is_signal_SRFake");


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
    Split_AllNtuples_ByCategory(v_samples, v_sample_groups, v_sel, v_years, make_nominal_samples, make_FakesMC_samples, make_FakesDATA_fullSample, v_TTrees, NPL_flag, nominal_tree_name, store_WCFit_inSMEFTsubsamples, split_WZ_byJetFlavour, hadd_subsamples_byGroup, hadd_fullSamples_byGroup, update_fullSMEFTSamples_withWCFit);

    return 0;
}
