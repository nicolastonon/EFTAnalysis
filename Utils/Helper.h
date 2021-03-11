#ifndef Helper_h
#define Helper_h

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

#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <sys/stat.h> // to be able to check file existence
#include <dirent.h> //list dir content
#include <numeric> //std::accumulate

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include <TObjArray.h>
#include <TObjString.h>

//Custom classes for EFT
#include "split_string.h"
#include "TH1EFT.h"

//--------------------------------------------
//-- HARDCODE THE BASE DIR. CONTAINING PRIVATE NTUPLES (defined in Helper.cxx, used by all codes)
extern TString NTUPLEDIR;

//--------------------------------------------
//-- Low-level helper functions
bool Check_File_Existence(const TString&, bool=false);
int MoveFile(TString, TString);
int CopyFile(TString, TString);
TString Convert_Number_To_TString(double, int=3);
float Convert_TString_To_Number(TString);
float Find_Number_In_TString(TString);
TString Convert_Sign_To_Word(TString);
std::pair<TString,TString> Break_Cuts_In_Two(TString);
bool Is_Event_Passing_Cut(TString, double);
double Compute_RangeScaled_Value(double, double, double, double, double);
double Compute_StdDevScaled_Value(double, double, double);
bool Get_Dir_Content(std::string, std::vector<TString>&);
TString Split_TString_Into_Keys(TString, TString);
int Count_nofHistos_inTFile(TString);
int Count_nofHistos_inTFile(TFile*);
float Round_Float_toDecimal(float, float, bool);
vector<TString> Break_TString(TString, TString);

//-- Basic analysis helper functions
void Fill_Last_Bin_TH1F(TH1F*, double); //Increment last bin of histo by specified weight
void Fill_First_Bin_TH1F(TH1F*, double); //Increment first bin of histo by specified weight
void Load_Canvas_Style();
TString Combine_Naming_Convention(TString);
void Extract_Ranking_Info(TString, TString);
void Get_Ranking_Vectors(TString, std::vector<TString>&, std::vector<double>&);
void Compare_Histograms(TString, TString, TString, TString);
float Rescale_Input_Variable(float, float, float);
void Get_WCFit(WCFit*&, vector<string>*, vector<float>*, const vector<float>&, float, float, float, int, int=25);
void Get_Mirror_Histogram(TH1F*&, TH1F*&, TH1F*&, bool);
void Get_TemplateSymm_Histogram(TH1F*&, TH1F*&, TH1F*&, bool);
void Inflate_Syst_inShapeTemplate(TH1F*&, TH1F*, float);

//-- Analysis-specific helper functions
bool Apply_CommandArgs_Choices(int, char**, std::vector<TString>&, TString&, TString&);
void Get_Samples_Colors(std::vector<int>&, std::vector<TColor*>&, std::vector<TString>, std::vector<TString>, int);
// void Set_Custom_ColorPalette(std::vector<TColor*>&, std::vector<int>&, std::vector<TString>); //Set custom color palette
bool Get_Variable_Range(TString, int&, float&, float&);
void Get_Template_Range(int&, float&, float&, TString, bool, bool, int, bool, int&, int&, int&, int&, vector<float>, bool=false);
TString Get_Variable_Name(TString);
TString Get_Category_Boolean_Name(TString, bool=false);
float Count_Total_Nof_Entries(TString, TString, std::vector<TString>, std::vector<TString>, std::vector<TString>, std::vector<TString>, std::vector<TString>, bool, bool);
TString Get_Modified_SystName(TString, TString, TString="");
void Get_Pointer_GENHisto(TH1F*&, TString);
vector<pair<TString,float>> Parse_EFTreweight_ID(TString);
float Get_x_jetCategory(float, float, int, int, int, int);
float Get_x_ZptCosCategory(float, float);
TString Get_MVAFile_InputPath(TString, TString, TString, bool, bool=true, bool=false, int=0, bool=false);
TString Get_HistoFile_InputPath(bool, TString, TString, TString, bool, TString, bool, int, bool, bool=false, bool=false, bool=true, bool=true);
bool Extract_Values_From_NNInfoFile(TString, vector<TString>&, vector<TString>&, TString&, TString&, int&, int&, vector<float>&, TString* NN_strategy=NULL);
TString Get_Region_Label(TString, TString);
void Fill_Variables_List(vector<TString>&, bool, TString, TString, bool, int, bool, TString, TString, vector<float>, vector<float>, bool, bool, bool=false);
vector<vector<float>> Get_nJets_SF(TString, TString, TString, vector<TString>);
float Apply_nJets_SF(vector<vector<float>>&, int, int, TString);
TString Get_Template_XaxisTitle(TString);
TString Get_EFToperator_label(TString);
TString Get_Unit_Variable(TString);
void Print_Yields_fromHistograms(TString, TString, vector<TString>, vector<TString>, vector<TString>);
//--------------------------------------------

//--------------------------------------------
//IN-PLACE IMPLEMENTATION
//Inline functions must be declared in header file
//Template functions: the compiler must be able to see the implementation in order to generate code for all specializations in your code

//Increment weight of first bin by 'weight'
inline void Fill_TH1F_UnderOverflow(TH1F* h, double value, double weight)
{
	if(value >= h->GetXaxis()->GetXmax() ) {h->Fill(h->GetXaxis()->GetXmax() - (h->GetXaxis()->GetBinWidth(1) / 2), weight);} //overflow in last bin
	else if(value <= h->GetXaxis()->GetXmin() ) {h->Fill(h->GetXaxis()->GetXmin() + (h->GetXaxis()->GetBinWidth(1) / 2), weight);} //underflow in first bin
	else {h->Fill(value, weight);}

	return;
};

inline void Fill_TH1EFT_UnderOverflow(TH1EFT* h, double value, float weight, WCFit& fit)
{
    if(value >= h->GetXaxis()->GetXmax() ) {h->Fill(h->GetXaxis()->GetXmax() - (h->GetXaxis()->GetBinWidth(1) / 2), weight, fit);} //overflow in last bin
    else if(value <= h->GetXaxis()->GetXmin() ) {h->Fill(h->GetXaxis()->GetXmin() + (h->GetXaxis()->GetBinWidth(1) / 2), weight, fit);} //underflow in first bin
    else {h->Fill(value, weight, fit);}
    return;
};

inline void Fill_TH1F_NoUnderOverflow(TH1F* h, double value, double weight)
{
    if(value < h->GetXaxis()->GetXmax() && value > h->GetXaxis()->GetXmin() ) {h->Fill(value, weight);}
    return;
};

//NB: got Combine error when setting yield/error way too small
template <class T> void Set_Histogram_FlatZero(T*& h, TString name="", bool printout=false)
{
    if(printout) {cout<<endl<<FRED("Histo "<<name<<" has integral = "<<h->Integral()<<" <= 0 ! Distribution set to ~>0 (flat), to avoid crashes in COMBINE !")<<endl;}

    for(int ibin=1; ibin<h->GetNbinsX()+1; ibin++)
    {
    	h->SetBinContent(ibin, pow(10, -5));
    	if(h->GetBinError(ibin) == 0) {h->SetBinError(ibin, pow(10, -6));}
    }

    return;
};

//-- Set empty or negative bins to 0, to avoid crashes in Combine (can't deal with <=0 PDFs)
template <class T> void Avoid_Histogram_EmptyOrNegativeBins(T*& h)
{
    for(int ibin=1; ibin<h->GetNbinsX()+1; ibin++)
    {
        //content <= 0 --> set content and error >~0
    	if(h->GetBinContent(ibin) <= 0) {h->SetBinContent(ibin, pow(10, -5)); h->SetBinError(ibin, pow(10, -6));}

        //error > content --> error = content //NB: not necessarily solid, be careful
        if(h->GetBinError(ibin) > h->GetBinContent(ibin)) {h->SetBinError(ibin, h->GetBinContent(ibin)-pow(10, -3));}
        // if(h->GetBinError(ibin) > h->GetBinContent(ibin)) {h->SetBinError(ibin, h->GetBinContent(ibin)/2.);}
        // if(h->GetBinError(ibin) > h->GetBinContent(ibin)/2.) {h->SetBinError(ibin, h->GetBinContent(ibin)/2.);}
    }

    return;
};

template <class T> void StoreEachHistoBinIndividually(TFile* f, T*& h, TString outname, int& nhistos, bool store_countExp=false)
{
    f->cd();
    if(h->GetNbinsX() == 1) {return;} //Dont need to split/merge anything if there's only 1 bin

    int firstBin = store_countExp? 0:1; //If want to also store histo as single bin (counting exp.), add corresponding 'bin 0' by convention
    for(int ibin=firstBin; ibin<h->GetNbinsX()+1; ibin++)
    {
        T* h_tmp = new T("", "", 1, 0, 1); //Create single-bin histo
        TString outname_tmp = "";

        if(!ibin) //ibin==0 --> Store entire histo content as a single bin (counting experiment)
        {
            Double_t integral=0, error=0;
            integral = h->IntegralAndError(0, h->GetNbinsX()+1, error); //Get integral and its error
            h_tmp->SetBinContent(1, integral);
            h_tmp->SetBinError(1, error);
            outname_tmp = "countExp_" + outname;
        }
        else //Store content of each histogram bin separately
        {
            h_tmp->SetBinContent(1, h->GetBinContent(ibin)); //Fill new single bin according to considered TH1EFT bin
            h_tmp->SetBinError(1, h->GetBinError(ibin));
            outname_tmp = "bin" + Convert_Number_To_TString(ibin) + "_" + outname;
        }

        if(!outname.Contains("NPL_MC") && !outname.Contains("data", TString::kIgnoreCase) && h_tmp->Integral() <= 0) {Set_Histogram_FlatZero(h_tmp, outname_tmp, false);} //If integral of histo is negative, set to 0 (else COMBINE crashes) -- must mean that norm is close to 0 anyway //Special case: NPL_MC is negative by design

        h_tmp->Write(outname_tmp);
        nhistos++;
        // cout<<"Wrote histo : "<<outname_tmp<<endl;

        delete h_tmp; h_tmp = NULL;
    }

    return;
};

/**
 * For theory systematics (only applied to signal samples), only want to consider the acceptance effect on normalization
 * <-> We apply a SF corresponding to (variation_beforeSel/nominal_beforeSel) <-> the normalization difference after selection should be purely due to differences in acceptance
 * NB: there is a small approximation, before the SFs were still obtained after some minimal preselection (applied at TopAnalysis level)
 */
template <class T> void Scale_THSyst_toBeforeSelection(T*& h, vector<float>& v_SWE_beforeSel, TString systname)
{
    int index_syst = -1;
    double nomYield_beforeSel = v_SWE_beforeSel[0]; //First element <-> nominal yield before event selection

    if(systname == "MEUp") {index_syst = 2-1;}
    else if(systname == "MEDown") {index_syst = 3-1;}

    else if(systname == "PDFUp") {index_syst = 9-1;}
    else if(systname == "PDFDown") {index_syst = 10-1;}

    else if(systname == "alphasUp") {index_syst = 11-1;}
    else if(systname == "alphasDown") {index_syst = 12-1;}

    else if(systname == "ISRUp") {index_syst = 13-1;}
    else if(systname == "ISRDown") {index_syst = 14-1;}

    else if(systname == "FSRUp") {index_syst = 15-1;}
    else if(systname == "FSRDown") {index_syst = 16-1;}

    if(index_syst <= 0 || v_SWE_beforeSel[index_syst] <= 0)
    {
        if(index_syst <= 0) {cout<<BOLD(FRED("[Scale_THSyst_toBeforeSelection] (Syst = "<<systname<<") ERROR: index_syst <= 0 !"))<<endl;}
        if(v_SWE_beforeSel[index_syst] <= 0) {cout<<BOLD(FRED("[Scale_THSyst_toBeforeSelection] (Syst = "<<systname<<") ERROR: v_SWE_beforeSel["<<index_syst<<"] = "<<v_SWE_beforeSel[index_syst]<<" !"))<<endl;}
        return;
    }

    // cout<<"index_syst "<<index_syst<<endl;
    // cout<<"[Scale_THSyst_toBeforeSelection / "<<systname<<"] Integral before: "<<h->Integral()<<endl;

    double SF = nomYield_beforeSel / v_SWE_beforeSel[index_syst];
    h->Scale(SF);

    // cout<<"[Scale_THSyst_toBeforeSelection / "<<systname<<"] h->Scale("<<SF<<")"<<endl;
    // cout<<"[Scale_THSyst_toBeforeSelection / "<<systname<<"] Integral after: "<<h->Integral()<<endl;

    return;
}

//--------------------------------------------

#endif
