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
// void Set_Histogram_FlatZero(TH1F*&, TString="", bool=false);
void Get_Mirror_Histogram(TH1F*&, TH1F*&, TH1F*&, bool);
void Get_TemplateSymm_Histogram(TH1F*&, TH1F*&, TH1F*&, bool);
void Inflate_Syst_inShapeTemplate(TH1F*&, TH1F*, float);

//-- Analysis-specific helper functions
bool Apply_CommandArgs_Choices(int, char **, std::vector<TString>&, TString&);
void Get_Samples_Colors(std::vector<int>&, std::vector<TColor*>&, std::vector<TString>, std::vector<TString>, int);
// void Set_Custom_ColorPalette(std::vector<TColor*>&, std::vector<int>&, std::vector<TString>); //Set custom color palette
bool Get_Variable_Range(TString, int&, double&, double&);
TString Get_Variable_Name(TString);
TString Get_Category_Boolean_Name(TString, bool=false);
float Count_Total_Nof_Entries(TString, TString, std::vector<TString>, std::vector<TString>, std::vector<TString>, std::vector<TString>, std::vector<TString>, bool, bool);
TString Get_Modified_SystName(TString, TString, TString="");
void Get_Pointer_GENHisto(TH1F*&, TString);
vector<pair<TString,float>> Parse_EFTreweight_ID(TString);
float Get_x_jetCategory(float, float, int, int, int, int);
float Get_x_ZptCosCategory(float, float);
TString Get_MVAFile_InputPath(TString, TString, TString, bool, bool=true, bool=false, int=0);
TString Get_HistoFile_InputPath(bool, TString, TString, TString, bool, TString, bool, int, bool=false);
bool Extract_Values_From_NNInfoFile(TString, vector<TString>&, vector<TString>&, TString&, TString&, int&, int&, TString* NN_strategy=NULL);
TString Get_Region_Label(TString, TString);
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


template <class T> void Set_Histogram_FlatZero(T*& h, TString name="", bool printout=false)
{
    if(printout)
    {
    	cout<<endl<<FRED("Histo "<<name<<" has integral = "<<h->Integral()<<" <= 0 ! Distribution set to ~>0 (flat), to avoid crashes in COMBINE !")<<endl;
    }

    for(int ibin=1; ibin<h->GetNbinsX()+1; ibin++)
    {
    	h->SetBinContent(ibin, pow(10, -9));

    	if(h->GetBinError(ibin) == 0) {h->SetBinError(ibin, 0.1);}
    }

    return;
};

template <class T> void StoreEachHistoBinIndividually(TFile* f, T*& h, TString outname)
{
    f->cd();
    for(int ibin=0; ibin<h->GetNbinsX()+1; ibin++)
    {
        T* h_tmp = new T("", "", 1, 0, 1);
        // TH1F* h_tmp = new TH1F("", "", 1, 0, 1);
        TString outname_tmp = "";

        if(!ibin) //ibin=0 --> Store entire histo content as single bin --> counting experiment
        {
            Double_t integral=0, error=0;
            integral = h->IntegralAndError(0, h->GetNbinsX()+1, error);
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

        if(h_tmp->Integral() <= 0) {Set_Histogram_FlatZero(h_tmp, outname_tmp, false);} //If integral of histo is negative, set to 0 (else COMBINE crashes) -- must mean that norm is close to 0 anyway

        h_tmp->Write(outname_tmp);
        // cout<<"Wrote histo : "<<outname_tmp<<endl;

        delete h_tmp; h_tmp = NULL;
    }

    return;
};
//--------------------------------------------

#endif
