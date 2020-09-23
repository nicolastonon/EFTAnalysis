/**
by Nicolas Tonon (DESY)

NB: the reason for separating the 'ROC_Plotter' (function declarations/definitions) and 'Compare_ROC_curves' codes is to allow function calls from other codes (i.e. without including twice a main() function)
**/

#ifndef ROC_h
#define ROC_h

#include "../Utils/Helper.h" //Helper functions

//--------------------------------------------

void Sum_Histograms(TH1F* &h, vector<TH1F*> v_histos);
int Get_Color(int index);
void Apply_Cosmetics(TCanvas* &c1);
double GetFloatPrecision(double value, double precision);
double Get_AUC_From_TMVAfile(TString filename, TString variable, bool use_TrainSample);
bool Create_Histogram_From_TMVA_Tree(bool is_bkg, TH1F*& h, TString filename, TString variable, bool use_TrainSample, TString cuts="1");
bool Retrieve_Histogram_From_TMVA_File(TH1F*& h, TString filename, TString variable, bool is_bkg, bool use_TrainSample);
bool Get_Histogram_From_KerasFile(TH1F*& h, TString filename, TString histo_name);
bool Get_Histogram_From_CustomFile(TH1F*& h, TString filename, TString histo_name);
bool Produce_Efficiency_TGraph(TGraph* &g, double& AUC, TH1F* h_sig, TH1F* h_bkg);
void Get_ROC_Curves(vector<TGraph*>& v_graph, vector<double>& v_AUC, vector<TString>& v_label, vector<TString> v_filepath, vector<TString> v_isTMVA_file, vector<TString> v_Filelabel, vector<bool> v_isTrainSample, TString region, vector<TString> v_processes, bool superimpose_allNodes_DNN, TString lumiYear, TString cuts="1", bool=false, TString="");
void Superimpose_ROC_Curves(vector<TGraph*> v_graph, vector<TString> v_label, vector<double> v_AUC);
void Make_Plot(vector<TString> v_filepath, vector<TString> v_Filelabel, vector<TString> v_isTMVA_file, vector<bool> v_isTrainSample, TString region, vector<TString> v_processes, bool superimpose_allNodes_DNN, TString lumiYear, TString cuts="1", bool=false, TString="");

//--------------------------------------------

#endif
