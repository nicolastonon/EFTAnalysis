#ifndef TopEFT_analysis_h
#define TopEFT_analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
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

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Timer.h"
#include "TMVA/Config.h"

#include <algorithm>
#include <cassert> 	//Can be used to terminate program if argument is not true. Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> //for mkdir

#include "Utils/Helper.h" //Helper functions
#include "Utils/TFModel.h" //Tensorflow functions

//Custom classes for EFT (see https://github.com/Andrew42/EFTGenReader/blob/maste)
#include "TH1EFT.h"
#include "WCPoint.h"
#include "WCFit.h"

using namespace std;

class TopEFT_analysis
{
	public :

	TopEFT_analysis(); //Default constructor
    TopEFT_analysis(vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<bool>, vector<TString>, TString, vector<TString>, bool, TString, TString, TString, bool);
	~TopEFT_analysis(); //Default destructor

//--- METHODS
	void Train_BDT(TString, bool); //Train BDT
    void Produce_Templates(TString, bool); //Produce templates
    void Draw_Templates(bool, TString, TString="", bool=true, bool=false); //Draw templates or input variables
    void Compare_TemplateShapes_Processes(TString, TString);

    void SetBranchAddress_SystVariationArray(TTree*, TString, vector<Double_t*>&, int); //Details in func comments
    void Merge_Templates_ByProcess(TString, TString, vector<TString>, bool=false);

    void Get_WCFit(WCFit*&, vector<string>*, vector<float>*, const vector<float>&, float);
    // void Fill_TH1EFT(TH1EFT*&, float, vector<string>*, vector<float>*, vector<float>, float);
    void Test_TH1EFT();

//--- MEMBERS
	bool stop_program;

	private :

//--- METHODS
    void Set_Luminosity(TString);

//--- MEMBERS
    TMVA::Reader *reader;
    TMVA::Reader *reader1; //1 reader for BDT xxx
    TMVA::Reader *reader2; //1 reader for BDT yyy
    TFModel* clfy1; //DNN classifier

    std::vector<TString> sample_list; //List of samples
    std::vector<TString> sample_groups; //List of "group naming", 1 per sample (many samples can have same naming)
    std::vector<TString> syst_list; //List of systematics stored as event weights
    std::vector<TString> systTree_list; //List of systematics stored as separate TTrees
	std::vector<TString> channel_list; //List of subchannels

	std::vector<TString> v_cut_name; //List of cut variable
	std::vector<TString> v_cut_def; //Cut definition
	std::vector<Float_t> v_cut_float; //Store variables needed for cut
	std::vector<Char_t> v_cut_char; //Categories are encoded into Char_t, because vector<bool> is 'broken' in C++
	std::vector<bool> v_cut_IsUsedForBDT; //true <-> the variable used for cut is also a input variable

	std::vector<TString> var_list; std::vector<Float_t> var_list_floats; //TMVA variables
	std::vector<TString> v_add_var_names; vector<Float_t> v_add_var_floats; //Additional vars only for CR plots

	std::vector<int> color_list;
	std::vector<TColor*> v_custom_colors;

	bool use_NeuralNetwork;
	TString classifier_name;
    TString DNN_inputLayerName, DNN_outputLayerName; int nNodes = 1; //DNN model params
    std::vector<TString> var_list_DNN; //Input features of DNN training may differ from those declared in 'analysis_main.cxx'
    std::vector<pair<float,float>> v_inputs_rescaling; //For now, can read rescaling params from DNN info file to rescale input features on the fly

    TString region; //Event category : "" / "tZq" / "ttZ" / "tWZ"
	TString categ_bool_name; //Name of boolean associated to category
	TString signal_process;
    TString dir_ntuples; //Path to base dir. containing Ntuples
    TString t_name;
	TString plot_extension;
    vector<TString> v_lumiYears;
    TString lumiName;
    double lumiValue;
	int nbins; //nbin_templates
	TString filename_suffix;
	bool show_pulls_ratio;
    bool is_blind;
    int nSampleGroups; //Nof sample groups (e.g. 'Rares',  ...)
    bool use_custom_colorPalette;

    //Systematics variations arrays //More details in comments of func Handle_SystVariationArray()
    double* array_PU;
    double* array_prefiringWeight;
    double* array_Btag;
    double* array_LepEff_mu;
    double* array_LepEff_el;
};

#endif
