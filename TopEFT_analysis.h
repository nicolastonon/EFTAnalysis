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

#include "ROCS/ROC_Plotter.h" //ROC functions

//Custom classes for EFT (see https://github.com/Andrew42/EFTGenReader/blob/maste)
#include "Utils/TH1EFT.h"
#include "Utils/WCPoint.h"
#include "Utils/WCFit.h"

using namespace std;

class TopEFT_analysis
{
	public :

	TopEFT_analysis(); //Default constructor
    TopEFT_analysis(vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<TString>, vector<bool>, vector<TString>, TString, vector<TString>, bool, TString, TString, TString, bool, TString, TString, vector<float>, vector<float>, bool, bool, int, bool, TString, bool, bool, bool, bool, bool, bool, bool);
	~TopEFT_analysis(); //Default destructor

//--- METHODS
	void Train_BDT(TString); //Train BDT
    void Produce_Templates(TString, bool, bool, bool, float, float, bool, bool); //Produce templates
    void Draw_Templates(bool, TString, bool, bool, TString="", bool=true, bool=false, TString="", bool=false, float* ymax_fixed1=NULL, float* ymax_fixed2=NULL); //Draw templates or input variables
    void Compare_TemplateShapes_Processes(TString, TString);

    void SetBranchAddress_SystVariationArray(TTree*, TString, vector<Double_t*>&, int); //Details in func comments
    void MergeSplit_Templates(bool, TString, vector<TString>, TString="", TString = "",bool=true);

    bool Get_VectorAllEvents_passMVACut(vector<int>&, vector<vector<float>>&, TString, TString, TTree*, TString, float, bool, bool, int, bool, int, TString="", bool=false, bool=false);

    void Make_PaperPlot_CommonRegions();
    void Make_PaperPlot_SignalRegions(TString);
    void Make_PaperPlot_ControlPlots();

    void Dump_Scores_allNNs();

//--- MEMBERS
	bool stop_program;

	private :

//--- METHODS
    void Set_Luminosity(TString);

//--- MEMBERS
    TString nominal_tree_name = "result"; //Name of the nominal tree to read in rootfiles

    TMVA::Reader* reader=NULL;
    //NB: if booking 2 BDTs, must make sure that they use the same input variables... or else, find some way to make it work in the code)
    TFModel* clfy1=NULL; //NN classifier
    TFModel* clfy2=NULL; //NN classifier

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

    std::vector<TString> var_list; std::vector<Float_t> var_list_floats; //Names + float storage of MVA input features
    std::vector<Float_t> var_list_floats_2; //float storage of MVA input features (additional MVA)
	std::vector<TString> v_add_var_names; vector<Float_t> v_add_var_floats; //Additional vars only for CR plots

    //CHANGED -- use vectors of pointers to read values (more flexible, can point several times to same address, etc.)
    std::vector<Float_t*> var_list_pfloats; //Names + float storage of MVA input features
    std::vector<Float_t*> var_list_pfloats_2; //float storage of MVA input features (additional MVA)
	std::vector<Float_t*> pv_add_var_floats; //Additional vars only for CR plots

	std::vector<int> color_list;
	std::vector<TColor*> v_custom_colors;

    bool use_NeuralNetwork;
    TString classifier_name;
    bool make_SMvsEFT_templates_plots; //When making/plotting templates only: false <-> consider SM vs SM templates; true <-> consider SM vs EFT templates (--> different input/output files)
    int categorization_strategy; //1 <-> define SRtZq/SRttZ with different jet multiplicities, apply dedicated binary classifiers; 2 <-> apply multi-classifier in merged SR; 0 <-> testing: read tmp MVA, no categ.
    bool use_specificMVA_eachYear; //true <-> look for year-specific MVA weight files
    TString NN_strategy; //'MVA_SM' <-> SM vs SM classifier; 'MVA_EFT' <-> fixed-point SM vs EFT classifier; 'MVA_param' <-> parametrized SM vs EFT classifier

    //-- Default NN
    TString NN_inputLayerName, NN_outputLayerName; int NN_nNodes = 1; //NN model params
    vector<TString> v_NN_nodeLabels; //Label(s) of the node(s), e.g. 'tZq'/'ttZ'/'Backgrounds'
    std::vector<TString> var_list_NN; //Input features of NN training may differ from those declared in 'analysis_main.cxx'
    // std::vector<pair<float,float>> v_inputs_rescaling; //For now, can read rescaling params from NN info file to rescale input features on the fly
    int NN_iMaxNode; //Index of the multiclass DNN node for which a given event has its max value
    std::vector<float> minmax_bounds; //Read the min/max prediction values evaluted by the NN (to adapt the template binning and avoid empty bins at boundaries)
    // std::vector<pair<float,float>> minmax_bounds; //Read the min/max prediction values evaluted by the NN (to adapt the template binning and avoid empty bins at boundaries) for all nodes

    //-- Additional NN -- in case we need to read 2 NNs in the same function (e.g. to classify both SM vs SM and SM vs EFT)
    TString NN2_inputLayerName, NN2_outputLayerName, NN2_strategy; int NN2_nNodes = 1; //NN model params
    vector<TString> v_NN2_nodeLabels; //Label(s) of the node(s), e.g. 'tZq'/'ttZ'/'Backgrounds'
    std::vector<TString> var_list_NN2; //Input features of NN training may differ from those declared in 'analysis_main.cxx'
    int NN2_iMaxNode; //Index of the multiclass DNN node for which a given event has its max value
    std::vector<int> v_varIndices_inMVA1; //If we read 2 MVAs at once, store indices of MVA1 variables which are also used by MVA2 --> Get values from MVA1 floats (can set address only once, not twice)
    std::vector<float> minmax_bounds2; //Read the min/max prediction values evaluted by the NN (to adapt the template binning and avoid empty bins at boundaries) for all nodes

    TString region; //Event category : "" / "tZq" / "ttZ" / "tWZ"
	TString categ_bool_name; //Name of boolean associated to category
	TString signal_process;
    TString dir_ntuples; //Path to base dir. containing Ntuples
    TString dir_ntuples_SR; //Path to dir. containing split Ntuples containing only SR events
    bool use_optimized_ntuples; //Set automatically in the constructor; true <-> special case: adapt code (+ class members) to deal with 'optimized' ntuples rather than the full list of ntuples (must have been produced beforehand with [Split_AllNtuples_ByCategory] code)
    TString plot_extension;
    vector<TString> v_lumiYears;
    TString lumiName;
    double lumiValue;
	int nbins; //nbin_templates
	TString filename_suffix;
	bool show_pulls_ratio;
    bool is_blind;
    int nSampleGroups; //Nof sample groups (e.g. 'Rares',  ...)
    bool use_DD_NPL;
    bool use_SMdiffAnalysis_strategy;
    bool make_fixedRegions_templates;

    //Systematics variations arrays //More details in comments of func Handle_SystVariationArray()
    double* array_PU;
    double* array_prefiringWeight;
    double* array_Btag;
    double* array_jetPileupID;
    double* array_fakeFactor;
    double* array_ME;
    double* array_alphaS;
    double* array_PDFtotal;
    double* array_partonShower;
    double* array_LepEffLoose_mu;
    double* array_LepEffTight_mu;
    double* array_LepEffLoose_el;
    double* array_LepEffTight_el;

    //-- For parameterized NN
    bool scanOperators_paramNN; //true <-> if considering a parametrized NN, multiple templates and plots will be created on a 1D or 2D grid of points (instead of a single point)
    TString operator_scan1, operator_scan2; //First and second EFT operators to scan
    int idx1_operator_scan1, idx1_operator_scan2;
    int idx2_operator_scan1, idx2_operator_scan2;
    vector<float> v_WCs_operator_scan1, v_WCs_operator_scan2; //Grid points for first and second scanned EFT operators

    vector<vector<float>> v_njets_SF_tZq; //May be filled with helper func 'Get_nJets_SF' to apply a shape uncertainty to PrivMC_tZq based on jet multiplicity disagreements w.r.t. central sample

    bool process_samples_byGroup; //true <-> read grouped samples (if already hadded together), else read individual samples and combine them when creating histograms if needed (default)

    bool split_EFTtemplates_perBin; //true <-> will also store separately each individual bin of SMvsEFT templates (--> for easy EFT parameterization in combine)

    bool use_paperStyle; //Enforce CMS paper styles

    bool use_NN_SRother; //Use NN-bkg node instead of mTW in SRother (testing)
};

#endif
