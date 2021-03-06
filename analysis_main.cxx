#include "TopEFT_analysis.h"

using namespace std;

int main(int argc, char **argv)
{
//---------------------------------------------------------------------------
//  #######  ########  ######## ####  #######  ##    ##  ######
// ##     ## ##     ##    ##     ##  ##     ## ###   ## ##    ##
// ##     ## ##     ##    ##      ##     ## ####  ## ##
// ##     ## ########     ##     ##  ##     ## ## ## ##  ######
// ##     ## ##           ##     ##  ##     ## ##  ####       ##
// ##     ## ##           ##     ##  ##     ## ##   ### ##    ##
//  #######  ##           ##    ####  #######  ##    ##  ######
//---------------------------------------------------------------------------

    //-- M A I N    A N A L  S I S    O P T I O N S --
    TString signal_process = "tZq"; //'tZq', 'ttZ', 'tWZ'
    TString region = ""; //Select a specific event category : '' (all preselected events) / 'tZq' / 'ttZ' / 'signal'
    bool use_systematics = true; //true <-> will compute/store systematics selected below
    bool split_JEC = true; //true <-> store all split JEC variations (in addition to total JEC -- much slower!)
    bool is_blind = false; //true <-> don't read/store data events
    bool make_fixedRegions_templates = false; //true <-> overrides some options, to enforce the creation of templates in SR/CR regions which are not expected to change (for now: ttZ 4l SR / WZ CR / ZZ CR / DY CR)

    //-- N T U P L E S --
    bool use_DD_NPL = true; //true <-> use data-driven fakes sample; otherwise use MC (ttbar+DY)
    bool include_PrivMC_samples = true; //true <-> process private SMEFT signal samples (necessary e.g. for limit-setting, but much slower)
    bool include_central_samples = false; //true <-> process central signal samples
    bool process_samples_byGroup = true; //true <-> read grouped samples (if already hadded together), else read individual samples and combine them when creating histograms if needed (default)

    //-- M V A    S T R A T E G Y --
    TString classifier_name = "NN"; //'BDT' or 'NN'
    bool use_specificMVA_eachYear = false; //true <-> look for year-specific MVA weight files

    bool make_SMvsEFT_templates_plots = true; //false <-> templates & plots are produced for SM scenario only (separate SM processes); true <-> consider SM vs EFT scenario (and apply beforehand the chosen categorization strategy)
        int categorization_strategy = 2; //1 <-> define SRtZq/SRttZ with different jet multiplicities, apply dedicated binary classifiers (events passing the cut fall into SRtZq/SRttZ, others into SRother); 2 <-> apply multi-classifier in merged SR (events fall into SRtZq/SRttZ/CR based on max node); 0 <-> testing: read tmp MVA, no categ. (retain all events, can't use multiple nodes simultaneously)
        float cut_value_tZq = 0.5, cut_value_ttZ = 0.3; //Hard-coded cut values to apply -- for templates (automatic) and plots (user-option)
        bool keep_aboveCut = true; //true <-> only keep events satisfying x>=cut
        bool also_applyCut_onMaxNodeValue = false; //true <-> for SM vs EFT strategy 2, don't only look for the max node, but also apply a cut on the corresponding node value (cut set here)

    bool scanOperators_paramNN = false; //true <-> if considering a parameterized NN, multiple templates and plots will be created on a 1D or 2D grid of points (instead of a single point)
        TString operator1 = ""; //First operator to scan (required)
        TString operator2 = ""; //Second operator to scan (optional)
        // vector<float> v_W = {-3, -2, -1.5, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.5, 2, 3}; //Grid points for first operator (required)
        vector<float> v_WCs_operator_scan1 = {-1.5, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.5}; //Grid points for first operator (required)
        vector<float> v_WCs_operator_scan2 = {}; //Grid points for second operator (optional)

    //-- T E M P L A T E S --
    bool use_NN_cpq3_SRttZ = true; //true <-> use a dedicated NN-cpq3-SRttZ (instead of the ttZ output node of NN-SM, previous default)
    bool use_NN_SRother = false; //true <-> use NN-bkg node instead of mTW in SRother
    bool split_EFTtemplates_perBin = true; //true <-> will also store separately each individual bin of SMvsEFT templates (--> for easy EFT parameterization in combine)
    bool split_analysis_by_channel = false; //true <-> will *also* produce templates/histos/plots for each subchannel (defined below)
    TString template_name = "NN_SM"; //Ex: 'BDT', 'NN', 'categ' (nbjet/njet bins), 'Zpt', 'ZptCos', 'NN_ctz', 'NN_all', 'NN_SM', ...

    //-- P L O T T I N G --
    bool show_pulls_ratio = false; //true <-> bottom pad shows pull; else shows data/mc ratio (w/ errors)
    TString plot_extension = ".png"; //Default file extension for plots //For 'paper plots', .pdf extension is enforced
    bool plot_onlyMaxNodeEvents = true; //For multiclass NN-SM template plots only (EFT strategy 2): true <-> only include events if they have their max output value in the corresponding node
    bool plot_onlyMVACutEvents = false; //For binary MVA-SM templates plots only: true <-> only include events which pass the specified tZq or ttZ cut values
    bool plot_EFTscan_eachPoint = true; //true <-> if making template plots for a parameterized NN, will make 1 plot per considered EFT point (if histograms are found)
    TString nominal_tree_name = "result"; //Name of the nominal tree to read in rootfiles
    bool use_paperStyle = true; //true <-> Add 'preliminary' label on plots


//-----------------------------------------------------------------------------------------
// ##       ##     ## ##     ## #### ##    ##  #######   ######  #### ######## ##    ##
// ##       ##     ## ###   ###  ##  ###   ## ##     ## ##    ##  ##     ##     ##  ##
// ##       ##     ## #### ####  ##  ####  ## ##     ## ##        ##     ##      ####
// ##       ##     ## ## ### ##  ##  ## ## ## ##     ##  ######   ##     ##       ##
// ##       ##     ## ##     ##  ##  ##  #### ##     ##       ##  ##     ##       ##
// ##       ##     ## ##     ##  ##  ##   ### ##     ## ##    ##  ##     ##       ##
// ########  #######  ##     ## #### ##    ##  #######   ######  ####    ##       ##
//-----------------------------------------------------------------------------------------
//-- Choose here what data you want to consider (separate ntuples, templates, ... per year)
//Naming convention enforced : 2016+2017 <-> "201617" ; etc.; 2016+2017+2018 <-> "Run2"
//NB : years must be placed in the right order !

	vector<TString> set_lumi_years;
    set_lumi_years.push_back("2016");
    set_lumi_years.push_back("2017");
    set_lumi_years.push_back("2018");


//-----------------------------------------------------------------------------------------
//   ######  ##     ## ########  ######
//  ##    ## ##     ##    ##    ##    ##
//  ##       ##     ##    ##    ##
//  ##       ##     ##    ##     ######
//  ##       ##     ##    ##          ##
//  ##    ## ##     ##    ##    ##    ##
//   ######   #######     ##     ######
//-----------------------------------------------------------------------------------------
//ex: set_v_cut_name.push_back("NBJets"); set_v_cut_def.push_back(">0 && <4"); set_v_cut_IsUsedForBDT.push_back(false);
//NB : if variable is to be used in BDT, can't ask it to take unique value (==) !

	vector<TString> set_v_cut_name;
	vector<TString> set_v_cut_def;
	vector<bool> set_v_cut_IsUsedForBDT;

    // set_v_cut_name.push_back("nJets");  set_v_cut_def.push_back("==3 || ==4"); set_v_cut_IsUsedForBDT.push_back(false);
    // set_v_cut_name.push_back("passedBJets");  set_v_cut_def.push_back("==1"); set_v_cut_IsUsedForBDT.push_back(false); //enforce final tZq 3l selection
    // set_v_cut_name.push_back("is_tzq_SR");  set_v_cut_def.push_back("==1"); set_v_cut_IsUsedForBDT.push_back(false);
    // set_v_cut_name.push_back("is_signal_SR");  set_v_cut_def.push_back("==1"); set_v_cut_IsUsedForBDT.push_back(false);


//---------------------------------------------------------------------------
//  ######  ##     ##    ###    ##    ## ##    ## ######## ##        ######
// ##    ## ##     ##   ## ##   ###   ## ###   ## ##       ##       ##    ##
// ##       ##     ##  ##   ##  ####  ## ####  ## ##       ##       ##
// ##       ######### ##     ## ## ## ## ## ## ## ######   ##        ######
// ##       ##     ## ######### ##  #### ##  #### ##       ##             ##
// ##    ## ##     ## ##     ## ##   ### ##   ### ##       ##       ##    ##
//  ######  ##     ## ##     ## ##    ## ##    ## ######## ########  ######
//---------------------------------------------------------------------------
    std::vector<TString > thechannellist;

    thechannellist.push_back(""); //KEEP ! (<-> no subcategorization, used for plots, etc.)

    if(split_analysis_by_channel)
    {
        thechannellist.push_back("uuu");
        thechannellist.push_back("uue");
        thechannellist.push_back("eeu");
        thechannellist.push_back("eee");
    }


//---------------------------------------------------------------------------
//  ######     ###    ##     ## ########  ##       ########  ######
// ##    ##   ## ##   ###   ### ##     ## ##       ##       ##    ##
// ##        ##   ##  #### #### ##     ## ##       ##       ##
//  ######  ##     ## ## ### ## ########  ##       ######    ######
//       ## ######### ##     ## ##        ##       ##             ##
// ##    ## ##     ## ##     ## ##        ##       ##       ##    ##
//  ######  ##     ## ##     ## ##        ######## ########  ######
//---------------------------------------------------------------------------

//-- List of sample names (as found in ./input_ntuples) //thesamplegroups <-> can merge multiple ntuples into same group (plotting)
    vector<TString> thesamplelist, thesamplegroups;

    //DATA (single sample, in first position)
    thesamplelist.push_back("DATA"); thesamplegroups.push_back("DATA");

    //Private MC production including EFT weights //Should be stored as separate sample groups (don't merge)
    if(include_PrivMC_samples)
    {
        thesamplelist.push_back("PrivMC_tZq"); thesamplegroups.push_back("PrivMC_tZq");
        thesamplelist.push_back("PrivMC_tWZ"); thesamplegroups.push_back("PrivMC_tWZ");
        thesamplelist.push_back("PrivMC_ttZ"); thesamplegroups.push_back("PrivMC_ttZ");
    }

    //Signals (central samples)
    if(include_central_samples)
    {
        thesamplelist.push_back("tZq"); thesamplegroups.push_back("tZq");
        thesamplelist.push_back("tWZ"); thesamplegroups.push_back("tWZ");
        thesamplelist.push_back("ttZ"); thesamplegroups.push_back("ttZ");
    }

    //t(t)X
    thesamplelist.push_back("ttZ_M1to10"); thesamplegroups.push_back("tX"); //Separate from ttZ because misses PDF weights, etcc.
    thesamplelist.push_back("tHq"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("tHW"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("ttH"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("ttW"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("ttZZ"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("ttWW"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("ttWZ"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("ttZH"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("ttWH"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("tttt"); thesamplegroups.push_back("tX");
    thesamplelist.push_back("ttHH"); thesamplegroups.push_back("tX");

    //WZ
    thesamplelist.push_back("WZ"); thesamplegroups.push_back("WZ");

    //VV(V)
    thesamplelist.push_back("ZZ4l"); thesamplegroups.push_back("VVV");
    thesamplelist.push_back("ZZZ"); thesamplegroups.push_back("VVV");
    thesamplelist.push_back("WZZ"); thesamplegroups.push_back("VVV");
    thesamplelist.push_back("WWW"); thesamplegroups.push_back("VVV");
    thesamplelist.push_back("WWZ"); thesamplegroups.push_back("VVV");

    //X+g
    thesamplelist.push_back("TTGamma_Dilep"); thesamplegroups.push_back("XG");
    thesamplelist.push_back("ZGToLLG_01J"); thesamplegroups.push_back("XG");
    // thesamplelist.push_back("tGJets"); thesamplegroups.push_back("XG"); //Useless ?
    // thesamplelist.push_back("WGToLNuG"); thesamplegroups.push_back("XG"); //Useless ?

    //NPL (Fakes)
    if(use_DD_NPL) //Data-driven
    {
        thesamplelist.push_back("NPL_DATA"); thesamplegroups.push_back("NPL");
        thesamplelist.push_back("NPL_MC"); thesamplegroups.push_back("NPL"); //Substract prompt MC contribution from AR
    }
    else //MC
    {
        thesamplelist.push_back("DY"); thesamplegroups.push_back("NPL");
        thesamplelist.push_back("TTba"); thesamplegroups.push_back("NPL");
    }


//---------------------------------------------------------------------------
// ########  ########  ########       ##     ##    ###    ########   ######
// ##     ## ##     ##    ##          ##     ##   ## ##   ##     ## ##    ##
// ##     ## ##     ##    ##          ##     ##  ##   ##  ##     ## ##
// ########  ##     ##    ##          ##     ## ##     ## ########   ######
// ##     ## ##     ##    ##           ##   ##  ######### ##   ##         ##
// ##     ## ##     ##    ##            ## ##   ##     ## ##    ##  ##    ##
// ########  ########     ##             ###    ##     ## ##     ##  ######
//---------------------------------------------------------------------------
//Variables used in BDT training (and evaluation)
//For NN, will read necessary input files in logfile, and include them automatically

    std::vector<TString > thevarlist;
    thevarlist.push_back("recoZ_Pt");


//---------------------------------------------------------------------------
//  #######  ######## ##     ## ######## ########       ##     ##    ###    ########   ######
// ##     ##    ##    ##     ## ##       ##     ##      ##     ##   ## ##   ##     ## ##    ##
// ##     ##    ##    ##     ## ##       ##     ##      ##     ##  ##   ##  ##     ## ##
// ##     ##    ##    ######### ######   ########       ##     ## ##     ## ########   ######
// ##     ##    ##    ##     ## ##       ##   ##         ##   ##  ######### ##   ##         ##
// ##     ##    ##    ##     ## ##       ##    ##         ## ##   ##     ## ##    ##  ##    ##
//  #######     ##    ##     ## ######## ##     ##         ###    ##     ## ##     ##  ######
//---------------------------------------------------------------------------
//Can add additionnal vars which are NOT used in TMVA NOR for cuts, only for CR plots !
//NOTE : Branch can be linked to only *one* variable via SetBranchAddress ; if additional variable is already present in other variable vector, it is removed from this vector !

    vector<TString> set_v_add_var_names;

    set_v_add_var_names.push_back("channel");
    set_v_add_var_names.push_back("njets");
    set_v_add_var_names.push_back("nbjets");
    set_v_add_var_names.push_back("metEt");

    set_v_add_var_names.push_back("recoZ_Pt");
    set_v_add_var_names.push_back("recoZ_Eta");
    set_v_add_var_names.push_back("recoZ_dPhill");
    set_v_add_var_names.push_back("mHT");
    set_v_add_var_names.push_back("mTW");
    set_v_add_var_names.push_back("Mass_3l");
    set_v_add_var_names.push_back("lAsymmetry");
    set_v_add_var_names.push_back("jPrimeAbsEta");
    set_v_add_var_names.push_back("maxEtaJet");
    set_v_add_var_names.push_back("maxDeepJet");
    set_v_add_var_names.push_back("cosThetaStarPolTop");
    set_v_add_var_names.push_back("cosThetaStarPolZ");
    set_v_add_var_names.push_back("maxDiJet_Pt");
    set_v_add_var_names.push_back("maxDiJet_M");
    set_v_add_var_names.push_back("jprime_Pt");
    set_v_add_var_names.push_back("recoLepTop_Pt");
    set_v_add_var_names.push_back("recoLepTop_Eta");
    set_v_add_var_names.push_back("recoLepTopLep_Pt");
    set_v_add_var_names.push_back("TopZsystem_M");
    set_v_add_var_names.push_back("mbjMax");

    set_v_add_var_names.push_back("dEta_tjprime");
    set_v_add_var_names.push_back("dR_tZ");
    set_v_add_var_names.push_back("dR_ZlW");
    set_v_add_var_names.push_back("dR_blW");
    set_v_add_var_names.push_back("dR_lWjprime");
    set_v_add_var_names.push_back("dEta_Zjprime");
    set_v_add_var_names.push_back("dR_tClosestLep");
    set_v_add_var_names.push_back("dEta_bjprime");
    set_v_add_var_names.push_back("dEta_lWjprime");
    set_v_add_var_names.push_back("dR_jprimeClosestLep");
    set_v_add_var_names.push_back("dR_bjprime");
    set_v_add_var_names.push_back("dR_Zjprime");

    set_v_add_var_names.push_back("lep1_pt");
    set_v_add_var_names.push_back("lep2_pt");
    set_v_add_var_names.push_back("jet1_pt");
    set_v_add_var_names.push_back("jet2_pt");

    set_v_add_var_names.push_back("lep3_pt");
    set_v_add_var_names.push_back("lep1_eta");
    set_v_add_var_names.push_back("lep2_eta");
    set_v_add_var_names.push_back("lep3_eta");
    set_v_add_var_names.push_back("lep1_phi");
    set_v_add_var_names.push_back("lep2_phi");
    set_v_add_var_names.push_back("lep3_phi");
    set_v_add_var_names.push_back("jet3_pt");
    set_v_add_var_names.push_back("jet1_eta");
    set_v_add_var_names.push_back("jet2_eta");
    set_v_add_var_names.push_back("jet3_eta");
    set_v_add_var_names.push_back("jet1_phi");
    set_v_add_var_names.push_back("jet2_phi");
    set_v_add_var_names.push_back("jet3_phi");
    set_v_add_var_names.push_back("jet1_DeepJet");
    set_v_add_var_names.push_back("jet2_DeepJet");
    set_v_add_var_names.push_back("jet3_DeepJet");


//---------------------------------------------------------------------------
//  ######  ##    ##  ######  ######## ######## ##     ##    ###    ######## ####  ######   ######
// ##    ##  ##  ##  ##    ##    ##    ##       ###   ###   ## ##      ##     ##  ##    ## ##    ##
// ##         ####   ##          ##    ##       #### ####  ##   ##     ##     ##  ##       ##
//  ######     ##     ######     ##    ######   ## ### ## ##     ##    ##     ##  ##        ######
//       ##    ##          ##    ##    ##       ##     ## #########    ##     ##  ##             ##
// ##    ##    ##    ##    ##    ##    ##       ##     ## ##     ##    ##     ##  ##    ## ##    ##
//  ######     ##     ######     ##    ######## ##     ## ##     ##    ##    ####  ######   ######
//---------------------------------------------------------------------------

    vector<TString> theSystWeights; //List of systematics implemented as event weights
    theSystWeights.push_back(""); //KEEP ! (<-> nominal event weight)

    vector<TString> theSystTree; //List of systematics implemented as separate TTrees
    theSystTree.push_back(""); //KEEP ! (<-> nominal TTree)

    if(use_systematics) //Define here the list of syst to run //Missing: JERC, leptonID, ME, PDFs, ...
    {
        //-- Implemented as separate TTrees

        //* Total JEV
        theSystTree.push_back("JESDown"); theSystTree.push_back("JESUp");
        theSystTree.push_back("JERDown"); theSystTree.push_back("JERUp");
        theSystTree.push_back("METDown"); theSystTree.push_back("METUp");

        //* Split JEC sources
        if(split_JEC)
        {
            theSystTree.push_back("AbsoluteStatDown"); theSystTree.push_back("AbsoluteStatUp");
            theSystTree.push_back("AbsoluteScaleDown"); theSystTree.push_back("AbsoluteScaleUp");
            theSystTree.push_back("AbsoluteMPFBiasDown"); theSystTree.push_back("AbsoluteMPFBiasUp");
            theSystTree.push_back("FragmentationDown"); theSystTree.push_back("FragmentationUp");
            theSystTree.push_back("SinglePionECALDown"); theSystTree.push_back("SinglePionECALUp");
            theSystTree.push_back("SinglePionHCALDown"); theSystTree.push_back("SinglePionHCALUp");
            theSystTree.push_back("FlavorQCDDown"); theSystTree.push_back("FlavorQCDUp");
            theSystTree.push_back("TimePtEtaDown"); theSystTree.push_back("TimePtEtaUp");
            theSystTree.push_back("RelativeJEREC1Down"); theSystTree.push_back("RelativeJEREC1Up");
            theSystTree.push_back("RelativeJEREC2Down"); theSystTree.push_back("RelativeJEREC2Up");
            theSystTree.push_back("RelativeJERHFDown"); theSystTree.push_back("RelativeJERHFUp");
            theSystTree.push_back("RelativePtBBDown"); theSystTree.push_back("RelativePtBBUp");
            theSystTree.push_back("RelativePtEC1Down"); theSystTree.push_back("RelativePtEC1Up");
            theSystTree.push_back("RelativePtEC2Down"); theSystTree.push_back("RelativePtEC2Up");
            theSystTree.push_back("RelativePtHFDown"); theSystTree.push_back("RelativePtHFUp");
            theSystTree.push_back("RelativeBalDown"); theSystTree.push_back("RelativeBalUp");
            theSystTree.push_back("RelativeSampleDown"); theSystTree.push_back("RelativeSampleUp");
            theSystTree.push_back("RelativeFSRDown"); theSystTree.push_back("RelativeFSRUp");
            theSystTree.push_back("RelativeStatFSRDown"); theSystTree.push_back("RelativeStatFSRUp");
            theSystTree.push_back("RelativeStatECDown"); theSystTree.push_back("RelativeStatECUp");
            theSystTree.push_back("RelativeStatHFDown"); theSystTree.push_back("RelativeStatHFUp");
            theSystTree.push_back("PileUpDataMCDown"); theSystTree.push_back("PileUpDataMCUp");
            theSystTree.push_back("PileUpPtRefDown"); theSystTree.push_back("PileUpPtRefUp");
            theSystTree.push_back("PileUpPtBBDown"); theSystTree.push_back("PileUpPtBBUp");
            theSystTree.push_back("PileUpPtEC1Down"); theSystTree.push_back("PileUpPtEC1Up");
            theSystTree.push_back("PileUpPtEC2Down"); theSystTree.push_back("PileUpPtEC2Up");
            theSystTree.push_back("PileUpPtHFDown"); theSystTree.push_back("PileUpPtHFUp");
        }

        //-- Implementend as event weights

        //* EXP
        theSystWeights.push_back("PUDown"); theSystWeights.push_back("PUUp");
        theSystWeights.push_back("prefireDown"); theSystWeights.push_back("prefireUp");
        theSystWeights.push_back("BtagHFDown"); theSystWeights.push_back("BtagHFUp");
        theSystWeights.push_back("BtagLFDown"); theSystWeights.push_back("BtagLFUp");
        theSystWeights.push_back("BtagHFstats1Down"); theSystWeights.push_back("BtagHFstats1Up");
        theSystWeights.push_back("BtagHFstats2Down"); theSystWeights.push_back("BtagHFstats2Up");
        theSystWeights.push_back("BtagLFstats1Down"); theSystWeights.push_back("BtagLFstats1Up");
        theSystWeights.push_back("BtagLFstats2Down"); theSystWeights.push_back("BtagLFstats2Up");
        theSystWeights.push_back("BtagCFerr1Down"); theSystWeights.push_back("BtagCFerr1Up");
        theSystWeights.push_back("BtagCFerr2Down"); theSystWeights.push_back("BtagCFerr2Up");
        theSystWeights.push_back("njets_tZqDown"); theSystWeights.push_back("njets_tZqUp");
        // theSystWeights.push_back("jetPUIDEffDown"); theSystWeights.push_back("jetPUIDEffUp");
        // theSystWeights.push_back("jetPUIDMTDown"); theSystWeights.push_back("jetPUIDMTUp");

        theSystWeights.push_back("FRm_normDown"); theSystWeights.push_back("FRm_normUp"); //FR from ttH: 3*2*2 sets of variations
        theSystWeights.push_back("FRm_ptDown"); theSystWeights.push_back("FRm_ptUp");
        theSystWeights.push_back("FRm_beDown"); theSystWeights.push_back("FRm_beUp");
        theSystWeights.push_back("FRe_normDown"); theSystWeights.push_back("FRe_normUp");
        theSystWeights.push_back("FRe_ptDown"); theSystWeights.push_back("FRe_ptUp");
        theSystWeights.push_back("FRe_beDown"); theSystWeights.push_back("FRe_beUp");

        theSystWeights.push_back("LepEff_muLooseDown"); theSystWeights.push_back("LepEff_muLooseUp");
        theSystWeights.push_back("LepEff_muTightDown"); theSystWeights.push_back("LepEff_muTightUp");
        theSystWeights.push_back("LepEff_elLooseDown"); theSystWeights.push_back("LepEff_elLooseUp");
        theSystWeights.push_back("LepEff_elTightDown"); theSystWeights.push_back("LepEff_elTightUp");

        //* THEORY
        theSystWeights.push_back("PDFDown"); theSystWeights.push_back("PDFUp"); //Signals only
        theSystWeights.push_back("alphasDown"); theSystWeights.push_back("alphasUp"); //Signals only
        theSystWeights.push_back("MEDown"); theSystWeights.push_back("MEUp"); //Signals only
        theSystWeights.push_back("ISRDown"); theSystWeights.push_back("ISRUp");
        theSystWeights.push_back("FSRDown"); theSystWeights.push_back("FSRUp");
    }


//---------------------------------------------------------------------------
// ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##        ######     ###    ##       ##        ######
// ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ##       ##    ##   ## ##   ##       ##       ##    ##
// ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ##       ##        ##   ##  ##       ##       ##
// ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##       ##       ##     ## ##       ##        ######
// ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##       ######### ##       ##             ##
// ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ###       ##    ## ##     ## ##       ##       ##    ##
// ##        #######  ##    ##  ######     ##    ####  #######  ##    ##        ######  ##     ## ######## ########  ######
//---------------------------------------------------------------------------

//*** CHOOSE HERE FROM BOOLEANS WHAT YOU WANT TO DO !

//-----------------    TRAINING
    bool train_BDT = false; //Train selected BDT in selected region (with events in training category)

//-----------------    TEMPLATES CREATION
    bool create_templates = false; //Create MVA templates

//-----------------    CONTROL HISTOGRAMS
    bool create_inputVar_histograms = false; //Create histograms of input variables, for plotting

//-----------------    PLOTS
    TString plotChannel = ""; //Can choose to plot particular subchannel //uu, ue, ee, ...

    bool draw_templates = false; //Plot templates of selected BDT, in selected region
        bool prefit = true; //true <-> plot prefit templates ; else postfit (requires combine output file)
        bool use_combine_file = false; //true <-> use MLF output file from Combine (can get postfit plots, total error, etc.)

    bool draw_input_vars = false; //Plot input variables
        bool draw_input_allChannels = false; //true <-> also draw for eachs split channel

    bool compare_template_shapes = false;

//-----------------    PAPER PLOTS

    bool make_paperPlot_commonRegions = false;
    bool make_paperPlot_signalRegions = false;
    bool make_paperPlot_controlPlots = false;

    bool make_animation = false;

//-----------------    OTHER

    bool run_test = false; //Run test function, if any

//-----------------













//--------------------------------------------
//--------------------------------------------
//--- Automated from here -- no need to modify
//--------------------------------------------
//--------------------------------------------

//--------------------------------------------
//    ###    ##     ## ########  #######  ##     ##    ###    ######## ####  ######
//   ## ##   ##     ##    ##    ##     ## ###   ###   ## ##      ##     ##  ##    ##
//  ##   ##  ##     ##    ##    ##     ## #### ####  ##   ##     ##     ##  ##
// ##     ## ##     ##    ##    ##     ## ## ### ## ##     ##    ##     ##  ##
// ######### ##     ##    ##    ##     ## ##     ## #########    ##     ##  ##
// ##     ## ##     ##    ##    ##     ## ##     ## ##     ##    ##     ##  ##    ##
// ##     ##  #######     ##     #######  ##     ## ##     ##    ##    ####  ######
//--------------------------------------------

//-- Apply choices given via command line, if any
	Apply_CommandArgs_Choices(argc, argv, set_lumi_years, region, template_name);

    // int nthreads = 4; ROOT::EnableImplicitMT(nthreads); //Enable multi-threading (I have 8 available threads)

    //#############################################
    //  CREATE INSTANCE OF CLASS & INITIALIZE
    //#############################################

    TopEFT_analysis* theAnalysis = new TopEFT_analysis(thesamplelist, thesamplegroups, theSystWeights, theSystTree, thechannellist, thevarlist, set_v_cut_name, set_v_cut_def, set_v_cut_IsUsedForBDT, set_v_add_var_names, plot_extension, set_lumi_years, show_pulls_ratio, region, signal_process, classifier_name, scanOperators_paramNN, operator1, operator2, v_WCs_operator_scan1, v_WCs_operator_scan2, make_SMvsEFT_templates_plots, is_blind, categorization_strategy, use_specificMVA_eachYear, nominal_tree_name, use_DD_NPL, make_fixedRegions_templates, process_samples_byGroup, split_EFTtemplates_perBin, use_paperStyle, use_NN_SRother, use_NN_cpq3_SRttZ);
    if(theAnalysis->stop_program) {cout<<"=== 'stop_program=true' ---> Exiting... "<<endl; return 1;}

    //#############################################
    // TRAINING
    //#############################################

    if(train_BDT) {theAnalysis->Train_BDT("");}

    //#############################################
    //  TEMPLATES CREATION
    //#############################################

    if(create_templates) {theAnalysis->Produce_Templates(template_name, false, plot_onlyMaxNodeEvents, plot_onlyMVACutEvents, cut_value_tZq, cut_value_ttZ, keep_aboveCut, also_applyCut_onMaxNodeValue);}

    //#############################################
    //  CONTROL HISTOGRAMS
    //#############################################

    if(create_inputVar_histograms) {theAnalysis->Produce_Templates(template_name, true, plot_onlyMaxNodeEvents, plot_onlyMVACutEvents, cut_value_tZq, cut_value_ttZ, keep_aboveCut, also_applyCut_onMaxNodeValue);}

    //#############################################
    //  DRAW PLOTS
    //#############################################

    //All channels
    if(draw_templates)
    {
        //-- Make plots for *each EFT point* considered in EFT scan
        if(plot_EFTscan_eachPoint && scanOperators_paramNN)
        {
            float* ymax_fixed1 = new float; float* ymax_fixed2 = new float; bool store_ymax_fixed = true;
            for(int iop1=0; iop1<v_WCs_operator_scan1.size(); iop1++)
            {
                for(int iop2=0; iop2<v_WCs_operator_scan2.size(); iop2++)
                {
                    TString EFTpoint = operator1+"_"+Convert_Number_To_TString(v_WCs_operator_scan1[iop1]);
                    if(operator2=="" && iop2>0) {break;} //Only loop on 1 operator
                    else if(operator2!="") {EFTpoint+= "_" + operator2 + "_" + Convert_Number_To_TString(v_WCs_operator_scan2[iop2]);}

                    theAnalysis->Draw_Templates(false, plotChannel, plot_onlyMaxNodeEvents, plot_onlyMVACutEvents, template_name, prefit, use_combine_file, EFTpoint, store_ymax_fixed, ymax_fixed1, ymax_fixed2); //chosen channel
                    store_ymax_fixed = false;
                }
            }
            delete ymax_fixed1; ymax_fixed1 = NULL; delete ymax_fixed2; ymax_fixed2 = NULL;
        }
        //-- Default: make plots corresponding to selected user-options
        else
        {
            theAnalysis->Draw_Templates(false, plotChannel, plot_onlyMaxNodeEvents, plot_onlyMVACutEvents, template_name, prefit, use_combine_file); //chosen channel

            if(plotChannel == "") //By default, also want to plot templates in subchannels
            {
                for(int ichan=1; ichan<thechannellist.size(); ichan++)
                {
                    theAnalysis->Draw_Templates(false, thechannellist[ichan], plot_onlyMaxNodeEvents, plot_onlyMVACutEvents, template_name, prefit, use_combine_file);
                }
            }
        }
    }

    if(draw_input_vars)
    {
        theAnalysis->Draw_Templates(true, plotChannel, plot_onlyMaxNodeEvents, plot_onlyMVACutEvents, "", prefit, use_combine_file);
        if(draw_input_allChannels)
        {
            for(int ichan=1; ichan<thechannellist.size(); ichan++)
            {
                theAnalysis->Draw_Templates(true, thechannellist[ichan], plot_onlyMaxNodeEvents, plot_onlyMVACutEvents, "", prefit, use_combine_file);
            }
        }
    }

    //#############################################
    //  PAPER PLOTS
    //#############################################

    if(make_paperPlot_commonRegions) {theAnalysis->Make_PaperPlot_CommonRegions();}
    if(make_paperPlot_signalRegions)
    {
        theAnalysis->Make_PaperPlot_SignalRegions("NN_ctz");
        theAnalysis->Make_PaperPlot_SignalRegions("NN_ctw");
        theAnalysis->Make_PaperPlot_SignalRegions("NN_cpq3");
        theAnalysis->Make_PaperPlot_SignalRegions("NN_5D");
        theAnalysis->Make_PaperPlot_SignalRegions("NN_SM");
    }
    if(make_paperPlot_controlPlots) {theAnalysis->Make_PaperPlot_ControlPlots();}


    //#############################################
    //  OTHER FUNCTIONS
    //#############################################

    if(make_animation)
    {
        TString outdir = "./plots/animations/";
        TString cmd = "rm " + outdir + "*png";
        system(cmd.Data());

        vector<TString> v_eftpoints;
        for(float ipt=0; ipt<3; ipt+=0.10)
        {
            TString pt = "rwgt_ctz_" + Convert_Number_To_TString(ipt);
            v_eftpoints.push_back(pt);
        }
        // v_eftpoints.push_back("rwgt_ctz_0");
        // v_eftpoints.push_back("rwgt_ctz_1");
        // v_eftpoints.push_back("rwgt_ctz_2");
        // v_eftpoints.push_back("rwgt_ctz_3");

        for(int ipt=0; ipt<v_eftpoints.size(); ipt++)
        {
            theAnalysis->Make_Animation_PhysicsBriefing(v_eftpoints[ipt]);
        }

        float delay = 20; //delay in 1/100th of seconds //Not below 20
        cmd = "convert -delay "+Convert_Number_To_TString(delay)+" -loop 0 ";
        for(int ipt=0; ipt<v_eftpoints.size(); ipt++)
        {
            cmd+= "./plots/animations/" + v_eftpoints[ipt] + ".png ";
        }
        int nrepeats = 300 / delay; //Number of times we want to repeat the last image (to stay longer on last frame)
        for(int i=0; i<nrepeats; i++)
        {
            cmd+= "./plots/animations/" + v_eftpoints[v_eftpoints.size()-1] + ".png "; //Repeat last frame
        }

        cmd+= outdir + "test.gif";
        cout<<"cmd = "<<cmd<<endl;
        system(cmd.Data());
    }

    if(compare_template_shapes) {theAnalysis->Compare_TemplateShapes_Processes(template_name, plotChannel);}

    if(run_test) //Whatever function that is currently under testing
    {
        theAnalysis->Dump_Scores_allNNs();
    }

    //#############################################
    //  FINALIZE
    //#############################################

    delete theAnalysis; theAnalysis = NULL;
}
