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

    //-- M A I N    A N A L Y S I S    O P T I O N S --
    TString signal_process = "tZq"; //'tZq' or 'ttZ'
    bool use_systematics = false; //true <-> will compute/store systematics selected below
    TString region = "wz"; //Select a specific event category : '' (all preselected events) / 'tZq' / 'ttZ' / 'signal'
    bool is_blind = false; //true <-> don't read/store data events

    //-- M V A    S T R A T E G Y --
    TString classifier_name = "NN"; //'BDT' or 'NN'
    bool use_specificMVA_eachYear = true; //true <-> look for year-specific MVA weight files

    bool make_SMvsEFT_templates_plots = true; //true <-> templates & plots are produced for SM scenario only (separate SM processes); else, consider SM vs EFT scenario (and apply beforehand the chosen categorization strategy)
        int categorization_strategy = 2; //1 <-> define SRtZq/SRttZ with different jet multiplicities, apply dedicated binary classifiers; 2 <-> apply multi-classifier in merged SR; 0 <-> testing: read tmp MVA, no categ.
        bool plot_onlyMaxNodeEvents = true; //For multiclass NN-SM template plots only: true <-> only include events if they have their max output value in the corresponding node
        float cut_value_tZq = 0.5, cut_value_ttZ = 0.5; //Hard-coded cut values to apply -- for templates (automatic) and plots (user-option)
        bool plot_onlyMVACutEvents = false; //For binary MVA-SM templates plots only: true <-> only include events which pass the specified tZq or ttZ cut values
        bool keep_aboveCut = true; //true <-> only keep events satisfying x>=cut

    bool scanOperators_paramNN = false; //true <-> if considering a parametrized NN, multiple templates and plots will be created on a 1D or 2D grid of points (instead of a single point)
        TString operator1 = "ctz"; //First operator to scan (required)
        TString operator2 = ""; //Second operator to scan (optional)
        vector<float> v_WCs_operator_scan1 = {-5,-4,-3,-2,-1,0,1,2,3,4,5}; //Grid points for first operator (required)
        vector<float> v_WCs_operator_scan2 = {}; //Grid points for second operator (optional)

    //-- T E M P L A T E S --
    bool split_analysis_by_channel = false; //true <-> will *also* produce templates/histos/plots for each subchannel (defined below)
    TString template_name = ""; //'BDT', 'NN', 'categ' (nbjet/njet bins), 'Zpt', ...

    //-- P L O T T I N G --
    bool show_pulls_ratio = false; //true <-> bottom pad shows pull; else shows data/mc ratio (w/ errors)
    TString plot_extension = ".png"; //extension of plots


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
    // set_lumi_years.push_back("2016");
    set_lumi_years.push_back("2017");
    // set_lumi_years.push_back("2018");


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

    //DATA (Single sample, in first position)
    thesamplelist.push_back("DATA"); thesamplegroups.push_back("DATA");

    //Private MC production including EFT weights
    // thesamplelist.push_back("PrivMC_tZq"); thesamplegroups.push_back("PrivMC_tZq");
    // thesamplelist.push_back("PrivMC_ttZ"); thesamplegroups.push_back("PrivMC_ttZ");

    //Signals (central samples)
    thesamplelist.push_back("tZq"); thesamplegroups.push_back("tZq");
    thesamplelist.push_back("ttZ"); thesamplegroups.push_back("ttZ");

    //t(t)X
    thesamplelist.push_back("tWZ"); thesamplegroups.push_back("tX");
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

    //VV(V)
    thesamplelist.push_back("ZZ4l"); thesamplegroups.push_back("VVV");
    thesamplelist.push_back("ggToZZTo4l"); thesamplegroups.push_back("VVV");
    thesamplelist.push_back("ZZZ"); thesamplegroups.push_back("VVV");
    thesamplelist.push_back("WZZ"); thesamplegroups.push_back("VVV");
    thesamplelist.push_back("WWW"); thesamplegroups.push_back("VVV");
    thesamplelist.push_back("WWZ"); thesamplegroups.push_back("VVV");

    //WZ
    thesamplelist.push_back("WZ"); thesamplegroups.push_back("WZ");

    //X+g
    thesamplelist.push_back("TTGamma_Dilep"); thesamplegroups.push_back("Xg");
    thesamplelist.push_back("tGJets"); thesamplegroups.push_back("Xg");
    thesamplelist.push_back("WGToLNuG"); thesamplegroups.push_back("Xg");
    thesamplelist.push_back("ZGToLLG_01J"); thesamplegroups.push_back("Xg");

    //NPL
    thesamplelist.push_back("DY"); thesamplegroups.push_back("NPL");
    thesamplelist.push_back("TTbar_DiLep"); thesamplegroups.push_back("NPL");


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
    thevarlist.push_back("mTW");
    thevarlist.push_back("mHT");
    thevarlist.push_back("Mass_3l");
    thevarlist.push_back("maxEtaJet");
    thevarlist.push_back("jPrimeAbsEta");
    thevarlist.push_back("lAsymmetry");
    // thevarlist.push_back("maxDelPhiLL");
    // thevarlist.push_back("maxDeepJet");
    // thevarlist.push_back("leptonCharge");
    // thevarlist.push_back("cosThetaStarPolTop");
    // thevarlist.push_back("cosThetaStarPolZ");
    // thevarlist.push_back("recoZ_Pt");
    // thevarlist.push_back("recoZ_Eta");
    // thevarlist.push_back("recoZ_M");
    // thevarlist.push_back("recoTopLep_Pt");
    // thevarlist.push_back("recoTopLep_Eta");
    // thevarlist.push_back("recoTopLep_M");
    // thevarlist.push_back("TopZsystem_Pt");
    // thevarlist.push_back("TopZsystem_M");
    // thevarlist.push_back("jprime_Pt");


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
    // set_v_add_var_names.push_back("nMediumBJets");

    set_v_add_var_names.push_back("channel");
    set_v_add_var_names.push_back("njets");
    set_v_add_var_names.push_back("nbjets");
    set_v_add_var_names.push_back("metEt");


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
        // theSystTree.push_back("JESDown"); theSystTree.push_back("JESUp");

        //-- Implementend as event weights
        theSystWeights.push_back("PUDown"); theSystWeights.push_back("PUUp");
        theSystWeights.push_back("prefiringWeightDown"); theSystWeights.push_back("prefiringWeightUp");
        theSystWeights.push_back("BtagHFDown"); theSystWeights.push_back("BtagHFUp");
        theSystWeights.push_back("BtagLFDown"); theSystWeights.push_back("BtagLFUp");
        theSystWeights.push_back("HFstats1Down"); theSystWeights.push_back("HFstats1Up");
        theSystWeights.push_back("HFstats2Down"); theSystWeights.push_back("HFstats2Up");
        theSystWeights.push_back("LFstats1Down"); theSystWeights.push_back("LFstats1Up");
        theSystWeights.push_back("LFstats2Down"); theSystWeights.push_back("LFstats2Up");
        theSystWeights.push_back("CFerr1Down"); theSystWeights.push_back("CFerr1Up");
        theSystWeights.push_back("CFerr2Down"); theSystWeights.push_back("CFerr2Up");

        // theSystWeights.push_back("LepEff_muLooseDown"); theSystWeights.push_back("LepEff_muLooseUp");
        // theSystWeights.push_back("LepEff_muTightDown"); theSystWeights.push_back("LepEff_muTightUp");
        // theSystWeights.push_back("LepEff_elLooseDown"); theSystWeights.push_back("LepEff_elLooseUp");
        // theSystWeights.push_back("LepEff_elTightDown"); theSystWeights.push_back("LepEff_elTightUp");
        // theSystWeights.push_back("PDFDown"); theSystWeights.push_back("PDFUp");
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

    bool draw_input_vars = true; //Plot input variables
        bool draw_input_allChannels = false; //true <-> also draw for eachs split channel

    bool compare_template_shapes = false;

//-----------------    OTHER

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

//Apply choices given via command line, if any
	Apply_CommandArgs_Choices(argc, argv, set_lumi_years, region);

    // int nthreads = 4; ROOT::EnableImplicitMT(nthreads); //Enable multi-threading (I have 8 available threads)

    //#############################################
    //  CREATE INSTANCE OF CLASS & INITIALIZE
    //#############################################

    TopEFT_analysis* theAnalysis = new TopEFT_analysis(thesamplelist, thesamplegroups, theSystWeights, theSystTree, thechannellist, thevarlist, set_v_cut_name, set_v_cut_def, set_v_cut_IsUsedForBDT, set_v_add_var_names, plot_extension, set_lumi_years, show_pulls_ratio, region, signal_process, classifier_name, scanOperators_paramNN, operator1, operator2, v_WCs_operator_scan1, v_WCs_operator_scan2, make_SMvsEFT_templates_plots, is_blind, categorization_strategy, use_specificMVA_eachYear);
    if(theAnalysis->stop_program) {return 1;}

    //#############################################
    // TRAINING
    //#############################################

    if(train_BDT) {theAnalysis->Train_BDT("");}

    //#############################################
    //  TEMPLATES CREATION
    //#############################################

    if(create_templates) {theAnalysis->Produce_Templates(template_name, false, plot_onlyMaxNodeEvents, plot_onlyMVACutEvents, cut_value_tZq, cut_value_ttZ, keep_aboveCut);}

    //#############################################
    //  CONTROL HISTOGRAMS
    //#############################################

    if(create_inputVar_histograms) {theAnalysis->Produce_Templates(template_name, true, plot_onlyMaxNodeEvents, plot_onlyMVACutEvents, cut_value_tZq, cut_value_ttZ, keep_aboveCut);}

    //#############################################
    //  DRAW PLOTS
    //#############################################

    //All channels
    if(draw_templates)
    {
        theAnalysis->Draw_Templates(false, plotChannel, template_name, prefit, use_combine_file); //chosen channel

        if(plotChannel == "") //By default, also want to plot templates in subchannels
        {
            for(int ichan=1; ichan<thechannellist.size(); ichan++)
            {
                theAnalysis->Draw_Templates(false, thechannellist[ichan], template_name, prefit, use_combine_file);
            }
        }
    }

    if(draw_input_vars)
    {
        theAnalysis->Draw_Templates(true, plotChannel);
        if(draw_input_allChannels)
        {
            for(int ichan=1; ichan<thechannellist.size(); ichan++)
            {
                theAnalysis->Draw_Templates(true, thechannellist[ichan]);
            }
        }
    }

    //#############################################
    //  OTHER FUNCTIONS
    //#############################################

    if(compare_template_shapes) {theAnalysis->Compare_TemplateShapes_Processes(template_name, plotChannel);}

    //#############################################
    //  FINALIZE
    //#############################################

    delete theAnalysis;
}
