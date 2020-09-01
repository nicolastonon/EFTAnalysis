//by Nicolas Tonon (DESY)

//--- LIST OF FUNCTIONS (for quick search) :
//--------------------------------------------
// Train_BDT

// Produce_Templates

// Draw_Templates

// Compare_TemplateShapes_Processes

// MergeSplit_Templates
// SetBranchAddress_SystVariationArray
//--------------------------------------------

#include "TopEFT_analysis.h"

#define MYDEBUG(msg) cout<<endl<<ITAL("-- DEBUG: " << __FILE__ << ":" << __LINE__ <<":")<<FRED(" " << msg  <<"")<<endl

using namespace std;

//---------------------------------------------------------------------------
// ####    ##    ##    ####    ########
//  ##     ###   ##     ##        ##
//  ##     ####  ##     ##        ##
//  ##     ## ## ##     ##        ##
//  ##     ##  ####     ##        ##
//  ##     ##   ###     ##        ##
// ####    ##    ##    ####       ##
//---------------------------------------------------------------------------

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

//Overloaded constructor
TopEFT_analysis::TopEFT_analysis(vector<TString> thesamplelist, vector<TString> thesamplegroups, vector<TString> thesystlist, vector<TString> thesystTreelist, vector<TString> thechannellist, vector<TString> thevarlist, vector<TString> set_v_cut_name, vector<TString> set_v_cut_def, vector<bool> set_v_cut_IsUsedForBDT, vector<TString> set_v_add_var_names, TString theplotextension, vector<TString> set_lumi_years, bool show_pulls, TString region, TString signal_process, TString classifier_name, bool scanOperators_paramNN, TString operator_scan1, TString operator_scan2, vector<float> v_WCs_operator_scan1, vector<float> v_WCs_operator_scan2, bool make_SMvsEFT_templates_plots, bool is_blind, int categorization_strategy, bool use_specificMVA_eachYear)
{
    //Canvas definition
    Load_Canvas_Style();

    TH1::SetDefaultSumw2();
    gStyle->SetErrorX(0.);

    nbins = 10; //default

	mkdir("outputs", 0777);
	mkdir("plots", 0777);

	stop_program = false;

	this->region = region;

    if(region == "")
    {
        this->region = "signal";
        cout<<FBLU("NB: you did not select a category; setting it to 'is_signal_SR' by default..."<<"")<<endl;
        usleep(1000000); //Pause for 1s (in microsec)
    }

    this->signal_process = signal_process;
    if(signal_process != "tZq" && signal_process != "ttZ" && signal_process != "tWZ") {cout<<BOLD(FRED("ERROR ! [signal_process] option not recognized ! Abort..."))<<endl; stop_program = true;}

    if(region == "tZq" && signal_process != "tZq") {cout<<BOLD(FRED("ERROR ! [region] and [signal_process] options are not compatible ! Abort..."))<<endl; stop_program = true;}
    else if(region == "ttZ" && signal_process != "ttZ") {cout<<BOLD(FRED("ERROR ! [region] and [signal_process] options are not compatible ! Abort..."))<<endl; stop_program = true;}

    this->make_SMvsEFT_templates_plots = make_SMvsEFT_templates_plots;

    this->categorization_strategy = categorization_strategy;

    this->use_specificMVA_eachYear = use_specificMVA_eachYear;

	this->is_blind = is_blind;

    sample_list.resize(thesamplelist.size());
    sample_groups.resize(thesamplelist.size());
	for(int i=0; i<thesamplelist.size(); i++)
	{
        sample_list[i] = thesamplelist[i];
        sample_groups[i] = thesamplegroups[i];
	}

    //Count nof different sample groups
    nSampleGroups=1; //minimum is 1 sample group
    for(int i=1; i<sample_groups.size(); i++)
    {
        if(sample_groups[i] != sample_groups[i-1]) {nSampleGroups++;}
    }
    if(region=="tWZ") {nSampleGroups++;} //In tWZ region, single out tWZ process from tX

//=== LUMI ===//
/*
# Luminosity conventions #
- The 'v_lumiYears' vector lists the years which are to be considered (its elements can be '2016', '2017' and '2018' - in correct order)
- Depending on the elements of 'v_lumiYears', we define a unique identifier 'lumiName' (e.g. "201617" for 2016+2017, or "Run2" for the sum of all 3 years)
- From there, we can compute the corresponding integrated luminosity and store it in the 'lumiValue' variable (for plots)

- Input ntuples are stored in separe "2016"/"2017"/"2018" sub-directories ; hence, to read them, we loop on the elements of v_lumiYears. Also specify in histo name the exact year it corresponds to.
- The 'lumiName' identifier is instead used for outputs (root files, histograms, ...) to differentiate them and chain them to subsequent codes (Combine, etc.)
 */

    //Lumi choice determines which ntuples are read
    v_lumiYears = set_lumi_years; //Vector containing the names of the years to process

    //Set a unique name to each combination of year(s)
    if(v_lumiYears.size() == 1) {lumiName = v_lumiYears[0];}
    else if(v_lumiYears.size() == 2)
    {
        if(v_lumiYears[0] == "2016" && v_lumiYears[1] == "2017") {lumiName = "201617";}
        else if(v_lumiYears[0] == "2016" && v_lumiYears[1] == "2018") {lumiName = "201618";}
        else if(v_lumiYears[0] == "2017" && v_lumiYears[1] == "2018") {lumiName = "201718";}
    }
    else if(v_lumiYears.size() == 3) {lumiName = "Run2";}
    else {cout<<BOLD(FRED("ERROR ! I don't understand the list of years to process, please stick to the conventions ! Abort..."))<<endl; stop_program = true;}
    // this->lumiYear = lumiYear;

    Set_Luminosity(lumiName); //Compute the corresponding integrated luminosity

    dir_ntuples = "./input_ntuples/";
    dir_ntuples_SR = dir_ntuples + "SR/";
	// cout<<"dir_ntuples : "<<dir_ntuples<<endl;

	//-- Get colors
    int color_scheme = 0; //Check color scheme definitions directly in Get_Samples_Colors()
    // int color_scheme = 2; //Check color scheme definitions directly in Get_Samples_Colors()
	color_list.resize(sample_list.size());
	Get_Samples_Colors(color_list, v_custom_colors, sample_list, sample_groups, color_scheme); //Read hard-coded sample colors

	this->classifier_name = classifier_name;

    t_name = "result";

	plot_extension = theplotextension;

	show_pulls_ratio = show_pulls;

	syst_list.resize(thesystlist.size());
	for(int i=0; i<thesystlist.size(); i++)
	{
		syst_list[i] = thesystlist[i];
	}
	if(syst_list.size() == 0 || syst_list[0] != "") {cout<<"ERROR : first element of 'syst_list' is not empty (<-> nominal event weight) ! If that's what you want, remove this protection !"<<endl; stop_program = true;}

	systTree_list.resize(thesystTreelist.size());
	for(int i=0; i<thesystTreelist.size(); i++)
	{
		systTree_list[i] = thesystTreelist[i];
	}
	if(systTree_list.size() == 0 || systTree_list[0] != "") {cout<<"ERROR : first element of 'systTree_list' is not empty (<-> nominal TTree) ! If that's what you want, remove this protection !"<<endl; stop_program = true;}

	channel_list.resize(thechannellist.size());
	for(int i=0; i<thechannellist.size(); i++)
	{
		channel_list[i] = thechannellist[i];
	}
	if(channel_list.size() == 0 || channel_list[0] != "") {cout<<"ERROR : first element of 'channel_list' is not empty (<-> no subcat.) or vector is empty ! If that's what you want, remove this protection !"<<endl; stop_program = true;}

	for(int i=0; i<set_v_cut_name.size(); i++) //Region cuts vars (e.g. NJets)
	{
		v_cut_name.push_back(set_v_cut_name[i]);
		v_cut_def.push_back(set_v_cut_def[i]);
		v_cut_IsUsedForBDT.push_back(set_v_cut_IsUsedForBDT[i]);
		v_cut_float.push_back(-999);
		v_cut_char.push_back(1); //flags set to true by default

		//NOTE : it is a problem if a variable is present in more than 1 list, because it will cause SetBranchAddress conflicts (only the last SetBranchAddress to a branch will work)
		//---> If a variable is present in 2 lists, erase it from other lists !
		for(int ivar=0; ivar<thevarlist.size(); ivar++)
		{
			if(thevarlist[ivar].BeginsWith("is") || thevarlist[ivar].BeginsWith("passed") )
			{
				cout<<BOLD(FRED("## Warning : categories should not been used as input/spectator variables ! Are you sure ? "))<<endl;
			}
			if(thevarlist[ivar] == set_v_cut_name[i])
			{
				cout<<FGRN("** Constructor")<<" : erased variable "<<thevarlist[ivar]<<" from vector thevarlist (possible conflict) !"<<endl;
				thevarlist.erase(thevarlist.begin() + ivar);
				ivar--; //modify index accordingly
			}

		}
		for(int ivar=0; ivar<set_v_add_var_names.size(); ivar++)
		{
			if(set_v_add_var_names[ivar].BeginsWith("is") || set_v_add_var_names[ivar].BeginsWith("passed") )
			{
				cout<<BOLD(FRED("## Warning : categories should not been used as input/spectator variables ! Are you sure ? "))<<endl;
			}
			if(set_v_add_var_names[ivar] == set_v_cut_name[i])
			{
				cout<<FGRN("** Constructor")<<" : erased variable "<<set_v_add_var_names[ivar]<<" from vector set_v_add_var_names (possible conflict) !"<<endl;
				set_v_add_var_names.erase(set_v_add_var_names.begin() + ivar);
				ivar--; //modify index accordingly
			}
		}

		// cout<<"Cuts : name = "<<v_cut_name[i]<<" / def = "<<v_cut_def[i]<<endl;
	}
	for(int i=0; i<thevarlist.size(); i++) //TMVA vars
	{
		var_list.push_back(thevarlist[i]);
		var_list_floats.push_back(-999);

		for(int ivar=0; ivar<set_v_add_var_names.size(); ivar++)
		{
			if(set_v_add_var_names[ivar] == thevarlist[i])
			{
				cout<<FGRN("** Constructor")<<" : erased variable "<<set_v_add_var_names[ivar]<<" from vector set_v_add_var_names (possible conflict) !"<<endl;
				set_v_add_var_names.erase(set_v_add_var_names.begin() + ivar);
				ivar--; //modify index accordingly
			}
		}
	}

	for(int i=0; i<set_v_add_var_names.size(); i++) //Additional vars, only for CR plots
	{
		v_add_var_names.push_back(set_v_add_var_names[i]);
		v_add_var_floats.push_back(-999);
	}

	//Make sure that the "==" sign is written properly, or rewrite it
	for(int ivar=0; ivar<v_cut_name.size(); ivar++)
	{
		if( v_cut_def[ivar].Contains("=") && !v_cut_def[ivar].Contains("!") && !v_cut_def[ivar].Contains("==") && !v_cut_def[ivar].Contains("<") && !v_cut_def[ivar].Contains(">") )
		{
			v_cut_def[ivar] = "==" + Convert_Number_To_TString(Find_Number_In_TString(v_cut_def[ivar]));

			cout<<endl<<BOLD(FBLU("##################################"))<<endl;
			cout<<"--- Changed cut on "<<v_cut_name[ivar]<<" to: "<<v_cut_def[ivar]<<" ---"<<endl;
			cout<<BOLD(FBLU("##################################"))<<endl<<endl;
		}
	}

    array_PU = NULL;
    array_prefiringWeight = NULL;
    array_Btag = NULL;
    // array_LepEff_mu = NULL;
    // array_LepEff_el = NULL;
    // array_LepEffLoose_mu = NULL;
    // array_LepEffLoose_el = NULL;
    // array_LepEffTight_mu = NULL;
    // array_LepEffTight_el = NULL;

	//Store the "cut name" that will be written as a suffix in the name of each output file
	this->filename_suffix = "";
	TString tmp = "";
	for(int ivar=0; ivar<v_cut_name.size(); ivar++)
	{
		if(v_cut_name[ivar].BeginsWith("is") || v_cut_name[ivar].BeginsWith("passed") ) {continue;} //No need to appear in filename

		if(v_cut_def[ivar] != "")
		{
            if(!v_cut_def[ivar].Contains("&&") && !v_cut_def[ivar].Contains("||") ) //Single condition
            {
                tmp+= "_" + v_cut_name[ivar] + Convert_Sign_To_Word(v_cut_def[ivar]) + Convert_Number_To_TString(Find_Number_In_TString(v_cut_def[ivar]));
            }
            else if(v_cut_def[ivar].Contains("&&")) //Double '&&' condition
            {
                TString cut1 = Break_Cuts_In_Two(v_cut_def[ivar]).first, cut2 = Break_Cuts_In_Two(v_cut_def[ivar]).second;
                tmp = "_" + v_cut_name[ivar] + Convert_Sign_To_Word(cut1) + Convert_Number_To_TString(Find_Number_In_TString(cut1));
                tmp+= Convert_Sign_To_Word(cut2) + Convert_Number_To_TString(Find_Number_In_TString(cut2));
            }
            else if(v_cut_def[ivar].Contains("||") )
            {
                TString cut1 = Break_Cuts_In_Two(v_cut_def[ivar]).first, cut2 = Break_Cuts_In_Two(v_cut_def[ivar]).second;
                tmp = "_" + v_cut_name[ivar] + Convert_Sign_To_Word(cut1) + Convert_Number_To_TString(Find_Number_In_TString(cut1));
                tmp+= "OR" + Convert_Sign_To_Word(cut2) + Convert_Number_To_TString(Find_Number_In_TString(cut2));
            }

			this->filename_suffix+= tmp;
		}
	}

    if(classifier_name == "NN")
    {
        //-- Parametrized NN options
        this->scanOperators_paramNN = scanOperators_paramNN;
        this->operator_scan1 = operator_scan1;
        this->operator_scan2 = operator_scan2;
        this->v_WCs_operator_scan1 = v_WCs_operator_scan1;
        this->v_WCs_operator_scan2 = v_WCs_operator_scan2;
        this->idx1_operator_scan1=-1; this->idx1_operator_scan2=-1; //Find indices of the Wilson Coefficients of the scanned EFT operators (provided as input features)
        this->idx2_operator_scan1=-1; this->idx2_operator_scan2=-1;
        if(!this->v_WCs_operator_scan2.size() || operator_scan2 == "") {this->v_WCs_operator_scan2.resize(1); this->v_WCs_operator_scan2[0]=-999;} //This vector must not be empty for 1D/2D EFT scan --> insert dummy value

        if(this->scanOperators_paramNN)
        {
            if(this->operator_scan1 == "" || !this->v_WCs_operator_scan1.size()) {cout<<BOLD(FRED("Warning: option [v_WCs_operator_scan1] is true, but you did not specify [operator_scan1] or [v_WCs_operator_scan1] ---> No MVA scan ! Abort !"))<<endl; stop_program = true;}
            if(!this->make_SMvsEFT_templates_plots) {cout<<BOLD(FRED("ERROR: option [scanOperators_paramNN=true] can only be used with option [make_SMvsEFT_templates_plots=true] ! Abort !"))<<endl; stop_program = true;}
            if( (this->operator_scan1 != "" && !this->v_WCs_operator_scan1.size()) || (this->operator_scan2 != "" && !this->v_WCs_operator_scan2.size()) ) {cout<<BOLD(FRED("ERROR: option [operator_scan1] or [operator_scan2] is non-empty, but the corresponding list of WC values is empty ! Abort !"))<<endl; stop_program = true;}
        }
    }

    cout<<endl<<endl<<BLINK(BOLD(FBLU("[Region : "<<region<<"]")))<<endl;
    cout<<endl<<BLINK(BOLD(FBLU("[Luminosity : "<<lumiName<<"]")))<<endl<<endl<<endl;

//--------------------------------------------

    usleep(1000000); //Pause for 1s (in microsec)
}

TopEFT_analysis::~TopEFT_analysis()
{
    // cout<<"~TopEFT_analysis()"<<endl;

    for(int icol=0; icol<v_custom_colors.size(); icol++)
    {
        if(v_custom_colors[icol] != 0) {delete v_custom_colors[icol];}
    }
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/**
 * Compute the luminosity re-scaling factor (MC),  to be used thoughout the code
 * @param desired_luminosity [Value of the desired lumi in fb-1]
 */
void TopEFT_analysis::Set_Luminosity(TString luminame)
{
    if(luminame == "2016") {lumiValue = 35.92;}
	else if(luminame == "2017") {lumiValue = 41.53;}
    else if(luminame == "2018") {lumiValue = 59.74;}
    else if(luminame == "201617") {lumiValue = 35.92+41.53;}
    else if(luminame == "201618") {lumiValue = 35.92+59.74;}
    else if(luminame == "201718") {lumiValue = 41.53+59.74;}
    else if(luminame == "Run2") {lumiValue = 35.92+41.53+59.74;}

	assert(lumiValue > 0 && "Using wrong lumi reference -- FIX IT !"); //Make sure we use a sensible lumi value

    // if(luminosity_rescale != 1)
	// {
	// 	cout<<endl<<BOLD(FBLU("##################################"))<<endl;
	// 	cout<<"--- Using luminosity scale factor : "<<desired_luminosity<<" / "<<lumiValue<<" = "<<luminosity_rescale<<" ! ---"<<endl;
	// 	cout<<BOLD(FBLU("##################################"))<<endl<<endl;
	// }

    return;
}













//---------------------------------------------------------------------------
// ########    ########        ###       ####    ##    ##    ####    ##    ##     ######
//    ##       ##     ##      ## ##       ##     ###   ##     ##     ###   ##    ##    ##
//    ##       ##     ##     ##   ##      ##     ####  ##     ##     ####  ##    ##
//    ##       ########     ##     ##     ##     ## ## ##     ##     ## ## ##    ##   ####
//    ##       ##   ##      #########     ##     ##  ####     ##     ##  ####    ##    ##
//    ##       ##    ##     ##     ##     ##     ##   ###     ##     ##   ###    ##    ##
//    ##       ##     ##    ##     ##    ####    ##    ##    ####    ##    ##     ######
//---------------------------------------------------------------------------


/**
 * Perform BDT training. Uses relevant samples, enforces cuts, etc. Modify training parameters/nof events/... in this function.
 * NB : Will sum all years given to member vector 'v_lumiYears' for the training ! If want separate training per year, run this func individually for each one!
 */
void TopEFT_analysis::Train_BDT(TString channel)
{
//--- OPTIONS --------------------------------
//--------------------------------------------
	bool use_relative_weights = false; //false <-> use abs(weight), much faster if there are many negative weights
    int nmaxEv = 10000; //max nof events for train or test
    float trainingEv_proportion = 0.7;
//--------------------------------------------
//--------------------------------------------

    cout<<endl<<BYEL("                          ")<<endl<<endl;
    cout<<FYEL("--- TRAINING ---")<<endl;
    cout<<endl<<BYEL("                          ")<<endl<<endl;

	if(use_relative_weights) {cout<<"-- Using "<<BOLD(FGRN("*RELATIVE weights*"))<<" --"<<endl<<endl<<endl;}
	else {cout<<"-- Using "<<BOLD(FGRN("*ABSOLUTE weights*"))<<" --"<<endl<<endl<<endl;}

    if(classifier_name != "BDT") {cout<<BOLD(FRED("ERROR : can only train BDTs within TMVA for now... Abort ! (classifier_name = "<<classifier_name<<")"))<<endl; return;}
    if(signal_process == "") {cout<<BOLD(FRED("ERROR: can only train a BDT if [signal_process] is specified !"))<<endl; return;}

	mkdir("weightsMVA", 0777);
    mkdir("weightsMVA/BDT", 0777);
    mkdir(("weightsMVA/BDT/"+lumiName).Data(), 0777);
    mkdir(("weightsMVA/BDT/"+lumiName+"/"+signal_process).Data(), 0777);

	usleep(1000000); //Pause for 1s (in microsec)

//--------------------------------
//  ####  #    # #####  ####
// #    # #    #   #   #
// #      #    #   #    ####
// #      #    #   #        #
// #    # #    #   #   #    #
//  ####   ####    #    ####
//--------------------------------------------

	//---Apply additional cuts on the signal and background samples (can be different)
	TCut mycuts = "";
	TCut mycutb = "";
	TString tmp = "";

	//--- CHOOSE TRAINING EVENTS <--> cut on corresponding category
	TString cat_tmp = "";
	cat_tmp = Get_Category_Boolean_Name(region);

	//Even if ask templates in the SR, need to use training (looser) category for training !
	// if(cat_tmp.Contains("_SR") )
	// {
	// 	int i = cat_tmp.Index("_SR"); //Find index of substrin g
	// 	cat_tmp.Remove(i); //Remove substring
	// }
    // tmp+= cat_tmp + "==1";

    tmp = "!TMath::IsNaN(mTW)"; //TMP fix... because trained on 0jet events ?

	//--- Define additionnal cuts
	for(int ivar=0; ivar<v_cut_name.size(); ivar++)
	{
		if(v_cut_def[ivar] != "")
		{
			if(tmp != "") {tmp+= " && ";}

			if(!v_cut_def[ivar].Contains("&&") && !v_cut_def[ivar].Contains("||")) {tmp+= v_cut_name[ivar] + v_cut_def[ivar];} //If cut contains only 1 condition
			else if(v_cut_def[ivar].Contains("&&") && v_cut_def[ivar].Contains("||")) {cout<<BOLD(FRED("ERROR ! Wrong cut definition !"))<<endl;}
			else if(v_cut_def[ivar].Contains("&&") )//If '&&' in the cut, break it in 2
			{
				tmp+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).first;
				tmp+= " && ";
				tmp+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).second;
			}
			else if(v_cut_def[ivar].Contains("||") )//If '||' in the cut, break it in 2
			{
				tmp+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).first;
				tmp+= " || ";
				tmp+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).second;
			}
		}
	}

	bool split_by_leptonChan = false;
	if(split_by_leptonChan && (channel != "all" && channel != ""))
	{
		if(channel == "uuu" || channel == "uu")	{mycuts = "channel==0"; mycutb = "channel==0";}
		else if(channel == "uue" || channel == "ue") {mycuts = "channel==1"; mycutb = "channel==1";}
		else if(channel == "eeu" || channel == "ee") {mycuts = "channel==2"; mycutb = "channel==2";}
		else if(channel == "eee") {mycuts = "channel==3"; mycutb = "channel==3";}
		else {cout << "WARNING : wrong channel name while training " << endl;}
	}

	cout<<"-- Will apply the following cut(s) : "<<BOLD(FGRN(""<<tmp<<""))<<endl<<endl<<endl<<endl;
	usleep(2000000); //Pause for 2s (in microsec)

	if(tmp != "") {mycuts+= tmp; mycutb+= tmp;}

	//--------------------------------------------
	//---------------------------------------------------------------
    // This loads the TMVA libraries
    TMVA::Tools::Instance();

	//Allows to bypass a protection in TMVA::Transplot_extensionionHandler, cf. description in source file:
	// if there are too many input variables, the creation of correlations plots blows up memory and basically kills the TMVA execution --> avoid above critical number (which can be user defined)
	(TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 300;

	TString weights_dir = "weightsMVA";

    //The TMVA::DataLoader object will hold the training and test data, and is later passed to the TMVA::Factory
	TMVA::DataLoader* myDataLoader = new TMVA::DataLoader(weights_dir); //If no TString given in arg, path of weightdir *in TTree* will be : default/weights/...

	//--- Could modify here the name of local dir. storing the BDT weights (default = "weights")
	//By setting it to "", weight files will be stored directly at the path given to myDataLoader
	//Complete path for weight files is : [path_given_toDataloader]/[fWeightFileDir]
	//Apparently, TMVAGui can't handle nested repos in path given to myDataLoader... so split path in 2 here
	TMVA::gConfig().GetIONames().fWeightFileDir = "BDT/"+lumiName+"/"+signal_process;
	// if(classifier_name == "NN" && NN_type == "TMVA") {TMVA::gConfig().GetIONames().fWeightFileDir = "NN/TMVA/"+lumiName;}

//--------------------------------------------
 // #    #   ##   #####  #   ##   #####  #      ######  ####
 // #    #  #  #  #    # #  #  #  #    # #      #      #
 // #    # #    # #    # # #    # #####  #      #####   ####
 // #    # ###### #####  # ###### #    # #      #           #
 //  #  #  #    # #   #  # #    # #    # #      #      #    #
 //   ##   #    # #    # # #    # #####  ###### ######  ####
//--------------------------------------------

	// Define the input variables that shall be used for the MVA training
	for(int i=0; i<var_list.size(); i++)
	{
		myDataLoader->AddVariable(var_list[i].Data(), 'F');
	}

	//Choose if the cut variables are used in BDT or not
	//Spectator vars are not used for training/evalution, but possible to check their correlations, etc.
    //NB: if a cut requires "var == x", all the selected events will be equal to x, so can't use it as discriminant variable !
	for(int i=0; i<v_cut_name.size(); i++)
	{
		if(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("==")) {myDataLoader->AddVariable(v_cut_name[i].Data(), 'F');}
	}

	double nEvents_sig = 0;
	double nEvents_bkg = 0;


//--------------------------------------------
//                       #
//  ####  #  ####       #  #####  #    #  ####      ####    ##   #    # #####  #      ######  ####
// #      # #    #     #   #    # #   #  #    #    #       #  #  ##  ## #    # #      #      #
//  ####  # #         #    #####  ####   #          ####  #    # # ## # #    # #      #####   ####
//      # # #  ###   #     #    # #  #   #  ###         # ###### #    # #####  #      #           #
// #    # # #    #  #      #    # #   #  #    #    #    # #    # #    # #      #      #      #    #
//  ####  #  ####  #       #####  #    #  ####      ####  #    # #    # #      ###### ######  ####
//--------------------------------------------

	//--- Only select few samples for training
	std::vector<TFile*> files_to_close;

    for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
    {
    	for(int isample=0; isample<sample_list.size(); isample++)
        {
            TString samplename_tmp = sample_list[isample];

            //-- Protections
            if(sample_list[isample] == "DATA") {continue;} //don't use data for training

            //Can hardcode here the backgrounds against which to train, instead of considering full list of samples
            // if(signal_process == "tZq")
            // {
            //     if(!samplename_tmp.Contains("tZq") && !samplename_tmp.EndsWith("ttZ") && !samplename_tmp.Contains("ttH") && samplename_tmp.Contains("ttW") && samplename_tmp.Contains("WZ") && samplename_tmp.Contains("ZZ4l") && samplename_tmp.Contains("TTbar_DiLep") ) {continue;}
            // }

            //-- TMP sample list
            if(samplename_tmp != "tZq" && samplename_tmp != "ttZ" && !samplename_tmp.Contains("ttH") && samplename_tmp.Contains("ttW") && samplename_tmp.Contains("WZ") && samplename_tmp.Contains("ZZ4l") && samplename_tmp.Contains("TTbar_DiLep") && samplename_tmp.Contains("DY") ) {continue;}

    		cout<<endl<<"-- Sample : "<<sample_list[isample]<<endl;

            // --- Register the training and test trees
            TString inputfile = dir_ntuples + v_lumiYears[iyear] + "/" + sample_list[isample] + ".root";


//--------------------------------------------
// ##### ##### #####  ###### ######  ####
//   #     #   #    # #      #      #
//   #     #   #    # #####  #####   ####
//   #     #   #####  #      #           #
//   #     #   #   #  #      #      #    #
//   #     #   #    # ###### ######  ####
//--------------------------------------------

    	    TFile *file_input = 0, *file_input_train = 0, *file_input_test = 0;
    		TTree *tree = 0, *tree_train = 0, *tree_test = 0;

            file_input = TFile::Open(inputfile);
            if(!file_input) {cout<<BOLD(FRED(<<inputfile<<" not found!"))<<endl; continue;}
            files_to_close.push_back(file_input); //store pointer to file, to close it after its events will have been read for training

            tree = (TTree*) file_input->Get(t_name);
            if(tree==0) {cout<<BOLD(FRED("ERROR :"))<<" file "<<inputfile<<" --> *tree = 0 !"<<endl; continue;}
            else if(tree->GetEntries()==0) {cout<<"Empty ! "<<endl; continue;}
            else {cout<<FMAG("=== Opened file : ")<<inputfile<<endl<<endl;}

            // global event weights per tree (see below for setting event-wise weights)
            Double_t signalWeight     = 1.0;
            Double_t backgroundWeight = 1.0;

        //-- Choose between absolute/relative weights for training
    		if(samplename_tmp.Contains(signal_process) )
    		{
                nEvents_sig+= tree->GetEntries(mycuts); myDataLoader->AddSignalTree(tree, signalWeight);

    			if(use_relative_weights)
    			{
                    // TString weightExp = "weight";
                    TString weightExp = "eventWeight*eventMCFactor";
    				myDataLoader->SetSignalWeightExpression(weightExp);
    				cout<<"Signal sample : "<<samplename_tmp<<" / Weight expression : "<<weightExp<<endl<<endl;
    			}
    			else
    			{
                    // TString weightExp = "fabs(weight)";
                    TString weightExp = "fabs(eventWeight*eventMCFactor)";
    				myDataLoader->SetSignalWeightExpression(weightExp);
    				cout<<"Signal sample : "<<samplename_tmp<<" / Weight expression : "<<weightExp<<endl<<endl;
    			}
    		}
    		else
    		{
                nEvents_bkg+= tree->GetEntries(mycutb); myDataLoader->AddBackgroundTree(tree, backgroundWeight);

                if(use_relative_weights)
    			{
                    // TString weightExp = "weight";
                    TString weightExp = "eventWeight*eventMCFactor";
    				myDataLoader->SetBackgroundWeightExpression(weightExp);
    				cout<<"Background sample : "<<samplename_tmp<<" / Weight expression : "<<weightExp<<endl;
    			}
    			else
    			{
                    // TString weightExp = "fabs(weight)";
                    TString weightExp = "fabs(eventWeight*eventMCFactor)";
    				myDataLoader->SetBackgroundWeightExpression(weightExp);
    				cout<<"Background sample : "<<samplename_tmp<<" / Weight expression : "<<weightExp<<endl;
    			}
    		}

        } //samples loop
    } //year loop


//--------------------------------
// #####  #####  ###### #####    ##   #####  ######    ##### #####  ###### ######  ####
// #    # #    # #      #    #  #  #  #    # #           #   #    # #      #      #
// #    # #    # #####  #    # #    # #    # #####       #   #    # #####  #####   ####
// #####  #####  #      #####  ###### #####  #           #   #####  #      #           #
// #      #   #  #      #      #    # #   #  #           #   #   #  #      #      #    #
// #      #    # ###### #      #    # #    # ######      #   #    # ###### ######  ####
//--------------------------------------------

	if(mycuts != mycutb) {cout<<__LINE__<<FRED("PROBLEM : cuts are different for signal and background ! If this is normal, modify code -- Abort")<<endl; delete myDataLoader; return;}

	// If nTraining_Events=nTesting_Events="0", half of the events in the tree are used for training, and the other half for testing
	//NB : when converting nEvents to TString, make sure to ask for sufficient precision !

	//-- Choose dataset splitting
	TString nTraining_Events_sig = "", nTraining_Events_bkg = "", nTesting_Events_sig = "", nTesting_Events_bkg = "";

    int nTrainEvSig = (nEvents_sig * trainingEv_proportion < nmaxEv) ? nEvents_sig * trainingEv_proportion : nmaxEv;
    int nTrainEvBkg = (nEvents_bkg * trainingEv_proportion < nmaxEv) ? nEvents_bkg * trainingEv_proportion : nmaxEv;
    int nTestEvSig = (nEvents_sig * (1-trainingEv_proportion) < nmaxEv) ? nEvents_sig * (1-trainingEv_proportion) : nmaxEv;
    int nTestEvBkg = (nEvents_bkg * (1-trainingEv_proportion) < nmaxEv) ? nEvents_bkg * (1-trainingEv_proportion) : nmaxEv;

    nTraining_Events_sig = Convert_Number_To_TString(nTrainEvSig, 12);
    nTraining_Events_bkg = Convert_Number_To_TString(nTrainEvBkg, 12);
    nTesting_Events_sig = Convert_Number_To_TString(nTestEvSig, 12);
    nTesting_Events_bkg = Convert_Number_To_TString(nTestEvBkg, 12);

	cout<<endl<<endl<<endl<<BOLD(FBLU("==================================="))<<endl;
	cout<<BOLD(FBLU("-- Requesting "<<nTraining_Events_sig<<" Training events [SIGNAL]"))<<endl;
	cout<<BOLD(FBLU("-- Requesting "<<nTesting_Events_sig<<" Testing events [SIGNAL]"))<<endl<<endl;
	cout<<BOLD(FBLU("-- Requesting "<<nTraining_Events_bkg<<" Training events [BACKGROUND]"))<<endl;
	cout<<BOLD(FBLU("-- Requesting "<<nTesting_Events_bkg<<" Testing events [BACKGROUND]"))<<endl;
	cout<<BOLD(FBLU("==================================="))<<endl<<endl<<endl<<endl;

    myDataLoader->PrepareTrainingAndTestTree(mycuts, mycutb, "nTrain_Signal="+nTraining_Events_sig+":nTrain_Background="+nTraining_Events_bkg+":nTest_Signal="+nTesting_Events_sig+":nTest_Background="+nTesting_Events_bkg+":SplitMode=Random:!V");

	//-- for quick testing -- few events
	// myDataLoader->PrepareTrainingAndTestTree(mycuts, mycutb, "nTrain_Signal=10:nTrain_Background=10:nTest_Signal=10:nTest_Background=10:SplitMode=Random:NormMode=NumEvents:!V");

    //Output rootfile containing TMVAGui infos, ROCS, ... for control
    TString basename = classifier_name;
    if(channel != "") {basename+= "_" + channel;}
    if(region != "") {basename+= "_" + region;}
    basename+= "_" + lumiName + this->filename_suffix;
    TString output_file_name = "outputs/" + basename + ".root";

	TFile* output_file = TFile::Open(output_file_name, "RECREATE");

    //The TMVA::Factory handles the training/testing/evaluation phases
	TMVA::Factory* factory = new TMVA::Factory(classifier_name, output_file, "V:!Silent:Color:DrawProgressBar:Correlations=True:AnalysisType=Classification");

	//So that the output weights are labelled differently
	TString method_title = signal_process;
	if(channel != "") {method_title = "_" + channel;}


//--------------------------------------------
//  ####  #####  ##### #  ####  #    #  ####     #    # ###### ##### #    #  ####  #####
// #    # #    #   #   # #    # ##   # #         ##  ## #        #   #    # #    # #    #
// #    # #    #   #   # #    # # #  #  ####     # ## # #####    #   ###### #    # #    #
// #    # #####    #   # #    # #  # #      #    #    # #        #   #    # #    # #    #
// #    # #        #   # #    # #   ## #    #    #    # #        #   #    # #    # #    #
//  ####  #        #   #  ####  #    #  ####     #    # ######   #   #    #  ####  #####
//--------------------------------------------

    TString method_options = "";

    //~ttH2017
    // method_options= "!H:!V:NTrees=200:BoostType=Grad:Shrinkage=0.10:!UseBaggedBoost:nCuts=200:MinNodeSize=5%:NNodesMax=5:MaxDepth=8:NegWeightTreatment=PairNegWeightsGlobal:CreateMVAPdfs:DoBoostMonitor=True";

    //tHq2017
    method_options = "!H:!V:NTrees=200:nCuts=40:MaxDepth=4:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:NegWeightTreatment=PairNegWeightsGlobal:CreateMVAPdfs";

    //Testing
    // method_options = "!H:!V:NTrees=800:nCuts=200:MaxDepth=4:MinNodeSize=5%:UseBaggedBoost=True:BaggedSampleFraction=0.5:BoostType=Grad:Shrinkage=0.10:NegWeightTreatment=PairNegWeightsGlobal";

    //Quick test
    // method_options = "!H:!V:NTrees=20:nCuts=5:MaxDepth=1:MinNodeSize=10%:UseBaggedBoost=True:BaggedSampleFraction=0.5:BoostType=Grad:Shrinkage=0.10";


//--------------------------------------------
 // ##### #####    ##   # #    #       ##### ######  ####  #####       ###### #    #   ##   #
 //   #   #    #  #  #  # ##   #         #   #      #        #         #      #    #  #  #  #
 //   #   #    # #    # # # #  #         #   #####   ####    #         #####  #    # #    # #
 //   #   #####  ###### # #  # #         #   #           #   #         #      #    # ###### #
 //   #   #   #  #    # # #   ##         #   #      #    #   #         #       #  #  #    # #
 //   #   #    # #    # # #    #         #   ######  ####    #         ######   ##   #    # ######
//--------------------------------------------

	if(classifier_name == "BDT") {factory->BookMethod(myDataLoader, TMVA::Types::kBDT, method_title, method_options);} //Book BDT -- pass dataLoader and options as arg

	output_file->cd();

    //-- Save TMVA ranking info
    mkdir("outputs/Rankings", 0777);
    mkdir(("outputs/Rankings/"+lumiName).Data(), 0777); //Dir. containing variable ranking infos
	TString ranking_file_path = "outputs/Rankings/"+lumiName+"/rank_"+basename+".txt";
    cout<<endl<<endl<<endl<<FBLU("NB : Temporarily redirecting standard output to file '"<<ranking_file_path<<"' in order to save Ranking Info !!")<<endl<<endl<<endl;
	std::ofstream out("ranking_info_tmp.txt"); //Temporary name
    out<<endl;
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to text file --> Ranking info will be saved !

    // Train MVAs using the set of training events
    factory->TrainAllMethods();

	std::cout.rdbuf(coutbuf); //reset to standard output again

	//-- NB : Test & Evaluation recap in the output files
    factory->TestAllMethods(); // ---- Evaluate all MVAs using the set of test events
    factory->EvaluateAllMethods(); // ----- Evaluate and compare performance of all configured MVAs

	//Could retrieve ROC graph directly
	// TMultiGraph* rocgraph = f.GetROCCurveAsMultiGraph("<datasetname>");

    // --------------------------------------------------------------
    // Save the output
    output_file->Close();
    std::cout << "==> Wrote root file: " << output_file_name << std::endl;
    std::cout << "==> TMVA is done!" << std::endl;

    MoveFile("./ranking_info_tmp.txt", ranking_file_path);
    Extract_Ranking_Info(ranking_file_path, channel); //Extract only ranking info from TMVA output

	for(unsigned int i=0; i<files_to_close.size(); i++) {files_to_close[i]->Close(); delete files_to_close[i];}

	delete myDataLoader; myDataLoader = NULL;
	delete factory; factory = NULL;
	output_file->Close(); output_file = NULL;

	return;
}
































//---------------------------------------------------------------------------
//  ######  ########  ########    ###    ######## ########       ######## ######## ##     ## ########  ##          ###    ######## ########  ######
// ##    ## ##     ## ##         ## ##      ##    ##                ##    ##       ###   ### ##     ## ##         ## ##      ##    ##       ##    ##
// ##       ##     ## ##        ##   ##     ##    ##                ##    ##       #### #### ##     ## ##        ##   ##     ##    ##       ##
// ##       ########  ######   ##     ##    ##    ######            ##    ######   ## ### ## ########  ##       ##     ##    ##    ######    ######
// ##       ##   ##   ##       #########    ##    ##                ##    ##       ##     ## ##        ##       #########    ##    ##             ##
// ##    ## ##    ##  ##       ##     ##    ##    ##                ##    ##       ##     ## ##        ##       ##     ##    ##    ##       ##    ##
//  ######  ##     ## ######## ##     ##    ##    ########          ##    ######## ##     ## ##        ######## ##     ##    ##    ########  ######
//---------------------------------------------------------------------------


void TopEFT_analysis::Produce_Templates(TString template_name, bool makeHisto_inputVars, bool plot_onlyMaxNodeEvents, bool plot_onlyMVACutEvents, float cut_value_tZq, float cut_value_ttZ, bool keep_aboveCut, bool also_applyCut_onMaxNodeValue)
{
//--- OPTIONS --------------------------------
//--------------------------------------------
    bool noSysts_inputVars = true; //true <-> don't compute syst weights for histos of input variables (not worth the CPU-time)
    int nentries_max = -1; //-1 <-> use all entries; N <-> use at most N entries per sample (for testing)
//--------------------------------------------
//--------------------------------------------

    if(template_name == "" && classifier_name != "BDT" && classifier_name != "NN") {cout<<BOLD(FRED("Error : classifier_name value ["<<classifier_name<<"] not supported !"))<<endl; return;}
    if(template_name=="") {template_name = classifier_name;}

    cout<<endl<<BYEL("                          ")<<endl<<endl;
	if(makeHisto_inputVars) {cout<<FYEL("--- Producing Input variables histograms ---")<<endl;}
	else {cout<<FYEL("--- Producing ["<<template_name<<"] Templates ---")<<endl;}
    cout<<endl<<BYEL("                          ")<<endl<<endl;


//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

	TH1::SetDefaultSumw2();

    TString restore_classifier_name = classifier_name;
    if(makeHisto_inputVars) {classifier_name = "";} //For naming conventions

    //-- Define category of interest
    TString cat_tmp = region;

    //Hard-coded here: determine whether we are producing EFT templates using a predefined strategy (cf. definitions in main)
    bool use_predefined_EFT_strategy = false;
    bool use_maxNode_events = false;
    bool apply_MVASM_cut = false;
    if(!makeHisto_inputVars && categorization_strategy > 0)
    {
        use_predefined_EFT_strategy = true;

        //By default, apply minimal event preselection and set tZq sa signal process
        //Otherwise, will only consider the region/signal of interest
        if(cat_tmp == "")
        {
            cat_tmp = "signal";
            signal_process = "tZq";
        }

        if(categorization_strategy == 2 && (make_SMvsEFT_templates_plots || (!make_SMvsEFT_templates_plots && plot_onlyMaxNodeEvents)) ) {use_maxNode_events = true;} //Cases for which we want to cut on a multiclass MVA-SM
        else if(categorization_strategy == 1 && (make_SMvsEFT_templates_plots || (!make_SMvsEFT_templates_plots && plot_onlyMVACutEvents)) ) {apply_MVASM_cut = true;}
    }

    TString dir_ntuples_tmp = dir_ntuples;
    if((region=="tZq" || region=="ttZ" || region=="signal") || use_predefined_EFT_strategy) {dir_ntuples_tmp = dir_ntuples_SR; cout<<endl<<FMAG("NB: will use SR sub-samples ! Make sure they are up-to-date !")<<endl<<endl;}
    if(dir_ntuples_tmp == dir_ntuples_SR && !Check_File_Existence((dir_ntuples_SR + v_lumiYears[0] + "/DATA.root")))
    {
        cout<<endl<<FMAG("Warning : SR sub-sample (produced with code input_ntuples/Split_FullSamples) "<<(dir_ntuples_SR + v_lumiYears[0] + "/DATA.root")<<" not found ! Assuming that split ntuples are not available --> Will use full ntuples ! (slower)")<<endl<<endl;
        dir_ntuples_tmp = dir_ntuples;
    }

    //Can fill these vectors via call to Get_VectorAllEvents_passMVACut() to determine beforehand, for each sample, whether each event passes or not a given MVA cut
    vector<int> v_isEventPassMVACut_multiclass; //For MVA-SM-multiclass
    vector<int> v_isEventPassMVACut_tZq; //For MVA-SM-tZq
    vector<int> v_isEventPassMVACut_ttZ; //For MVA-SM-ttZ

    vector<TString> restore_syst_list = syst_list;
    vector<TString> restore_systTree_list = systTree_list;
    //Don't make systematics shifted his3tos for input vars (too long)
    //Force removal of systematics ; restore values at end of this func
    if(makeHisto_inputVars && noSysts_inputVars)
    {
        syst_list.resize(1);
        syst_list[0] = "";
        systTree_list.resize(1);
        systTree_list[0] = "";
    }
    else
    {
        //Must allocate memory to these arrays, which will hold the values of the systematics variations weights
        //NB : the sizes of the arrays are hardcoded, depends on Potato implementation (could also set same large size for all)
        array_PU = new double[2];
        array_prefiringWeight = new double[2];
        array_Btag = new double[16];
        // array_LepEff_mu = new double[4];
        // array_LepEff_el = new double[4];
        // array_LepEffLoose_mu = new double[2];
        // array_LepEffLoose_el = new double[2];
        // array_LepEffTight_mu = new double[2];
        // array_LepEffTight_el = new double[2];
    }

	//-- Create output file
    TString output_file_name = "";
    TString MVA_type = "";
	if(makeHisto_inputVars) {output_file_name = "outputs/ControlHistograms_" + region + "_" + lumiName + filename_suffix +".root";}
    else
    {
        if(use_predefined_EFT_strategy)
        {
            if(make_SMvsEFT_templates_plots)
            {
                MVA_type+= "EFT";
                if(categorization_strategy > 0) {MVA_type+= Convert_Number_To_TString(categorization_strategy);}
                if(this->scanOperators_paramNN) {MVA_type+= "param";}
            }
            else {MVA_type+= "SM";}
        }
        output_file_name = "outputs/Templates_" + template_name + "_" + MVA_type + "_" + region + "_" + lumiName + filename_suffix + ".root";
    }
	TFile* file_output = TFile::Open(output_file_name, "RECREATE");

    if(!makeHisto_inputVars)
    {
        if(classifier_name == "BDT") //BDT //NB: only consider 1 single list of variables (although training weight file may differ between years, read later)
        {
            //NB : TMVA requires floats, and nothing else, to ensure reproducibility of results (training done with floats) => Need to recast e.g. doubles as floats //See: https://sourceforge.net/p/tmva/mailman/message/836453/
            reader = new TMVA::Reader("!Color:!Silent");

            // Name & adress of local variables which carry the updated input values during the event loop
            // NB: the variable names MUST corresponds in name and type to those given in the weight file(s) used -- same order
            // NB: if booking 2 BDTs, must make sure that they use the same input variables... or else, find some way to make it work in the code)
            for(int i=0; i<var_list.size(); i++)
            {
                //cout<<"Added variable "<<var_list[i]<<endl;
                reader->AddVariable(var_list[i].Data(), &var_list_floats[i]);
            }
            for(int i=0; i<v_cut_name.size(); i++)
            {
                if(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("==")) {reader->AddVariable(v_cut_name[i].Data(), &v_cut_float[i]);}
            }
        }
    }

	//Input TFile and TTree, called for each sample
	TFile* file_input = NULL;
	TTree* tree = NULL;

    //Define ranges of jet/bjets multiplciities -- for 'categ' templates only
    int nbjets_min=1, nbjets_max=2, njets_min=2, njets_max=6;

	//Template binning
    nbins = 20;
	double xmin = -1, xmax = 1; //BDT: [-1,1]
    if(template_name == "NN") {xmin = 0;} //NN: [0,1]
    else if(template_name == "categ") {nbins = (nbjets_max-nbjets_min+1)*(njets_max-njets_min+1); xmin = 0; xmax = nbins;} //1 bin per sub-category
    else if(template_name == "Zpt") {nbins = 5; xmin = 0; xmax = 500;} //1 bin per sub-category

	//Want to plot ALL selected variables
	vector<TString> total_var_list;
	if(makeHisto_inputVars)
	{
		for(int i=0; i<var_list.size(); i++)
		{
			total_var_list.push_back(var_list[i]);
		}
		for(int i=0; i<v_add_var_names.size(); i++)
		{
			total_var_list.push_back(v_add_var_names[i]);
		}
	}
	else
	{
        if(template_name != "NN") //For NN, depends on content of info file, read directly within main loops below
        {
            total_var_list.push_back(template_name);
        }
	}
    vector<float> total_var_floats(total_var_list.size()); //NB : can not read/cut on BDT... (would conflict with input var floats ! Can not set address twice)


// #    #      ##      #    #    #
// ##  ##     #  #     #    ##   #
// # ## #    #    #    #    # #  #
// #    #    ######    #    #  # #
// #    #    #    #    #    #   ##
// #    #    #    #    #    #    #

// #       ####   ####  #####   ####
// #      #    # #    # #    # #
// #      #    # #    # #    #  ####
// #      #    # #    # #####       #
// #      #    # #    # #      #    #
// ######  ####   ####  #       ####

	// cout<<endl<<ITAL("-- Ntuples directory : "<<dir_ntuples_tmp<<"")<<endl<<endl;

    // float tmp_compare = 0;
    float total_nentries_toProcess = Count_Total_Nof_Entries(dir_ntuples_tmp, t_name, sample_list, systTree_list, v_cut_name, v_cut_def, v_lumiYears, makeHisto_inputVars, noSysts_inputVars);

    cout<<endl<<FBLU(OVERLINE("                           "))<<endl;
    cout<<FBLU(BOLD("Will process "<<std::setprecision(12)<<total_nentries_toProcess<<" entries..."))<<endl;
    cout<<FBLU(UNDL("                           "))<<endl<<endl<<endl;

    //-- Draw progress bar
    bool draw_progress_bar = true;
    if(total_nentries_toProcess < 200000) {draw_progress_bar = false;}
    Int_t ibar = 0; //event counter
    TMVA::Timer timer(total_nentries_toProcess, "", true);
    TMVA::gConfig().SetDrawProgressBar(1);
    TMVA::gConfig().SetUseColor(1);

    bool MVA_alreadyLoaded = false; //If reading same MVA for multiple years, need to load weights only once
    TString BDT_method_name = "BDT"; //Need to store BDT method name for evaluation

    //-- YEARS LOOP
    for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
    {
        cout<<endl<<UNDL(FMAG("=== YEAR : "<<v_lumiYears[iyear]<<""))<<endl<<endl;


// #    # #    #   ##
// ##  ## #    #  #  #
// # ## # #    # #    #
// #    # #    # ######
// #    #  #  #  #    #
// #    #   ##   #    #

        //-- BDT
        if(!makeHisto_inputVars && (template_name == "BDT" || template_name == "NN") && (use_specificMVA_eachYear || !MVA_alreadyLoaded)) //Read MVA input file
        {
            MVA_alreadyLoaded = true;
            TString MVA_input_path = Get_MVAFile_InputPath(template_name, signal_process, v_lumiYears[iyear], use_specificMVA_eachYear, make_SMvsEFT_templates_plots, false, categorization_strategy);
            if(MVA_input_path == "") {return;} //MVA file not found

            if(template_name == "BDT") //Book TMVA reader from .xml weight file
            {
                reader->BookMVA(BDT_method_name, MVA_input_path);
            }
            else if(template_name == "NN") //Read NN info and weights
            {
                //-- Get path of NN info text file --> Read list of input variables (and more)
                var_list_NN.resize(0); NN_iMaxNode = -1; NN_strategy = ""; NN_inputLayerName = ""; NN_outputLayerName = ""; NN_nNodes = -1;
                TString NNinfo_input_path = Get_MVAFile_InputPath(template_name, signal_process, v_lumiYears[iyear], use_specificMVA_eachYear, make_SMvsEFT_templates_plots, true, categorization_strategy);
                if(NNinfo_input_path == "") {return;} //MVA file not found

                if(Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN, v_NN_nodeLabels, NN_inputLayerName, NN_outputLayerName, NN_iMaxNode, NN_nNodes, &NN_strategy)) {clfy1 = new TFModel(MVA_input_path.Data(), var_list_NN.size(), NN_inputLayerName.Data(), NN_nNodes, NN_outputLayerName.Data());} //Load neural network model

                if(this->scanOperators_paramNN && NN_strategy != "MVA_param") {cout<<BOLD(FRED("ERROR ! You have set [scanOperators_paramNN=true] to scan over EFT operators, but you are reading a non-parametrized MVA ! Please fix this ! Aborting..."))<<endl; return;}
                else if(this->scanOperators_paramNN && template_name != "NN") {cout<<BOLD(FRED("ERROR ! You have set [scanOperators_paramNN=true] to scan over EFT operators, but you have selected 'template_name != NN' (MVA scan only makes sense for parametrized NN) ! Please fix this ! Aborting..."))<<endl; return;}

                var_list = var_list_NN; //Use NN input features
                var_list_floats.resize(var_list.size());
                total_var_list.clear();

                if(use_predefined_EFT_strategy)
                {
                    //-- Special case: if making [SM vs EFT] templates using strategy 1 or 2, need to access 2 MVA-EFTs (1 for tZq + 1 for ttZ)
                    if(cat_tmp != "tZq") //Book a second MVA (EFT-ttZ)
                    {
                        TString MVA_input_path = Get_MVAFile_InputPath(template_name, "ttZ", v_lumiYears[iyear], use_specificMVA_eachYear, true, false, categorization_strategy);
                        if(MVA_input_path == "") {return;} //MVA file not found

                        var_list_NN2.resize(0); NN2_iMaxNode = -1; NN2_strategy = ""; NN2_inputLayerName = ""; NN2_outputLayerName = ""; NN2_nNodes = -1;
                        TString NNinfo_input_path = Get_MVAFile_InputPath(template_name, "ttZ", v_lumiYears[iyear], use_specificMVA_eachYear, true, true, categorization_strategy);
                        if(NNinfo_input_path == "") {return;} //MVA file not found

                        //-- Read info file for second NN //NB: it is assumed that the NN strategy (e.g. parametrized NN) is the same as for the first NN !
                        if(Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN2, v_NN2_nodeLabels, NN2_inputLayerName, NN2_outputLayerName, NN2_iMaxNode, NN2_nNodes)) {clfy2 = new TFModel(MVA_input_path.Data(), var_list_NN2.size(), NN2_inputLayerName.Data(), NN2_nNodes, NN2_outputLayerName.Data());} //Load neural network model

                        var_list_floats_2.resize(var_list_NN2.size());

                        //-- Can set variables addresses only once --> Get correspondance between list of variables of MVA1 and MVA2
                        v_varIndices_inMVA1.resize(var_list_NN2.size());
                        for(int i=0; i<var_list_NN2.size(); i++)
                        {
                            v_varIndices_inMVA1[i] = -1;
                            for(int j=0; j<var_list_NN.size(); j++)
                            {
                                if(var_list_NN2[i] == var_list_NN[j]) {v_varIndices_inMVA1[i] = j;} //Store index of corresponding MVA1 feature
                            }
                        }
                    }
                }
                else
                {
                    if(NN_nNodes==1) {total_var_list.push_back(template_name);} //Binary
                    else //Multiclass --> Different label for each output node
                    {
                        for(int inode=0; inode<NN_nNodes; inode++) {total_var_list.push_back(template_name + Convert_Number_To_TString(inode));}
                    }
                }

                //-- Extend variable names (1 variable per EFT point)
                if(this->scanOperators_paramNN)
                {
                    //-- Find feature indices of input WCs
                    for(int ivar=0; ivar<var_list_NN.size(); ivar++)
                    {
                        if(var_list_NN[ivar]==operator_scan1) {idx1_operator_scan1 = ivar;}
                        if(var_list_NN[ivar]==operator_scan2) {idx1_operator_scan2 = ivar;}
                    }
                    if(clfy2)
                    {
                        for(int ivar=0; ivar<var_list_NN2.size(); ivar++)
                        {
                            if(var_list_NN2[ivar]==operator_scan1) {idx2_operator_scan1 = ivar;}
                            if(var_list_NN2[ivar]==operator_scan2) {idx2_operator_scan2 = ivar;}
                        }
                    }
                } //Param. MVA

            } //NN
        } //Read MVA

        //-- If using predefined SM vs SM categ. strategy, consider the variable in 3 hard-coded regions
        if(!makeHisto_inputVars && use_predefined_EFT_strategy)
        {
            vector<TString> total_var_list_tmp(total_var_list); //Tmp copy of variable list
            total_var_list.clear();
            for(int ivar=0; ivar<total_var_list_tmp.size(); ivar++)
            {
                if(cat_tmp == "signal" || cat_tmp == "tZq") total_var_list.push_back(total_var_list_tmp[ivar] + "_" + MVA_type + "_SRtZq");
                if(cat_tmp == "signal" || cat_tmp == "ttZ") total_var_list.push_back(total_var_list_tmp[ivar] + "_" + MVA_type + "_SRttZ");
                if(cat_tmp == "signal" && !this->scanOperators_paramNN) //Only include CR for non-parametrized MVA (simpler)
                {
                    total_var_list.push_back(total_var_list_tmp[ivar] + "_" + MVA_type + "_CR");
                }
            }

            //-- EFT scan: need to consider 1 set of variable *per EFT point*
            if(this->scanOperators_paramNN)
            {
                vector<TString> total_var_list_tmp(total_var_list); //Tmp copy of variable list
                total_var_list.clear(); //Reset the actual list
                for(int ivar=0; ivar<total_var_list_tmp.size(); ivar++) //For each variable initially found in the list, duplicate it once for each considered EFT point
                {
                    for(int iop1=0; iop1<v_WCs_operator_scan1.size(); iop1++)
                    {
                        TString opname1 = operator_scan1 + "_" + v_WCs_operator_scan1[iop1]; //Default: scan 1 operator
                        TString opname2 = ""; //Optional: scan 2nd operator
                        for(int iop2=0; iop2<v_WCs_operator_scan2.size(); iop2++)
                        {
                            if(operator_scan2=="" && iop2>0) {break;}
                            else if(operator_scan2!="") {opname2 = "_" + operator_scan2 + "_" + v_WCs_operator_scan2[iop2];}
                            total_var_list.push_back(total_var_list_tmp[ivar] + "_" + opname1 + opname2);
                            // cout<<"Added variable: "<<total_var_list[total_var_list.size()-1]<<endl;
                        }
                    }
                }
            } //MVA EFT scan
        } //Template + predefined strategy

        total_var_floats.resize(total_var_list.size());


//  ####    ##   #    # #####  #      ######    #       ####   ####  #####
// #       #  #  ##  ## #    # #      #         #      #    # #    # #    #
//  ####  #    # # ## # #    # #      #####     #      #    # #    # #    #
//      # ###### #    # #####  #      #         #      #    # #    # #####
// #    # #    # #    # #      #      #         #      #    # #    # #
//  ####  #    # #    # #      ###### ######    ######  ####   ####  #

    	//-- SAMPLE LOOP
    	for(int isample=0; isample<sample_list.size(); isample++)
    	{
    		// cout<<endl<<endl<<UNDL(FBLU("Sample : "<<sample_list[isample]<<""))<<endl;

    		//Open input TFile
    		TString inputfile = dir_ntuples_tmp + v_lumiYears[iyear] + "/" + sample_list[isample] + ".root";

            v_isEventPassMVACut_tZq.clear(); v_isEventPassMVACut_ttZ.clear(); v_isEventPassMVACut_multiclass.clear();
            if(use_predefined_EFT_strategy)
            {
                if(categorization_strategy == 1)
                {
                    if(apply_MVASM_cut)
                    {
                        if(cat_tmp == "signal" || cat_tmp == "tZq") {Get_VectorAllEvents_passMVACut(v_isEventPassMVACut_tZq, "tZq", classifier_name, "result", inputfile, v_lumiYears[iyear], cut_value_tZq, keep_aboveCut, use_specificMVA_eachYear, categorization_strategy, false, nentries_max, "is_tzq_SR", also_applyCut_onMaxNodeValue);}
                        if(cat_tmp == "signal" || cat_tmp == "ttZ") {Get_VectorAllEvents_passMVACut(v_isEventPassMVACut_ttZ, "ttZ", classifier_name, "result", inputfile, v_lumiYears[iyear], cut_value_ttZ, keep_aboveCut, use_specificMVA_eachYear, categorization_strategy, false, nentries_max, "is_ttz_SR", also_applyCut_onMaxNodeValue);}
                    }
                }
                else if(categorization_strategy == 2 && use_maxNode_events) {Get_VectorAllEvents_passMVACut(v_isEventPassMVACut_multiclass, "Multiclass", classifier_name, "result", inputfile, v_lumiYears[iyear], -1., keep_aboveCut, use_specificMVA_eachYear, categorization_strategy, false, nentries_max, "is_signal_SR", also_applyCut_onMaxNodeValue);}
            }

            //-- Obsolete ? Apply ad-hoc scale factor to private ttZ sample so that SM yield matches that of central ttZ sample
            // float SF_SMEFT_ttZ = 1.;
            // if(sample_list[isample].Contains("PrivMC_ttZ"))
            // {
            //     if(v_lumiYears[iyear] == "2016") {SF_SMEFT_ttZ = 1.;}
            //     else if(v_lumiYears[iyear] == "2017") {SF_SMEFT_ttZ = 0.316;}
            //     else if(v_lumiYears[iyear] == "2018") {SF_SMEFT_ttZ = 1.;}
            // }

    		// cout<<"inputfile "<<inputfile<<endl;
    		if(!Check_File_Existence(inputfile))
    		{
    			cout<<endl<<"File "<<inputfile<<FRED(" not found!")<<endl;
    			continue;
    		}

    		file_input = TFile::Open(inputfile, "READ");
            cout<<endl<<FBLU("Opened file "<<inputfile<<" ...")<<endl<<endl;

            bool isPrivMC = false;
            if(sample_list[isample].Contains("PrivMC") && !sample_list[isample].Contains("_c")) {isPrivMC = true;}

            vector<float> v_SWE; //Store Sums of Weights (SWE) for all reweight points -- for private MC samples only
            if(isPrivMC) //Not available yet
            {
                TString hSWE_name = "EFT_SumWeights";
                if(!file_input->GetListOfKeys()->Contains(hSWE_name)) {cout<<"ERROR ! Histogram "<<hSWE_name<<" containing the sums of weights not found... Abort !"<<endl; return;}

                //Read and store sums of weights (SWE)
                TH1F* h_SWE = (TH1F*) file_input->Get(hSWE_name);
                for(int ibin=0; ibin<100; ibin++)
                // for(int ibin=0; ibin<h_SWE->GetNbinsX(); ibin++)
                {
                    v_SWE.push_back(h_SWE->GetBinContent(ibin+1)); //1 SWE stored for each stored weight
                    // cout<<"v_SWE["<<ibin<<"] = "<<v_SWE[ibin]<<endl;
                }
                delete h_SWE;
            }

    		//-- Loop on TTrees : first empty element corresponds to nominal TTree ; additional TTrees may correspond to JES/JER TTrees (defined in main)
    		//NB : only nominal TTree contains systematic weights ; others only contain the nominal weight (but variables have different values)
    		for(int itree=0; itree<systTree_list.size(); itree++)
    		{
                if(sample_list[isample] == "DATA" && systTree_list[itree] != "") {continue;}

    			tree = 0;
                TString tmp = systTree_list[itree];
    			if(tmp == "") {tmp = t_name;}

                tree = (TTree*) file_input->Get(tmp);

    			if(!tree)
    			{
    				cout<<BOLD(FRED("ERROR : tree '"<<tmp<<"' not found in file : "<<inputfile<<" ! Skip !"))<<endl;
    				continue; //Skip sample
    			}


//   ##   #####  #####  #####  ######  ####   ####  ######  ####
//  #  #  #    # #    # #    # #      #      #      #      #
// #    # #    # #    # #    # #####   ####   ####  #####   ####
// ###### #    # #    # #####  #           #      # #           #
// #    # #    # #    # #   #  #      #    # #    # #      #    #
// #    # #####  #####  #    # ######  ####   ####  ######  ####

    			//Disactivate all un-necessary branches ; below, activate only needed ones
    			tree->SetBranchStatus("*", 0); //disable all branches by default, speed up
    			// tree->SetBranchStatus("xxx", 1);

                float njets, nbjets; //Needed for 'categ' templates
                float* mTW = 0; //Needed for CR EFT template
    			if(makeHisto_inputVars)
    			{
    				for(int i=0; i<total_var_list.size(); i++)
    				{
                        tree->SetBranchStatus(total_var_list[i], 1);
                        tree->SetBranchAddress(total_var_list[i], &total_var_floats[i]);
    				}
    			}
    			else //Book input variables in same order as for MVA training
    			{
                    if(template_name=="BDT" || template_name=="NN") //Activate input features needed for MVA evaluation (same as used for training)
                    {
                        int idx_mTW = -1; //Must keep track of this variable's index, and whether its address is already set
                        for(int i=0; i<var_list.size(); i++)
                        {
                            if(!isample) {cout<<DIM("MVA 1 -- "<<i<<" -- Activate variable '"<<var_list[i]<<"'")<<endl;}
                            if(var_list[i] == "ctz" || var_list[i] == "ctw" || var_list[i] == "cpq3" || var_list[i] == "cpqm" || var_list[i] == "cpt") {continue;} //WC input values are arbitrary, there is no address to set !
                            tree->SetBranchStatus(var_list[i], 1);
                            tree->SetBranchAddress(var_list[i], &var_list_floats[i]);
                            if(var_list[i] == "mTW") {idx_mTW = i;}
                        }

                        if(clfy2) //Also set branch addresses for second MVA
                        {
                            cout<<endl<<endl;
                            for(int i=0; i<var_list_NN2.size(); i++)
                            {
                                if(!isample) {cout<<DIM("MVA 2 "<<i<<" -- Activate variable '"<<var_list_NN2[i]<<"'")<<endl;}

                                if(var_list_NN2[i] == "ctz" || var_list_NN2[i] == "ctw" || var_list_NN2[i] == "cpq3" || var_list_NN2[i] == "cpqm" || var_list_NN2[i] == "cpt") {continue;} //WC input values are arbitrary, there is no address to set !

                                int index_sameVar_in_NN1_list = -1; //Check whether same input variable was already set in first MVA (can not set branch addresses twice !)
                                for(int j=0; j<var_list_NN.size(); j++)
                                {
                                    if(var_list_NN2[i] == var_list_NN[j]) {index_sameVar_in_NN1_list = j; break;} //Found same variable in var_list_NN
                                }

                                if(index_sameVar_in_NN1_list < 0) //Variable only found in 2nd MVA
                                {
                                    tree->SetBranchStatus(var_list_NN2[i], 1);
                                    tree->SetBranchAddress(var_list_NN2[i], &var_list_floats_2[i]);
                                }
                            }
                        } //2nd MVA

                        //-- mTW address
                        tree->SetBranchStatus("mTW", 1);
                        if(idx_mTW >= 0) {mTW = &var_list_floats[idx_mTW];} //Point to float already storing the value
                        else {tree->SetBranchAddress("mTW", mTW);} //Set address

                    } //MVA templates
                    else if(template_name=="categ") //Need to read jet/bjet multiplicities of each event
                    {
                        tree->SetBranchStatus("njets", 1);
                        tree->SetBranchAddress("njets", &njets);
                        tree->SetBranchStatus("nbjets", 1);
                        tree->SetBranchAddress("nbjets", &nbjets);
                    }
                    else if(template_name=="Zpt") //Need to read Zpt variable
                    {
                        tree->SetBranchStatus("recoZ_Pt", 1);
                        tree->SetBranchAddress("recoZ_Pt", &total_var_floats[0]);
                    }
    			}

    			for(int i=0; i<v_cut_name.size(); i++)
    			{
                    // cout<<"v_cut_name[i] "<<v_cut_name[i]<<endl;

    				tree->SetBranchStatus(v_cut_name[i], 1);
                    if(v_cut_name[i].BeginsWith("is") || v_cut_name[i].BeginsWith("passed") ) //Categories are encoded into Char_t, not float
    				{
    					tree->SetBranchAddress(v_cut_name[i], &v_cut_char[i]);
    				}
    				else //All others are floats
    				{
    					tree->SetBranchAddress(v_cut_name[i], &v_cut_float[i]);
    				}
    			}

    			//--- Cut on relevant event selection stored as flag (Char_t)
    			Char_t is_goodCategory = true; //Cut on event category flag
                if(cat_tmp != "")
                {
                    TString cat_name = Get_Category_Boolean_Name(cat_tmp);

                    bool alreadyCutOnFlag = false; //Check if address of category flag is already set and cut on (don't cut again)
                    for(int icut=0; icut<v_cut_name.size(); icut++)
                    {
                        if(cat_name == v_cut_name[icut]) {alreadyCutOnFlag = true;}
                    }

                    if(!alreadyCutOnFlag)
                    {
                        tree->SetBranchStatus(cat_name, 1);
                        tree->SetBranchAddress(cat_name, &is_goodCategory);
                        // cout<<"Will apply cut on flag ["<<cat_name<<"] !"<<endl;
                    }
                }

                //-- For strategy 1, need to cut on tZq / ttZ SR flags // Hard-coded
                Char_t is_tZqSR, is_ttZSR;
                if(use_predefined_EFT_strategy && categorization_strategy == 1)
                {
                    tree->SetBranchStatus(Get_Category_Boolean_Name("tZq"), 1); tree->SetBranchAddress(Get_Category_Boolean_Name("tZq"), &is_tZqSR);
                    tree->SetBranchStatus(Get_Category_Boolean_Name("ttZ"), 1); tree->SetBranchAddress(Get_Category_Boolean_Name("ttZ"), &is_ttZSR);
                }

    			//--- Cut on relevant categorization (lepton flavour, btagging, charge)
    			float channel;
    			tree->SetBranchStatus("channel", 1);
    			tree->SetBranchAddress("channel", &channel);

                //--- Event weights
                double eventWeight; //Corresponds to 'nominalMEWeight_' in TopAnalysis code
                float eventMCFactor; //Sample-specific SF (xsec*lumi/SWE)
                Float_t weightMENominal;
                tree->SetBranchStatus("eventWeight", 1);
    			tree->SetBranchAddress("eventWeight", &eventWeight);
                tree->SetBranchStatus("eventMCFactor", 1);
    			tree->SetBranchAddress("eventMCFactor", &eventMCFactor);
                tree->SetBranchStatus("weightMENominal", 1);
    			tree->SetBranchAddress("weightMENominal", &weightMENominal);

                //-- EFT reweights
                vector<float>* v_wgts = 0;
                vector<string>* v_ids = 0;
                if(isPrivMC)
                {
                    v_wgts = new vector<float>;
                    v_ids = new vector<string>;

                    tree->SetBranchStatus("mc_EFTweights", 1);
                    tree->SetBranchStatus("mc_EFTweightIDs", 1);
                    tree->SetBranchAddress("mc_EFTweights", &v_wgts);
                    tree->SetBranchAddress("mc_EFTweightIDs", &v_ids);
                }

                //Reserve 1 float for each systematic weight (also for nominal to keep ordering, but not used)
    			vector<Double_t*> v_double_systWeights(syst_list.size(), NULL);
    			for(int isyst=0; isyst<syst_list.size(); isyst++)
    			{
    				//-- Protections : not all syst weights apply to all samples, etc.
    				if(sample_list[isample] == "DATA") {break;}
    				if(systTree_list[itree] != "") {break;} //Syst event weights only stored in nominal TTree

                    //Set proper branch address (hard-coded mapping)
                    SetBranchAddress_SystVariationArray(tree, syst_list[isyst], v_double_systWeights, isyst);
    			}

    			//Reserve memory for 1 TH1F per category, per systematic, per variable //v3 <-> vec of vec of vec
    			vector<vector<vector<TH1F*>>> v3_histo_chan_syst_var(channel_list.size());

                //Idem for TH1EFT
                vector<vector<vector<TH1EFT*>>> v3_TH1EFT_chan_syst_var(channel_list.size());

    			for(int ichan=0; ichan<channel_list.size(); ichan++)
    			{
    				if((channel_list.size() > 1 && channel_list[ichan] == "") || sample_list[isample] == "DATA" || systTree_list[itree] != "") {v3_histo_chan_syst_var[ichan].resize(1);} //Cases for which we only need to store the nominal weight

                    //Reserve memory for TH1Fs
    				if(sample_list[isample] == "DATA" || systTree_list[itree] != "") //1 single weight
    				{
                        v3_histo_chan_syst_var[ichan].resize(1); //Cases for which we only need to store the nominal weight
                        v3_TH1EFT_chan_syst_var[ichan].resize(1);
    				}
    				else //Subcategories -> 1 histo for nominal + 1 histo per systematic
    				{
                        v3_histo_chan_syst_var[ichan].resize(syst_list.size());
                        v3_TH1EFT_chan_syst_var[ichan].resize(syst_list.size());
    				}

    				//Init TH1Fs
    				for(int isyst=0; isyst<v3_histo_chan_syst_var[ichan].size(); isyst++)
    				{
                        v3_histo_chan_syst_var[ichan][isyst].resize(total_var_list.size());
                        v3_TH1EFT_chan_syst_var[ichan][isyst].resize(total_var_list.size());

    					for(int ivar=0; ivar<total_var_list.size(); ivar++)
    					{
    						if(makeHisto_inputVars && !Get_Variable_Range(total_var_list[ivar], nbins, xmin, xmax)) {cout<<FRED("Unknown variable name : "<<total_var_list[ivar]<<"! (include it in function Get_Variable_Range() in Helper.cxx)")<<endl; continue;} //Get binning for this input variable

                            //-- Hard-coded template binnings
                            int nbins_tmp = nbins; double xmin_tmp = xmin, xmax_tmp = xmax;
                            if(total_var_list[ivar].Contains("mTW")) {nbins_tmp=20; xmin_tmp=0; xmax_tmp=200;}

                            v3_histo_chan_syst_var[ichan][isyst][ivar] = new TH1F("", "", nbins_tmp, xmin_tmp, xmax_tmp);

                            if(isPrivMC) {v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = new TH1EFT("", "", nbins_tmp, xmin_tmp, xmax_tmp);}
                            else {v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = NULL;}
    					}
    				} //syst
    			} //chan

                WCFit* eft_fit = NULL; //1 WCFit object per PrivMC event

                //For private EFT samples, get and store index of SM reweight
                int idx_sm = -1;
                if(isPrivMC)
                {
                    tree->GetEntry(0); //Read 1 entry

                    for(int iwgt=0; iwgt<v_ids->size(); iwgt++)
                    {
                        cout<<"ToLower(v_ids->at(iwgt) "<<ToLower(v_ids->at(iwgt))<<endl; //FIXME
                        if(ToLower(v_ids->at(iwgt)).Contains("_sm") ) {idx_sm = iwgt; break;} //SM weight found
                    }

                    if(idx_sm == -1) {cout<<BOLD(FRED("Error: SM reweight not found in private sample ! Abort ! "))<<endl; return;}
                }


// ###### #    # ###### #    # #####    #       ####   ####  #####
// #      #    # #      ##   #   #      #      #    # #    # #    #
// #####  #    # #####  # #  #   #      #      #    # #    # #    #
// #      #    # #      #  # #   #      #      #    # #    # #####
// #       #  #  #      #   ##   #      #      #    # #    # #
// ######   ##   ###### #    #   #      ######  ####   ####  #

    			// cout<<"* Tree '"<<systTree_list[itree]<<"' :"<<endl;

    			// int nentries = 1000;
    			int nentries = tree->GetEntries();
                if(nentries_max > 0 && nentries > nentries_max)
                {
                    nentries = nentries_max;
                    cout<<FRED("Warning: you have set the limit nentries_max = "<<nentries_max<<"")<<endl;
                }

                // if(!draw_progress_bar) {cout<<endl<< "--- "<<sample_list[isample]<<" : Processing: " << tree->GetEntries() << " events" << std::endl;}
                cout<<endl<< "--- "<<sample_list[isample]<<" : Processing " << nentries << " events" << std::endl;
                cout<<DIM("(path: "<<inputfile<<")")<<endl;

    			for(int ientry=0; ientry<nentries; ientry++)
    			{
    				// cout<<FGRN("ientry "<<ientry<<"")<<endl;

                    if(isPrivMC && ientry%5000==0) {cout<<DIM("Entry "<<ientry<<"")<<endl;} //Very slow, print progress

    				std::fill(var_list_floats.begin(), var_list_floats.end(), 0); //Reset vectors reading inputs to 0

    				tree->GetEntry(ientry);

    				if(isnan(eventWeight*eventMCFactor) || isinf(eventWeight*eventMCFactor))
    				{
    					cout<<BOLD(FRED("* Found event with eventWeight*eventMCFactor = "<<eventWeight*eventMCFactor<<" ; remove it..."))<<endl; continue;
    				}

//---- APPLY CUTS HERE -------------------------

                    //--- Cut on category value
                    if(cat_tmp != "" && !is_goodCategory) {continue;}

    				bool pass_all_cuts = true;
    				for(int icut=0; icut<v_cut_name.size(); icut++)
    				{
    					if(v_cut_def[icut] == "") {continue;}

    					//Categories are encoded into Char_t. Convert them to float for code automation
    					if(v_cut_name[icut].BeginsWith("is") || v_cut_name[icut].BeginsWith("passed") ) {v_cut_float[icut] = (float) v_cut_char[icut];}
    					// cout<<"Cut : name="<<v_cut_name[icut]<<" / def="<<v_cut_def[icut]<<" / value="<<v_cut_float[icut]<<" / pass ? "<<Is_Event_Passing_Cut(v_cut_def[icut], v_cut_float[icut])<<endl;
    					if(!Is_Event_Passing_Cut(v_cut_def[icut], v_cut_float[icut])) {pass_all_cuts = false; break;}
    				}
    				if(!pass_all_cuts) {continue;}

//--------------------------------------------

                    ibar++;
                    if(draw_progress_bar && ibar%50000==0) {timer.DrawProgressBar(ibar, ""); cout<<ibar<<" / "<<total_nentries_toProcess<<endl; }

                    // cout<<"eventWeight "<<eventWeight<<endl;
                    // cout<<"eventMCFactor "<<eventMCFactor<<endl;
                    // for(int ivar=0; ivar<var_list_floats.size(); ivar++) {cout<<"var_list_floats "<<ivar<<" = "<<var_list_floats[ivar]<<endl;} //Debug

    				//-- Compute MVA value
                    std::vector<float> clfy1_outputs, clfy2_outputs;
                    if(!makeHisto_inputVars)
                    {
                        if(clfy2) //Also read 2nd MVA
                        {
                            for(int i=0; i<var_list_NN2.size(); i++)
                            {
                                if(v_varIndices_inMVA1[i] >= 0) {var_list_floats_2[i] = var_list_floats[v_varIndices_inMVA1[i]];} //Read value from MVA1 float vector
                            }
                        }

                        if(template_name == "BDT") {total_var_floats[0] = reader->EvaluateMVA(BDT_method_name);} //BDT output value

                        else if(template_name == "NN" && !use_predefined_EFT_strategy) //NN output value //Default
                        {
                            clfy1_outputs = clfy1->evaluate(var_list_floats); //Evaluate output node(s) value(s)
                            float mva_tmp = -1;
                            NN_iMaxNode = -1;
                            for(int inode=0; inode<NN_nNodes; inode++)
                            {
                                total_var_floats[inode] = clfy1_outputs[inode];
                                if(clfy1_outputs[inode] > mva_tmp) {mva_tmp = clfy1_outputs[inode]; NN_iMaxNode = inode;}
                            }

                            if(clfy2) //Also read 2nd MVA
                            {
                                NN2_iMaxNode = -1;
                                clfy2_outputs = clfy2->evaluate(var_list_floats_2); //Evaluate output node(s) value(s)
                                float mva_tmp = -1;
                                for(int inode=0; inode<NN2_nNodes; inode++)
                                {
                                    if(clfy2_outputs[inode] > mva_tmp) {mva_tmp = clfy2_outputs[inode]; NN2_iMaxNode = inode;}
                                }
                            } //2nd NN
                        } //NN

                        else if(template_name == "categ") //Arbitrary binning depending on jets/bjets multiplicities
                        {
                            total_var_floats[0] = Get_x_jetCategory(njets, nbjets, nbjets_min, nbjets_max, njets_min, njets_max);
                            // cout<<"njets "<<njets<<" / nbjets "<<nbjets<<" --> categ "<<total_var_floats[0]<<endl;
                        }
                        //NB: no need to set anything more for Zpt template, its address is already set
                    }

                    double weight_tmp = eventWeight * eventMCFactor; //Fill histo with this weight ; manipulate differently depending on syst
                    float w_SMpoint = 0;
                    if(isPrivMC)
                    {
                        w_SMpoint = weight_tmp * v_wgts->at(idx_sm) / (weightMENominal * v_SWE[idx_sm]);

                        eft_fit = new WCFit("myfit"); //Init object, filled in function
                        Get_WCFit(eft_fit, v_ids, v_wgts, v_SWE, weight_tmp, weightMENominal, w_SMpoint, idx_sm); //Fill w/ nominal weight ?

                        //Obsolete? -- Apply ad-hoc scale factor to private ttZ sample so that SM yield matches that of central ttZ sample
                        // if(sample_list[isample] == "PrivMC_ttZ") {eft_fit->scale(SF_SMEFT_ttZ);}
                    }

    				//-- Can choose to fill event only into histogram corresponding to relevant leptonic channel
    				for(int ichan=0; ichan<channel_list.size(); ichan++)
    				{
                        //-- Fill histos for sub-channels with corresponding events
                        if(channel_list[ichan] == "uuu" && channel != 0) {continue;}
                        if(channel_list[ichan] == "uue" && channel != 1) {continue;}
                        if(channel_list[ichan] == "eeu" && channel != 2) {continue;}
                        if(channel_list[ichan] == "eee" && channel != 3) {continue;}

                        bool cat_alreadyFound = false; //For pre-defined strategies, fill orthogonal categories
                        for(int ivar=0; ivar<total_var_list.size(); ivar++)
                        {
                            if(use_predefined_EFT_strategy) //If event does not pass the required cut, don't fill the corresponding template
                            {
                                int idx_EFT_var = ivar; //Trick: by default, for predefined EFT strategies, there are 3 variables (hard-coded); but if we scan operators, there are 3 * N1 * N2 variables --> Determine which 'original variable' (tzq, ttz, cr) they correspond to
                                if(scanOperators_paramNN && this->scanOperators_paramNN)
                                {
                                    if(ivar < (v_WCs_operator_scan1.size() * v_WCs_operator_scan2.size())) {idx_EFT_var = 0;}
                                    else {idx_EFT_var = 1;}
                                }

                                if(categorization_strategy == 1) //Strategy 1
                                {
                                    if(cat_alreadyFound && !this->scanOperators_paramNN) {continue;}
                                    else if(idx_EFT_var == 0 && (!is_tZqSR || (apply_MVASM_cut && !v_isEventPassMVACut_tZq[ientry]))) {continue;} //Event did not pass the required MVA-SM-tZq cut
                                    else if(idx_EFT_var == 1 && (!is_ttZSR || (apply_MVASM_cut && !v_isEventPassMVACut_ttZ[ientry]))) {continue;} //Event did not pass the required MVA-SM-ttZ cut
                                    cat_alreadyFound = true;
                                }
                                else //Strategy 2
                                {
                                    // cout<<"ivar "<<ivar<<" / v_isEventPassMVACut_multiclass[ientry] "<<v_isEventPassMVACut_multiclass[ientry]<<endl;
                                    if(use_maxNode_events && v_isEventPassMVACut_multiclass[ientry] != idx_EFT_var) {continue;} //Current node is not max node
                                }

                                //-- If scanning EFT operators, need to set the values of input WCs accordingly
                                if(scanOperators_paramNN && this->scanOperators_paramNN)
                                {
                                    for(int iop1=0; iop1<v_WCs_operator_scan1.size(); iop1++)
                                    {
                                        for(int iop2=0; iop2<v_WCs_operator_scan2.size(); iop2++)
                                        {
                                            if(operator_scan2 == "" && iop2>0) {break;} //1D scan only
                                            if(ivar == (iop1*v_WCs_operator_scan2.size()+iop2+idx_EFT_var*v_WCs_operator_scan1.size()*v_WCs_operator_scan2.size()))
                                            {
                                                if(idx_EFT_var == 0) //tZq MVA
                                                {
                                                    var_list_floats[idx1_operator_scan1] = v_WCs_operator_scan1[iop1];
                                                    var_list_floats[idx1_operator_scan2] = v_WCs_operator_scan2[iop2];
                                                    // cout<<iop1<<" var "<<idx1_operator_scan1<<" ==> "<<v_WCs_operator_scan1[iop1]<<endl;
                                                    // cout<<iop2<<" var "<<idx1_operator_scan2<<" ==> "<<v_WCs_operator_scan2[iop2]<<endl;
                                                }
                                                else if(idx_EFT_var == 1) //ttZ MVA
                                                {
                                                    var_list_floats_2[idx2_operator_scan1] = v_WCs_operator_scan1[iop1];
                                                    var_list_floats_2[idx2_operator_scan2] = v_WCs_operator_scan2[iop2];
                                                }
                                            }
                                        }
                                    }

                                    //-- Re-evaluate output nodes values *at this particular EFT point*
                                    if(idx_EFT_var == 0)
                                    {
                                        std::vector<float> clfy1_outputs = clfy1->evaluate(var_list_floats);
                                        total_var_floats[ivar] = clfy1_outputs[0];
                                    }
                                    else if(idx_EFT_var == 1)
                                    {
                                        std::vector<float> clfy2_outputs = clfy2->evaluate(var_list_floats_2);
                                        total_var_floats[ivar] = clfy2_outputs[0];
                                    }
                                }
                                else if(template_name == "NN" && make_SMvsEFT_templates_plots) //SM vs EFT --> use EFT MVA or mTW
                                {
                                    if(idx_EFT_var==0) {std::vector<float> clfy1_outputs = clfy1->evaluate(var_list_floats); total_var_floats[ivar] = clfy1_outputs[0];} //MVA-EFT-tZq
                                    else if(idx_EFT_var==1) {std::vector<float> clfy2_outputs = clfy2->evaluate(var_list_floats_2); total_var_floats[ivar] = clfy2_outputs[0];} //MVA-EFT-ttZ
                                    else if(idx_EFT_var==2) {total_var_floats[ivar] = *mTW;}
                                    else {cout<<"ERROR: wrong index !"<<endl; continue;}
                                }
                                else if(template_name == "NN") {total_var_floats[ivar] = clfy1_outputs[idx_EFT_var];} //SM vs SM --> use multiclass MVA nodes
                                else if(template_name == "Zpt")
                                {
                                    total_var_floats[ivar] = total_var_floats[0]; //Branch address of Zpt is set to 'total_var_floats[0]' by default
                                }
                                else if(template_name == "categ") //Arbitrary binning depending on jets/bjets multiplicities
                                {
                                    total_var_floats[ivar] = Get_x_jetCategory(njets, nbjets, nbjets_min, nbjets_max, njets_min, njets_max);
                                    // cout<<"njets "<<njets<<" / nbjets "<<nbjets<<" --> categ "<<total_var_floats[0]<<endl;
                                }
                                // cout<<"ivar "<<ivar<<" / total_var_floats[ivar] "<<total_var_floats[ivar]<<endl;
                            }
                            else if(makeHisto_inputVars) //Some features have their address set twice (but only works once) --> need to read the proper variable
                            {
                                if(total_var_list[ivar] == "channel") {total_var_floats[ivar] = channel;}
                            }

        					for(int isyst=0; isyst<syst_list.size(); isyst++)
        					{
                                //-- Protections : not all syst weights apply to all samples, etc.
                                if(sample_list[isample] == "DATA" && syst_list[isyst] != "") {break;}
                                else if(systTree_list[itree] != "") {break;} //Syst event weights only stored in nominal TTree
                                else if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefiring") ) {continue;} //no prefire in 2018

        						// cout<<"-- sample "<<sample_list[isample]<<" / channel "<<channel_list[ichan]<<" / syst "<<syst_list[isyst]<<endl;

        						weight_tmp = eventWeight*eventMCFactor; //Nominal weight (no syst) //Re-initialize for each syst

                                //No prefiring uncertainty for 2018 samples ; still write the histo (=nominal) to avoid troubles down the way... ?
                                // if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefiring")) {*(v_double_systWeights[isyst]) = 1;}

                                // cout<<"nominal : "<<weight_tmp<<endl;
                                if(syst_list[isyst] != "") {weight_tmp*= *(v_double_systWeights[isyst]);} //Syst weights were already divided by nominal weight
                                // cout<<"syst : "<<weight_tmp<<endl;

        						if(isnan(weight_tmp) || isinf(weight_tmp))
        						{
        							cout<<BOLD(FRED("* Found event with syst. weight = "<<weight_tmp<<" ; remove it..."))<<endl;
        							cout<<"(sample "<<sample_list[isample]<<" / channel "<<channel_list[ichan]<<" / syst "<<syst_list[isyst]<<")"<<endl;
        							continue;
        						}

                                //-- Fill histogram
                                if(isPrivMC)
                                {
                                    Fill_TH1F_UnderOverflow(v3_histo_chan_syst_var[ichan][isyst][ivar], total_var_floats[ivar], w_SMpoint);
                                    Fill_TH1EFT_UnderOverflow(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar], total_var_floats[ivar], w_SMpoint, *eft_fit);
                                }
                                else {Fill_TH1F_UnderOverflow(v3_histo_chan_syst_var[ichan][isyst][ivar], total_var_floats[ivar], weight_tmp);}

                            } //syst loop
    					} //var loop

    					if(channel_list[ichan] != "") {break;} //subcategories are orthogonal ; if already found, can break subcat. loop
    				} //subcat/chan loop

                    if(eft_fit) {delete eft_fit; eft_fit = NULL;}
    			} //TTree entries loop
//--------------------------------------------


// #    # #####  # ##### ######
// #    # #    # #   #   #
// #    # #    # #   #   #####
// # ## # #####  #   #   #
// ##  ## #   #  #   #   #
// #    # #    # #   #   ######

    			// --- Write histograms
    			TString samplename = sample_list[isample];
    			if(sample_list[isample] == "DATA") {samplename = "data_obs";}

    			for(int ichan=0; ichan<channel_list.size(); ichan++)
    			{
    				// cout<<"channel "<<channel_list[ichan]<<endl;

    				for(int isyst=0; isyst<syst_list.size(); isyst++)
    				{
                        //-- Protections : not all syst weights apply to all samples, etc.
                        if(sample_list[isample] == "DATA" && syst_list[isyst] != "") {break;}
                        else if(systTree_list[itree] != "") {break;} //Syst event weights only stored in nominal TTree
                        else if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefiring") ) {continue;} //no prefire in 2018

    					// cout<<"isyst "<<isyst<<endl;

    					for(int ivar=0; ivar<total_var_list.size(); ivar++)
    					{
    						TString output_histo_name;

                            output_histo_name = total_var_list[ivar];
							if(channel_list[ichan] != "") {output_histo_name+= "_" + channel_list[ichan];}
                            output_histo_name+= "_" + v_lumiYears[iyear] + "__" + samplename;
							if(syst_list[isyst] != "") {output_histo_name+= "__" + Get_Modified_SystName(syst_list[isyst], v_lumiYears[iyear]);}
							else if(systTree_list[itree] != "") {output_histo_name+= "__" + systTree_list[itree];}

                            //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                            output_histo_name.ReplaceAll('-', 'm');

    						file_output->cd();

                            if(isPrivMC) //Private SMEFT samples, nominal -- Only used in Combine to extract yield parametrizations; also used in this code to plot SMEFT samples, etc.
                            {
                                if(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Integral() <= 0) {Set_Histogram_FlatZero((TH1F*&) v3_TH1EFT_chan_syst_var[ichan][isyst][ivar], output_histo_name, false);} //If integral of histo is negative, set to 0 (else COMBINE crashes) -- must mean that norm is close to 0 anyway

                                if(!makeHisto_inputVars && this->scanOperators_paramNN) //ParametrizedNN: rescale TH1EFT according to current WC value
                                {
                                    TString wcpt_name = "";
                                    for(int iop1=0; iop1<v_WCs_operator_scan1.size(); iop1++)
                                    {
                                        for(int iop2=0; iop2<v_WCs_operator_scan2.size(); iop2++)
                                        {
                                            if(operator_scan2=="" && iop2>0) {break;}
                                            if((iop1*v_WCs_operator_scan2.size()+iop2) == ivar)
                                            {
                                                wcpt_name = "rwgt_"+operator_scan1+"_"+v_WCs_operator_scan1[iop1];
                                                if(operator_scan2 != "") {wcpt_name+= "_"+operator_scan2+"_"+v_WCs_operator_scan2[iop2];}
                                                // cout<<"wcpt_name "<<wcpt_name<<endl;
                                            }
                                        }
                                    }
                                    WCPoint wcp = WCPoint((string) wcpt_name, 1.);
                                    v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Scale(wcp); //FIXME //-- actually, should keep SM norm ? (rescaling done in Combine)
                                }

                                //Obsolete? -- Apply ad-hoc scale factor to private ttZ sample so that SM yield matches that of central ttZ sample
                                // if(sample_list[isample] == "PrivMC_ttZ") {v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Scale(SF_SMEFT_ttZ);}

                                v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Write(output_histo_name);
                                // v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Write("TH1EFT_"+output_histo_name); //Store with specific prefix

                                //-- MVA-EFT: need to store each histogram bin separately so that they can be scaled independently in Combine
                                if(!makeHisto_inputVars && make_SMvsEFT_templates_plots) {StoreEachHistoBinIndividually(file_output, v3_TH1EFT_chan_syst_var[ichan][isyst][ivar], output_histo_name);}

                                /*
                                bool debug = false;
                                if(debug) //Printout WC values and extrapolated integrals for each benchmark point
                                {
                                    WCFit fit = v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->GetSumFit(); //Get summed fit (over all bins)

                                    cout<<"//-------------------------------------------- "<<endl;
                                    cout<<"chan "<<channel_list[ichan]<<" / syst "<<syst_list[isyst]<<" / var"<<var_list[ivar]<<endl;
                                    fit.dump(); //Printout all names and coefficients of WC pairs
                                    for (uint i=0; i < fit.points.size(); i++)
                                    {
                                        WCPoint wc_pt = fit.points.at(i);
                                        double fit_val = fit.evalPoint(&wc_pt);
                                        wc_pt.dump(); //Printout names and values of all WCs for this point
                                        std::cout << "===> " << std::setw(3) << i << ": " << std::setw(12) << std::setw(12) << fit_val << std::endl; //Printout : i / true weight / evaluated weight / diff
                                    }
                                    cout<<endl<<endl<<endl;
                                }
                                */
                            }
                            else
                            {
                                if(v3_histo_chan_syst_var[ichan][isyst][ivar]->Integral() <= 0) {Set_Histogram_FlatZero(v3_histo_chan_syst_var[ichan][isyst][ivar], output_histo_name, false);} //If integral of histo is negative, set to 0 (else COMBINE crashes) -- must mean that norm is close to 0 anyway

                                v3_histo_chan_syst_var[ichan][isyst][ivar]->Write(output_histo_name);
                                // cout<<"Wrote histo : "<<output_histo_name<<endl;

                                //-- MVA-EFT: need to store each histogram bin separately so that they can be scaled independently in Combine
                                if(!makeHisto_inputVars && make_SMvsEFT_templates_plots) {StoreEachHistoBinIndividually(file_output, v3_histo_chan_syst_var[ichan][isyst][ivar], output_histo_name);}
                            }

    						delete v3_histo_chan_syst_var[ichan][isyst][ivar]; v3_histo_chan_syst_var[ichan][isyst][ivar] = NULL;
                            if(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]) {delete v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]; v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = NULL;}
    					} //var loop
    				} //syst loop
    			} //chan loop

    			// cout<<"Done with "<<sample_list[isample]<<" sample"<<endl;

                if(isPrivMC) {delete v_wgts; delete v_ids;}
    			tree->ResetBranchAddresses(); //Detach tree from local variables (safe)
    			delete tree; tree = NULL;
    		} //end tree loop

    		file_input->Close(); file_input = NULL;
    	} //end sample loop

        if(!makeHisto_inputVars)
        {
            if(classifier_name == "BDT")  {delete reader; reader = NULL;}
            else
            {
                if(clfy1 && use_specificMVA_eachYear) {delete clfy1; clfy1 = NULL; MVA_alreadyLoaded = false;} //Only delete if not going to be reused for another iteration
                //-- NB: do not delete clfy2 //Causes crash (because both it uses the same TF session as clfy1, hence double-free... ?)
                // if(clfy2) {delete clfy2; clfy2 = NULL;}
            }
        }

    } //years loop


//  ####  #       ####   ####  ######
// #    # #      #    # #      #
// #      #      #    #  ####  #####
// #      #      #    #      # #
// #    # #      #    # #    # #
//  ####  ######  ####   ####  ######

	cout<<endl<<FYEL("==> Created root file: ")<<file_output->GetName()<<endl;
	cout<<FYEL("containing the "<<classifier_name<<" templates as histograms for : all samples / all channels")<<endl<<endl;

	file_output->Close(); file_output = NULL;

    //-- Can verify that the total nof processed entries computed from Count_Total_Nof_Entries() was effectively the nof processed entries
    // cout<<"total_nentries_toProcess --> "<<total_nentries_toProcess<<endl;
    // cout<<"tmp_compare --> "<<tmp_compare<<endl;

    //Restore potentially modified variables
    classifier_name = restore_classifier_name;
    syst_list = restore_syst_list;
    systTree_list = restore_systTree_list;

    if(!makeHisto_inputVars || !noSysts_inputVars) //free memory
    {
        if(array_PU) {delete array_PU; array_PU = NULL;}
        if(array_prefiringWeight) {delete array_prefiringWeight; array_prefiringWeight = NULL;}
        if(array_Btag) {delete array_Btag; array_Btag = NULL;}
        // if(array_LepEff_mu) {delete array_LepEff_mu; array_LepEff_mu = NULL; }
        // if(array_LepEff_el) {delete array_LepEff_el; array_LepEff_el = NULL;}
        // if(array_LepEffLoose_mu) {delete array_LepEffLoose_mu; array_LepEffLoose_mu = NULL;}
        // if(array_LepEffLoose_el) {delete array_LepEffLoose_el; array_LepEffLoose_el = NULL;}
        // if(array_LepEffTight_mu) {delete array_LepEffTight_mu; array_LepEffTight_mu = NULL;}
        // if(array_LepEffTight_el) {delete array_LepEffTight_el; array_LepEffTight_el = NULL;}
    }

    if(clfy1) {delete clfy1; clfy1 = NULL;}


// #    # ###### #####   ####  ######
// ##  ## #      #    # #    # #
// # ## # #####  #    # #      #####
// #    # #      #####  #  ### #
// #    # #      #   #  #    # #
// #    # ###### #    #  ####  ######

    if(!makeHisto_inputVars) //For COMBINE fit, want to directly merge contributions from different processes into single histograms
    {
        MergeSplit_Templates(output_file_name, total_var_list, template_name, region, true);
    }

	return;
}
































//-----------------------------------------------------------------------------------------
// ########  ########     ###    ##      ##
// ##     ## ##     ##   ## ##   ##  ##  ##
// ##     ## ##     ##  ##   ##  ##  ##  ##
// ##     ## ########  ##     ## ##  ##  ##
// ##     ## ##   ##   ######### ##  ##  ##
// ##     ## ##    ##  ##     ## ##  ##  ##
// ########  ##     ## ##     ##  ###  ###

// ######## ######## ##     ## ########  ##          ###    ######## ########  ######
//    ##    ##       ###   ### ##     ## ##         ## ##      ##    ##       ##    ##
//    ##    ##       #### #### ##     ## ##        ##   ##     ##    ##       ##
//    ##    ######   ## ### ## ########  ##       ##     ##    ##    ######    ######
//    ##    ##       ##     ## ##        ##       #########    ##    ##             ##
//    ##    ##       ##     ## ##        ##       ##     ##    ##    ##       ##    ##
//    ##    ######## ##     ## ##        ######## ##     ##    ##    ########  ######
//-----------------------------------------------------------------------------------------

//Create plots from histograms representing templates (used for signal extraction) or input variables (for validation).
//NB: sums histograms from all years included in 'v_lumiYears' ! If you want to plot for a single year, activate only that specific year !
/**
 * Create plots from histograms representing templates (used for signal extraction) or input variables (for validation)
 * NB: sums histograms from all years included in 'v_lumiYears' ! If you want to plot for a single year, activate only that specific year !
 * @param drawInputVars    [true <-> plot input features / selected kinematic variables; else plot templates]
 * @param channel          [If specified, make plots for a specific (lepton) sub-channel]
 * @param template_name    [Specify the name of the templates to plot]
 * @param prefit           [true <-> plot prefit templates; else postfit]
 * @param use_combine_file [true <-> look for Combine FitDiagnostics output rootfile (needed e.g. to make postfit plots); else use my own template file]
 * @param EFTpoint         [If specified, look for/plot histogram corresponding to this specific EFT point (--> plot different points in EFT scan)]
 * @param store_ymax_fixed [true <-> will store the plot's ymax value into ymax_fixed1 or ymax_fixed2, so that a common y-scale is used between different scan points]
 */
void TopEFT_analysis::Draw_Templates(bool drawInputVars, TString channel, TString template_name, bool prefit, bool use_combine_file, TString EFTpoint, bool store_ymax_fixed, float* ymax_fixed1, float* ymax_fixed2)
{
//--- OPTIONS --------------------------------
//--------------------------------------------
    bool draw_errors = true; //true <-> superimpose error bands on plot/ratio plot

	bool draw_logarithm = false; //true <-> plot y-axis in log scale

    bool superimpose_GENhisto = false; //true <-> superimpose corresponding GEN-level EFT histogram, for shape comparison... //Experimental

    bool superimpose_EFThist = true; //true <-> superimpose shape of EFT hists
        bool normalize_EFThist = false; //true <-> normalize EFT hists (arbitrary)
        vector<TString> v_EFT_samples;//Names of the private EFT samples to superimpose
        v_EFT_samples.push_back("PrivMC_tZq");
        v_EFT_samples.push_back("PrivMC_ttZ");
        vector<TString> v_EFT_points; //Names of the EFT points at which to reweight the histos //Must follow naming convention used for private generation
        v_EFT_points.push_back("rwgt_ctz_0");
        v_EFT_points.push_back("rwgt_ctz_3");
        // v_EFT_points.push_back("rwgt_ctz_0_ctw_0_cpqm_0_cpq3_0_cpt_15");

    // bool doNot_stack_signal = false; //Obsolete
//--------------------------------------------
//--------------------------------------------

    if(template_name == "" && classifier_name != "BDT" && classifier_name != "NN") {cout<<BOLD(FRED("Error : classifier_name value ["<<classifier_name<<"] not supported !"))<<endl; return;}
    if(template_name=="") {template_name = classifier_name;}

    //-- For parametrized NN, can produce plots for *each EFT point* considered in the scan (adapt histo name accordingly, etc.)
    if(EFTpoint != "")
    {
        v_EFT_points.resize(0);
        v_EFT_points.push_back("rwgt_"+EFTpoint);
        EFTpoint.ReplaceAll('-', 'm'); //Can't have minus sign in histogram name
    }

    cout<<endl<<BYEL("                          ")<<endl<<endl;
	if(drawInputVars) {cout<<FYEL("--- Producing Input Variables Plots / channel : "<<channel<<" ---")<<endl;}
    else {cout<<FYEL("--- Producing "<<template_name<<" Template Plots ---")<<endl;}
    cout<<endl<<BYEL("                          ")<<endl<<endl;

	if(drawInputVars)
	{
		classifier_name = ""; //For naming conventions
		use_combine_file = false;
		if(drawInputVars && !prefit) {cout<<"Error ! Can not draw postfit input vars yet !"<<endl; return;}
	}

//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

	//Can use 2 different files :
	//- the files containing the template histograms, produced with this code (-> only prefit plots)
	//- or, better, the file produced by Combine from the templates : contains the prefit distributions with total errors, and the postfit distribution
	//If want postfit plots, must use the Combine file. If want prefit plots, can use both of them (NB : errors will be different)

    // TString cat_tmp = (region=="") ? "SR" : region+"Cat";
    TString cat_tmp = region;

	if(!prefit) {use_combine_file = true;}
	if(drawInputVars && use_combine_file)
	{
		cout<<"-- Setting 'use_combine_file = false' !"<<endl;
		use_combine_file = false;
	}

	//-- Want to plot ALL selected variables
	vector<TString> total_var_list;
	if(drawInputVars)
	{
		for(int i=0; i<var_list.size(); i++)
		{
			total_var_list.push_back(var_list[i]);
		}
		for(int i=0; i<v_add_var_names.size(); i++)
		{
			total_var_list.push_back(v_add_var_names[i]);
		}
	}
	else
	{
        if(template_name == "NN" && !this->categorization_strategy)
        {
            //-- Read the relevant NN info file, just to know if wer are dealing with multiclass NN templates... !
            TString NNinfo_input_path = Get_MVAFile_InputPath(template_name, signal_process, v_lumiYears[0], use_specificMVA_eachYear, make_SMvsEFT_templates_plots, true, this->categorization_strategy);
            Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN, v_NN_nodeLabels, NN_inputLayerName, NN_outputLayerName, NN_iMaxNode, NN_nNodes);

            for(int inode=0; inode<NN_nNodes; inode++) {total_var_list.push_back("NN" + (NN_nNodes == 1? "" : Convert_Number_To_TString(inode)));} //1 template per output node
        }
        else {total_var_list.push_back(template_name);} //Single BDT template

        //-- Add suffix to variables (SM/EFT, category, ...)
        if(!drawInputVars && this->categorization_strategy > 0) //Hard-coded
        {
            TString MVA_type = "SM";
            if(make_SMvsEFT_templates_plots)
            {
                MVA_type = "EFT" + Convert_Number_To_TString(categorization_strategy);
                if(EFTpoint != "") {MVA_type+= "param";}
            }

            vector<TString> total_var_list_tmp(total_var_list); //Tmp copy of variable list
            total_var_list.clear();
            for(int ivar=0; ivar<total_var_list_tmp.size(); ivar++)
            {
                total_var_list.push_back(total_var_list_tmp[ivar] + "_" + MVA_type + "_SRtZq");
                total_var_list.push_back(total_var_list_tmp[ivar] + "_" + MVA_type + "_SRttZ");
                if(EFTpoint == "") //No CR for parametrized MVA scan
                {
                    total_var_list.push_back(total_var_list_tmp[ivar] + "_" + MVA_type + "_CR");
                }
            }
        } //Use categ. strategy
	} //Templates

    //-- Read input file (may be year-dependent)
    TString inputFile_path = Get_HistoFile_InputPath(!drawInputVars, template_name, signal_process, cat_tmp, lumiName, use_combine_file, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, (EFTpoint!=""));
    if(inputFile_path == "") {cat_tmp = signal_process; inputFile_path = Get_HistoFile_InputPath(!drawInputVars, template_name, cat_tmp, signal_process, lumiName, use_combine_file, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, (EFTpoint!=""));} //Retry with 'signal_process' as 'region' argument
    if(inputFile_path == "") {return;}
    TFile* file_input = TFile::Open(inputFile_path, "READ");


// #       ####   ####  #####   ####
// #      #    # #    # #    # #
// #      #    # #    # #    #  ####
// #      #    # #    # #####       #
// #      #    # #    # #      #    #
// ######  ####   ####  #       ####

	for(int ivar=0; ivar<total_var_list.size(); ivar++)
	{
		if(drawInputVars) {cout<<endl<<FBLU("* Variable : "<<total_var_list[ivar]<<" ")<<endl<<endl;}

		TH1F* h_tmp = NULL; //Tmp storing histo
		TH1F* h_tzq = NULL; //Store tZq shape
		TH1F* h_ttz = NULL; //Store ttZ shape
		TH1F* h_sum_data = NULL; //Will store data histogram
		vector<TH1F*> v_MC_histo; //Will store all MC histograms (1 TH1F* per MC sample)

		TGraphAsymmErrors* g_data = NULL; //If using Combine file, data are stored in TGAE
		TGraphAsymmErrors* g_tmp = NULL; //Tmp storing graph

		vector<TString> MC_samples_legend; //List the MC samples to mention in legend

		//-- Init error vectors
		double x, y, errory_low, errory_high;

		vector<double> v_eyl, v_eyh, v_exl, v_exh, v_x, v_y; //Contain the systematic errors (used to create the TGraphError)
		int nofbins=-1;

		vector<float> v_yield_sig, v_yield_bkg;

        TString data_histo_name = "";

        vector<vector<TH1EFT*>> v2_th1eft(v_EFT_points.size()); //Store TH1EFT objects (inner: samples, outer: EFT points)
        vector<vector<TString>> v2_th1eft_labels(v_EFT_points.size()); //Corresponding names for legend (inner: samples, outer: EFT points)

        //-- All histos are for given lumiYears and sub-channels --> Need to sum them all for plots
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
        	//Need to rescale signal to fitted signal strength manually, and add its error in quadrature in each bin (verify)
        	double sigStrength = 0;
        	double sigStrength_Error = 0;
        	if(!prefit)
        	{
        		TTree* t_postfit = (TTree*) file_input->Get("tree_fit_sb");
        		t_postfit->SetBranchAddress("r", &sigStrength);
        		t_postfit->SetBranchAddress("rErr", &sigStrength_Error);
        		t_postfit->GetEntry(0); //Only 1 entry = fit results
        		delete t_postfit; t_postfit = NULL;
        	}

    		for(int ichan=0; ichan<channel_list.size(); ichan++)
    		{
    			//If using my own template file, there is already a "summed channels" version of the histograms
    			if(channel_list[ichan] != channel)
    			{
    				if(use_combine_file) {if(channel != "") {continue;} } //In combine file, to get inclusive plot, must sum all subcategories
    				else {continue;}
    			}

    			//Combine file : histos stored in subdirs -- define dir name
    			TString dir_hist = "";
    			if(prefit) {dir_hist = "shapes_prefit/";}
    			else {dir_hist = "shapes_fit_s/";}
                dir_hist+= template_name;
                // dir_hist+= classifier_name + template_name;
                if(cat_tmp != "") {dir_hist+= "_" + cat_tmp;}
    			if(channel_list[ichan] != "") {dir_hist+= "_" + channel_list[ichan];} //for combine file
                dir_hist+= "_" + v_lumiYears[iyear] + "/";
    			if(use_combine_file && !file_input->GetDirectory(dir_hist)) {cout<<FRED("Directory "<<dir_hist<<" not found ! Skip !")<<endl; continue;}


// #    #  ####
// ##  ## #    #
// # ## # #
// #    # #
// #    # #    #
// #    #  ####

    			//--- Retrieve all MC samples
    			int nof_skipped_samples = 0; //Get sample index right

    			vector<bool> v_isSkippedSample(sample_list.size()); //Get sample index right (some samples are skipped)

    			for(int isample = 0; isample < sample_list.size(); isample++)
    			{
    				int index_MC_sample = isample - nof_skipped_samples; //Sample index, but not counting data/skipped sample

    				//In Combine, some individual contributions are merged as "Rares"/"EWK", etc.
    				//If using Combine file, change the names of the samples we look for, and look only once for histogram of each "group"
    				TString samplename = sample_list[isample];
    				if(use_combine_file)
    				{
    					if(isample > 0 && sample_groups[isample] == sample_groups[isample-1]) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //if same group as previous sample, skip it
    					else {samplename = sample_groups[isample];}
    				}

    				//Protections, special cases
    				if(sample_list[isample].Contains("DATA") ) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}
                    if(sample_list[isample].Contains("PrivMC"))  {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //Only central samples get stacked, not private samples

                    if(iyear==0 && !samplename.Contains("DATA")) //Add sample name to list (used for legend)
                    {
                        if(v_MC_histo.size() <=  index_MC_sample) {MC_samples_legend.push_back(samplename);}
                        // cout<<"v_MC_histo.size() "<<v_MC_histo.size()<<" / index_MC_sample "<<index_MC_sample<<" / MC_samples_legend.size() "<<MC_samples_legend.size()<<endl;
                    }
                    if(v_isSkippedSample[isample] == true) {continue;} //Skip this sample

    				// cout<<endl<<UNDL(FBLU("-- Sample : "<<sample_list[isample]<<" : "))<<endl;

    				h_tmp = 0;

    				TString histo_name = "";
    				if(use_combine_file) {histo_name = dir_hist + samplename;}
    				else
    				{
    					histo_name = total_var_list[ivar];
	                    if(channel != "") {histo_name+= "_" + channel;}
                        if(EFTpoint != "") {histo_name+= "_" + EFTpoint;}
    					histo_name+= + "_" + v_lumiYears[iyear];
                        histo_name+= + "__" + samplename;
    				}

    				if(use_combine_file && !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(samplename) ) {cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found ! Skip..."))<<endl; v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}
    				else if(!use_combine_file && !file_input->GetListOfKeys()->Contains(histo_name) ) {cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found ! Skip..."))<<endl; v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}

    				h_tmp = (TH1F*) file_input->Get(histo_name);
    				// cout<<"histo_name "<<histo_name<<endl;
    	            // cout<<"h_tmp->Integral() = "<<h_tmp->Integral()<<endl;

    				if(draw_errors)
    				{
    					//Initialize error vectors (only once at start)
    					if(nofbins == -1) //if not yet init, get histo parameters
    					{
    						nofbins = h_tmp->GetNbinsX();
    						for(int ibin=0; ibin<nofbins; ibin++)
    						{
    							v_eyl.push_back(0); v_eyh.push_back(0);
    							v_exl.push_back(h_tmp->GetXaxis()->GetBinWidth(ibin+1) / 2); v_exh.push_back(h_tmp->GetXaxis()->GetBinWidth(ibin+1) / 2);
    							v_x.push_back( (h_tmp->GetXaxis()->GetBinLowEdge(nofbins+1) - h_tmp->GetXaxis()->GetBinLowEdge(1) ) * ((ibin+1 - 0.5)/nofbins) + h_tmp->GetXaxis()->GetBinLowEdge(1));
    							v_y.push_back(0);
    						}
    					}

    					//Increment errors
    					for(int ibin=0; ibin<nofbins; ibin++) //Start at bin 1
    					{
    						// NOTE : for postfit, the bin error accounts for all systematics !
    						//If using Combine output file (from MLF), bin error contains total error. Else if using template file directly, just stat. error
    						v_eyl[ibin]+= pow(h_tmp->GetBinError(ibin+1), 2);
    						v_eyh[ibin]+= pow(h_tmp->GetBinError(ibin+1), 2);

    						v_y[ibin]+= h_tmp->GetBinContent(ibin+1); //This vector is used to know where to draw the error zone on plot (= on top of stack)

    						// if(ibin != 4) {continue;} //cout only 1 bin
    						// cout<<"x = "<<v_x[ibin]<<endl;    cout<<", y = "<<v_y[ibin]<<endl;    cout<<", eyl = "<<v_eyl[ibin]<<endl;    cout<<", eyh = "<<v_eyh[ibin]<<endl; //cout<<", exl = "<<v_exl[ibin]<<endl;    cout<<", exh = "<<v_exh[ibin]<<endl;
    					} //loop on bins

    					//-- Draw all errors
    					//--------------------------------------------
    					if(!use_combine_file) //In Combine file, already accounted in binError
    					{
    						for(int itree=0; itree<systTree_list.size(); itree++)
    						{
    							for(int isyst=0; isyst<syst_list.size(); isyst++)
    							{
    								//-- Protections : not all syst weights apply to all samples, etc.
    								if(syst_list[isyst] != "" && systTree_list[itree] != "") {break;} //JES,JER,... -> read first element only
                                    else if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefiring") ) {continue;} //no prefire in 2018

    								// cout<<"sample "<<sample_list[isample]<<" / channel "<<channel_list[ichan]<<" / syst "<<syst_list[isyst]<<endl;

    								TH1F* histo_syst = 0; //Store the "systematic histograms"

    								TString histo_name_syst = histo_name + "__";
                                    if(syst_list[isyst] != "") {histo_name_syst+= Get_Modified_SystName(syst_list[isyst], v_lumiYears[iyear]);}
                                    else {histo_name_syst+= systTree_list[itree];}

    								if(!file_input->GetListOfKeys()->Contains(histo_name_syst)) {continue;} //No error messages if systematics histos not found

    								histo_syst = (TH1F*) file_input->Get(histo_name_syst);

    								//Add up here the different errors (quadratically), for each bin separately
    								for(int ibin=0; ibin<nofbins; ibin++)
    								{
    									if(histo_syst->GetBinContent(ibin+1) == 0) {continue;} //Some syst may be null, don't compute diff

    									double tmp = 0;

    									//For each systematic, compute (shifted-nominal), check the sign, and add quadratically to the corresponding bin error
    									tmp = histo_syst->GetBinContent(ibin+1) - h_tmp->GetBinContent(ibin+1);

    									if(tmp>0) {v_eyh[ibin]+= pow(tmp,2);}
    									else if(tmp<0) {v_eyl[ibin]+= pow(tmp,2);}

    									if(ibin > 0) {continue;} //cout only first bin
    									// cout<<"//--------------------------------------------"<<endl;
    									// cout<<"Sample "<<sample_list[isample]<<" / Syst "<<syst_list[isyst]<< " / chan "<<channel_list[ichan]<<endl;
    									// cout<<"x = "<<v_x[ibin]<<endl;    cout<<", y = "<<v_y[ibin]<<endl;    cout<<", eyl = "<<v_eyl[ibin]<<endl;    cout<<", eyh = "<<v_eyh[ibin]<<endl; //cout<<", exl = "<<v_exl[ibin]<<endl;    cout<<", exh = "<<v_exh[ibin]<<endl;
    									// cout<<"(nominal value = "<<h_tmp->GetBinContent(ibin+1)<<" - shifted value = "<<histo_syst->GetBinContent(ibin+1)<<") = "<<h_tmp->GetBinContent(ibin+1)-histo_syst->GetBinContent(ibin+1)<<endl;
    								}

    								delete histo_syst;
    							} //end syst loop
    						} //systTree_list loop
    					} //use combine file?
    				} //draw errors?

                    //-- obsolete ?
    				// if(!samplename.Contains("DATA") )
    				// {
    				// 	if(v_MC_histo.size() <=  index_MC_sample) {MC_samples_legend.push_back(samplename);}
    				// }


 //  ####   ####  #       ####  #####   ####
 // #    # #    # #      #    # #    # #
 // #      #    # #      #    # #    #  ####
 // #      #    # #      #    # #####       #
 // #    # #    # #      #    # #   #  #    #
 //  ####   ####  ######  ####  #    #  ####

    				//Use color vector filled in main()
    				h_tmp->SetFillStyle(1001);
    				if(samplename == "Fakes") {h_tmp->SetFillStyle(3005);}
    		        else if(samplename == "QFlip" ) {h_tmp->SetFillStyle(3006);}
    		        // else if(samplename.Contains("TTbar") || samplename.Contains("TTJet") )
    				// {
    					// h_tmp->SetFillStyle(3005);
    				// 	if(samplename.Contains("Semi") ) {h_tmp->SetLineWidth(0);}
    				// }

    				h_tmp->SetFillColor(color_list[isample]);
    				h_tmp->SetLineColor(kBlack);
                    h_tmp->SetLineColor(color_list[isample]);
                    // cout<<"sample_list[isample] "<<sample_list[isample];
                    // cout<<" => color_list[isample] "<<color_list[isample]<<endl;

    				//Check color of previous *used* sample (up to 6, to account for potentially skipped samples)
    				for(int k=1; k<6; k++)
    				{
    					if(isample - k >= 0)
    					{
    						if(v_isSkippedSample[isample-k]) {continue;} //If previous sample was skipped, don't check its color
    						else if(color_list[isample] == color_list[isample-k]) {h_tmp->SetLineColor(color_list[isample]); break;} //If previous sample had same color, don't draw line
    					}
    					else {break;}
    				}

                    if(region == "tWZ" && samplename == "tWZ") {h_tmp->SetFillColor(kBlue);} //Different color

    				// v_MC_histo.push_back((TH1F*) h_tmp->Clone());
    				if(v_MC_histo.size() <=  index_MC_sample) {v_MC_histo.push_back((TH1F*) h_tmp->Clone());}
    				else if(!v_MC_histo[index_MC_sample]) {v_MC_histo[index_MC_sample] = (TH1F*) h_tmp->Clone();} //For FakeEle and FakeMu
    				else {v_MC_histo[index_MC_sample]->Add((TH1F*) h_tmp->Clone());}

                    v_MC_histo[index_MC_sample]->SetDirectory(0); //Dis-associate histo from TFile

    				if(channel_list[ichan] == "")
    				{
    					if(v_yield_sig.size() == 0) {v_yield_sig.resize(h_tmp->GetNbinsX()); v_yield_bkg.resize(h_tmp->GetNbinsX());}

    					for(int ibin=0; ibin<h_tmp->GetNbinsX(); ibin++)
    					{
    						if(sample_list[isample].Contains(signal_process) ) {v_yield_sig[ibin]+= h_tmp->GetBinContent(ibin+1);}
    						else {v_yield_bkg[ibin]+= h_tmp->GetBinContent(ibin+1);}
    					}
    				}

    				// cout<<"sample : "<<sample_list[isample]<<" / color = "<<color_list[isample]<<" fillstyle = "<<h_tmp->GetFillStyle()<<endl;
    				// cout<<"index_MC_sample "<<index_MC_sample<<endl;
    				// cout<<"v_MC_histo.size() "<<v_MC_histo.size()<<endl;
    				// cout<<"MC_samples_legend.size() "<<MC_samples_legend.size()<<endl<<endl;

    				delete h_tmp; h_tmp = NULL;
    			} //end sample loop
    		} //channel loop
        // } //years loop


// #####    ##   #####   ##
// #    #  #  #    #    #  #
// #    # #    #   #   #    #
// #    # ######   #   ######
// #    # #    #   #   #    #
// #####  #    #   #   #    #

    		if(use_combine_file) {data_histo_name = "data";}
    		else
    		{
                data_histo_name = total_var_list[ivar];
                if(channel != "") {data_histo_name+= "_" + channel;}
                if(EFTpoint != "") {data_histo_name+= "_" + EFTpoint;}
    			data_histo_name+= "_" + v_lumiYears[iyear];
                data_histo_name+= "__data_obs";
    		}

    		if(use_combine_file)
    		{
    			for(int ichan=0; ichan<channel_list.size(); ichan++) //With combine file, to get inclusive plot, must sum all sub-categories
    			{
    				if(channel != "" && channel_list[ichan] != channel) {continue;}

    				//Combine file : histos stored in subdirs -- define dir name
    				TString dir_hist = "";
    				if(prefit) {dir_hist = "shapes_prefit/";}
    				else {dir_hist = "shapes_fit_s/";}
                    if(cat_tmp != "") {dir_hist+= "_" + cat_tmp;}
    				if(channel_list[ichan] != "") {dir_hist+= "_" + channel_list[ichan] + "/";} //for combine file
                    dir_hist+= v_lumiYears[iyear] + "/";
    				if(!file_input->GetDirectory(dir_hist)) {cout<<ITAL("Directory "<<dir_hist<<" not found ! Skip !")<<endl; continue;}

    				if(!file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains("data")) {cout<<FRED(""<<dir_hist<<"data : not found ! Skip...")<<endl; continue;}

    				data_histo_name = dir_hist + "data";
    				// cout<<"data_histo_name "<<data_histo_name<<endl;
    				g_tmp = (TGraphAsymmErrors*) file_input->Get(data_histo_name); //stored as TGraph

    				//Remove X-axis error bars, not needed for plot
    				for(int ipt=0; ipt<g_tmp->GetN(); ipt++)
    				{
    					g_tmp->SetPointEXhigh(ipt, 0);
    					g_tmp->SetPointEXlow(ipt, 0);
    				}

    				if(!g_data) {g_data = (TGraphAsymmErrors*) g_tmp->Clone();}
    				else //Need to sum TGraphs content by hand //not anymore, 1 channel only !
    				{
    					double x_tmp,y_tmp,errory_low_tmp,errory_high_tmp;
    					for(int ipt=0; ipt<g_data->GetN(); ipt++)
    					{
    						g_data->GetPoint(ipt, x, y);
    						errory_low = g_data->GetErrorYlow(ipt);
    						errory_high = g_data->GetErrorYhigh(ipt);

    						g_tmp->GetPoint(ipt, x_tmp, y_tmp);
    						errory_low_tmp = g_tmp->GetErrorYlow(ipt);
    						errory_high_tmp = g_tmp->GetErrorYhigh(ipt);

    						double new_error_low = sqrt(errory_low*errory_low+errory_low_tmp*errory_low_tmp);
    						double new_error_high = sqrt(errory_high_tmp*errory_high_tmp+errory_high_tmp*errory_high_tmp);
    						g_data->SetPoint(ipt, x, y+y_tmp);
    						g_data->SetPointError(ipt,0,0, new_error_low, new_error_high); //ok to add errors in quadrature ?

    						// cout<<"ipt "<<ipt<<" / x1 "<<x<<" / y1 "<<y<<" / error1 "<<errory_low<<", "<<errory_high<<endl;
    						// cout<<"ipt "<<ipt<<" / x2 "<<x_tmp<<" / y2 "<<y_tmp<<" / error2 "<<errory_low_tmp<<", "<<errory_high_tmp<<endl;
    						// cout<<"=> y1+y2 = "<<y+y_tmp<<" / error = "<<new_error_low<<", "<<new_error_high<<endl;
    					}
    				}
    			} //chan loop
    		}
    		else //If using template file made from this code
    		{
    			if(!file_input->GetListOfKeys()->Contains(data_histo_name)) {cout<<data_histo_name<<" : not found"<<endl;}
    			else
    			{
    				h_tmp = (TH1F*) file_input->Get(data_histo_name);
    				if(h_sum_data == NULL) {h_sum_data = (TH1F*) h_tmp->Clone();}
    				else {h_sum_data->Add((TH1F*) h_tmp->Clone());} //not needed anymore (1 channel only)
                    h_sum_data->SetDirectory(0); //Dis-associate from TFile
    				delete h_tmp; h_tmp = NULL;
    			}
    		}


// #####  #####  # #    #         ####    ##   #    # #####  #      ######
// #    # #    # # #    #        #       #  #  ##  ## #    # #      #
// #    # #    # # #    #         ####  #    # # ## # #    # #      #####
// #####  #####  # #    # ###         # ###### #    # #####  #      #
// #      #   #  #  #  #  ###    #    # #    # #    # #      #      #
// #      #    # #   ##   ###     ####  #    # #    # #      ###### ######

            if(superimpose_EFThist)
            {
                int icolor = 2; //Use different colors for each histo
                for(int ipoint=0; ipoint<v_EFT_points.size(); ipoint++)
                {
                    for(int isample=0; isample<v_EFT_samples.size(); isample++)
                    {
                        // cout<<"isample "<<isample<<" / ipoint "<<ipoint<<endl;

                        TH1EFT* th1eft_tmp = 0;

                        TString histo_name = total_var_list[ivar];
	                    if(channel != "") {histo_name+= "_" + channel;}
                        if(EFTpoint != "") {histo_name+= "_" + EFTpoint;}
    					histo_name+= "_" + v_lumiYears[iyear] + "__" + v_EFT_samples[isample];
                        // cout<<"histo_name "<<histo_name<<endl;

                        if(!file_input->GetListOfKeys()->Contains(histo_name) ) {cout<<ITAL("TH1EFT object '"<<histo_name<<"' : not found ! Skip...")<<endl; continue;}
                        th1eft_tmp = (TH1EFT*) file_input->Get(histo_name);
                        // cout<<"th1eft_tmp->Integral() "<<th1eft_tmp->Integral()<<endl;

                        //Rescale TH1EFT accordingly to current reweight //Pay attention to operator exact names !
                        WCPoint wcp = WCPoint((string) v_EFT_points[ipoint], 1.);
                        th1eft_tmp->Scale(wcp);
                        // cout<<"th1eft_tmp->Integral() "<<th1eft_tmp->Integral()<<endl;

                        // th1eft_tmp->SetLineColor(kRed);
                        th1eft_tmp->SetLineColor(icolor);
                        icolor++;
                        th1eft_tmp->SetLineWidth(5);
                        th1eft_tmp->SetLineStyle(2);

                        if(iyear == 0 || v2_th1eft[ipoint].size() <= isample || v2_th1eft[ipoint][isample] == 0) //New point or sample --> Add new elements
                        {
                            v2_th1eft[ipoint].push_back(th1eft_tmp);

                            //-- Get/save corresponding label
                            vector<pair<TString,float>> v = Parse_EFTreweight_ID(v_EFT_points[ipoint]);
                            TString EFTpointlabel = "";
                            for(int i=0; i<v.size(); i++)
                            {
                                if(v[i].second != 0)
                                {
                                    if(EFTpointlabel != "") {EFTpointlabel+= ",";}
                                    EFTpointlabel+= v[i].first + "=" + v[i].second;
                                }
                            }
                            if(EFTpointlabel == "") {EFTpointlabel = "SM";}
                            // cout<<"EFTpointlabel "<<EFTpointlabel<<endl;

                            std::vector<std::string> words;
                            split_string((string) v_EFT_samples[isample], words, "_");
                            TString leg_name = words.at(1) + " ("+EFTpointlabel+")";
                            // TString leg_name = "#splitline{"+words.at(1)+"}{("+EFTpointlabel+")}"; //Process name + EFT point //Split over 2 lines
                            v2_th1eft_labels[ipoint].push_back(leg_name);
                        }
                        else {v2_th1eft[ipoint][isample]->Add(th1eft_tmp);} //Sum different years together

                        v2_th1eft[ipoint][isample]->SetDirectory(0); //Dis-associate from TFile
                    } //sample loop
                } //EFT points loop
            } //Superimpose private SMEFT histograms
        } //years loop

//--------------------------------------------

		bool data_notEmpty = true;
		if(use_combine_file && !g_data) {cout<<endl<<BOLD(FRED("--- Empty data TGraph !"))<<endl<<endl; data_notEmpty = false;}
		if(!use_combine_file && !h_sum_data) {cout<<endl<<BOLD(FRED("--- Empty data histogram "<<data_histo_name<<" !"))<<endl<<endl; data_notEmpty = false;}

		//Make sure there are no negative bins
		if(data_notEmpty)
		{
			if(use_combine_file)
			{
				for(int ipt=0; ipt<g_data->GetN(); ipt++)
				{
					g_data->GetPoint(ipt, x, y);
					if(y<0) {g_data->SetPoint(ipt, x, 0); g_data->SetPointError(ipt,0,0,0,0);}
                    if(is_blind) {g_data->SetPoint(ipt, x, 0); g_data->SetPointError(ipt,0,0,0,0);}
                }
			}
			else
			{
				for(int ibin = 1; ibin<h_sum_data->GetNbinsX()+1; ibin++)
				{
					if(h_sum_data->GetBinContent(ibin) < 0) {h_sum_data->SetBinContent(ibin, 0);}
                    if(is_blind) {h_sum_data->SetBinContent(ibin, 0);}
                }
			}
			for(int k=0; k<v_MC_histo.size(); k++)
			{
				if(!v_MC_histo[k]) {continue;} //Fakes templates can be null
				for(int ibin=0; ibin<v_MC_histo[k]->GetNbinsX(); ibin++)
				{
					if(v_MC_histo[k]->GetBinContent(ibin) < 0) {v_MC_histo[k]->SetBinContent(ibin, 0);}
				}
			}
		}


// # #    # #####  ###### #    #
// # ##   # #    # #       #  #
// # # #  # #    # #####    ##
// # #  # # #    # #        ##
// # #   ## #    # #       #  #
// # #    # #####  ###### #    #

        //-- Get indices of particular samples, sum the others into 1 single histo (used for ratio subplot)
    	TH1F* histo_total_MC = 0; //Sum of all MC samples

    	//Indices of important samples, for specific treatment
    	int index_tZq_sample = -9;
        int index_ttZ_sample = -9;
        int index_tWZ_sample = -9;
    	int index_NPL_sample = -9;

    	// cout<<"v_MC_histo.size() "<<v_MC_histo.size()<<endl;
    	// cout<<"MC_samples_legend.size() "<<MC_samples_legend.size()<<endl;

    	//Merge all the MC nominal histograms (contained in v_MC_histo)
    	for(int i=0; i<v_MC_histo.size(); i++)
    	{
    		if(!v_MC_histo[i]) {continue;} //Fakes templates may be null

    		// cout<<"MC_samples_legend[i] "<<MC_samples_legend[i]<<endl;

    		if(MC_samples_legend[i].Contains("tZq") )
    		{
    			if(index_tZq_sample<0) {index_tZq_sample = i;}
    			if(!h_tzq) {h_tzq = (TH1F*) v_MC_histo[i]->Clone();}
    			else {h_tzq->Add((TH1F*) v_MC_histo[i]->Clone());}
    			// if(doNot_stack_signal) continue; //don't stack
    		}
    		else if(MC_samples_legend[i].EndsWith("ttZ") )
    		{
    			if(index_ttZ_sample<0) {index_ttZ_sample = i;}
    			if(!h_ttz) {h_ttz = (TH1F*) v_MC_histo[i]->Clone();}
    			else {h_ttz->Add((TH1F*) v_MC_histo[i]->Clone());}
    			// if(doNot_stack_signal) continue; //don't stack
    		}
            else if(MC_samples_legend[i] == "tWZ") {index_tWZ_sample = i;}

    		// cout<<"Adding sample "<<MC_samples_legend[i]<<" to histo_total_MC"<<endl;

    		if(!histo_total_MC) {histo_total_MC = (TH1F*) v_MC_histo[i]->Clone();}
    		else {histo_total_MC->Add((TH1F*) v_MC_histo[i]->Clone());}
    	}

        //If histo_total_MC is null, variable was not found, skip it
        if(!histo_total_MC)
        {
            cout<<FRED("Error ! Variable '"<<total_var_list[ivar]<<"' not found ! Skip it...")<<endl;
            continue;
        }


// ##### #    #  ####  #####   ##    ####  #    #
//   #   #    # #        #    #  #  #    # #   #
//   #   ######  ####    #   #    # #      ####
//   #   #    #      #   #   ###### #      #  #
//   #   #    # #    #   #   #    # #    # #   #
//   #   #    #  ####    #   #    #  ####  #    #

		THStack* stack_MC = new THStack;

		//Add legend entries -- iterate backwards, so that last histo stacked is on top of legend
		//Also add MC histograms to the THStack
		for(int i=v_MC_histo.size()-1; i>=0; i--)
		{
			if(!v_MC_histo[i]) {continue;} //Some templates may be null

            if(region=="tWZ" && MC_samples_legend[i] == "tWZ") {continue;} //Put tWZ on top in that region

			stack_MC->Add(v_MC_histo[i]);
			// cout<<"Stacking sample "<<MC_samples_legend[i]<<" / integral "<<v_MC_histo[i]->Integral()<<endl;
		}
        if(region=="tWZ") {stack_MC->Add(v_MC_histo[index_tWZ_sample]);} //Put tWZ on top in that region


 // #####  #####  # #    #         ####    ##   #    # #####  #      ######
 // #    # #    # # #    #        #       #  #  ##  ## #    # #      #
 // #    # #    # # #    #         ####  #    # # ## # #    # #      #####
 // #####  #####  # #    # ###         # ###### #    # #####  #      #
 // #      #   #  #  #  #  ###    #    # #    # #    # #      #      #
 // #      #    # #   ##   ###     ####  #    # #    # #      ###### ######

        if(superimpose_EFThist)
        {
            for(int ipoint=0; ipoint<v_EFT_points.size(); ipoint++)
            {
                for(int isample=0; isample<v_EFT_samples.size(); isample++)
                {
                    if(normalize_EFThist) {v2_th1eft[ipoint][isample]->Scale(((TH1*) stack_MC->GetStack()->Last())->Integral()/(2*v2_th1eft[ipoint][isample]->Integral()));} //Normalize to half-integral of stack (arbitrary)  //NB: access stack integral via summed object '((TH1*) stack_MC->GetStack()->Last())'
                }
            }
        }


// ####### #
//    #    #       ######  ####  ###### #    # #####
//    #    #       #      #    # #      ##   # #    #
//    #    #       #####  #      #####  # #  # #    #
//    #    #       #      #  ### #      #  # # #    #
//    #    #       #      #    # #      #   ## #    #
//    #    ####### ######  ####  ###### #    # #####

        if(superimpose_EFThist) {nSampleGroups+= v_EFT_points.size() * v_EFT_samples.size();}

        int n_columns = ceil(nSampleGroups/2.) > 6 ? 6 : ceil(nSampleGroups/2.);
		TLegend* qw = 0;
        float x_left = 0.94-n_columns*0.10; //each column allocated same x-space
        // qw = new TLegend(x_left,0.88,0.99,0.99);
        // qw = new TLegend(x_left,0.75,0.95,0.89); //y-space for 2 rows
        if(superimpose_EFThist)
        {
            qw = new TLegend(x_left-0.05,0.72,0.95,0.87); //y-space for 2 rows
            qw->SetTextSize(0.03);
        }
        else
        {
            qw = new TLegend(x_left,0.78,0.94,0.87); //y-space for 2 rows
            qw->SetTextSize(0.04);
        }
        qw->SetNColumns(n_columns);
        qw->SetBorderSize(0);
        qw->SetFillStyle(0); //transparent
        // cout<<"x_left "<<x_left<<endl;
        // cout<<"ceil(nSampleGroups/2.) "<<ceil(nSampleGroups/2.)<<endl;

		//--Data on top of legend
        if(!is_blind && use_combine_file && g_data != 0) {qw->AddEntry(g_data, "Data" , "ep");}
        else if(!use_combine_file && h_sum_data != 0) {qw->AddEntry(h_sum_data, "Data" , "ep");}
        else {cout<<__LINE__<<BOLD(FRED(" : null data !"))<<endl;}

		for(int i=0; i<v_MC_histo.size(); i++)
		{
            // cout<<"MC_samples_legend[i] "<<MC_samples_legend[i]<<endl;
			if(!v_MC_histo[i]) {continue;} //Fakes templates can be null

			//Merged with other samples in legend -- don't add these
			// if(MC_samples_legend[i].Contains("xxx")
	        // || MC_samples_legend[i] == "xxx"
	        // ) {continue;}

            //Decide here which sample are mentioned in the legend
            if(MC_samples_legend[i].Contains("tZq")) {qw->AddEntry(v_MC_histo[i], "tZq", "f");}
            else if(MC_samples_legend[i].EndsWith("ttZ") ) {qw->AddEntry(v_MC_histo[i], "t#bar{t}Z", "f");}
            else if(region=="tWZ" && MC_samples_legend[i] == "tWZ") {qw->AddEntry(v_MC_histo[i], "tWZ", "f");}
            else if(MC_samples_legend[i] == "ttW") {qw->AddEntry(v_MC_histo[i], "t(#bar{t})X", "f");}
            else if(MC_samples_legend[i] == "WZ") {qw->AddEntry(v_MC_histo[i], "WZ", "f");}
            else if(MC_samples_legend[i] == "WWZ") {qw->AddEntry(v_MC_histo[i], "VV(V)", "f");}
            else if(MC_samples_legend[i] == "TTGamma_Dilep") {qw->AddEntry(v_MC_histo[i], "X+#gamma", "f");}
            else if(MC_samples_legend[i] == "TTbar_DiLep") {qw->AddEntry(v_MC_histo[i], "NPL", "f");}
            // else if(MC_samples_legend[i] == "ZZ4l") {qw->AddEntry(v_MC_histo[i], "ZZ", "f");}
            // else if(MC_samples_legend[i] == "tHq") {qw->AddEntry(v_MC_histo[i], "tX", "f");}
            // else if(MC_samples_legend[i] == "WZ") {qw->AddEntry(v_MC_histo[i], "VV(V)", "f");}
            // else if(MC_samples_legend[i] == "ttW") {qw->AddEntry(v_MC_histo[i], "t#bar{t}X", "f");}
            // else if(MC_samples_legend[i] == "GammaConv") {qw->AddEntry(v_MC_histo[i], "#gamma-conv.", "f");}
            // else if(MC_samples_legend[i] == "DY" ) {qw->AddEntry(v_MC_histo[i], "V+jets", "f");}
            // else if(MC_samples_legend[i] == "TTbar_DiLep" || MC_samples_legend[i] == "TTbar") {qw->AddEntry(v_MC_histo[i], "t#bar{t}", "f");}
		}

        for(int isample=0; isample<v_EFT_samples.size(); isample++)
        {
            for(int ipoint=0; ipoint<v_EFT_points.size(); ipoint++)
            {
                qw->AddEntry(v2_th1eft[ipoint][isample], v2_th1eft_labels[ipoint][isample], "L");
            }
        }


// #####  #####    ##   #    #
// #    # #    #  #  #  #    #
// #    # #    # #    # #    #
// #    # #####  ###### # ## #
// #    # #   #  #    # ##  ##
// #####  #    # #    # #    #

		//Canvas definition
		Load_Canvas_Style();
		TCanvas* c1 = new TCanvas("c1","c1", 1000, 800);
		// TCanvas* c1 = new TCanvas("c1","c1", 600, 800);
		c1->SetTopMargin(0.1);
		c1->SetBottomMargin(0.25);

		if(draw_logarithm) {c1->SetLogy();}

		//Draw stack
		stack_MC->Draw("hist");

		//Draw data
		if(data_notEmpty)
		{
			if(use_combine_file)
			{
				g_data->SetMarkerStyle(20);
				g_data->Draw("e0psame");
			}
			else
			{
				h_sum_data->SetMarkerStyle(20);
				h_sum_data->SetMinimum(0.) ;
				h_sum_data->Draw("e0psame");
			}
		}

        for(int i=0; i<v2_th1eft.size(); i++)
        {
            for(int j=0; j<v2_th1eft[i].size(); j++)
            {
                if(v2_th1eft[i][j]) {v2_th1eft[i][j]->Draw("same hist");}
            }
        }

		//Superimpose shape of signal
		// if(doNot_stack_signal)
		// {
		// 	if(h_tzq != 0) {h_tzq->Draw("same hist");}
		// 	if(h_ttz != 0) {h_ttz->Draw("same hist");}
		// }

        qw->Draw("same"); //Draw legend

        //-- Possible option to manipulate each legend entry individually //Not needed
        // TList* list_leg_entries = qw->GetListOfPrimitives();
        // TIter next(list_leg_entries);
        // TObject *obj_tmp;
        // TLegendEntry *leg_entry_tmp;
        // int i = 0;
        // while((obj_tmp = next()))
        // {
        //     leg_entry_tmp = (TLegendEntry*) obj_tmp;
        //     i++;
        //     if(((TString) leg_entry_tmp->GetLabel()).Contains(" (")) {leg_entry_tmp->SetTextSize(qw->GetTextSize()-0.005);}
        // }


// #   # #    #   ##   #    #
//  # #  ##  ##  #  #   #  #
//   #   # ## # #    #   ##
//   #   #    # ######   ##
//   #   #    # #    #  #  #
//   #   #    # #    # #    #

        c1->Update();

        //Set Yaxis maximum & minimum
        double ymax_EFT = -1; //Keep track of y-axis maximum value for EFT histograms
        if(superimpose_EFThist && !normalize_EFThist) //Check whether max y-value is from EFT histo
        {
            for(int ipoint=0; ipoint<v_EFT_points.size(); ipoint++)
            {
                for(int isample=0; isample<v_EFT_samples.size(); isample++)
                {
                    if(v2_th1eft[ipoint][isample]->GetMaximum() > ymax_EFT) {ymax_EFT = v2_th1eft[ipoint][isample]->GetMaximum();}
                }
            }
        }

        double ymax = 0;
        if(use_combine_file)
		{
			Long64_t locmax = TMath::LocMax(g_data->GetN(), g_data->GetY()); //the corresponding x value can be obtained with double xmax = gr->GetX()[locmax];
            if(data_notEmpty) {ymax = g_data->GetY()[locmax];} //Data ymax
            if(ymax < stack_MC->GetMaximum()) {ymax = stack_MC->GetMaximum();} //MC ymax
		}
		else if(!use_combine_file)
		{
            if(data_notEmpty) {ymax = h_sum_data->GetMaximum();} //Data ymax
            if(ymax < stack_MC->GetMaximum()) {ymax = stack_MC->GetMaximum();} //MC ymax
            if(ymax < ymax_EFT) {ymax = ymax_EFT;} //EFT ymax
		}
        ymax*= 1.4;
        stack_MC->SetMaximum(ymax);
        // if(ymax > qw->GetY1()) {ymax = qw->GetY1();} //Avoid overlap with TLegend
        // cout<<"qw->GetY1() "<<qw->GetY1()<<endl;
        // cout<<"c1->GetUymax() "<<c1->GetUymax()<<endl;

		stack_MC->SetMinimum(0.0001); //Remove '0' label

		if(draw_logarithm) //?
		{
			stack_MC->SetMinimum(0.5);
			stack_MC->SetMaximum(stack_MC->GetMaximum()*6);
		}

        //-- Trick so that different plots made for different EFT scan points all share the same y-axis range: save ymax for the first point (according to boolean in arg.), then read the store value to set ymax for next points
        if(EFTpoint != "")
        {
            float* pointer_tmp = ymax_fixed1; if(ivar>0) pointer_tmp = ymax_fixed2;
            if(store_ymax_fixed) {*pointer_tmp = ymax;}
            else {stack_MC->SetMaximum(*pointer_tmp);}
            // cout<<"*pointer_tmp "<<*pointer_tmp<<endl;
        }

        c1->Modified();


// ###### ###### #####     ####  ###### #    #
// #      #        #      #    # #      ##   #
// #####  #####    #      #      #####  # #  #
// #      #        #      #  ### #      #  # #
// #      #        #      #    # #      #   ##
// ###### #        #       ####  ###### #    #

        TH1F* h_GEN = 0;
        if(superimpose_GENhisto && drawInputVars) {Get_Pointer_GENHisto(h_GEN, var_list[ivar]);}

        if(h_GEN)
        {
            h_GEN->Scale(h_sum_data->Integral()/(2*h_GEN->Integral())); //Arbitrary norm (half of data)
            h_GEN->SetLineColor(1);
            h_GEN->SetLineWidth(2);

            h_GEN->Draw("hist SAME");
        }


// ###### #####  #####   ####  #####   ####      ####  #####   ##    ####  #    #
// #      #    # #    # #    # #    # #         #        #    #  #  #    # #   #
// #####  #    # #    # #    # #    #  ####      ####    #   #    # #      ####
// #      #####  #####  #    # #####       #         #   #   ###### #      #  #
// #      #   #  #   #  #    # #   #  #    #    #    #   #   #    # #    # #   #
// ###### #    # #    #  ####  #    #  ####      ####    #   #    #  ####  #    #

		//-- Compute sqrt of quadratic errors
		if(draw_errors)
		{
			for(int ibin=0; ibin<nofbins; ibin++)
			{
				v_eyh[ibin] = pow(v_eyh[ibin], 0.5);
				v_eyl[ibin] = pow(v_eyl[ibin], 0.5);

				// if(ibin > 0) {continue;} //cout only first bin
				// cout<<"x = "<<v_x[ibin]<<endl;    cout<<", y = "<<v_y[ibin]<<endl;    cout<<", eyl = "<<v_eyl[ibin]<<endl;    cout<<", eyh = "<<v_eyh[ibin]<<endl; //cout<<", exl = "<<v_exl[ibin]<<endl;    cout<<", exh = "<<v_exh[ibin]<<endl;
			}
		}

		//Use pointers to vectors : need to give the adress of first element (all other elements can then be accessed iteratively)
		double* eyl = &v_eyl[0];
		double* eyh = &v_eyh[0];
		double* exl = &v_exl[0];
		double* exh = &v_exh[0];
		double* xx = &v_x[0];
		double* yy = &v_y[0];

		//Create TGraphAsymmErrors with the error vectors / (x,y) coordinates --> Can superimpose it on plot
		TGraphAsymmErrors* gr_error = 0;

		gr_error = new TGraphAsymmErrors(nofbins,xx,yy,exl,exh,eyl,eyh);
		gr_error->SetFillStyle(3002);
		gr_error->SetFillColor(kBlack);
		gr_error->Draw("e2 same"); //Superimposes the uncertainties on stack


// #####    ##   ##### #  ####
// #    #  #  #    #   # #    #
// #    # #    #   #   # #    #
// #####  ######   #   # #    #
// #   #  #    #   #   # #    #
// #    # #    #   #   #  ####

// #####  #       ####  #####
// #    # #      #    #   #
// #    # #      #    #   #
// #####  #      #    #   #
// #      #      #    #   #
// #      ######  ####    #

		//-- create subpad to plot ratio
		TPad *pad_ratio = new TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 1.0);
		pad_ratio->SetTopMargin(0.75);
		pad_ratio->SetFillColor(0);
		pad_ratio->SetFillStyle(0);
		pad_ratio->SetGridy(1);
		pad_ratio->Draw();
		pad_ratio->cd(0);

		if(use_combine_file && data_notEmpty) //Copy the content of the data graph into a TH1F (NB : symmetric errors...?)
		{
			if(!v_MC_histo[0]) {cout<<__LINE__<<FRED("Error : v_MC_histo[0] is null ! Abort")<<endl; return;}

	        h_sum_data = (TH1F*) v_MC_histo[0]->Clone(); //To clone binning of the MC histos
			h_sum_data->SetFillColor(kBlack);
			h_sum_data->SetLineColor(kBlack);
			// cout<<"h_sum_data->GetNbinsX() "<<h_sum_data->GetNbinsX()<<endl;

	        for(int ipt=0; ipt<g_data->GetN(); ipt++)
	        {
	            g_data->GetPoint(ipt, x, y);
	            double error = g_data->GetErrorY(ipt);

	            h_sum_data->SetBinContent(ipt+1, y);
	            h_sum_data->SetBinError(ipt+1, error);
	        }
		}

		TH1F* histo_ratio_data = 0;
		if(data_notEmpty)
		{
			histo_ratio_data = (TH1F*) h_sum_data->Clone();

			//debug printout
			// cout<<"h_sum_data->GetBinContent(5) "<<h_sum_data->GetBinContent(5)<<endl;
			// cout<<"h_sum_data->GetBinError(5) "<<h_sum_data->GetBinError(5)<<endl;
			// cout<<"histo_total_MC->GetBinContent(5) "<<histo_total_MC->GetBinContent(5)<<endl;
			// cout<<"histo_total_MC->GetBinError(5) "<<histo_total_MC->GetBinError(5)<<endl;

			if(!show_pulls_ratio)
			{
				//To get error bars correct in ratio plot, must only account for errors from data, not MC ! (MC error shown as gray bad)
				for(int ibin=1; ibin<histo_total_MC->GetNbinsX()+1; ibin++)
				{
					histo_total_MC->SetBinError(ibin, 0);
				}

				histo_ratio_data->Divide(histo_total_MC);
			} //Ratio
		 	else //--- Compute pull distrib
			{
				for(int ibin=1; ibin<histo_ratio_data->GetNbinsX()+1; ibin++)
				{
					//Add error on signal strength (since we rescale signal manually)
					// double bin_error_mu = v_MC_histo.at(index_tZq_sample)->GetBinError(ibin) * sig_strength_err;
					// cout<<"bin_error_mu = "<<bin_error_mu<<endl;

					double bin_error_mu = 0; //No sig strength uncert. for prefit ! //-- postfit -> ?

					//Quadratic sum of systs, stat error, and sig strength error
					double bin_error = pow(pow(histo_total_MC->GetBinError(ibin), 2) + pow(histo_ratio_data->GetBinError(ibin), 2) + pow(bin_error_mu, 2), 0.5);

					// if(ibin==1) {cout<<"Data = "<<histo_ratio_data->GetBinContent(1)<<" / Total MC = "<<histo_total_MC->GetBinContent(1)<<" / error = "<<bin_error<<endl;}

					if(!histo_total_MC->GetBinError(ibin)) {histo_ratio_data->SetBinContent(ibin,-99);} //Don't draw null markers
					else{histo_ratio_data->SetBinContent(ibin, (histo_ratio_data->GetBinContent(ibin) - histo_total_MC->GetBinContent(ibin)) / bin_error );}
				}
			}

			//debug printout
			// cout<<"histo_ratio_data->GetBinContent(5) "<<histo_ratio_data->GetBinContent(5)<<endl;
			// cout<<"histo_ratio_data->GetBinError(5) "<<histo_ratio_data->GetBinError(5)<<endl;

			//Don't draw null data
			for(int ibin=1; ibin<histo_ratio_data->GetNbinsX()+1; ibin++)
			{
				// cout<<"histo_ratio_data["<<ibin<<"] = "<<histo_ratio_data->GetBinContent(ibin)<<endl;

				if(std::isnan(histo_ratio_data->GetBinContent(ibin)) || std::isinf(histo_ratio_data->GetBinContent(ibin)) || histo_ratio_data->GetBinContent(ibin) == 0) {histo_ratio_data->SetBinContent(ibin, -99);}
			}
		}
		else {histo_ratio_data = (TH1F*) histo_total_MC->Clone();}

		if(show_pulls_ratio) {histo_ratio_data->GetYaxis()->SetTitle("Pulls");}
		else {histo_ratio_data->GetYaxis()->SetTitle("Data/MC");}
		histo_ratio_data->GetYaxis()->SetTickLength(0.);
        histo_ratio_data->GetYaxis()->SetTitleOffset(1.15);
        // histo_ratio_data->GetYaxis()->SetTitleOffset(1.2);
        histo_ratio_data->GetYaxis()->SetLabelSize(0.04);
        // histo_ratio_data->GetYaxis()->SetLabelSize(0.048);
		histo_ratio_data->GetXaxis()->SetLabelFont(42);
		histo_ratio_data->GetYaxis()->SetLabelFont(42);
		histo_ratio_data->GetXaxis()->SetTitleFont(42);
		histo_ratio_data->GetYaxis()->SetTitleFont(42);
        histo_ratio_data->GetYaxis()->SetNdivisions(503); //grid draw on primary tick marks only
		histo_ratio_data->GetXaxis()->SetNdivisions(505);
		histo_ratio_data->GetYaxis()->SetTitleSize(0.06);
		histo_ratio_data->GetXaxis()->SetTickLength(0.04);
		histo_ratio_data->SetMarkerStyle(20);
		histo_ratio_data->SetMarkerSize(1.2); //changed from 1.4

		//NB : when using SetMaximum(), points above threshold are simply not drawn
		//So for ratio plot, even if point is above max but error compatible with 1, point/error bar not represented!
		if(show_pulls_ratio)
		{
			histo_ratio_data->SetMinimum(-2.99);
			histo_ratio_data->SetMaximum(2.99);
		}
		else
		{
            histo_ratio_data->SetMinimum(0.4);
            histo_ratio_data->SetMaximum(1.6);
            // histo_ratio_data->SetMaximum(2.2); //503 divisions, label size 0.048
            // histo_ratio_data->SetMinimum(-0.1);
		}

		if(drawInputVars) {histo_ratio_data->GetXaxis()->SetTitle(Get_Variable_Name(total_var_list[ivar]));}
		else
		{
            // histo_ratio_data->GetXaxis()->SetTitle(classifier_name+" (vs "+template_name + ")");
            histo_ratio_data->GetXaxis()->SetTitle(total_var_list[ivar]);

            //Hardcode NN output nodes names...?
            if(total_var_list[ivar] == "NN") {histo_ratio_data->GetXaxis()->SetTitle("NN output");}
            else if(total_var_list[ivar] == "NN0") {histo_ratio_data->GetXaxis()->SetTitle("NN (tZq node)");}
            else if(total_var_list[ivar] == "NN1" && NN_nNodes == 3) {histo_ratio_data->GetXaxis()->SetTitle("NN (ttZ node)");}
            else if(total_var_list[ivar] == "NN2" && NN_nNodes == 3) {histo_ratio_data->GetXaxis()->SetTitle("NN (Bkgs node)");}

			if(template_name == "categ") //Vertical text X labels (categories names)
			{
                //Hard-coded -- should automate labels
                {
                    const char *labels[10]  = {"1bj,2j","1bj,3j","1bj,4j","1bj,5j","1bj,6j","2bj,2j","2bj,3j","2bj,4j","2bj,5j","2bj,6j"};
                    for(int i=1;i<=10;i++) {histo_ratio_data->GetXaxis()->SetBinLabel(i,labels[i-1]);}
                }

                histo_ratio_data->GetXaxis()->SetTitle("");
                // histo_ratio_data->GetXaxis()->SetTitle("Categ.");
				histo_ratio_data->GetXaxis()->SetLabelSize(0.06);
				histo_ratio_data->GetXaxis()->SetLabelOffset(0.02);
                histo_ratio_data->LabelsOption("v", "X"); //X labels vertical
			}
		}

		pad_ratio->cd(0);
		if(show_pulls_ratio) {histo_ratio_data->Draw("HIST P");} //Draw ratio points
		else {histo_ratio_data->Draw("E1X0 P");} //Draw ratio points ; E1 : perpendicular lines at end ; X0 : suppress x errors

// ###### #####  #####   ####  #####   ####     #####    ##   ##### #  ####
// #      #    # #    # #    # #    # #         #    #  #  #    #   # #    #
// #####  #    # #    # #    # #    #  ####     #    # #    #   #   # #    #
// #      #####  #####  #    # #####       #    #####  ######   #   # #    #
// #      #   #  #   #  #    # #   #  #    #    #   #  #    #   #   # #    #
// ###### #    # #    #  ####  #    #  ####     #    # #    #   #   #  ####

		TGraphAsymmErrors* gr_ratio_error = 0;
		if(draw_errors)
		{
			//Copy previous TGraphAsymmErrors, then modify it -> error TGraph for ratio plot
			TGraphAsymmErrors *thegraph_tmp;
			double *theErrorX_h;
			double *theErrorY_h;
			double *theErrorX_l;
			double *theErrorY_l;
			double *theY;
			double *theX;

			thegraph_tmp = (TGraphAsymmErrors*) gr_error->Clone();
			theErrorX_h = thegraph_tmp->GetEXhigh();
			theErrorY_h = thegraph_tmp->GetEYhigh();
			theErrorX_l = thegraph_tmp->GetEXlow();
			theErrorY_l = thegraph_tmp->GetEYlow();
			theY        = thegraph_tmp->GetY() ;
			theX        = thegraph_tmp->GetX() ;

			//Divide error --> ratio
			for(int i=0; i<thegraph_tmp->GetN(); i++)
			{
				theErrorY_l[i] = theErrorY_l[i]/theY[i];
				theErrorY_h[i] = theErrorY_h[i]/theY[i];
				theY[i]=1; //To center the filled area around "1"
			}

			gr_ratio_error = new TGraphAsymmErrors(thegraph_tmp->GetN(), theX , theY ,  theErrorX_l, theErrorX_h, theErrorY_l, theErrorY_h);
			gr_ratio_error->SetFillStyle(3002);
			gr_ratio_error->SetFillColor(kBlack);
			// gr_ratio_error->SetFillColor(kCyan);

			pad_ratio->cd(0);
			if(!show_pulls_ratio) {gr_ratio_error->Draw("e2 same");} //Draw error bands in ratio plot

            //-- Add sub-legend here ? (stat or stat+syst)
		} //draw errors


//  ####   ####   ####  #    # ###### ##### #  ####   ####
// #    # #    # #      ##  ## #        #   # #    # #
// #      #    #  ####  # ## # #####    #   # #       ####
// #      #    #      # #    # #        #   # #           #
// #    # #    # #    # #    # #        #   # #    # #    #
//  ####   ####   ####  #    # ######   #   #  ####   ####

		//-- Draw ratio y-lines manually
		TH1F *h_line1 = 0;
		TH1F *h_line2 = 0;
		if(data_notEmpty)
		{
			h_line1 = new TH1F("","",this->nbins, h_sum_data->GetXaxis()->GetXmin(), h_sum_data->GetXaxis()->GetXmax());
			h_line2 = new TH1F("","",this->nbins, h_sum_data->GetXaxis()->GetXmin(), h_sum_data->GetXaxis()->GetXmax());
			// TH1F *h_line3 = new TH1F("","",this->nbins, h_sum_data->GetXaxis()->GetXmin(), h_sum_data->GetXaxis()->GetXmax());
			for(int ibin=1; ibin<this->nbins +1; ibin++)
			{
				if(show_pulls_ratio)
				{
					h_line1->SetBinContent(ibin, -1);
					h_line2->SetBinContent(ibin, 1);
				}
				else
				{
                    h_line1->SetBinContent(ibin, 0.75);
					h_line2->SetBinContent(ibin, 1.25);
                    // h_line1->SetBinContent(ibin, 0.5);
					// h_line2->SetBinContent(ibin, 1.5);
				}
			}
			h_line1->SetLineStyle(6);
			h_line2->SetLineStyle(6);
			h_line1->Draw("hist same");
			h_line2->Draw("hist same");
		}

		double xmax_stack = stack_MC->GetXaxis()->GetXmax();
		double xmin_stack = stack_MC->GetXaxis()->GetXmin();
		TString Y_label = "Events";
		if(data_notEmpty)
		{
			double xmax_data = h_sum_data->GetXaxis()->GetXmax();
			double xmin_data = h_sum_data->GetXaxis()->GetXmin();
			Y_label = "Events / " + Convert_Number_To_TString( (xmax_data - xmin_data) / h_sum_data->GetNbinsX(), 2); //Automatically get the Y label depending on binning
		}

		if(stack_MC!= 0)
		{
			stack_MC->GetXaxis()->SetLabelFont(42);
			stack_MC->GetYaxis()->SetLabelFont(42);
			stack_MC->GetYaxis()->SetTitleFont(42);
			stack_MC->GetYaxis()->SetTitleSize(0.06);
			stack_MC->GetYaxis()->SetTickLength(0.04);
			stack_MC->GetXaxis()->SetLabelSize(0.0);
			stack_MC->GetYaxis()->SetLabelSize(0.048);
			stack_MC->GetXaxis()->SetNdivisions(505);
			stack_MC->GetYaxis()->SetNdivisions(506);
            stack_MC->GetYaxis()->SetTitleOffset(1.28);
            // stack_MC->GetYaxis()->SetTitleOffset(1.2);
			stack_MC->GetYaxis()->SetTitle(Y_label);
		}

	//----------------
	// CAPTIONS //
	//----------------
	// -- using https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines

        bool draw_cms_prelim_label = true;

		float l = c1->GetLeftMargin();
		float t = c1->GetTopMargin();

		TString cmsText = "CMS";
		TLatex latex;
		latex.SetNDC();
		latex.SetTextColor(kBlack);
		latex.SetTextFont(61);
		latex.SetTextAlign(11);
		latex.SetTextSize(0.06);
		if(draw_cms_prelim_label) {latex.DrawLatex(l + 0.01, 0.92, cmsText);}

		TString extraText = "Preliminary";
		latex.SetTextFont(52);
		latex.SetTextSize(0.05);
		if(draw_cms_prelim_label) {latex.DrawLatex(l + 0.12, 0.92, extraText);}

		float lumi = lumiValue;
		TString lumi_ts = Convert_Number_To_TString(lumi);
		lumi_ts += " fb^{-1} (13 TeV)";
		latex.SetTextFont(42);
		latex.SetTextAlign(31);
		latex.SetTextSize(0.04);
        latex.DrawLatex(0.96, 0.92,lumi_ts);

		//------------------
		//-- channel info
		TLatex text2 ;
		text2.SetNDC();
		text2.SetTextAlign(13);
		text2.SetTextSize(0.045);
		text2.SetTextFont(42);

		TString info_data = "l^{#pm}l^{#pm}l^{#pm}";
        if (channel=="eee")    info_data = "eee";
		else if (channel=="eeu")  info_data = "ee#mu";
		else if (channel=="uue")  info_data = "#mu#mu e";
		else if (channel=="uuu") info_data = "#mu#mu #mu";

		// if(h_sum_data->GetBinContent(h_sum_data->GetNbinsX() ) > h_sum_data->GetBinContent(1) ) {text2.DrawLatex(0.55,0.87,info_data);}
		// else {text2.DrawLatex(0.20,0.87,info_data);}
		if(!superimpose_EFThist) {text2.DrawLatex(0.23,0.86,info_data);}


// #    # #####  # ##### ######     ####  #    # ##### #####  #    # #####
// #    # #    # #   #   #         #    # #    #   #   #    # #    #   #
// #    # #    # #   #   #####     #    # #    #   #   #    # #    #   #
// # ## # #####  #   #   #         #    # #    #   #   #####  #    #   #
// ##  ## #   #  #   #   #         #    # #    #   #   #      #    #   #
// #    # #    # #   #   ######     ####   ####    #   #       ####    #

        TString outdir = "";
		if(drawInputVars)
		{
            outdir = "plots/input_vars/";
            mkdir(outdir.Data(), 0777);
            outdir+= cat_tmp;
            mkdir(outdir.Data(), 0777);
            outdir+= "/" + lumiName;
            mkdir(outdir.Data(), 0777);
		}
		else
		{
            outdir = "plots/templates/";
			mkdir(outdir.Data(), 0777);
            outdir+= cat_tmp;
            mkdir(outdir.Data(), 0777);
            outdir+= "/"+lumiName;
            mkdir(outdir.Data(), 0777);
            if(prefit) {outdir+= "/prefit";}
            else {outdir+= "/postfit";}
            mkdir(outdir.Data(), 0777);
            if(categorization_strategy>0 && make_SMvsEFT_templates_plots) {outdir+= "/strategy" + Convert_Number_To_TString(categorization_strategy); if(EFTpoint!="") {outdir+= "param";}}
            mkdir(outdir.Data(), 0777);
		}
        outdir+= "/";

		//Output
		TString output_plot_name;
		if(drawInputVars)
		{
			output_plot_name = outdir + total_var_list[ivar];
		}
		else
		{
            output_plot_name = outdir + total_var_list[ivar] + (EFTpoint==""? "":"_"+EFTpoint) + "_template_" + signal_process;
		}
		if(channel != "") {output_plot_name+= "_" + channel;}
		output_plot_name+= this->filename_suffix;
		if(draw_logarithm) {output_plot_name+= "_log";}
		output_plot_name+= this->plot_extension;

		c1->SaveAs(output_plot_name);

		if(data_notEmpty)
		{
			delete h_line1; h_line1 = NULL;
			delete h_line2; h_line2 = NULL;
		}
        if(h_GEN) {delete h_GEN;}
		if(draw_errors) {delete gr_error; delete gr_ratio_error; gr_error = NULL; gr_ratio_error = NULL;}
		delete pad_ratio; pad_ratio = NULL;

		delete c1; c1 = NULL;
		delete qw; qw = NULL;
		delete stack_MC; stack_MC = NULL;
        for(int i=0; i<v2_th1eft.size(); i++)
        {
            for(int j=0; j<v2_th1eft[i].size(); j++)
            {
                if(v2_th1eft[i][j]) {delete v2_th1eft[i][j];}
            }
        }
		if(use_combine_file) {delete g_data; g_data = NULL;}
	} //Var loop

	file_input->Close();

	return;
}





























//--------------------------------------------
//  ######   #######  ##     ## ########     ###    ########  ########
// ##    ## ##     ## ###   ### ##     ##   ## ##   ##     ## ##
// ##       ##     ## #### #### ##     ##  ##   ##  ##     ## ##
// ##       ##     ## ## ### ## ########  ##     ## ########  ######
// ##       ##     ## ##     ## ##        ######### ##   ##   ##
// ##    ## ##     ## ##     ## ##        ##     ## ##    ##  ##
//  ######   #######  ##     ## ##        ##     ## ##     ## ########

//  ######  ##     ##    ###    ########  ########  ######
// ##    ## ##     ##   ## ##   ##     ## ##       ##    ##
// ##       ##     ##  ##   ##  ##     ## ##       ##
//  ######  ######### ##     ## ########  ######    ######
//       ## ##     ## ######### ##        ##             ##
// ##    ## ##     ## ##     ## ##        ##       ##    ##
//  ######  ##     ## ##     ## ##        ########  ######
//--------------------------------------------

//Compare shapes of several (hard-coded) samples for some variables
//Reads pre-existing histograms, which must have been previously produced with Draw_Templates(drawInputVars=True)
//== Add ratio pad ==
void TopEFT_analysis::Compare_TemplateShapes_Processes(TString template_name, TString channel)
{
//--- OPTIONS --------------------------------
//--------------------------------------------
	bool drawInputVars = true;

	bool normalize = true;

    TString type = "ttZ"; //'' or 'tZq' or 'ttZ' --> Compare corresponding private/central samples

    vector<TString> total_var_list;
    total_var_list.push_back("mTW");
    total_var_list.push_back("mHT");
    total_var_list.push_back("Mass_3l");
    total_var_list.push_back("maxEtaJet");
    total_var_list.push_back("jPrimeAbsEta");
    total_var_list.push_back("channel");
    total_var_list.push_back("njets");
    total_var_list.push_back("nbjets");
    total_var_list.push_back("metEt");
    // total_var_list.push_back("lAsymmetry");
    // total_var_list.push_back("maxDelPhiLL");
    // total_var_list.push_back("maxDeepCSV");
    // total_var_list.push_back("deepCSV_2nd");
    // total_var_list.push_back("leptonCharge");
    // total_var_list.push_back("cosThetaStarPolTop");
    // total_var_list.push_back("cosThetaStarPolZ");
    // total_var_list.push_back("recoZ_Pt");
    // total_var_list.push_back("recoZ_Eta");
    // total_var_list.push_back("recoZ_M");
    // total_var_list.push_back("recoLepTop_Pt");
    // total_var_list.push_back("recoLepTop_Eta");
    // total_var_list.push_back("recoLepTop_M");
    // total_var_list.push_back("TopZsystem_Pt");
    // total_var_list.push_back("TopZsystem_M");

    TString theyear = "2017"; //2016,2017,2018

    //-- Hardcode samples here... or could filter the main sample list
	vector<TString> v_samples; vector<TString> v_groups; vector<int> v_colors;
    v_samples.push_back("tZq"); v_groups.push_back("tZq (Central)"); v_colors.push_back(kRed);
    v_samples.push_back("PrivMC_tZq"); v_groups.push_back("tZq (Private)"); v_colors.push_back(kBlue);
    // v_samples.push_back("ttZ"); v_groups.push_back("ttZ (Central)"); v_colors.push_back(kRed);
    // v_samples.push_back("PrivMC_ttZ"); v_groups.push_back("ttZ (Private)"); v_colors.push_back(kBlue);

    if(type == "tZq")
    {
        v_samples.resize(0); v_groups.resize(0); v_colors.resize(0);
        v_samples.push_back("tZq"); v_groups.push_back("tZq (Central)"); v_colors.push_back(kRed);
        v_samples.push_back("PrivMC_tZq"); v_groups.push_back("tZq (Private)"); v_colors.push_back(kBlue);
    }
    else if(type == "ttZ")
    {
        v_samples.resize(0); v_groups.resize(0); v_colors.resize(0);
        v_samples.push_back("ttZ"); v_groups.push_back("ttZ (Central)"); v_colors.push_back(kRed);
        v_samples.push_back("PrivMC_ttZ"); v_groups.push_back("ttZ (Private)"); v_colors.push_back(kBlue);
    }

    vector<TString> v_syst;
    v_syst.push_back("");
    // v_syst.push_back("JESUp");

//--------------------------------------------
//--------------------------------------------

    if(!drawInputVars)
    {
        total_var_list.clear();
        if(classifier_name == "BDT") {total_var_list.push_back(classifier_name);}
        else
        {
            if(NN_nNodes==1) {total_var_list.push_back(classifier_name);}
            else
            {
                for(int inode=0; inode<NN_nNodes; inode++) {total_var_list.push_back(classifier_name + Convert_Number_To_TString(inode));}
            }
        }
    }

    cout<<endl<<BYEL("                          ")<<endl<<endl;
    if(drawInputVars) {cout<<FYEL("--- Producing Input Vars Plots / channel : "<<channel<<" ---")<<endl;}
    else {cout<<FYEL("--- Producing Template Plots / channel : "<<channel<<" ---")<<endl;}
    cout<<endl<<BYEL("                          ")<<endl<<endl;

//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

    // TString cat_tmp = (region=="") ? "SR" : region+"Cat";
    TString cat_tmp = region;

	//Get input TFile
	TString input_name;
	if(drawInputVars)
	{
        input_name = "outputs/ControlHistograms_" + cat_tmp + "_" + lumiName + filename_suffix +".root";
	}
	else //OBSOLETE
	{
        input_name = "outputs/Templates_" + classifier_name + template_name + "_" + cat_tmp + "_" + lumiName + filename_suffix + ".root";
    }

	if(!Check_File_Existence(input_name))
	{
		cout<<FRED("File "<<input_name<<" (prefit templates) not found ! Did you specify the region/background ? Abort")<<endl;
		return;
	}
	cout<<endl<<endl<<endl;

	//Input file containing histos
	TFile* file_input = 0;
	file_input = TFile::Open(input_name);
    cout<<DIM("Opening file: "<<input_name<<"")<<endl;

	//TH1F* to retrieve distributions
	TH1F* h_tmp = 0; //Tmp storing histo

	vector<vector<vector<TH1F*>>> v3_histos_var_sample_syst(total_var_list.size()); //store histos, for each var/sample/syst


// #       ####   ####  #####   ####
// #      #    # #    # #    # #
// #      #    # #    # #    #  ####
// #      #    # #    # #####       #
// #      #    # #    # #      #    #
// ######  ####   ####  #       ####

    for(int ivar=0; ivar<total_var_list.size(); ivar++)
    {
        v3_histos_var_sample_syst[ivar].resize(v_samples.size());

    	//Combine output : all histos are given for subcategories --> Need to sum them all
    	for(int ichan=0; ichan<channel_list.size(); ichan++)
    	{
    		if(channel_list[ichan] != channel) {continue;}

    		for(int isample = 0; isample<v_samples.size(); isample++)
    		{
    			// cout<<endl<<UNDL(FBLU("-- Sample : "<<v_samples[isample]<<" : "))<<endl;

    			v3_histos_var_sample_syst[ivar][isample].resize(v_syst.size());

                TFile* f;
                f = file_input;

    			TString samplename = v_samples[isample];

                for(int isyst=0; isyst<v_syst.size(); isyst++)
                {
    				if(v_samples[isample].Contains("Fake") && !v_syst[isyst].Contains("Clos") && !v_syst[isyst].Contains("FR") && v_syst[isyst] != "") {continue;}

    				// cout<<"syst "<<v_syst[isyst]<<endl;

                    h_tmp = 0;

                    TString histo_name = total_var_list[ivar];
                    if(channel_list[ichan] != "") {histo_name+= "_" + channel_list[ichan];}
                    // histo_name+= "_" + v_lumiYears[iyear];
                    histo_name+= "_" + theyear;
                    histo_name+= "__" + samplename;

        			if(!f->GetListOfKeys()->Contains(histo_name) ) {cout<<ITAL("Histogram '"<<histo_name<<"' : not found ! Skip...")<<endl; continue;}

        			h_tmp = (TH1F*) f->Get(histo_name);
    				h_tmp->SetDirectory(0); //Dis-associate from TFile
        			// cout<<"histo_name "<<histo_name<<endl;
                    // cout<<"h_tmp->Integral() = "<<h_tmp->Integral()<<endl;

     //  ####   ####  #       ####  #####   ####
     // #    # #    # #      #    # #    # #
     // #      #    # #      #    # #    #  ####
     // #      #    # #      #    # #####       #
     // #    # #    # #      #    # #   #  #    #
     //  ####   ####  ######  ####  #    #  ####

        			//Use color vector filled in main()
        			// h_tmp->SetFillStyle(1001);
    				h_tmp->SetFillColor(kWhite);

    				h_tmp->SetLineColor(v_colors[isample]);
                    if(samplename.Contains("PrivMC")) {h_tmp->SetLineStyle(2);}

    				//HARDCODED
    				if(v_syst[isyst] == "JESUp") {h_tmp->SetLineColor(kRed);}
    				else if(v_syst[isyst] == "JESDown") {h_tmp->SetLineColor(kBlue);}

    				// h_tmp->SetLineColor(v_colors[isample]+isyst);
    				// cout<<"v_colors[isample] "<<v_colors[isample]<<endl;

        			h_tmp->SetLineWidth(3);

    				h_tmp->SetMaximum(h_tmp->GetMaximum()*1.5);
    				if(normalize) {h_tmp->SetMaximum(0.5);}

                    if(v_syst[isyst] != "") {h_tmp->SetLineStyle(2);}

        			if(normalize) {h_tmp->Scale(1./h_tmp->Integral() );}

                    v3_histos_var_sample_syst[ivar][isample][isyst] = (TH1F*) h_tmp->Clone();

    				// cout<<"v3_histos_var_sample_syst[ivar]["<<isample<<"]["<<isyst<<"]->Integral() "<<v3_histos_var_sample_syst[ivar][isample][isyst]->Integral()<<endl;

        			delete h_tmp; h_tmp = 0;
                } //end syst loop

    			// f->Close();
    		} //end sample loop
    	} //subcat loop

    } //var loop

// #####  #####    ##   #    #
// #    # #    #  #  #  #    #
// #    # #    # #    # #    #
// #    # #####  ###### # ## #
// #    # #   #  #    # ##  ##
// #####  #    # #    # #    #

    for(int ivar=0; ivar<v3_histos_var_sample_syst.size(); ivar++)
    {
        //Canvas definition
        Load_Canvas_Style();
        gStyle->SetOptTitle(1);
        TCanvas* c = new TCanvas("", "", 1000, 800);
        c->SetTopMargin(0.1);
        // c->SetRightMargin(0.1);
        c->SetBottomMargin(0.25);
        c->cd();

        // c->SetLogy();

        TLegend* qw;
        qw = new TLegend(0.75,.70,1.,1.);

    	for(int isample=0; isample<v3_histos_var_sample_syst[ivar].size(); isample++)
    	{
    		TString systlist = "";

    		for(int isyst=0; isyst<v_syst.size(); isyst++)
    		{
    			if(v_samples[isample].Contains("Fake") && !v_syst[isyst].Contains("Clos") && !v_syst[isyst].Contains("FR") && v_syst[isyst] != "") {continue;}

    			if(v3_histos_var_sample_syst[ivar][isample][isyst] == 0) {cout<<"Null histo ! Skip"<<endl; continue;}

    			if(isample == 0)
    			{
    				if(v_syst[isyst] != "")
    				{
    					systlist+= " / " + v_syst[isyst];
    				}
                }

    			if(isample < v3_histos_var_sample_syst[ivar].size() - 1)
    			{
    				if(v_groups[isample] == v_groups[isample+1])
    				{
    					v3_histos_var_sample_syst[ivar][isample+1][isyst]->Add(v3_histos_var_sample_syst[ivar][isample][isyst]); //Merge with next sample
    					continue; //Will draw merged histo, not this single one
    				}
    			}

                //-- CHANGED
                // v3_histos_var_sample_syst[ivar][isample][isyst]->GetXaxis()->SetTitle(Get_Variable_Name(total_var_list[ivar]));
                v3_histos_var_sample_syst[ivar][isample][isyst]->GetXaxis()->SetLabelSize(0.);

                if(normalize) {v3_histos_var_sample_syst[ivar][isample][isyst]->GetYaxis()->SetTitle("Normalized");}
                else {v3_histos_var_sample_syst[ivar][isample][isyst]->GetYaxis()->SetTitle("Events");}

                //-- CHANGED
                // if(normalize) {v3_histos_var_sample_syst[ivar][isample][isyst]->SetMaximum(0.5);}
                v3_histos_var_sample_syst[ivar][isample][isyst]->SetMaximum(v3_histos_var_sample_syst[ivar][isample][isyst]->GetMaximum()*1.5);
                v3_histos_var_sample_syst[ivar][isample][isyst]->SetMinimum(0.00001);

    			v3_histos_var_sample_syst[ivar][isample][isyst]->Draw("hist E same");

    			if(v_syst[isyst] == "")
    			{
                    if(v_groups[isample] == "tZq") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "tZq", "L");}
                    else if(v_groups[isample] == "ttZ" ) {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "t#bar{t}Z", "L");}
                    else if(v_groups[isample] == "ttW") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "t#bar{t}X", "L");}
                    else if(v_groups[isample] == "tHq") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "tX", "L");}
        			else if(v_groups[isample] == "WZ") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "VV(V)", "L");}
        			else if(v_groups[isample] == "Fakes") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "Non-prompt", "L");}
        			else if(v_groups[isample] == "QFlip") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "Flip", "L");}
                    else if(v_groups[isample] == "GammaConv") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "#gamma-conv.", "L");}
                    else if(v_groups[isample] == "TTGamma_Dilep") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "X+#gamma", "L");}
                    else if(v_groups[isample] == "DY" ) {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "V+jets", "L");}
                    else if(v_groups[isample] == "TTbar_DiLep" ) {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "t#bar{t}", "L");}
                    else {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], v_groups[isample], "L");}
    			}

    			//HARDCODED
    			if(v_syst[isyst] == "JESUp") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "JES Up", "L");}
    			if(v_syst[isyst] == "JESDown") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "JES Down", "L");}
    		} //syst loop
    	} //sample loop

        qw->Draw("same");


// #####    ##   ##### #  ####
// #    #  #  #    #   # #    #
// #    # #    #   #   # #    #
// #####  ######   #   # #    #
// #   #  #    #   #   # #    #
// #    # #    #   #   #  ####

// #####  #       ####  #####
// #    # #      #    #   #
// #    # #      #    #   #
// #####  #      #    #   #
// #      #      #    #   #
// #      ######  ####    #

		//-- create subpad to plot ratio
		TPad *pad_ratio = new TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 1.0);
		pad_ratio->SetTopMargin(0.75);
		pad_ratio->SetFillColor(0);
		pad_ratio->SetFillStyle(0);
		pad_ratio->SetGridy(1);
		pad_ratio->Draw();
		pad_ratio->cd(0);

        if(v3_histos_var_sample_syst[ivar].size()<2) {return;} //Can't make ratio
        TH1F* histo_ratio = (TH1F*) v3_histos_var_sample_syst[ivar][1][0]->Clone();
        TH1F* histo_ratio_denominator = (TH1F*) v3_histos_var_sample_syst[ivar][0][0]->Clone();

        histo_ratio->Divide(histo_ratio_denominator);

		histo_ratio->GetYaxis()->SetTitle("Ratio");
		histo_ratio->GetYaxis()->SetTickLength(0.);
        histo_ratio->GetYaxis()->SetTitleOffset(1.15);
        histo_ratio->GetYaxis()->SetLabelSize(0.04);
        histo_ratio->GetXaxis()->SetLabelSize(0.04);
		histo_ratio->GetXaxis()->SetLabelFont(42);
		histo_ratio->GetYaxis()->SetLabelFont(42);
		histo_ratio->GetXaxis()->SetTitleFont(42);
		histo_ratio->GetYaxis()->SetTitleFont(42);
        histo_ratio->GetYaxis()->SetNdivisions(507); //grid drawn on primary tick marks only
		histo_ratio->GetXaxis()->SetNdivisions(505);
		histo_ratio->GetYaxis()->SetTitleSize(0.06);
		histo_ratio->GetXaxis()->SetTickLength(0.04);

		//NB : when using SetMaximum(), points above threshold are simply not drawn
		//So for ratio plot, even if point is above max but error compatible with 1, point/error bar not represented!
        histo_ratio->SetMinimum(0.6);
        histo_ratio->SetMaximum(1.4);

		if(drawInputVars) {histo_ratio->GetXaxis()->SetTitle(Get_Variable_Name(total_var_list[ivar]));}
		else
		{
            histo_ratio->GetXaxis()->SetTitle(total_var_list[ivar]);

            //Hardcode NN output nodes names...?
            if(total_var_list[ivar] == "NN") {histo_ratio->GetXaxis()->SetTitle("NN output");}
            else if(total_var_list[ivar] == "NN0") {histo_ratio->GetXaxis()->SetTitle("NN (tZq node)");}
            else if(total_var_list[ivar] == "NN1" && NN_nNodes == 3) {histo_ratio->GetXaxis()->SetTitle("NN (ttZ node)");}
            else if(total_var_list[ivar] == "NN2" && NN_nNodes == 3) {histo_ratio->GetXaxis()->SetTitle("NN (Bkgs node)");}

			if(template_name == "categ") //Vertical text X labels (categories names)
			{
                //Hard-coded -- should automate labels
                {
                    const char *labels[10]  = {"1bj,2j","1bj,3j","1bj,4j","1bj,5j","1bj,6j","2bj,2j","2bj,3j","2bj,4j","2bj,5j","2bj,6j"};
                    for(int i=1;i<=10;i++) {histo_ratio->GetXaxis()->SetBinLabel(i,labels[i-1]);}
                }

                histo_ratio->GetXaxis()->SetTitle("");
                // histo_ratio->GetXaxis()->SetTitle("Categ.");
				histo_ratio->GetXaxis()->SetLabelSize(0.06);
				histo_ratio->GetXaxis()->SetLabelOffset(0.02);
                histo_ratio->LabelsOption("v", "X"); //X labels vertical
			}
		}

		pad_ratio->cd(0);
        histo_ratio->Draw("hist E same");

        histo_ratio_denominator->Divide(histo_ratio_denominator);
        histo_ratio_denominator->Draw("E same");
        // for(int ibin=1; ibin<histo_ratio_denominator->GetNbinsX()+1; ibin++)
        // {
        //     cout<<"bin "<<ibin<<" / "<<histo_ratio_denominator->GetBinContent(ibin)<<" / "<<histo_ratio_denominator->GetBinError(ibin)<<" / "<<histo_ratio_denominator->GetBinError(ibin)/histo_ratio_denominator->GetBinContent(ibin)<<endl;
        // }


//  ####   ####   ####  #    # ###### ##### #  ####   ####
// #    # #    # #      ##  ## #        #   # #    # #
// #      #    #  ####  # ## # #####    #   # #       ####
// #      #    #      # #    # #        #   # #           #
// #    # #    # #    # #    # #        #   # #    # #    #
//  ####   ####   ####  #    # ######   #   #  ####   ####

//----------------
// CAPTIONS //
//----------------
// -- using https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines

    	float l = c->GetLeftMargin();
    	float t = c->GetTopMargin();

    	TString cmsText = "CMS";
    	TLatex latex;
    	latex.SetNDC();
    	latex.SetTextAngle(0);
    	latex.SetTextColor(kBlack);
    	latex.SetTextFont(61);
    	latex.SetTextAlign(11);
    	latex.SetTextSize(0.06);
    	// latex.DrawLatex(l + 0.01, 0.92, cmsText);

    	TString extraText = "Preliminary";
    	latex.SetTextFont(52);
    	latex.SetTextSize(0.05);
    	// if(draw_preliminary_label)
    	{
    		// latex.DrawLatex(l + 0.12, 0.92, extraText);
    	}

    	float lumi = lumiValue;
    	TString lumi_ts = Convert_Number_To_TString(lumi);
    	lumi_ts += " fb^{-1} (13 TeV)";

    	latex.SetTextFont(42);
    	latex.SetTextAlign(31);
    	latex.SetTextSize(0.04);
    	// latex.DrawLatex(0.96, 0.92,lumi_ts);

    	//------------------
    	//-- channel info
    	TLatex text2 ;
    	text2.SetNDC();
    	text2.SetTextAlign(13);
    	text2.SetTextSize(0.045);
    	text2.SetTextFont(42);

        TString info_data = "l^{#pm}l^{#pm}l^{#pm}";
        if (channel=="eee")    info_data = "eee";
        else if (channel=="eeu")  info_data = "ee#mu";
        else if (channel=="uue")  info_data = "#mu#mu e";
        else if (channel=="uuu") info_data = "#mu#mu #mu";

    	// if(h_sum_data->GetBinContent(h_sum_data->GetNbinsX() ) > h_sum_data->GetBinContent(1) ) {text2.DrawLatex(0.55,0.87,info_data);}
    	// else {text2.DrawLatex(0.20,0.87,info_data);}
    	text2.DrawLatex(0.23,0.86,info_data);


// #    # #####  # ##### ######     ####  #    # ##### #####  #    # #####
// #    # #    # #   #   #         #    # #    #   #   #    # #    #   #
// #    # #    # #   #   #####     #    # #    #   #   #    # #    #   #
// # ## # #####  #   #   #         #    # #    #   #   #####  #    #   #
// ##  ## #   #  #   #   #         #    # #    #   #   #      #    #   #
// #    # #    # #   #   ######     ####   ####    #   #       ####    #

        TString outdir = "plots/templates_shapes/";
        mkdir(outdir.Data(), 0777);
        if(type != "") {outdir+=  type + "/"; mkdir(outdir.Data(), 0777);}

    	//Output
    	TString output_plot_name = outdir + total_var_list[ivar] + "_" + region +"_templatesShapes";
    	if(channel != "") {output_plot_name+= "_" + channel;}
    	output_plot_name+= this->filename_suffix + this->plot_extension;

    	c->SaveAs(output_plot_name);

        delete c; c = NULL;
        delete qw; qw = NULL;
        delete histo_ratio; histo_ratio = NULL;
        delete histo_ratio_denominator; histo_ratio_denominator = NULL;
    } //var loop

    for(int ivar=0; ivar<v3_histos_var_sample_syst.size(); ivar++)
    {
    	for(int isample=0; isample<v3_histos_var_sample_syst[ivar].size(); isample++)
    	{
    		for(int isyst=0; isyst<v_syst.size(); isyst++) {delete v3_histos_var_sample_syst[ivar][isample][isyst];}
        }
    }

	return;
}

















//--------------------------------------------
//  ######  ######## ########     ######  ##    ##  ######  ########
// ##    ## ##          ##       ##    ##  ##  ##  ##    ##    ##
// ##       ##          ##       ##         ####   ##          ##
//  ######  ######      ##        ######     ##     ######     ##
//       ## ##          ##             ##    ##          ##    ##
// ##    ## ##          ##       ##    ##    ##    ##    ##    ##    ###
//  ######  ########    ##        ######     ##     ######     ##    ###

//    ###    ########  ########  ########  ########  ######   ######  ########  ######
//   ## ##   ##     ## ##     ## ##     ## ##       ##    ## ##    ## ##       ##    ##
//  ##   ##  ##     ## ##     ## ##     ## ##       ##       ##       ##       ##
// ##     ## ##     ## ##     ## ########  ######    ######   ######  ######    ######
// ######### ##     ## ##     ## ##   ##   ##             ##       ## ##             ##
// ##     ## ##     ## ##     ## ##    ##  ##       ##    ## ##    ## ##       ##    ##
// ##     ## ########  ########  ##     ## ########  ######   ######  ########  ######
//--------------------------------------------

//Problem : in Potato code, systematics variations weights are encoded into arrays. Hence, can not directly set branch addresses via SetBranchAddress('name', &var)
//==> My solution : for each of these arrays, hard-code 1 class member array with corresponding size. This member array 'reads' the array stored in the ntuple.
//==> To make things simpler, I then use a vector<double*> so that each element will effectively read the value of a given systematic variation. From there, can proceed as usual...
//Use this function to hard-code which vector element corresponds to which systematics, and set the address
//NB : some inconsistent indices : down variation may be element 0 or 1... hence, correspondances must be hard-coded !
void TopEFT_analysis::SetBranchAddress_SystVariationArray(TTree* t, TString systname, vector<Double_t*> &v_doubles, int isyst)
{
    TString array_name = ""; //Name of the systematic array as stored in the ntuple
    double* address_memberArray = NULL; //Address to the class member which will hold the array values (hardcoded)
    int index=-1; //index of the specific array member which holds the value of the syst of interest (hardcoded)

    if(systname == "") {return;}
    else if(systname.BeginsWith("PU"))
    {
        address_memberArray = array_PU;
        array_name = "varWeightPU";
        if(systname.EndsWith("Down")) {index = 0;}
        else if(systname.EndsWith("Up")) {index = 1;}
    }
    else if(systname.BeginsWith("prefiringWeight"))
    {
        address_memberArray = array_prefiringWeight;
        array_name = "varWeightPrefire";
        if(systname.EndsWith("Down")) {index = 0;}
        else if(systname.EndsWith("Up")) {index = 1;}
    }
    else if(systname.BeginsWith("Btag"))
    {
        address_memberArray = array_Btag;
        array_name = "btagEventWeightVar";

        if(systname.EndsWith("HFUp")) {index = 0;}
        else if(systname.EndsWith("HFDown")) {index = 1;}
        else if(systname.EndsWith("LFUp")) {index = 2;}
        else if(systname.EndsWith("LFDown")) {index = 3;}
        else if(systname.EndsWith("HFstats1Up")) {index = 4;}
        else if(systname.EndsWith("HFstats1Down")) {index = 5;}
        else if(systname.EndsWith("HFstats2Up")) {index = 6;}
        else if(systname.EndsWith("HFstats2Down")) {index = 7;}
        else if(systname.EndsWith("LFstats1Up")) {index = 8;}
        else if(systname.EndsWith("LFstats1Down")) {index = 9;}
        else if(systname.EndsWith("LFstats2Up")) {index = 10;}
        else if(systname.EndsWith("LFstats2Down")) {index = 11;}
        else if(systname.EndsWith("CFerr1Down")) {index = 12;}
        else if(systname.EndsWith("CFerr1Up")) {index = 13;}
        else if(systname.EndsWith("CFerr2Down")) {index = 14;}
        else if(systname.EndsWith("CFerr2Up")) {index = 15;}
    }
    else if(systname.BeginsWith("LepEff_mu"))
    {
        if(systname.Contains("Loose")) {address_memberArray = array_LepEffLoose_mu; array_name = "varWeightMuonLoose";}
        else {address_memberArray = array_LepEffTight_mu; array_name = "varWeightMuonTight";}
        if(systname.EndsWith("Down")) {index = 0;}
        else if(systname.EndsWith("Up")) {index = 1;}
    }
    else if(systname.BeginsWith("LepEff_el"))
    {
        if(systname.Contains("Loose")) {address_memberArray = array_LepEffLoose_mu; array_name = "varWeightElectronLoose";}
        else {address_memberArray = array_LepEffTight_mu; array_name = "varWeightElectronTight";}
        if(systname.EndsWith("Down")) {index = 0;}
        else if(systname.EndsWith("Up")) {index = 1;}
    }
    /*
    else if(systname.BeginsWith("LepEff_mu"))
    {
        address_memberArray = array_LepEff_mu;
        array_name = "varWeightMuon";
        if(systname.EndsWith("LooseDown")) {index = 0;}
        else if(systname.EndsWith("LooseUp")) {index = 1;}
        else if(systname.EndsWith("TightDown")) {index = 2;}
        else if(systname.EndsWith("TightUp")) {index = 3;}
    }
    else if(systname.BeginsWith("LepEff_el"))
    {
        address_memberArray = array_LepEff_el;
        array_name = "varWeightElectron";
        if(systname.EndsWith("LooseDown")) {index = 0;}
        else if(systname.EndsWith("LooseUp")) {index = 1;}
        else if(systname.EndsWith("TightDown")) {index = 2;}
        else if(systname.EndsWith("TightUp")) {index = 3;}
    }
    */
    else{cout<<FRED("ERROR ! Systematic '"<<systname<<"' not included in function SetBranchAddress_SystVariation() from Helper.cxx ! Can *not* compute it !")<<endl; return;}

    t->SetBranchStatus(array_name, 1);
    t->SetBranchAddress(array_name, address_memberArray);
    v_doubles[isyst] = &address_memberArray[index];

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

// ######## ######## ##     ## ########  ##          ###    ######## ########  ######
//    ##    ##       ###   ### ##     ## ##         ## ##      ##    ##       ##    ##
//    ##    ##       #### #### ##     ## ##        ##   ##     ##    ##       ##
//    ##    ######   ## ### ## ########  ##       ##     ##    ##    ######    ######
//    ##    ##       ##     ## ##        ##       #########    ##    ##             ##
//    ##    ##       ##     ## ##        ##       ##     ##    ##    ##       ##    ##
//    ##    ######## ##     ## ##        ######## ##     ##    ##    ########  ######
//--------------------------------------------

/**
 * The main code is producing/plotting template histograms for each sample separately
 * But we may want to group some processes together for Combine fit (e.g. "Rares", "EWK", ...)
 * ===> In addition to storing individual histos, also merge+store histos for the relevant process groups
 * NB: here the order of loops is important because we sum histograms recursively, and the 'sample_list' loop must be the most nested one !
 */
void TopEFT_analysis::MergeSplit_Templates(TString filename, vector<TString> total_var_list, TString template_name, TString category, bool force_normTemplate_positive)
{
	cout<<endl<<FYEL("==> Merging some templates in file : ")<<filename<<endl;

	if(!Check_File_Existence(filename) ) {cout<<endl<<FRED("File "<<filename<<" not found! Abort template merging !")<<endl; return;}
	TFile* f = TFile::Open(filename, "UPDATE");

	//NB :here the order of loops is important because we sum histograms recursively ! The 'sample_list' loop *must be the most nested one* !
    for(int ivar=0; ivar<total_var_list.size(); ivar++) //There may be more than 1 template (e.g. several NN output nodes)
    {
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
            for(int ichan=0; ichan<channel_list.size(); ichan++)
        	{
                // cout<<"channel_list[ichan] "<<channel_list[ichan]<<endl;

        		for(int itree=0; itree<systTree_list.size(); itree++)
        		{
        			// if(systTree_list[itree] != "" && channel_list.size() > 1 && channel_list[ichan] == "") {continue;}

        			for(int isyst=0; isyst<syst_list.size(); isyst++)
        			{
        				// if(((channel_list.size() > 1 && channel_list[ichan] == "") || systTree_list[itree] != "") && syst_list[isyst] != "") {continue;}
        				if(systTree_list[itree] != "" && syst_list[isyst] != "") {break;}

        				TH1F* h_merging = NULL; //Tmp merged histogram

                        //-- For EFT strategies, need to store all histogram bins separately
                        int n_singleBins = 0; //Default: will only merge: a) entire histograms, b) histograms stored as single bin //Gets updated below depending on strategy
                        for(int ibin=-1; ibin<n_singleBins+1; ibin++) //Convention (depending on NN_strategy): bin==-1 <-> merge full histogram; bin==0 <-> merge histo stored as single bin (counting exp.); bin>0 <-> merge corresponding template bin (individually)
                        {
                            // cout<<"ibin "<<ibin<<" ("<<"n_singleBins "<<n_singleBins<<")"<<endl;
            				for(int isample=0; isample<sample_list.size(); isample++)
            				{
            					//-- Protections : not all syst weights apply to all samples, etc.
                                if(sample_list[isample] == "DATA" && (systTree_list[itree] != "" || syst_list[isyst] != "")) {continue;} //nominal data only
                                if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefiring") ) {continue;} //no prefire in 2018

            					// cout<<endl<<"Syst "<<syst_list[isyst]<<systTree_list[itree]<<" / chan "<<channel_list[ichan]<<" / sample "<<sample_list[isample]<<endl;

            					//Check if this sample needs to be merged, i.e. if the samples before/after belong to the same "group of samples"
            					bool merge_this_sample = false;
            					if(!isample && sample_groups.size() > 1 && sample_groups[isample+1] == sample_groups[isample]) {merge_this_sample = true;}
            					else if(isample == sample_list.size()-1 && sample_groups[isample-1] == sample_groups[isample]) {merge_this_sample = true;}
            					else if(isample > 0 && isample < sample_list.size()-1 && (sample_groups[isample+1] == sample_groups[isample] || sample_groups[isample-1] == sample_groups[isample])) {merge_this_sample = true;}

            					// cout<<"merge_this_sample "<<merge_this_sample<<endl;
            					if(!merge_this_sample) {continue;} //Only care about samples to merge : others are already stored in file

            					TString samplename = sample_list[isample];
            					if(samplename == "DATA") {samplename = "data_obs";}

                                TString histoname = ""; //Default: merge full histograms
                                if(ibin==0) {histoname+= "countExp_";} //Merge full histos stored as single bins (counting exp.)
                                else if(n_singleBins > 0 && ibin>0) {histoname+= (TString) "bin"+Form("%d",ibin)+"_";} //Merge bin per bin
                                histoname+= total_var_list[ivar];
                                // if(category != "") {histoname+= "_" + category;}
            					if(channel_list[ichan] != "") {histoname+= "_" + channel_list[ichan];}
                                histoname+= "_" + v_lumiYears[iyear] + "__" + samplename;
            					if(syst_list[isyst] != "") {histoname+= "__" + Get_Modified_SystName(syst_list[isyst], v_lumiYears[iyear]);}
            					else if(systTree_list[itree] != "") {histoname+= "__" + systTree_list[itree];}
            					// cout<<"histoname "<<histoname<<endl;

                                //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                                histoname.ReplaceAll('-', 'm');

            					if(!f->GetListOfKeys()->Contains(histoname) && systTree_list[itree] == "" && syst_list[isyst] == "")
            					{
            						cout<<FRED("Histo "<<histoname<<" not found in file "<<filename<<" !")<<endl;
            					 	continue;
            					}

            					TH1F* h_tmp = (TH1F*) f->Get(histoname); //Get individual histograms
            					// cout<<"h_tmp->Integral() = "<<h_tmp->Integral()<<endl;

                                if(make_SMvsEFT_templates_plots) //For EFT templates, will also consider inidividual-bin histos
                                {
                                    if(n_singleBins==0 && ibin==-1) //Only read the binning once, from full histograms
                                    {
                                        n_singleBins = h_tmp->GetNbinsX(); //Update the limit for the for-loop *within the loop* (--> first read full histo to infer the correct binning)
                                    }
                                }

            					int factor = +1; //Addition
            					// if(sample_list[isample] == "Fakes_MC") {factor = -1;} //Substraction of 'MC Fakes' (prompt contribution to fakes)

            					if(h_tmp != 0)
            					{
            						if(!h_merging) {h_merging = (TH1F*) h_tmp->Clone();}
            						else {h_merging->Add(h_tmp, factor);}
            					}
            					else {cout<<"h_tmp null !"<<endl;}

            					// cout<<"h_merging->Integral() = "<<h_merging->Integral()<<endl;

            					delete h_tmp; h_tmp = NULL;
            					if(!h_merging) {cout<<"Syst "<<syst_list[isyst]<<systTree_list[itree]<<" / chan "<<channel_list[ichan]<<" / sample "<<sample_list[isample]<<endl; cout<<"h_merging is null ! Fix this first"<<endl; return;}

            					//Check if next sample will be merged with this one, or else if must write the histogram
            					if(isample < sample_list.size()-1 && sample_groups[isample+1] == sample_groups[isample]) {continue;}
            					else
            					{
                                    TString histoname_new = "";
                                    if(ibin==0) {histoname_new+= "countExp_";} //Merge full histos stored as single bins (counting exp.)
                                    else if(n_singleBins > 0 && ibin>0) {histoname_new+= (TString) "bin"+Form("%d",ibin)+"_";} //Merge bin per bin
                                    histoname_new+= total_var_list[ivar];
                                    // if(category != "") {histoname_new+= "_" + category;}
                                    if(channel_list[ichan] != "") {histoname_new+="_"  + channel_list[ichan];}
                                    histoname_new+= "_" + v_lumiYears[iyear] + "__" + sample_groups[isample];
                                    if(syst_list[isyst] != "") {histoname_new+= "__" + Get_Modified_SystName(syst_list[isyst], v_lumiYears[iyear]);}
                                    else if(systTree_list[itree] != "") {histoname_new+= "__" + systTree_list[itree];}

                                    //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                                    histoname_new.ReplaceAll('-', 'm');


            						if(force_normTemplate_positive)
            						{
            							//If integral of histo is negative, set to 0 (else COMBINE crashes) -- must mean that norm is close to 0 anyway
            							if(h_merging->Integral() <= 0)
            							{
            								// cout<<endl<<"While merging processes by groups ('Rares'/...) :"<<endl<<FRED(" h_merging->Integral() = "<<h_merging->Integral()<<" (<= 0) ! Distribution set to ~>0 (flat), to avoid crashes in COMBINE !")<<endl;
            								Set_Histogram_FlatZero(h_merging, histoname_new, false);
            								cout<<"(Syst "<<syst_list[isyst]<<systTree_list[itree]<<" / chan "<<channel_list[ichan]<<" / sample "<<sample_list[isample]<<")"<<endl;
            							}
            						}
            						// cout<<"h_merging->Integral() = "<<h_merging->Integral()<<endl;

            						f->cd();
            						h_merging->Write(histoname_new, TObject::kOverwrite);

            						// cout<<"-- Writing merged histo "<<histoname_new<<" with integral "<<h_merging->Integral()<<endl;

            						delete h_merging; h_merging = NULL;
            					} //write histo
                            } //sample loop
                        } //bin loop

        			} //syst loop
        		} //tree loop
        	} //channel loop
        } //years loop
    } //vars loop

    cout<<DIM("(Total: "<<Count_nofHistos_inTFile(f)<<" histograms)")<<endl;

	f->Close();

    // cout<<endl<<FYEL("Updated file: " )<<filename<<endl<<endl<<endl;
    cout<<FYEL("... Done")<<endl<<endl<<endl;

	return;
}




















//--------------------------------------------
// ##     ## ##     ##    ###        ######  ##     ## ########
// ###   ### ##     ##   ## ##      ##    ## ##     ##    ##
// #### #### ##     ##  ##   ##     ##       ##     ##    ##
// ## ### ## ##     ## ##     ##    ##       ##     ##    ##
// ##     ##  ##   ##  #########    ##       ##     ##    ##
// ##     ##   ## ##   ##     ##    ##    ## ##     ##    ##
// ##     ##    ###    ##     ##     ######   #######     ##
//--------------------------------------------

/**
 * Read a MVA file, loop over the input file (path given in arg.), and stores in a vector (passed in arg.) either 1) if the event passes the required MVA cut, or 2) which of the multiclass MVA node has the max value
 * Intended use: call this function to fill a per-sample vector already containing the information on a MVA cut for all events
 */
void TopEFT_analysis::Get_VectorAllEvents_passMVACut(vector<int>& v, TString signal, TString classifier_name, TString tree_name, TString input_file_path, TString year, float cut_value, bool keep_aboveCut, bool use_specificMVA_eachYear, int categorization_strategy, bool MVA_EFT, int nentries_max, TString event_cat, bool also_applyCut_onMaxNodeValue)
{
    cout<<endl<<FYEL("=== Filling vector with specific MVA information (pass/fail cut or max. node)... ")<<endl;

    if(classifier_name != "BDT" && classifier_name != "NN") {cout<<BOLD(FRED("ERROR: wrong [classifier_name] option !"))<<endl; return;}

    int nevents_passingCut = 0;
    TString BDT_method_name = "BDT";
    vector<TString> var_list_tmp;
    vector<float> var_floats_tmp;

    TString MVA_input_path = Get_MVAFile_InputPath(classifier_name, signal, year, use_specificMVA_eachYear, MVA_EFT, false, categorization_strategy);
    if(MVA_input_path == "") {return;} //MVA file not found

    //Use local variables to avoid conflict with class members (different MVAs)
    vector<TString> var_list_NN; TString NN_inputLayerName = ""; TString NN_outputLayerName = ""; int NN_iMaxNode = -1; int NN_nNodes = -1;

    if(classifier_name == "BDT") //BDT
    {
        var_list_tmp = var_list;

        reader_tmp = new TMVA::Reader("!Color:!Silent");
        for(int i=0; i<var_list_tmp.size(); i++) {reader_tmp->AddVariable(var_list_tmp[i].Data(), &var_floats_tmp[i]);}
        for(int i=0; i<v_cut_name.size(); i++)
        {
            if(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("==")) {reader_tmp->AddVariable(v_cut_name[i].Data(), &v_cut_float[i]);}
        }

        // BDT_method_name = "BDT_"+signal+"_"+year + " method";
        reader_tmp->BookMVA(BDT_method_name, MVA_input_path);
    }
    else //NN
    {
        //-- Get path of NN info text file --> Read list of input variables (and more) //NB: only use local variables to avoid conflicts
        TString NNinfo_input_path = Get_MVAFile_InputPath(classifier_name, signal, year, use_specificMVA_eachYear, MVA_EFT, true, categorization_strategy);
        if(NNinfo_input_path == "") {return;} //MVA file not found

        if(Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN, v_NN_nodeLabels, NN_inputLayerName, NN_outputLayerName, NN_iMaxNode, NN_nNodes)) {clfy_tmp = new TFModel(MVA_input_path.Data(), var_list_NN.size(), NN_inputLayerName.Data(), NN_nNodes, NN_outputLayerName.Data());} //Load neural network model

        var_list_tmp = var_list_NN;
        var_floats_tmp.resize(var_list_tmp.size());
    }

	if(!Check_File_Existence(input_file_path) ) {cout<<BOLD(FRED("ERROR: "<<input_file_path<<" not found!"))<<endl; return;}
	TFile* file_input = TFile::Open(input_file_path, "READ");
	TTree* tree = (TTree*) file_input->Get(tree_name);

    //-- Set addresses of input features
	for(int ivar=0; ivar<var_list_tmp.size(); ivar++)
    {
        if(var_list_tmp[ivar] == "ctz" || var_list_tmp[ivar] == "ctw" || var_list_tmp[ivar] == "cpq3" || var_list_tmp[ivar] == "cpqm" || var_list_tmp[ivar] == "cpt") {continue;} //WC input values are arbitrary, there is no address to set !
        tree->SetBranchAddress(var_list_tmp[ivar], &var_floats_tmp[ivar]);
    }

    //-- May cut on an 'event category flag' whose name is given as argument (<-> no need to evaluate MVA for events which do not enter the region of interest)
    Char_t is_goodCategory = true;
    if(event_cat != "") {tree->SetBranchAddress(event_cat, &is_goodCategory);}

	int nentries = tree->GetEntries();
    if(nentries_max > 0 && nentries > nentries_max) {nentries = nentries_max;}
    v.clear(); v.resize(nentries); std::fill(v.begin(),v.end(),0); //Fill with '0' for all tree entries

	for(int ientry=0; ientry<nentries; ientry++)
	{
		if(ientry && ientry%50000==0) {cout<<DIM(" --- "<<ientry<<" / "<<nentries<<"")<<endl;}
		tree->GetEntry(ientry);

        if(!is_goodCategory) {continue;}

        float mva_output = 0.;

        if(classifier_name == "BDT") {mva_output = reader_tmp->EvaluateMVA(BDT_method_name);}
        else //NN
        {
            // for(int ivar=0; ivar<var_floats_tmp.size(); ivar++) {cout<<"var_floats_tmp[ivar] "<<var_floats_tmp[ivar]<<endl;} //Debug

            std::vector<float> clfy_outputs = clfy_tmp->evaluate(var_floats_tmp); //Evaluate output node(s) value(s)
            NN_iMaxNode = -1;
            for(int inode=0; inode<NN_nNodes; inode++)
            {
                if(clfy_outputs[inode] > mva_output) {mva_output = clfy_outputs[inode]; NN_iMaxNode = inode;}
                // cout<<"clfy_outputs[inode] "<<clfy_outputs[inode]<<" / mva_output "<<mva_output<<" / NN_iMaxNode "<<NN_iMaxNode<<endl;
            }

            if(NN_nNodes > 1) //Multiclass --> Store max. node information
            {
                if(!also_applyCut_onMaxNodeValue || ((keep_aboveCut && mva_output >= cut_value) || (!keep_aboveCut && mva_output < cut_value))) {v[ientry] = NN_iMaxNode;}
                continue;
            }
        }

        //Binary classifier --> Determine/store whether the event passes the MVA cut or not
        bool pass_cut = false;
        if((keep_aboveCut && mva_output >= cut_value) || (!keep_aboveCut && mva_output < cut_value)) {pass_cut = true; nevents_passingCut++;}
        v[ientry] = pass_cut;

	} //loop on entries

	if(v.size() != nentries) {cout<<BOLD(FRED("Wrong number of entries in BDT cut vector ! Check it please !"))<<endl; return;}
	file_input->Close();

	// cout<<FMAG("---- Vector containing BDTfakeSR cut results is filled !")<<endl;
	if(NN_nNodes <= 1) {cout<<DIM("(Number of events passing the MVA cut : "<<nevents_passingCut<<")")<<endl;} //Meaningless for multiclass NN

    if(reader_tmp) {delete reader_tmp; reader_tmp = NULL;}
    if(clfy_tmp) {delete clfy_tmp; clfy_tmp = NULL;}

    cout<<FYEL("... Done ! ===")<<endl<<endl;

	return;
}























//--------------------------------------------
// ######## ########  ######  ######## #### ##    ##  ######
//    ##    ##       ##    ##    ##     ##  ###   ## ##    ##
//    ##    ##       ##          ##     ##  ####  ## ##
//    ##    ######    ######     ##     ##  ## ## ## ##   ####
//    ##    ##             ##    ##     ##  ##  #### ##    ##
//    ##    ##       ##    ##    ##     ##  ##   ### ##    ##
//    ##    ########  ######     ##    #### ##    ##  ######
//--------------------------------------------
