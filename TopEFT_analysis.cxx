//-- by Nicolas Tonon (DESY) --//

//--- LIST OF FUNCTIONS (for quick search) :
//--------------------------------------------
// Train_BDT

// Produce_Templates

// Draw_Templates

// Compare_TemplateShapes_Processes

// MergeSplit_Templates
// SetBranchAddress_SystVariationArray
// Get_VectorAllEvent

// Make_PaperPlot_CommonRegions
// Make_PaperPlot_SignalRegions
// Make_PaperPlot_ControlPlots
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
TopEFT_analysis::TopEFT_analysis(vector<TString> thesamplelist, vector<TString> thesamplegroups, vector<TString> thesystlist, vector<TString> thesystTreelist, vector<TString> thechannellist, vector<TString> thevarlist, vector<TString> set_v_cut_name, vector<TString> set_v_cut_def, vector<bool> set_v_cut_IsUsedForBDT, vector<TString> set_v_add_var_names, TString theplotextension, vector<TString> set_lumi_years, bool show_pulls, TString region, TString signal_process, TString classifier_name, bool scanOperators_paramNN, TString operator_scan1, TString operator_scan2, vector<float> v_WCs_operator_scan1, vector<float> v_WCs_operator_scan2, bool make_SMvsEFT_templates_plots, bool is_blind, int categorization_strategy, bool use_specificMVA_eachYear, TString nominal_tree_name, bool use_DD_NPL, bool make_fixedRegions_templates, bool process_samples_byGroup, bool split_EFTtemplates_perBin, bool use_paperStyle, bool use_NN_SRother, bool use_NN_cpq3_SRttZ)
{
    //Canvas definition
    Load_Canvas_Style();

    TH1::SetDefaultSumw2();
    gStyle->SetErrorX(0.);

    nbins = 10; //default

	mkdir("outputs", 0777);
	mkdir("plots", 0777);

    if(nominal_tree_name != "") {this->nominal_tree_name = nominal_tree_name;}

	stop_program = false;

	this->region = region;

    this->split_EFTtemplates_perBin = split_EFTtemplates_perBin;

    this->use_paperStyle = use_paperStyle;
    this->use_NN_SRother = use_NN_SRother;
    this->use_NN_cpq3_SRttZ = use_NN_cpq3_SRttZ;

    if(region == "" && !make_fixedRegions_templates)
    {
        this->region = "signal";
        // this->make_SMvsEFT_templates_plots = true;
        cout<<FBLU("NB: you did not select a category; setting it to 'is_signal_SR' by default..."<<"")<<endl;
        usleep(1000000); //Pause for 1s (in microsec)
    }

    this->signal_process = signal_process;
    if(signal_process != "tZq" && signal_process != "ttZ" && signal_process != "tWZ") {cout<<BOLD(FRED("ERROR ! [signal_process] option not recognized ! Abort..."))<<endl; stop_program = true;}

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
    if(region=="twz") {nSampleGroups++;} //In tWZ region, single out tWZ process from tX

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

//-------------

    dir_ntuples = NTUPLEDIR; //Defined in Utils/Helper
    TString dir_ntuples_tmp = dir_ntuples;

    this->use_optimized_ntuples = false; //Obsolete
    /* //Obsolete
    dir_ntuples_SR = dir_ntuples + "SR/";
    if(region=="tZq" || region=="ttZ" || region=="signal") {dir_ntuples = dir_ntuples_SR; cout<<endl<<FMAG("NB: will use SR sub-samples ! Make sure they are up-to-date !")<<endl<<endl;}
    if(dir_ntuples == dir_ntuples_SR && !Check_File_Existence((dir_ntuples_SR + v_lumiYears[0] + "/DATA.root")))
    {
        cout<<endl<<FMAG("Warning : SR sub-sample (produced with code input_ntuples/Split_FullSamples) "<<(dir_ntuples_SR + v_lumiYears[0] + "/DATA.root")<<" not found ! Assuming that split ntuples are not available --> Will use full ntuples ! (slower)")<<endl<<endl;
        dir_ntuples = NTUPLEDIR;
    }
    // cout<<"dir_ntuples : "<<dir_ntuples<<endl;

    //-- To considerably speed-up the usual AN workflow, use [Split_AllNtuples_ByCategory] code to 1) split full ntuples by sub-categories, 2) merge all samples which belong to same sample groups
    //---> if found, will only use these 'optimized' ntuples <-> update the member lists accordingly
    if(dir_ntuples == dir_ntuples_SR && Check_File_Existence((dir_ntuples_SR + v_lumiYears[0] + "/NPL.root")) ) //Hard-coded check: if 'xx/SR/xx/NPL.root' sample is present, it must mean that we have split the ntuples by categ. and merged them by groups (here: NPL group)
    {
        this->use_optimized_ntuples = true;
        cout<<BOLD(FMAG("=== Optimized sub-ntuples found ! Will use ntuples split by categories & merged by groups (++ performance) ==="))<<endl;
        cout<<DIM("--> Set sample_list = sample_groups")<<endl<<endl;

        //Update member lists
        sample_list.clear();
        for(int igroup=0; igroup<sample_groups.size(); igroup++)
        {
            if(igroup > 0 && sample_groups[igroup]==sample_groups[igroup-1]) {continue;} //Group already considered --> Skip
            sample_list.push_back(sample_groups[igroup]); //Add sample group
        }
        sample_groups = sample_list;
    }
    */

    this->process_samples_byGroup = process_samples_byGroup;
    if(this->process_samples_byGroup) //Run on group ntuples to reduce the number of files to read
    {
        bool found_allGroupNtuples = true;
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
            for(int igroup=0; igroup<sample_groups.size(); igroup++)
            {
                if(igroup > 0 && sample_groups[igroup] == sample_groups[igroup-1]) {continue;}
                if(!Check_File_Existence(dir_ntuples + v_lumiYears[iyear] + "/" + sample_groups[igroup] + ".root")) //If some group ntuple is missing, will read all individual ntuples instead (default)
                {
                    cout<<FRED("WARNING: you have set [process_samples_byGroup=true] but I could not find group ntuple ["<<(dir_ntuples + v_lumiYears[iyear] + "/" + sample_groups[igroup] + ".root")<<"]... !")<<endl;
                    // usleep(100000); //Pause

                    // cout<<BOLD(FRED("WARNING: you have set [process_samples_byGroup=true] but I could not find group ntuple ["<<(dir_ntuples + v_lumiYears[iyear] + "/" + sample_groups[igroup] + ".root")<<"]... Setting [process_samples_byGroup=false] and processing individual ntuples instead !"))<<endl;
                    // found_allGroupNtuples = false;
                }
            }
        }
        if(found_allGroupNtuples)
        {
           sample_groups.erase( unique( sample_groups.begin(), sample_groups.end() ), sample_groups.end() ); //Remove all duplicates in 'sample_groups' (<-> only keep 1 occurence for each group)
           sample_list = sample_groups; //Samples will now refer to sample groups

           // for(int igroup=0; igroup<sample_groups.size(); igroup++) {cout<<"sample_groups[igroup] "<<sample_groups[igroup]<<endl;}
           // for(int isample=0; isample<sample_list.size(); isample++) {cout<<"sample_list[isample] "<<sample_list[isample]<<endl;}
        }
    }

	//-- Get colors
    int color_scheme = 0; //Check color scheme definitions directly in Get_Samples_Colors()
	color_list.resize(sample_list.size());
	Get_Samples_Colors(color_list, v_custom_colors, sample_list, sample_groups, color_scheme); //Read hard-coded sample colors

	this->classifier_name = classifier_name;
    if(classifier_name != "BDT" && classifier_name != "NN") {cout<<BOLD(FRED("Warning : classifier_name value ["<<classifier_name<<"] not supported !"))<<endl;}

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
		// var_list_floats.push_back(-999);
        var_list_pfloats.push_back(NULL);

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
    array_jetPileupID = NULL;
    array_fakeFactor = NULL;
    array_ME = NULL;
    array_alphaS = NULL;
    array_PDFtotal = NULL;
    array_partonShower = NULL;
    array_LepEffLoose_mu = NULL;
    array_LepEffTight_mu = NULL;
    array_LepEffLoose_el = NULL;
    array_LepEffTight_el = NULL;

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
            if(make_fixedRegions_templates) {cout<<BOLD(FRED("ERROR: options [scanOperators_paramNN=true] / [make_fixedRegions_templates=true] are incompatible ! Abort !"))<<endl; stop_program = true;}
        }
    }

    //NPL background
    this->use_DD_NPL = use_DD_NPL;
    if(!Check_File_Existence((dir_ntuples + v_lumiYears[0] + "/NPL.root")) && !Check_File_Existence((dir_ntuples + v_lumiYears[0] + "/NPL_DATA.root")))
    {
        cout<<endl<<BOLD(FMAG("ERROR : option [use_DD_NPL=true] but could not find corresponding file (either "<<dir_ntuples + v_lumiYears[0] + "/NPL.root, or " + dir_ntuples + v_lumiYears[0] + "/NPL_DATA.root"<<") ! Expect a crash..."))<<endl<<endl;
        // stop_program = true;
    }

    this->make_fixedRegions_templates = make_fixedRegions_templates;
    if(make_fixedRegions_templates)
    {
        if(region != "") {cout<<FRED("Warning: you have set option [make_fixedRegions_templates=true] but set a value to the [region] option ! This option will be overriden (hardcoded region selection) ")<<endl; region = "";}
        if(v_cut_name.size()) {cout<<FRED("Warning: you have set option [make_fixedRegions_templates=true] but set some manual cuts in the main()... is this intended ? Make sure that you don't cut on the same category flags used by this hardcoded strategy !")<<endl;}
    }

    //-- Hard-coded: for njet-PrivMC_tZq shape systematics, need to read pre-existing histos
    for(int isyst=0; isyst<syst_list.size(); isyst++)
    {
        if(syst_list[isyst] == "njets_tZqDown")
        {
            v_njets_SF_tZq = Get_nJets_SF("njets", "tZq", "PrivMC_tZq", v_lumiYears);
            if(v_njets_SF_tZq[0].size() == 0) {cout<<BOLD(FMAG("Warning: Get_nJets_SF() failed (missing histogram input file ? Did you include both central/private samples ?) --> Removing this systematic from the list !"))<<endl; syst_list.erase(syst_list.begin() + isyst); syst_list.erase(syst_list.begin() + isyst);} //Erase down/up variations for this syst
        }
    }
    // for(int isyst=0; isyst<syst_list.size(); isyst++) {cout<<"syst_list[isyst] "<<syst_list[isyst]<<endl;}

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
    if(luminame == "2016") {lumiValue = 36.33;} //Preliminary was: 35.92
	else if(luminame == "2017") {lumiValue = 41.53;}
    else if(luminame == "2018") {lumiValue = 59.74;}
    else if(luminame == "201617") {lumiValue = 36.33+41.53;}
    else if(luminame == "201618") {lumiValue = 36.33+59.74;}
    else if(luminame == "201718") {lumiValue = 41.53+59.74;}
    else if(luminame == "Run2") {lumiValue = 36.33+41.53+59.74;}

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
    int nmaxEv = 100000; //max nof events for train or test; -1 <-> use all available events //NB: full Run 2 is too much
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
	TCut mycut = "";
	TString tmp = "";

	//--- CHOOSE TRAINING EVENTS <--> cut on corresponding category
	TString cat_tmp = "";
	cat_tmp = Get_Category_Boolean_Name(region); //NB: does not work for NPL samples (different flag)... would need workaround

	//Even if ask templates in the SR, need to use training (looser) category for training !
	// if(cat_tmp.Contains("_SR") )
	// {
	// 	int i = cat_tmp.Index("_SR"); //Find index of substrin g
	// 	cat_tmp.Remove(i); //Remove substring
	// }
    // tmp+= cat_tmp + "==1";

    if(cat_tmp == "") {cat_tmp = "1";}

    tmp = cat_tmp;
    // tmp = cat_tmp + " && !TMath::IsNaN(mTW) && !std::isinf(abs(TopZsystem_M))"; //TMP fix... because trained on 0jet events ?

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
		if(channel == "uuu" || channel == "uu")	{mycut = "channel==0";}
		else if(channel == "uue" || channel == "ue") {mycut = "channel==1";}
		else if(channel == "eeu" || channel == "ee") {mycut = "channel==2";}
		else if(channel == "eee") {mycut = "channel==3";}
		else {cout << "WARNING : wrong channel name while training " << endl;}
	}

	cout<<"-- Will apply the following cut(s) : "<<BOLD(FGRN(""<<tmp<<""))<<endl<<endl<<endl<<endl;
	usleep(2000000); //Pause for 2s (in microsec)

	if(tmp != "") {mycut+= tmp;}

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

            TString mycut_string_tmp = mycut.GetTitle();
            if(samplename_tmp.Contains("NPL", TString::kIgnoreCase) || samplename_tmp.Contains("DY", TString::kIgnoreCase) || samplename_tmp.Contains("ttbar", TString::kIgnoreCase)) {mycut_string_tmp = mycut_string_tmp.ReplaceAll("_SR", "_SRFake"); mycut_string_tmp = mycut_string_tmp.ReplaceAll("_CR", "_CRFake");} //Modify flag for NPL samples only
            TCut mycut_tmp = mycut_string_tmp.Data();

            //-- Protections
            if(sample_list[isample] == "DATA") {continue;} //Don't use data for training

            //Can hardcode here the backgrounds against which to train, instead of considering full list of samples
            // if(signal_process == "tZq")
            // {
            //     if(!samplename_tmp.Contains("tZq") && !samplename_tmp.EndsWith("ttZ") && !samplename_tmp.Contains("ttH") && samplename_tmp.Contains("ttW") && samplename_tmp.Contains("WZ") && samplename_tmp.Contains("ZZ4l") && samplename_tmp.Contains("TTbar_DiLep") ) {continue;}
            // }

            //-- TMP sample list
            if(samplename_tmp != "tZq" && samplename_tmp != "ttZ" && samplename_tmp != "ttH" && samplename_tmp != "ttW" && samplename_tmp != "WZ" && samplename_tmp != "ZZ4l" && samplename_tmp != "TTbar_DiLep" ) {continue;}

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

    	    TFile *file_input = NULL, *file_input_train = NULL, *file_input_test = NULL;
    		TTree *tree = NULL, *tree_train = NULL, *tree_test = NULL;

            file_input = TFile::Open(inputfile);
            if(!file_input) {cout<<BOLD(FRED(<<inputfile<<" not found!"))<<endl; continue;}
            files_to_close.push_back(file_input); //store pointer to file, to close it after its events will have been read for training

            tree = (TTree*) file_input->Get(nominal_tree_name);
            if(tree==0) {cout<<BOLD(FRED("ERROR :"))<<" file "<<inputfile<<" --> *tree = 0 !"<<endl; continue;}
            else if(tree->GetEntries()==0) {cout<<"Empty ! "<<endl; continue;}
            else {cout<<FMAG("=== Opened file : ")<<inputfile<<endl<<endl;}

            // global event weights per tree (see below for setting event-wise weights)
            Double_t signalWeight     = 1.0;
            Double_t backgroundWeight = 1.0;

        //-- Choose between absolute/relative weights for training
    		if(samplename_tmp.Contains(signal_process) )
    		{
                nEvents_sig+= tree->GetEntries(mycut); myDataLoader->AddSignalTree(tree, signalWeight);

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
                nEvents_bkg+= tree->GetEntries(mycut); myDataLoader->AddBackgroundTree(tree, backgroundWeight);

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

	// if(mycuts != mycutb) {cout<<__LINE__<<FRED("PROBLEM : cuts are different for signal and background ! If this is normal, modify code -- Abort")<<endl; delete myDataLoader; return;}

	// If nTraining_Events=nTesting_Events="0", half of the events in the tree are used for training, and the other half for testing
	//NB : when converting nEvents to TString, make sure to ask for sufficient precision !

	//-- Choose dataset splitting
	TString nTraining_Events_sig = "", nTraining_Events_bkg = "", nTesting_Events_sig = "", nTesting_Events_bkg = "";

    int nTrainEvSig = ((nmaxEv == -1 || nEvents_sig * trainingEv_proportion     < nmaxEv)) ? nEvents_sig * trainingEv_proportion     : nmaxEv;
    int nTrainEvBkg = ((nmaxEv == -1 || nEvents_bkg * trainingEv_proportion     < nmaxEv)) ? nEvents_bkg * trainingEv_proportion     : nmaxEv;
    int nTestEvSig  = ((nmaxEv == -1 || nEvents_sig * (1-trainingEv_proportion) < nmaxEv)) ? nEvents_sig * (1-trainingEv_proportion) : nmaxEv;
    int nTestEvBkg  = ((nmaxEv == -1 || nEvents_bkg * (1-trainingEv_proportion) < nmaxEv)) ? nEvents_bkg * (1-trainingEv_proportion) : nmaxEv;

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

    myDataLoader->PrepareTrainingAndTestTree(mycut, mycut, "nTrain_Signal="+nTraining_Events_sig+":nTrain_Background="+nTraining_Events_bkg+":nTest_Signal="+nTesting_Events_sig+":nTest_Background="+nTesting_Events_bkg+":SplitMode=Random:EqualTrainSample:!V"); //EqualTrainSample = ?

	//-- for quick testing -- few events
	// myDataLoader->PrepareTrainingAndTestTree(mycut, mycut, "nTrain_Signal=10:nTrain_Background=10:nTest_Signal=10:nTest_Background=10:SplitMode=Random:NormMode=NumEvents:!V");

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
    // method_options = "!H:!V:NTrees=200:nCuts=40:MaxDepth=4:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:NegWeightTreatment=PairNegWeightsGlobal:CreateMVAPdfs";

    //Testing
    // method_options = "!H:!V:NTrees=200:nCuts=40:MaxDepth=4:MinNodeSize=5%:UseBaggedBoost=True:BaggedSampleFraction=0.5:BoostType=Grad:Shrinkage=0.10";

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
    int nentries_max = -1; //-1 <-> use all entries; N <-> use at most N entries per sample (for testing)

    bool noSysts_inputVars = true; //true <-> don't compute syst weights for histos of input variables (not worth the CPU-time)
    bool storeWCFit_PrivSamples = true; //true <-> default; false <-> for PrivMC samples, don't store WCFit objects (only e.g. the SM point) --> Much faster, but could not do reweighting (for debugging only)
    bool include_EFTparam_inControlHist = true; //true <-> use TH1EFTs for SMEFT signals also for control histograms; false <-> use regular TH1s (much faster, but then can not superimpose EFT scenarios, etc.)
    TString make_controlPlots_inSRsubregion = ""; //"" (default) <-> make control histograms in the current region; "SRtZq/SRttZ/SRother" --> will apply cut on NN-SM (like for templates) so that we only consider the subset of event entering the corresponding subregion
    bool store_emptyPrivMC_inCR = true; //true (default) <-> don't fill templates for private SMEFT samples in CRs (negligible, speedup), but still store empty hists for combine; false <-> fill templates (and systs) as usual

    //-- Hardcoded
    bool make_njetstZqSF_file = false; //true <-> create hardcoded file containing njets histograms for tZq and PrivMC_tZq (to compute njets_tZq shape systematic)
    bool make_histos_forControlPlotsPaper = false; //true <-> hardcoded to make only the necessary histograms needed to produce the control plots for the paper
//--------------------------------------------
//--------------------------------------------

    if(template_name == "" && classifier_name != "BDT" && classifier_name != "NN") {cout<<BOLD(FRED("Error : classifier_name value ["<<classifier_name<<"] not supported !"))<<endl; return;}
    if(template_name=="") {template_name = classifier_name;}
    if(!makeHisto_inputVars && !make_SMvsEFT_templates_plots && categorization_strategy>0 && template_name!="BDT" && !template_name.Contains("NN")) {cout<<BOLD(FRED("Error : template_name value ["<<template_name<<"] not supported for this strategy !"))<<endl; return;} //Special case: if want to produce SM vs SM classifier plots, make sure we consider NN templates (for BDT, must specify in main)
    if(this->make_fixedRegions_templates && makeHisto_inputVars) {cout<<FRED("Warning: option [make_fixedRegions_templates] is incompatible with control plots for input variables... Setting it to false !")<<endl; this->make_fixedRegions_templates = false;}
    if(template_name.Contains("NN") && template_name!="NN" && template_name!="NN_ctz" && template_name!="NN_ctw" && template_name!="NN_cpq3" && template_name!="NN_5D" && template_name!="NN_SM") {cout<<BOLD(FRED("Error : template_name value ["<<template_name<<"] not recognized !"))<<endl; return;}

    cout<<endl<<BYEL("                          ")<<endl<<endl;
	if(makeHisto_inputVars) {cout<<FYEL("--- Producing Input variables histograms ---")<<endl;}
    else if(this->make_fixedRegions_templates) {cout<<FYEL("--- Producing [otherRegions] Templates ---")<<endl;}
	else {cout<<FYEL("--- Producing "<<template_name<<" Templates ---")<<endl;}
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

    //-- Hard-coded: determine whether we are producing EFT templates using a 'predefined strategy' (cf. definitions in main --> different strategies to apply automatic cuts)
    bool use_predefined_EFT_strategy = false;
    bool use_maxNode_events = false;
    bool apply_MVASM_cut = false;
    TString MVA_type = "";
    if(categorization_strategy > 0 && ((!makeHisto_inputVars && !this->make_fixedRegions_templates) || (makeHisto_inputVars && make_controlPlots_inSRsubregion != ""))) //Only for templates, an for SR-3l (--> used to subdivide into SRtZq/SRttZ/SRother subregions)
    {
        use_predefined_EFT_strategy = true;

        //By default, apply minimal event preselection and set tZq sa signal process
        //Otherwise, will only consider the region/signal of interest
        if(cat_tmp == "")
        {
            cat_tmp = "signal";
            signal_process = "tZq";
        }

        if(make_SMvsEFT_templates_plots)
        {
            MVA_type = "EFT" + Convert_Number_To_TString(this->categorization_strategy);
            if(this->scanOperators_paramNN) {MVA_type+= "param";}
        }
        else {MVA_type = "SM";}

        if(this->categorization_strategy == 2 && (make_SMvsEFT_templates_plots || (!make_SMvsEFT_templates_plots && plot_onlyMaxNodeEvents)) ) {use_maxNode_events = true;} //Cases for which we want to cut on a multiclass MVA-SM
        else if(this->categorization_strategy == 1 && (make_SMvsEFT_templates_plots || (!make_SMvsEFT_templates_plots && plot_onlyMVACutEvents)) ) {apply_MVASM_cut = true;}
    }

    //-- Can fill these vectors via call to Get_VectorAllEvents_passMVACut() to determine beforehand, for each sample, whether each event passes or not a given MVA cut
    vector<int> v_isEventPassMVACut_multiclass; //For MVA-SM-multiclass
    vector<int> v_isEventPassMVACut_tZq; //For MVA-SM-tZq
    vector<int> v_isEventPassMVACut_ttZ; //For MVA-SM-ttZ
    vector<vector<float>> v_nodes_multiclass, v_nodes_tZq, v_nodes_ttZ; //New: also store output values

    //-- Don't make systematics shifted histos for input vars (too long)
    //Force removal of systematics ; restore values at end of this func
    if(make_histos_forControlPlotsPaper) {noSysts_inputVars = false; make_histos_forControlPlotsPaper = true;} //Must include systematics
    vector<TString> restore_syst_list = syst_list;
    vector<TString> restore_systTree_list = systTree_list;
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
        array_jetPileupID = new double[4];
        // array_fakeFactor = new double[2]; //David's FR
        array_fakeFactor = new double[12]; //ttH FR
        array_ME = new double[6];
        array_alphaS = new double[2];
        array_PDFtotal = new double[2];
        array_partonShower = new double[4];
        array_LepEffLoose_mu = new double[2];
        array_LepEffTight_mu = new double[2];
        array_LepEffLoose_el = new double[2];
        array_LepEffTight_el = new double[2];
    }

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
                // reader->AddVariable(var_list[i].Data(), &var_list_floats[i]);
                var_list_pfloats[i] = new float; reader->AddVariable(var_list[i].Data(), var_list_pfloats[i]);
            }
            for(int i=0; i<v_cut_name.size(); i++)
            {
                if(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("==")) {reader->AddVariable(v_cut_name[i].Data(), &v_cut_float[i]);}
            }
        }
    }

    //-- Default template binning (modified in 'Get_Template_Range()')
    nbins = 15;
    float xmin = -1, xmax = 1; //BDT: [-1,1]

    //-- Define ranges of jet/bjets multiplicities -- for 'categ' templates only (modified in 'Get_Template_Range()')
    int nbjets_min = 1, nbjets_max=2, njets_min=2, njets_max=6;

    if(!makeHisto_inputVars) {make_njetstZqSF_file = false;}

    //-- Create output file //NB: the same conventions must be enforced in function 'Get_HistoFile_InputPath' (to automatically find/read this output file later)
    TString output_file_name = "";
    if(makeHisto_inputVars) {output_file_name = (TString) "outputs/ControlHistograms" + (region == ""? "":"_"+region) + "_" + lumiName + filename_suffix +".root";}
    else {output_file_name = (TString) "outputs/Templates" + (this->make_fixedRegions_templates? "_otherRegions":(template_name == ""? "":"_"+template_name)) + (MVA_type == ""? "":"_"+MVA_type) + ((use_predefined_EFT_strategy || region == "")? "":"_"+region) + "_" + lumiName + filename_suffix + ".root";}
    // else {output_file_name = (TString) "outputs/Templates" + (this->make_fixedRegions_templates? "_otherRegions":(template_name == ""? "":"_"+template_name)) + (MVA_type == "" || !make_SMvsEFT_templates_plots? "":"_"+MVA_type) + ((use_predefined_EFT_strategy || region == "")? "":"_"+region) + "_" + lumiName + filename_suffix + ".root";}
    if(make_njetstZqSF_file) {output_file_name = "./outputs/njets_tZq.root";}
    if(make_histos_forControlPlotsPaper) {output_file_name = "./outputs/ControlPlotsPaper.root";}
    TFile* file_output = TFile::Open(output_file_name, "RECREATE");

    int nhistos = 0; //Count nof histos written to the output file

    vector<vector<vector<vector<int>>>> v4_nEntries_year_proc_var_bin(v_lumiYears.size()); //Count nof actual entries (events) entering each bin of each nominal histo (info not stored in TH1Fs)

    if(make_controlPlots_inSRsubregion != "" && make_controlPlots_inSRsubregion != "SRtZq" && make_controlPlots_inSRsubregion != "SRttZ" && make_controlPlots_inSRsubregion != "SRother") {cout<<FRED("Error: wrong [make_controlPlots_inSRsubregion] value !")<<endl; make_controlPlots_inSRsubregion = "";}


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

	// cout<<endl<<ITAL("-- Ntuples directory : "<<dir_ntuples<<"")<<endl<<endl;

    // float tmp_compare = 0;
    float total_nentries_toRead = Count_Total_Nof_Entries(dir_ntuples, nominal_tree_name, sample_list, systTree_list, v_cut_name, v_cut_def, v_lumiYears, makeHisto_inputVars, noSysts_inputVars);

    cout<<endl<<FBLU(OVERLINE("                           "))<<endl;
    cout<<FBLU(BOLD("...Will READ (not necessarily process) a total of "<<std::setprecision(12)<<total_nentries_toRead<<" entries..."))<<endl;
    cout<<FBLU(UNDL("                           "))<<endl<<endl<<endl;

    //-- Draw progress bar
    bool draw_progress_bar = false;
    if(total_nentries_toRead < 200000) {draw_progress_bar = false;}
    Int_t ibar = 0; //event counter
    TMVA::Timer timer(total_nentries_toRead, "", true);
    TMVA::gConfig().SetDrawProgressBar(1);
    TMVA::gConfig().SetUseColor(1);

    bool MVA_alreadyLoaded = false; //If reading same MVA for multiple years, need to load weights only once
    TString BDT_method_name = "BDT"; //Need to store BDT method name for evaluation

    //-- YEARS LOOP
    vector<TString> total_var_list; vector<float> total_var_floats; vector<float*> total_var_pfloats;
    for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
    {
        cout<<endl<<UNDL(FMAG("=== YEAR : "<<v_lumiYears[iyear]<<""))<<endl<<endl;

        v4_nEntries_year_proc_var_bin[iyear].resize(sample_list.size());


// #    # #    #   ##
// ##  ## #    #  #  #
// # ## # #    # #    #
// #    # #    # ######
// #    #  #  #  #    #
// #    #   ##   #    #

        int inode_clfy2 = 0; //When reading 2 NNs, may either need to read first (NN-EFTs) or second (NN-SM) node of clfy2 //Obsolete?

        //-- BDT
        if(!makeHisto_inputVars && (template_name == "BDT" || template_name.Contains("NN")) && (use_specificMVA_eachYear || !MVA_alreadyLoaded)) //Read MVA input file
        {
            MVA_alreadyLoaded = true;
            TString MVA_input_path = Get_MVAFile_InputPath(template_name, signal_process, v_lumiYears[iyear], use_specificMVA_eachYear, make_SMvsEFT_templates_plots, false, categorization_strategy, this->scanOperators_paramNN, this->use_NN_cpq3_SRttZ);
            if(MVA_input_path == "") {return;} //MVA file not found

            if(template_name == "BDT") //Book TMVA reader from .xml weight file
            {
                reader->BookMVA(BDT_method_name, MVA_input_path);
            }
            else if(template_name.Contains("NN")) //Read NN info and weights
            {
                //-- Get path of NN info text file --> Read list of input variables (and more)
                var_list_NN.resize(0); NN_iMaxNode = -1; NN_strategy = ""; NN_inputLayerName = ""; NN_outputLayerName = ""; NN_nNodes = -1; minmax_bounds.clear();
                TString NNinfo_input_path = Get_MVAFile_InputPath(template_name, signal_process, v_lumiYears[iyear], use_specificMVA_eachYear, make_SMvsEFT_templates_plots, true, categorization_strategy, this->scanOperators_paramNN, this->use_NN_cpq3_SRttZ);
                if(NNinfo_input_path == "") {return;} //MVA file not found

                if(Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN, v_NN_nodeLabels, NN_inputLayerName, NN_outputLayerName, NN_iMaxNode, NN_nNodes, minmax_bounds, &NN_strategy)) {clfy1 = new TFModel(MVA_input_path.Data(), var_list_NN.size(), NN_inputLayerName.Data(), NN_nNodes, NN_outputLayerName.Data());} //Load neural network model
                else {return;} //Error: missing NN infos

                if(this->scanOperators_paramNN && NN_strategy != "MVA_param") {cout<<BOLD(FRED("ERROR ! You have set [scanOperators_paramNN=true] to scan over EFT operators, but you are reading a non-parametrized MVA ! Please fix this ! Aborting..."))<<endl; return;}
                else if(this->scanOperators_paramNN && !template_name.Contains("NN")) {cout<<BOLD(FRED("ERROR ! You have set [scanOperators_paramNN=true] to scan over EFT operators, but you have selected 'template_name != NN' (MVA scan only makes sense for parametrized NN) ! Please fix this ! Aborting..."))<<endl; return;}

                var_list = var_list_NN; //Use NN input features
                // var_list_floats.resize(var_list.size());
                var_list_pfloats.resize(var_list.size()); for(int ivar=0; ivar<var_list_pfloats.size(); ivar++) {var_list_pfloats[ivar] = new float;} //Allocate necessary memory

                //-- If making NN-EFT templates using strategy 1 or 2, need to access 2 MVA-EFTs (1 for tZq + 1 for ttZ)
                // if(this->make_SMvsEFT_templates_plots && use_predefined_EFT_strategy >= 1 && cat_tmp != "tZq" && template_name != "NN_SM")
                if(this->make_SMvsEFT_templates_plots && use_predefined_EFT_strategy >= 1 && cat_tmp != "tZq" && template_name != "NN_SM" && (use_NN_cpq3_SRttZ || template_name != "NN_cpq3")) //Don't need clfy2 for cpq3 by default
                {
                    TString MVA_input_path = Get_MVAFile_InputPath(template_name, "ttZ", v_lumiYears[iyear], use_specificMVA_eachYear, true, false, categorization_strategy, this->scanOperators_paramNN, this->use_NN_cpq3_SRttZ);
                    if(MVA_input_path == "") {return;} //MVA file not found

                    TString NNinfo_input_path = Get_MVAFile_InputPath(template_name, "ttZ", v_lumiYears[iyear], use_specificMVA_eachYear, true, true, categorization_strategy, this->scanOperators_paramNN, this->use_NN_cpq3_SRttZ);
                    if(NNinfo_input_path == "") {return;} //MVA file not found

                    //-- Read info file for second NN //NB: it is assumed that the NN strategy (e.g. parametrized NN) is the same as for the first NN !
                    var_list_NN2.resize(0); NN2_iMaxNode = -1; NN2_strategy = ""; NN2_inputLayerName = ""; NN2_outputLayerName = ""; NN2_nNodes = -1; minmax_bounds2.clear();
                    if(Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN2, v_NN2_nodeLabels, NN2_inputLayerName, NN2_outputLayerName, NN2_iMaxNode, NN2_nNodes, minmax_bounds2)) {clfy2 = new TFModel(MVA_input_path.Data(), var_list_NN2.size(), NN2_inputLayerName.Data(), NN2_nNodes, NN2_outputLayerName.Data());} //Load neural network model
                    else {return;} //Error: missing NN infos

                    // var_list_floats_2.resize(var_list_NN2.size());
                    var_list_pfloats_2.resize(var_list_NN2.size()); for(int ivar=0; ivar<var_list_pfloats_2.size(); ivar++) {var_list_pfloats_2[ivar] = new float; *var_list_pfloats_2[ivar] = 0;} //Reserve necessary memory

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

                    // if(template_name == "NN_cpq3") {inode_clfy2 = 1;} //Special case: ttZ/cpq3 --> Read NN-SM (ttZ node) in SRttZ //Obsolete ?
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


 // #    #   ##   #####  #   ##   #####  #      ######  ####
 // #    #  #  #  #    # #  #  #  #    # #      #      #
 // #    # #    # #    # # #    # #####  #      #####   ####
 // #    # ###### #####  # ###### #    # #      #           #
 //  #  #  #    # #   #  # #    # #    # #      #      #    #
 //   ##   #    # #    # # #    # #####  ###### ######  ####

        //-- Define the list of variables there, because may need to first read a MVA info file (year-dependent) //Define list only for first year, assume it's identical for all years
        //-- NB: for NN templates, list of variables may be hard-coded here based on the info found in the MVA file --> Make sure the correct MVA files are found
        if(!iyear)
        {
        	//-- Vector of variables --> Loop
        	if(makeHisto_inputVars) //Input variables
        	{
        		for(int i=0; i<var_list.size(); i++) {total_var_list.push_back(var_list[i]);}
        		for(int i=0; i<v_add_var_names.size(); i++) {total_var_list.push_back(v_add_var_names[i]);}
        	}
        	else //Templates
        	{
                //-- Could use this helper function to fill the (hardcoded) list of variables to consider, based on user-options
                Fill_Variables_List(total_var_list, use_predefined_EFT_strategy, template_name, this->region, this->scanOperators_paramNN, this->NN_nNodes, this->make_SMvsEFT_templates_plots, operator_scan1, operator_scan2, v_WCs_operator_scan1, v_WCs_operator_scan2, this->make_fixedRegions_templates);
        	} //Templates

            // total_var_floats.resize(total_var_list.size());
            total_var_pfloats.resize(total_var_list.size()); for(int ivar=0; ivar<total_var_pfloats.size(); ivar++) {total_var_pfloats[ivar] = new float;} //Reserve necessary memory
        }
        // for(int ivar=0; ivar<total_var_list.size(); ivar++) {cout<<"total_var_list[ivar] "<<total_var_list[ivar]<<endl;}

        if(make_njetstZqSF_file) //Hardcoded
        {
            sample_list.clear();
            sample_list.push_back("tZq");
            sample_list.push_back("PrivMC_tZq");
            total_var_list.clear();
            total_var_list.push_back("njets");
        }
        if(make_histos_forControlPlotsPaper)
        {
            //-- Hardcoded list of variables
            total_var_list.clear();
            total_var_list.push_back("njets");
            total_var_list.push_back("nbjets");
            total_var_list.push_back("jPrimeAbsEta");
            total_var_list.push_back("jprime_Pt");
            total_var_list.push_back("metEt");
            total_var_list.push_back("mTW");
            total_var_list.push_back("lAsymmetry");
            total_var_list.push_back("recoZ_Pt");
            total_var_list.push_back("recoZ_dPhill");
            total_var_list.push_back("dR_tZ");
            total_var_list.push_back("maxDeepJet");
        }


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

            bool isFake = false; //Identify NPL samples (different flags)
            if(sample_list[isample].Contains("NPL", TString::kIgnoreCase) || sample_list[isample].Contains("DY", TString::kIgnoreCase) || sample_list[isample].Contains("ttbar", TString::kIgnoreCase)) {isFake = true;}

            //-- For current year/sample, store the nominal integral(s) per channel and variable (for nominal TTree only) --> Can use it to rescale some systematics later (normalize to nominal)
            vector<vector<float>> v_storeNominalIntegral_chan_var(channel_list.size());
            for(int ichan=0; ichan<v_storeNominalIntegral_chan_var.size(); ichan++) {v_storeNominalIntegral_chan_var[ichan].resize(total_var_list.size());}

            //-- Store nominal bin contents/errors for all channel/variable/histo bin --> for systematics, can then adjust the bin content to avoid that either nominal or variation is set to ~0, but not the other (causes combine failures)
            //-- For current year/sample, store the nominal integral(s) per channel and variable (for nominal TTree only) --> Can use it to rescale some systematics later (normalize to nominal)
            vector<vector<vector<float>>> v_storeNominalBinContents_chan_var_bin(channel_list.size()), v_storeNominalBinErrors_chan_var_bin(channel_list.size());
            for(int ichan=0; ichan<v_storeNominalBinContents_chan_var_bin.size(); ichan++) {v_storeNominalBinContents_chan_var_bin[ichan].resize(total_var_list.size()); v_storeNominalBinErrors_chan_var_bin[ichan].resize(total_var_list.size());}

    		//-- Open input TFile
    		TString inputfile = dir_ntuples + v_lumiYears[iyear] + "/" + sample_list[isample] + ".root";

            //-- Apply ad-hoc scale factor to private ttZ sample so that SM yield matches that of central ttZ sample
            //-- NB: should correct the xsec at potato level instead !
            // float SF_SMEFT_ttZ = 0.86;

    		// cout<<"inputfile "<<inputfile<<endl;
    		if(!Check_File_Existence(inputfile))
    		{
    			cout<<endl<<"File "<<inputfile<<FRED(" not found!")<<endl;
    			continue;
    		}

            //-- Open input file a first time only to read and store some information (sums of weights, etc.)
            TFile* file_input = TFile::Open(inputfile, "READ");
            cout<<endl<<FBLU("Opened file "<<inputfile<<" ...")<<endl;

            //-- Read/store SWE at different EFT points (but only actually need the SM one for normalization, since this is the ref xsec value we use always)
            bool isPrivMC = false;
            vector<float> v_SWE_EFT; //Store Sums of Weights (SWE) for all reweight points -- for private MC samples only
            if(sample_list[isample].Contains("PrivMC") && !sample_list[isample].Contains("_c")) {isPrivMC = true;}
            if(isPrivMC)
            {
                TString hSWE_name = "EFT_SumWeights";

                if(!file_input->GetListOfKeys()->Contains(hSWE_name)) {cout<<"ERROR ! Histogram "<<hSWE_name<<" containing the sums of weights not found... Abort !"<<endl; return;}

                //-- Read and store sums of weights (SWE)
                TH1F* h_SWE = (TH1F*) file_input->Get(hSWE_name);
                for(int ibin=0; ibin<h_SWE->GetNbinsX(); ibin++)
                {
                    v_SWE_EFT.push_back(h_SWE->GetBinContent(ibin+1)); //1 SWE stored for each stored weight
                    // cout<<"v_SWE_EFT["<<ibin<<"] = "<<v_SWE_EFT[ibin]<<endl;
                }
                delete h_SWE; h_SWE = NULL;
            }

            //-- Read/store SWEs for different TH systematics before preselection, for renormalization (only keep acceptance effects)
            TH1F* h_sumWeights_beforeSel = NULL; //Read histogram containing the SFs (yield_variation/yield_nominal) *before event selection* --> Only used for signal samples, to retain only acceptance effects for TH systematics
            vector<float> v_SWE_beforeSel; //Store Sums of Weights (SWE) for nominal and theory systematics *before applying any event preselection*
            if(isPrivMC || sample_list[isample] == "tZq" || sample_list[isample] == "ttZ" || sample_list[isample] == "tWZ")
            {
                TString h_sumWeights_beforeSel_name = "h_sumWeights_beforeSel";
                if(!file_input->GetListOfKeys()->Contains(h_sumWeights_beforeSel_name)) {cout<<endl<<FRED("Histogram "<<h_sumWeights_beforeSel_name<<" not found in file "<<inputfile<<" !")<<endl;}

                //-- Read and store sums of weights (SWE)
                h_sumWeights_beforeSel = (TH1F*) file_input->Get(h_sumWeights_beforeSel_name);
                for(int ibin=0; ibin<h_sumWeights_beforeSel->GetNbinsX(); ibin++)
                {
                    v_SWE_beforeSel.push_back(h_sumWeights_beforeSel->GetBinContent(ibin+1));
                    // cout<<"v_SWE_beforeSel["<<ibin<<"] = "<<v_SWE_beforeSel[ibin]<<endl;
                }
                delete h_sumWeights_beforeSel; h_sumWeights_beforeSel = NULL;
            }

            //-- Close the input file to delete associated objects cleany
            file_input->Close(); file_input = NULL;

    		//-- Loop on TTrees : first empty element corresponds to nominal TTree ; additional TTrees may correspond to JES/JER TTrees (defined in main)
    		//NB : only nominal TTree contains systematic weights ; others only contain the nominal weight (but variables have different values)
            for(int itree=0; itree<systTree_list.size(); itree++)
    		{
                if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name && (sample_list[isample] == "DATA" || sample_list[isample] == "DY" || sample_list[isample].Contains("NPL") || sample_list[isample].Contains("TTbar") ) ) {continue;} //Special cases for which only nominal tree must be considered

                TString treename_tmp = systTree_list[itree];
    			if(treename_tmp == "") {treename_tmp = nominal_tree_name;}
                cout<<DIM("-- Tree "<<treename_tmp<<" --")<<endl;

                //-- Re-open the input file for each individual TTree (easier to cleanly delete all associated objects, prevents conflicts)
                file_input = TFile::Open(inputfile, "READ");

                //-- Read TTree
                TTree* tree = (TTree*) file_input->Get(treename_tmp);
                // tree->SetDirectory(0); //Need to dis-associated the TTree from the TFile, else may get (rare) crashes

    			if(!tree || !tree->GetEntries())
    			{
    				cout<<BOLD(FRED("ERROR : tree '"<<treename_tmp<<"' not found in file : "<<inputfile<<" ! Skip !"))<<endl;
                    if(!itree) {break;} //If we did not find the nominal TTree, no need to look for the others
                    continue; //Skip TTree
    			}

                //-- If making templates 'predefined strategy', first need to evaluate the response of the NN_SM for all events in the tree
                //NB: passing the TTree as argument (make sure to properly reset addresses, etc. afterwards to avoid conflicts)
                v_isEventPassMVACut_tZq.clear(); v_isEventPassMVACut_ttZ.clear(); v_isEventPassMVACut_multiclass.clear();
                if(!makeHisto_inputVars && use_predefined_EFT_strategy)
                {
                    if(this->categorization_strategy == 1)
                    {
                        if(apply_MVASM_cut)
                        {
                            if(cat_tmp == "signal" || cat_tmp == "tZq") {if(!Get_VectorAllEvents_passMVACut(v_isEventPassMVACut_tZq, v_nodes_tZq, "tZq", classifier_name, tree, v_lumiYears[iyear], cut_value_tZq, keep_aboveCut, use_specificMVA_eachYear, categorization_strategy, false, nentries_max, "tZq", also_applyCut_onMaxNodeValue, isFake)) {continue;}}
                            if(cat_tmp == "signal" || cat_tmp == "ttZ") {if(!Get_VectorAllEvents_passMVACut(v_isEventPassMVACut_ttZ, v_nodes_ttZ, "ttZ", classifier_name, tree, v_lumiYears[iyear], cut_value_ttZ, keep_aboveCut, use_specificMVA_eachYear, categorization_strategy, false, nentries_max, "ttZ", also_applyCut_onMaxNodeValue, isFake)) {continue;}}
                        }
                    }
                    else if(this->categorization_strategy == 2 && use_maxNode_events) {if(!Get_VectorAllEvents_passMVACut(v_isEventPassMVACut_multiclass, v_nodes_multiclass, "Multiclass", classifier_name, tree, v_lumiYears[iyear], -1., keep_aboveCut, use_specificMVA_eachYear, categorization_strategy, false, nentries_max, "signal", also_applyCut_onMaxNodeValue, isFake)) {continue;}}
                }
                else if(makeHisto_inputVars && make_controlPlots_inSRsubregion != "" && this->categorization_strategy == 2 && use_maxNode_events) //New: can also make control plots in a specific subregion
                {
                    if(!Get_VectorAllEvents_passMVACut(v_isEventPassMVACut_multiclass, v_nodes_multiclass, "Multiclass", classifier_name, tree, v_lumiYears[iyear], -1., keep_aboveCut, use_specificMVA_eachYear, categorization_strategy, false, nentries_max, "signal", also_applyCut_onMaxNodeValue, isFake)) {continue;}
                }
                tree->ResetBranchAddresses(); //Make sure to reset all branch addresses before continuing

                //-- SMEFT SAMPLES: need to compute/read EFT parameterization, so that we can propagate the EFT parameterization later in combine
                bool EFTparameterization_alreadyStored = false; //Check whether the SMEFT parameterization was already stored previously (--> only need to read it)
                TTree* tree_EFTparameterization = NULL; //The per-event SMEFT parameterization may have been stored already (using [Split_FullSamples] code) --> If this is the case, can simply read it (much faster)
                WCFit* eft_fit = NULL;
                if(isPrivMC && (systTree_list[itree] == "" || systTree_list[itree] == nominal_tree_name) && (include_EFTparam_inControlHist || !makeHisto_inputVars))
                {
                    eft_fit = new WCFit("myfit"); //Need to initialize the object for all TH1EFT objects, even if not used

                    //-- Check whether the SMEFT parameterization was already stored previously (--> only need to read it)
                    TString tree_EFTparameterization_name = "EFTparameterization";
                    // if(systTree_list[itree] == "TotalDown") {tree_EFTparameterization_name+= "_JESDown";}
                    if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name) {tree_EFTparameterization_name+= "_" + systTree_list[itree];}
                    // cout<<"tree_EFTparameterization_name "<<tree_EFTparameterization_name<<endl;

                    if(file_input->GetListOfKeys()->Contains(tree_EFTparameterization_name))
                    // if(systTree_list[itree] == "" && file_input->GetListOfKeys()->Contains(tree_EFTparameterization_name))
                    {
                        EFTparameterization_alreadyStored = true;
                        tree_EFTparameterization = (TTree*) file_input->Get(tree_EFTparameterization_name);
                        tree_EFTparameterization->SetBranchAddress("eft_fit", &eft_fit);
                        cout<<DIM("Reading EFT parameterization directly from TTree '"<<tree_EFTparameterization_name<<"'...")<<endl;
                    }
                    else if(!sample_list[isample].Contains("TOP19001")) {cout<<endl<<FRED("Warning: EFT parameterization not found in file "<<inputfile<<", will need to compute it on-the-fly (slow)...")<<endl<<endl;}
                }

                v4_nEntries_year_proc_var_bin[iyear][isample].resize(total_var_list.size());


//   ##   #####  #####  #####  ######  ####   ####  ######  ####
//  #  #  #    # #    # #    # #      #      #      #      #
// #    # #    # #    # #    # #####   ####   ####  #####   ####
// ###### #    # #    # #####  #           #      # #           #
// #    # #    # #    # #   #  #      #    # #    # #      #    #
// #    # #####  #####  #    # ######  ####   ####  ######  ####

    			//Disactivate all un-necessary branches ; below, activate only needed ones
    			tree->SetBranchStatus("*", 0); //disable all branches by default, speed up

                //-- Always keep track of these vars (needed for different templates, syst, etc.)
                float* njets = new float;
                float* mTW = new float;
                float* channel = new float;
                tree->SetBranchStatus("mTW", 1); tree->SetBranchAddress("mTW", mTW);
                tree->SetBranchStatus("njets", 1); tree->SetBranchAddress("njets", njets);
                tree->SetBranchStatus("channel", 1); tree->SetBranchAddress("channel", channel);

                float nbjets; //Needed for 'categ' templates

                if(makeHisto_inputVars)
    			{
    				for(int i=0; i<total_var_list.size(); i++)
    				{
                        if(total_var_list[i] == "mTW") {total_var_pfloats[i] = mTW; continue;}
                        if(total_var_list[i] == "njets") {total_var_pfloats[i] = njets; continue;}
                        if(total_var_list[i] == "channel") {total_var_pfloats[i] = channel; continue;}

                        tree->SetBranchStatus(total_var_list[i], 1);
                        // tree->SetBranchAddress(total_var_list[i], &total_var_floats[i]);
                        tree->SetBranchAddress(total_var_list[i], total_var_pfloats[i]);
    				}

    			}
    			else //Templates
    			{
                    if(template_name=="BDT" || template_name.Contains("NN")) //Book input variables in same order as for MVA training //Activate input features needed for MVA evaluation (same as used for training)
                    {
                        for(int i=0; i<var_list.size(); i++)
                        {
                            if(!isample) {cout<<DIM("MVA 1 -- "<<i<<" -- Activate variable '"<<var_list[i]<<"'")<<endl;}
                            if(var_list[i] == "ctz" || var_list[i] == "ctw" || var_list[i] == "cpq3" || var_list[i] == "cpqm" || var_list[i] == "cpt") {continue;} //WC input values are arbitrary, there is no address to set !

                            if(var_list[i] == "mTW") {var_list_pfloats[i] = mTW; continue;}
                            if(var_list[i] == "njets") {var_list_pfloats[i] = njets; continue;}
                            if(var_list[i] == "channel") {var_list_pfloats[i] = channel; continue;}

                            // cout<<"Activate branch "<<var_list[i]<<"("<<var_list_pfloats[i]<<")"<<endl;
                            tree->SetBranchStatus(var_list[i], 1);
                            tree->SetBranchAddress(var_list[i], var_list_pfloats[i]);
                        }

                        if(clfy2) //Also set branch addresses for second MVA
                        {
                            cout<<endl<<endl;
                            for(int i=0; i<var_list_NN2.size(); i++)
                            {
                                if(!isample) {cout<<DIM("MVA 2 "<<i<<" -- Activate variable '"<<var_list_NN2[i]<<"'")<<endl;}

                                if(var_list_NN2[i] == "ctz" || var_list_NN2[i] == "ctw" || var_list_NN2[i] == "cpq3" || var_list_NN2[i] == "cpqm" || var_list_NN2[i] == "cpt") {continue;} //WC input values are arbitrary, there is no address to set !
                                if(var_list_NN2[i] == "mTW") {var_list_pfloats_2[i] = mTW; continue;}
                                if(var_list_NN2[i] == "njets") {var_list_pfloats_2[i] = njets; continue;}
                                if(var_list_NN2[i] == "channel") {var_list_pfloats_2[i] = channel; continue;}

                                int index_sameVar_in_NN1_list = -1; //Check whether same input variable was already set in first MVA (can not set branch addresses twice !)
                                for(int j=0; j<var_list_NN.size(); j++)
                                {
                                    if(var_list_NN2[i] == var_list_NN[j])
                                    {
                                        // cout<<"Var "<<var_list_NN2[i]<<" already found in clfy1..."<<endl;
                                        var_list_pfloats_2[i] = var_list_pfloats[j];
                                        index_sameVar_in_NN1_list = j; break; //Found same variable in var_list_NN
                                    }
                                }

                                if(index_sameVar_in_NN1_list < 0) //Variable only found in 2nd MVA
                                {
                                    // cout<<"Activate branch "<<var_list_NN2[i]<<"("<<var_list_pfloats_2[i]<<")"<<endl;
                                    tree->SetBranchStatus(var_list_NN2[i], 1);
                                    tree->SetBranchAddress(var_list_NN2[i], var_list_pfloats_2[i]);
                                    // tree->SetBranchAddress(var_list_NN2[i], &var_list_floats_2[i]);
                                }
                            }
                        } //2nd MVA
                    } //MVA templates
                    else if(template_name=="categ") //Need to read jet/bjet multiplicities of each event
                    {
                        //NB: njets branch read by default
                        tree->SetBranchStatus("nbjets", 1);
                        tree->SetBranchAddress("nbjets", &nbjets);
                    }
                    else if(template_name=="Zpt") //Need to read Zpt variable //NB: no conflict with any MVA input, etc.
                    {
                        tree->SetBranchStatus("recoZ_Pt", 1);
                        // tree->SetBranchAddress("recoZ_Pt", &total_var_floats[0]);
                        tree->SetBranchAddress("recoZ_Pt", total_var_pfloats[0]);
                    }
                    else if(template_name=="ZptCos") //Need to read Zpt and cosThetaStarPolZ variables
                    {
                        tree->SetBranchStatus("recoZ_Pt", 1);
                        // tree->SetBranchAddress("recoZ_Pt", &total_var_floats[0]);
                        tree->SetBranchAddress("recoZ_Pt", total_var_pfloats[0]);
                        tree->SetBranchStatus("cosThetaStarPolZ", 1);
                        tree->SetBranchAddress("cosThetaStarPolZ", total_var_pfloats[1]);
                        // tree->SetBranchAddress("cosThetaStarPolZ", &total_var_floats[1]);
                    }
                    else if(template_name=="mTW") {total_var_pfloats[0] = mTW;} //Set address to dedicated mTW variable
                } //Templates

    			for(int i=0; i<v_cut_name.size(); i++)
    			{
                    // cout<<"v_cut_name[i] "<<v_cut_name[i]<<endl;

    				tree->SetBranchStatus(v_cut_name[i], 1);
                    if(v_cut_name[i].BeginsWith("is") || v_cut_name[i].BeginsWith("passed") ) //Categories are encoded into Char_t, not float
    				{
                        if(!isFake) {tree->SetBranchAddress(v_cut_name[i], &v_cut_char[i]);}
                        else {tree->SetBranchStatus(v_cut_name[i]+"Fake", 1); tree->SetBranchAddress(v_cut_name[i]+"Fake", &v_cut_char[i]);} //Fake <-> separate flag
    				}
    				else //All others are floats
    				{
    					tree->SetBranchAddress(v_cut_name[i], &v_cut_float[i]);
    				}
    			}

    			//--- Cut on relevant event selection stored as flag (Char_t)
                //NB: don't 'call' the flag directly, use function Get_Category_Boolean_Name() instead (need different flags for NPL samples)
    			Char_t is_goodCategory = true; //Cut on event category flag
                vector<Char_t> v_is_goodCategory(4, 1); //Needed for 'make_fixedRegions_templates' scenario (multiple regions)
                if(this->make_fixedRegions_templates)
                {
                    //-- For now: consider 3 regions (ttZ 4l, CR WZ, CR ZZ, CR DY)
                    TString cat_name = Get_Category_Boolean_Name("ttz4l", isFake);
                    tree->SetBranchStatus(cat_name, 1);
                    tree->SetBranchAddress(cat_name, &v_is_goodCategory[0]);
                    cat_name = Get_Category_Boolean_Name("wz", isFake);
                    tree->SetBranchStatus(cat_name, 1);
                    tree->SetBranchAddress(cat_name, &v_is_goodCategory[1]);
                    cat_name = Get_Category_Boolean_Name("zz", isFake);
                    tree->SetBranchStatus(cat_name, 1);
                    tree->SetBranchAddress(cat_name, &v_is_goodCategory[2]);
                    cat_name = Get_Category_Boolean_Name("dy", isFake);
                    tree->SetBranchStatus(cat_name, 1);
                    tree->SetBranchAddress(cat_name, &v_is_goodCategory[3]);
                }
                else if(cat_tmp != "")
                {
                    TString cat_name = Get_Category_Boolean_Name(cat_tmp, isFake);

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
                if(!makeHisto_inputVars && use_predefined_EFT_strategy && categorization_strategy == 1)
                {
                    tree->SetBranchStatus(Get_Category_Boolean_Name("tZq", isFake), 1); tree->SetBranchAddress(Get_Category_Boolean_Name("tZq", isFake), &is_tZqSR);
                    tree->SetBranchStatus(Get_Category_Boolean_Name("ttZ", isFake), 1); tree->SetBranchAddress(Get_Category_Boolean_Name("ttZ", isFake), &is_ttZSR);
                }

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
                vector<float>* v_wgts = NULL;
                vector<string>* v_ids = NULL;
                float SMweight = 0.; //Read SM weight directly (hardcoded in potato)
                if(isPrivMC)
                {
                    tree->SetBranchStatus("SMweight", 1);
                    tree->SetBranchAddress("SMweight", &SMweight);

                    //CHANGED -- only need reweighting features (<-> read all reweights) for nominal TTree ! For syst TTrees, only consider SM-like TH1Fs
                    if(systTree_list[itree] == "" || systTree_list[itree] == nominal_tree_name)
                    {
                        v_wgts = new vector<float>;
                        tree->SetBranchStatus("mc_EFTweights", 1);
                        tree->SetBranchAddress("mc_EFTweights", &v_wgts);

                        v_ids = new vector<string>;
                        tree->SetBranchStatus("mc_EFTweightIDs", 1);
                        tree->SetBranchAddress("mc_EFTweightIDs", &v_ids);
                    }
                }

                //-- Reserve 1 float for each systematic weight (also for nominal to keep ordering, but not used)
    			vector<Double_t*> v_double_systWeights(syst_list.size(), NULL);
    			for(int isyst=0; isyst<syst_list.size(); isyst++)
    			{
    				//-- Protections : not all syst weights apply to all samples, etc.
    				if(sample_list[isample] == "DATA") {break;}
    				else if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name) {break;} //Syst event weights only stored in nominal TTree
                    else if((sample_list[isample].Contains("NPL") && syst_list[isyst] != "" && !syst_list[isyst].BeginsWith("FR")) || (!sample_list[isample].Contains("NPL") && syst_list[isyst].BeginsWith("FR"))) {continue;} //NPL <-> only fakes sytematics; all others <-> no fakes systematics
                    else if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefir") ) {continue;} //no prefire in 2018
                    else if((syst_list[isyst].BeginsWith("PDF") || syst_list[isyst].BeginsWith("ME") || syst_list[isyst].BeginsWith("alpha") || syst_list[isyst].BeginsWith("ISR") || syst_list[isyst].BeginsWith("FSR")) && !sample_list[isample].Contains("PrivMC") && sample_list[isample] != "tZq" && sample_list[isample] != "ttZ") {continue;} //NB: no TH weights available in central tWZ V12
                    else if(syst_list[isyst].Contains("njets_tZq") && !sample_list[isample].Contains("PrivMC_tZq")) {continue;} //Only applies to LO tZq

                    //Set proper branch address (hard-coded mapping)
                    SetBranchAddress_SystVariationArray(tree, syst_list[isyst], v_double_systWeights, isyst);
    			}


 //                                    #
 // # #    # # #####    ##### #    #  ##
 // # ##   # #   #        #   #    # # #
 // # # #  # #   #        #   ######   #
 // # #  # # #   #        #   #    #   #
 // # #   ## #   #        #   #    #   #
 // # #    # #   #        #   #    # #####

    			//-- Reserve memory for 1 TH1F per category, per systematic, per variable //v3 <-> vec of vec of vec
                //-- Idem for TH1EFT
    			vector<vector<vector<TH1F*>>> v3_histo_chan_syst_var(channel_list.size());
                vector<vector<vector<TH1EFT*>>> v3_TH1EFT_chan_syst_var(channel_list.size());

    			for(int ichan=0; ichan<channel_list.size(); ichan++)
    			{
                    //-- Cases for which we only need to store the nominal weight
    				if((channel_list.size() > 1 && channel_list[ichan] == "") || sample_list[isample] == "DATA"  || (systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name))
                    {
                        v3_histo_chan_syst_var[ichan].resize(1);
                        v3_TH1EFT_chan_syst_var[ichan].resize(1);
                    }

                    //-- Reserve memory for TH1Fs/TH1EFTs
    				if(sample_list[isample] == "DATA" || (systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name)) //1 single weight
    				{
                        v3_histo_chan_syst_var[ichan].resize(1); //Cases for which we only need to store the nominal weight
                        v3_TH1EFT_chan_syst_var[ichan].resize(1);
    				}
    				else //Need 1 histo for nominal + 1 histo per systematic
    				{
                        v3_histo_chan_syst_var[ichan].resize(syst_list.size());
                        v3_TH1EFT_chan_syst_var[ichan].resize(syst_list.size()); //Actually only need TH1EFTs for nominal histos... but keep same structure as TH1Fs for symmetry
    				}

    				//-- Init TH1Fs/TH1EFTs
    				for(int isyst=0; isyst<v3_histo_chan_syst_var[ichan].size(); isyst++)
    				{
                        // cout<<"syst_list[isyst] "<<syst_list[isyst]<<endl;

                        v3_histo_chan_syst_var[ichan][isyst].resize(total_var_list.size());
                        v3_TH1EFT_chan_syst_var[ichan][isyst].resize(total_var_list.size());

    					for(int ivar=0; ivar<total_var_list.size(); ivar++)
    					{
                            // cout<<"total_var_list[ivar] "<<total_var_list[ivar]<<endl;

                            //-- Init to null
                            v3_histo_chan_syst_var[ichan][isyst][ivar] = NULL;
                            v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = NULL;

    						if(makeHisto_inputVars) //Kinematic variables
                            {
                                if(!Get_Variable_Range(total_var_list[ivar], nbins, xmin, xmax)) {cout<<FRED("Unknown variable name : "<<total_var_list[ivar]<<"! (include it in function Get_Variable_Range() in Helper.cxx)")<<endl; continue;} //Get binning for this input variable
                                if(region.Contains("ttz4l", TString::kIgnoreCase) && !total_var_list[ivar].Contains("njets") && !total_var_list[ivar].Contains("nbjets") && !total_var_list[ivar].Contains("channel")) {nbins = (int) (nbins / 2);} //Need tighter binning in low-stat regions

                                if(make_histos_forControlPlotsPaper) //Hardcoded ranges
                                {
                                    nbins = 10;

                                    if(total_var_list[ivar] == "njets") {xmin = 2; xmax = 8; nbins = 6;}
                                    else if(total_var_list[ivar] == "nbjets") {xmin = 1; xmax = 4; nbins = 3;}
                                    else if(total_var_list[ivar] == "jPrimeAbsEta" || total_var_list[ivar] == "jprime_Pt" || total_var_list[ivar] == "dR_tZ" || total_var_list[ivar] == "maxDeepJet") {nbins = 8;}
                                    // else if(total_var_list[ivar] == "metEt" || total_var_list[ivar] == "mTW" || total_var_list[ivar] == "lAsymmetry" | total_var_list[ivar] == "jprime_Pt") {nbins = 15;}
                                }
                            }
                            else //Templates
                            {
                                if(total_var_list[ivar].Contains("SRtZq")) {Get_Template_Range(nbins, xmin, xmax, total_var_list[ivar], this->make_SMvsEFT_templates_plots, this->categorization_strategy, plot_onlyMaxNodeEvents, nbjets_min, nbjets_max, njets_min, njets_max, minmax_bounds, use_NN_SRother, use_NN_cpq3_SRttZ);}
                                else {Get_Template_Range(nbins, xmin, xmax, total_var_list[ivar], this->make_SMvsEFT_templates_plots, this->categorization_strategy, plot_onlyMaxNodeEvents, nbjets_min, nbjets_max, njets_min, njets_max, minmax_bounds2, use_NN_SRother, use_NN_cpq3_SRttZ);}
                            }
                            // cout<<"nbins "<<nbins<<" / xmin "<<xmin<<" / xmax "<<xmax<<endl;

                            // v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = new TH1EFT("", "", nbins, xmin, xmax);
                            // v3_histo_chan_syst_var[ichan][isyst][ivar] = new TH1F("", "", nbins, xmin, xmax);
                            if(isPrivMC && (systTree_list[itree] == "" || systTree_list[itree] == nominal_tree_name) && syst_list[isyst] == "" && (!this->make_fixedRegions_templates || total_var_list[ivar].Contains("ttZ4l")) && (include_EFTparam_inControlHist || !makeHisto_inputVars)) //Only fill/store TH1EFTs for nominal private SMEFT samples in SRs //For all the rest, can use regular TH1s
                            {
                                v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = new TH1EFT("", "", nbins, xmin, xmax);
                            }
                            else //Don't need TH1EFTs in this case, use regular TH1s
                            {
                                v3_histo_chan_syst_var[ichan][isyst][ivar] = new TH1F("", "", nbins, xmin, xmax);
                            }

                            v4_nEntries_year_proc_var_bin[iyear][isample][ivar].resize(nbins);
    					} //var
    				} //syst
    			} //chan

                //-- For private EFT samples, get and store index of SM reweight
                int idx_sm = 0; //Hardcoded: SM index = 0 !
                int nweights = 25; //For my SMEFT samples, only need >= 21 EFT weights for parameterization
                //NB: obsolete, now reading SMweight directly via hardcoded variable //Keep this more generic approach
                /*
                if(isPrivMC)
                {
                    // if(sample_list[isample].Contains("TOP19001")) {nweights = 160;} //For TOP19001 samples, need many more weights for parameterization (hard-coded) //Actually, don't parameterize these samples, too slow !

                    tree->GetEntry(0); //Read 1 entry

                    for(int iwgt=0; iwgt<v_ids->size(); iwgt++)
                    {
                        if(ToLower(v_ids->at(iwgt)).Contains("_sm") ) {idx_sm = iwgt; break;} //SM weight found
                        if(((TString) v_ids->at(iwgt)).Contains("EFTrwgt183_") ) {idx_sm = iwgt; break;} //SM weight found //TOP19001 sample
                    }

                    if(EFTparameterization_alreadyStored) {tree->SetBranchStatus("mc_EFTweightIDs", 0);} //If simply reading EFT parameterization from file, don't need to read names of EFT points for each event; however, need to read vector of weights to get the proper SM event weights
                    if(idx_sm == -1) {cout<<BOLD(FRED("Error: SM reweight not found in private sample ! Abort ! "))<<endl; return;}
                    else {cout<<DIM("idx_sm = "<<idx_sm<<")")<<endl;}
                }*/


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
                cout<<"--- Reading " << nentries << " events" <<endl;
                // cout<< "--- "<<sample_list[isample]<<" : Processing " << nentries << " events" << std::endl;
                cout<<DIM("(path: "<<inputfile<<")")<<endl;

    			for(int ientry=0; ientry<nentries; ientry++)
    			{
    				// cout<<FGRN("ientry "<<ientry<<"")<<endl;

                    if(isPrivMC && ientry%5000==0) {cout<<DIM("Entry "<<ientry<<"")<<endl;} //Very slow, print progress

                    // std::fill(var_list_floats.begin(), var_list_floats.end(), 0); //Reset vectors reading inputs to 0

    				tree->GetEntry(ientry);

    				if(std::isnan(eventWeight*eventMCFactor) || std::isinf(eventWeight*eventMCFactor))
    				{
    					cout<<BOLD(FRED("* Found event with eventWeight*eventMCFactor = "<<eventWeight*eventMCFactor<<" ; remove it..."))<<endl; continue;
    				}

//---- APPLY EVENT CUTS HERE -------------------------

                    //--- Cut on category flag
                    if(!is_goodCategory) {continue;}
                    if(this->make_fixedRegions_templates)
                    {
                        if(std::accumulate(v_is_goodCategory.begin(), v_is_goodCategory.end(), 0) != 1) {continue;} //Check whether at least 1 cat. flag is true
                    }

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

                    if(makeHisto_inputVars && make_controlPlots_inSRsubregion != "")
                    {
                        if(make_controlPlots_inSRsubregion == "SRtZq" && v_isEventPassMVACut_multiclass[ientry] != 0) {continue;}
                        else if(make_controlPlots_inSRsubregion == "SRttZ" && v_isEventPassMVACut_multiclass[ientry] != 1) {continue;}
                        else if(make_controlPlots_inSRsubregion == "SRother" && v_isEventPassMVACut_multiclass[ientry] != 2) {continue;}
                    }

//--------------------------------------------

                    ibar++;
                    if(draw_progress_bar && ibar%50000==0) {timer.DrawProgressBar(ibar, ""); cout<<DIM("Processed "<<ibar<<" / "<<total_nentries_toRead<<"")<<endl; }

                    // cout<<"//-------------------------------------------- "<<endl;
                    // cout<<"eventWeight "<<eventWeight<<endl;
                    // cout<<"eventMCFactor "<<eventMCFactor<<endl;
                    // for(int ivar=0; ivar<var_list_floats.size(); ivar++) {cout<<"var_list_floats "<<ivar<<" = "<<var_list_floats[ivar]<<endl;} //Debug
                    // for(int ivar=0; ivar<var_list_floats.size(); ivar++) {cout<<"var_list_floats "<<ivar<<" "<<var_list[ivar]<<" / "<<var_list_pfloats[ivar]<<" = "<<*var_list_pfloats[ivar]<<endl;} //Debug
                    // cout<<"njets "<<*njets<<endl;
                    // cout<<"mTW "<<*mTW<<endl;

    				//-- S i n g l e    t e m p l a t e     v a l u e (<-> only consider 'template_name')
                    std::vector<float> clfy1_outputs, clfy2_outputs; //Compute NN responses only once per event
                    if(!makeHisto_inputVars) //Templates
                    {
                        if(template_name == "BDT") {*total_var_pfloats[0] = reader->EvaluateMVA(BDT_method_name);} //BDT output value
                        else if(template_name.Contains("NN") && !use_predefined_EFT_strategy) //NN output value //Default
                        {
                            // clfy1_outputs = clfy1->evaluate(var_list_floats); //Evaluate output node(s) value(s)
                            clfy1_outputs = clfy1->evaluate(var_list_pfloats); //Evaluate output node(s) value(s)
                            float mva_tmp = -1;
                            NN_iMaxNode = -1;
                            for(int inode=0; inode<NN_nNodes; inode++)
                            {
                                // total_var_floats[inode] = clfy1_outputs[inode];
                                // if(clfy1_outputs[inode] > mva_tmp) {mva_tmp = clfy1_outputs[inode]; NN_iMaxNode = inode;} //Store max. node info
                                *total_var_pfloats[inode] = clfy1_outputs[inode];
                                if(clfy1_outputs[inode] > mva_tmp) {mva_tmp = clfy1_outputs[inode]; NN_iMaxNode = inode;} //Store max. node info
                            }

                            if(clfy2) //Also read 2nd MVA
                            {
                                NN2_iMaxNode = -1;
                                // clfy2_outputs = clfy2->evaluate(var_list_floats_2); //Evaluate output node(s) value(s)
                                clfy2_outputs = clfy2->evaluate(var_list_pfloats_2); //Evaluate output node(s) value(s)
                                float mva_tmp = -1;
                                for(int inode=0; inode<NN2_nNodes; inode++)
                                {
                                    if(clfy2_outputs[inode_clfy2] > mva_tmp) {mva_tmp = clfy2_outputs[inode_clfy2]; NN2_iMaxNode = inode;}
                                }
                            } //2nd NN
                        } //NN && !use_predefined_EFT_strategy

                        else if(template_name == "categ") //Arbitrary binning depending on jets/bjets multiplicities
                        {
                            // total_var_floats[0] = Get_x_jetCategory(njets, nbjets, nbjets_min, nbjets_max, njets_min, njets_max);
                            *total_var_pfloats[0] = Get_x_jetCategory(*njets, nbjets, nbjets_min, nbjets_max, njets_min, njets_max);
                            // cout<<"njets "<<njets<<" / nbjets "<<nbjets<<" --> categ "<<total_var_floats[0]<<endl;
                        }
                        else if(template_name == "countExp")
                        {
                            *total_var_pfloats[0] = 0.5;
                        }
                        else if(template_name == "channel") //Hard-coded from main SM analysis
                        {
                            *total_var_pfloats[0] = (*channel==1 || *channel==3)? 0.5:1.5;
                            // total_var_floats[0] = (*channel==1 || *channel==3)? 0.5:1.5;
                        }
                        else if(template_name == "ZptCos")
                        {
                            *total_var_pfloats[0] = Get_x_ZptCosCategory(*total_var_pfloats[0], *total_var_pfloats[1]);
                            // total_var_floats[0] = Get_x_ZptCosCategory(total_var_floats[0], total_var_floats[1]);
                        }
                        //NB: no need to set anything more for Zpt/mTW/... templates: variable address already set before event loop
                    } //Templates

                    if(v_lumiYears[iyear] =="2016") {eventMCFactor*= 36.33/35.92;} //Quick fix: account here for the updated 2016 luminosity (current NTuples produced with preliminary value) //FIXME

                    double weight_tmp = eventWeight * eventMCFactor; //Fill histo with this weight ; manipulate differently depending on syst
                    // cout<<"eventWeight "<<eventWeight<<" / eventMCFactor "<<eventMCFactor<<" / weight_tmp "<<weight_tmp<<endl;

                    //-- S M E F T   w e i g h t
                    float w_SMpoint = 0.;
                    if(isPrivMC)
                    {
                        w_SMpoint = weight_tmp * SMweight / (weightMENominal * v_SWE_EFT[idx_sm]); //Hardcoded SM index
                        // w_SMpoint = weight_tmp * v_wgts->at(idx_sm) / (weightMENominal * v_SWE_EFT[idx_sm]); //Keep also more generic definition

                        // cout<<"v_wgts->at(idx_sm) "<<v_wgts->at(idx_sm)<<" / weightMENominal "<<weightMENominal<<" / v_SWE_EFT[idx_sm] "<<v_SWE_EFT[idx_sm]<<endl;

                        //-- Tmp fix: wrong eventMCFactor and wrong SWEs
                        // if(sample_list[isample] == "PrivMC_ttZ_TOP19001") {w_SMpoint*= 2.482*20;}
                        // else if(sample_list[isample] == "PrivMC_tZq_TOP19001") {w_SMpoint*= 3.087*20;}

                        if(storeWCFit_PrivSamples && (systTree_list[itree] == "" || systTree_list[itree] == nominal_tree_name) && (include_EFTparam_inControlHist || !makeHisto_inputVars)) //Need EFT parameterization //NB: not looping on shapeSsysts/variables yet
                        {
                            if(!EFTparameterization_alreadyStored && !sample_list[isample].Contains("TOP19001")) //EFT parameterization not stored, need to determine it //NB: don't compute parameterization for TOP19001 samples, way too slow
                            {
                                Get_WCFit(eft_fit, v_ids, v_wgts, v_SWE_EFT, weight_tmp, weightMENominal, w_SMpoint, 0, nweights); //Hardcoded SM index=0

                                //-- Apply ad-hoc scale factor to private sample so that SM yield matches that of central sample
                                // if(sample_list[isample] == "PrivMC_ttZ") {eft_fit->scale(SF_SMEFT_ttZ);}
                            }
                            else {tree_EFTparameterization->GetEntry(ientry);} //Else, eft_fit is read automatically from dedicated TTree

                            //-- Debug printout
                            // cout<<"first WCFit name "<<eft_fit->getNames()[0]<<endl;
                        }
                    } //Private SMEFT samples

                    //-- F i l l     a l l      v a r i a b l e s (<-> consider full list of variables 'total_var_list')

    				//-- Can choose to fill event only into histogram corresponding to relevant leptonic channel
    				for(int ichan=0; ichan<channel_list.size(); ichan++)
    				{
                        //-- Fill histos for sub-channels with corresponding events
                        if(channel_list[ichan] == "uuu" && *channel != 0) {continue;}
                        if(channel_list[ichan] == "uue" && *channel != 1) {continue;}
                        if(channel_list[ichan] == "eeu" && *channel != 2) {continue;}
                        if(channel_list[ichan] == "eee" && *channel != 3) {continue;}

                        bool cat_alreadyFound = false; //For pre-defined strategies, fill orthogonal categories
                        for(int ivar=0; ivar<total_var_list.size(); ivar++)
                        {
                            //-- TEMPLATES
                            //NB: only need to fill additional variables here in 2 cases: a) if 'use_predefined_EFT_strategy=true'; b) if 'make_fixedRegions_templates=true' --> predefined set of multiple templates ! (otherwise, considering single template, which was filled above)
                            if(!makeHisto_inputVars)
                            {
                                if(store_emptyPrivMC_inCR && isPrivMC && this->make_fixedRegions_templates && !total_var_list[ivar].Contains("ttZ4l") ) {continue;} //Special case: neglect completely EFT signals in CRs (negligible, speedup) -- also systs ! //NB: we will still create/store the corresponding (nominal) histograms, because combine needs some signal (even empty) in all regions //EDIT: actually not, simply use to add option '--X-allow-no-background' in fits ?

                                //NB: case [template_name == "NN/BDT" && !use_predefined_EFT_strategy] already taken care of above
                                if(use_predefined_EFT_strategy) //Apply 'predefined strategy' --> if event does not pass the required cut, don't fill the corresponding template
                                {
                                    int idx_EFT_var = ivar; //Trick: by default, for predefined EFT strategies, there are 3 variables (hard-coded); but if we scan operators, there are 3 * N1 * N2 variables --> Determine which 'original variable' (tzq, ttz, cr) they correspond to
                                    if(scanOperators_paramNN && this->scanOperators_paramNN)
                                    {
                                        if(ivar < (v_WCs_operator_scan1.size() * v_WCs_operator_scan2.size())) {idx_EFT_var = 0;}
                                        else if(ivar < 2 * (v_WCs_operator_scan1.size() * v_WCs_operator_scan2.size())) {idx_EFT_var = 1;}
                                        else {idx_EFT_var = 2;} //CR template
                                    }

                                    if(this->categorization_strategy == 1) //Strategy 1 //Cut on flags & MVA-SM
                                    {
                                        if(cat_alreadyFound && !this->scanOperators_paramNN) {continue;} //Don't consider same event twice
                                        else if(idx_EFT_var == 0 && (!is_tZqSR || (apply_MVASM_cut && !v_isEventPassMVACut_tZq[ientry]))) {continue;} //Event did not pass the required MVA-SM-tZq cut
                                        else if(idx_EFT_var == 1 && (!is_ttZSR || (apply_MVASM_cut && !v_isEventPassMVACut_ttZ[ientry]))) {continue;} //Event did not pass the required MVA-SM-ttZ cut
                                        //else: idx_EFT_var == 2 (SRother) --> contains all events which did not enter SRtZq/SRttZ

                                        cat_alreadyFound = true; //Did not hit 'continue' so far <-> event category has been found
                                    }
                                    else //Strategy 2 //Cut only on MVA-SM
                                    {
                                        // cout<<"ivar "<<ivar<<" / v_isEventPassMVACut_multiclass[ientry] "<<v_isEventPassMVACut_multiclass[ientry]<<endl;
                                        if(use_maxNode_events && v_isEventPassMVACut_multiclass[ientry] != idx_EFT_var) {continue;} //Current node is not max node

                                        //-- Could hardcode different cuts on NN-SM nodes here
                                        // if(idx_EFT_var == 0) {if(v_isEventPassMVACut_multiclass[ientry] != 0 || v_nodes_multiclass[ientry][2]>0.2) {continue;}}
                                        // else if(idx_EFT_var == 1) {if(v_isEventPassMVACut_multiclass[ientry] != 1 || v_nodes_multiclass[ientry][2]>0.2) {continue;}}
                                        // else if(idx_EFT_var == 2) {if(v_isEventPassMVACut_multiclass[ientry] != 2 && v_nodes_multiclass[ientry][2]<0.2) {continue;}}
                                    }

                                    //-- If scanning EFT operators, need to set the values of input WCs accordingly for param. NN templates
                                    if(scanOperators_paramNN && this->scanOperators_paramNN && idx_EFT_var < 2)
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
                                                        *var_list_pfloats[idx1_operator_scan1] = v_WCs_operator_scan1[iop1];
                                                        if(idx1_operator_scan2>=0) {*var_list_pfloats[idx1_operator_scan2] = v_WCs_operator_scan2[iop2];}
                                                        // var_list_floats[idx1_operator_scan1] = v_WCs_operator_scan1[iop1];
                                                        // var_list_floats[idx1_operator_scan2] = v_WCs_operator_scan2[iop2];
                                                        // cout<<iop1<<" var "<<idx1_operator_scan1<<" ==> "<<v_WCs_operator_scan1[iop1]<<endl;
                                                        // cout<<iop2<<" var "<<idx1_operator_scan2<<" ==> "<<v_WCs_operator_scan2[iop2]<<endl;
                                                    }
                                                    else if(idx_EFT_var == 1) //ttZ MVA
                                                    {
                                                        *var_list_pfloats_2[idx2_operator_scan1] = v_WCs_operator_scan1[iop1];
                                                        if(idx2_operator_scan2>=0) {*var_list_pfloats_2[idx2_operator_scan2] = v_WCs_operator_scan2[iop2];}
                                                        // var_list_floats_2[idx2_operator_scan1] = v_WCs_operator_scan1[iop1];
                                                        // var_list_floats_2[idx2_operator_scan2] = v_WCs_operator_scan2[iop2];
                                                    }
                                                }
                                            }
                                        }

                                        //-- Re-evaluate output nodes values *at this particular EFT point*
                                        if(idx_EFT_var == 0)
                                        {
                                            // clfy1_outputs = clfy1->evaluate(var_list_floats);
                                            // total_var_floats[ivar] = clfy1_outputs[0];
                                            clfy1_outputs = clfy1->evaluate(var_list_pfloats);
                                            *total_var_pfloats[ivar] = clfy1_outputs[0];
                                        }
                                        else if(idx_EFT_var == 1)
                                        {
                                            // clfy2_outputs = clfy2->evaluate(var_list_floats_2);
                                            // total_var_floats[ivar] = clfy2_outputs[0];
                                            clfy2_outputs = clfy2->evaluate(var_list_pfloats_2);
                                            *total_var_pfloats[ivar] = clfy2_outputs[inode_clfy2];
                                        }
                                    } //MVA EFT scan
                                    else if(template_name.Contains("NN"))
                                    {
                                        if(make_SMvsEFT_templates_plots && template_name != "NN_SM") //SM vs EFT --> use EFT MVA or mTW)
                                        {
                                            if(idx_EFT_var==0) {*total_var_pfloats[ivar] = clfy1->evaluate(var_list_pfloats)[0];} //MVA-EFT-tZq
                                            else if(idx_EFT_var==1 && template_name.Contains("cpq3") && !use_NN_cpq3_SRttZ) {*total_var_pfloats[ivar] = v_nodes_multiclass[ientry][1];} //Can choose to use MVA-SM ttZ node for cpq3 operator in SRttZ
                                            else if(idx_EFT_var==1) {*total_var_pfloats[ivar] = clfy2->evaluate(var_list_pfloats_2)[inode_clfy2];} //MVA-EFT-ttZ
                                            else if(idx_EFT_var==2 && use_NN_SRother) {*total_var_pfloats[ivar] = v_nodes_multiclass[ientry][2];} //SRother --> Use NN-bkg node
                                            else if(idx_EFT_var==2) {*total_var_pfloats[ivar] = *mTW;} //SRother --> mTW
                                            else {cout<<"ERROR: wrong index !"<<endl; continue;}
                                            // cout<<"idx_EFT_var "<<idx_EFT_var<<" / pfloat = "<<*total_var_pfloats[ivar]<<endl;
                                        }
                                        else //SM vs SM --> use multiclass MVA nodes
                                        {
                                            if(idx_EFT_var <= 1) {*total_var_pfloats[ivar] = clfy1->evaluate(var_list_pfloats)[idx_EFT_var];} //SRtZq/SRttZ --> Use NN-SM nodes
                                            else if(use_NN_SRother) {*total_var_pfloats[ivar] = v_nodes_multiclass[ientry][2];} //SRother --> Use NN-bkg node
                                            else {*total_var_pfloats[ivar] = *mTW;} //SRother --> mTW

                                            //-- Use all NN-SM nodes ? //Testing, to remove
                                            // *total_var_pfloats[ivar] = clfy1->evaluate(var_list_pfloats)[idx_EFT_var];
                                            // cout<<"idx_EFT_var "<<idx_EFT_var<<" / clfy1->evaluate(var_list_pfloats)[idx_EFT_var] = "<<clfy1->evaluate(var_list_pfloats)[idx_EFT_var]<<endl;
                                        }
                                    }
                                    else //E.g. Zpt template, ...--> Need to fill all variables (with same value) //NB: copy value, not address, else will get deleted twice and segfault !
                                    {
                                        // if(idx_EFT_var==1) {total_var_floats[ivar] = total_var_floats[0];} //Force use of first variable in both SRs
                                        // else if(idx_EFT_var==2) {total_var_floats[ivar] = *mTW;} //Force use of mTW in CR
                                        if(idx_EFT_var==1) {*total_var_pfloats[ivar] = *total_var_pfloats[0];} //Force use of first variable in both SRs
                                        else if(idx_EFT_var==2 && use_NN_SRother) {*total_var_pfloats[ivar] = v_nodes_multiclass[ientry][2];} //SRother --> Use NN-bkg node
                                        else if(idx_EFT_var==2) {*total_var_pfloats[ivar] = *mTW;} //SRother --> mTW
                                    }
                                } //predefined EFT strategy
                                else if(this->make_fixedRegions_templates)
                                {
                                    //Hardcoded: fill variables with corresponding variables, provided the event satisfies the relevant flag
                                    // if(ivar == 0 && v_is_goodCategory[ivar]) {total_var_floats[ivar] = 0.5;} //countExp ttZ 4l SR
                                    // else if(ivar == 1 && v_is_goodCategory[ivar]) {total_var_floats[ivar] = *mTW;} //mTW WZ CR
                                    // else if(ivar == 2 && v_is_goodCategory[ivar]) {total_var_floats[ivar] = 0.5;} //countExp ZZ CR
                                    // else if(ivar == 3 && v_is_goodCategory[ivar]) {total_var_floats[ivar] = 0.5;} //countExp DY CR
                                    if(ivar == 0 && v_is_goodCategory[ivar]) {*total_var_pfloats[ivar] = 0.5;} //countExp ttZ 4l SR
                                    else if(ivar == 1 && v_is_goodCategory[ivar]) {*total_var_pfloats[ivar] = *mTW;} //mTW WZ CR //NB: copy value, not pointer, else will get deleted twice !
                                    else if(ivar == 2 && v_is_goodCategory[ivar]) {*total_var_pfloats[ivar] = 0.5;} //countExp ZZ CR
                                    else if(ivar == 3 && v_is_goodCategory[ivar]) {*total_var_pfloats[ivar] = 0.5;} //countExp DY CR
                                    else {continue;} //Categories are non-orthogonal, so make sure we don't fill unwanted entries
                                }
                            } //Templates //NB: ALL OTHER CASES SHOULD HAVE BEEN TAKEN CARE OF ABOVE THIS
                            // cout<<total_var_list[ivar]<<" = "<<*total_var_pfloats[ivar]<<endl; //Debug printout

                            //-- S y s t  e m a t i c    w e i g h t s
        					for(int isyst=0; isyst<syst_list.size(); isyst++)
        					{
                                //-- Protections : not all syst weights apply to all samples, etc.
                                if((sample_list[isample] == "DATA") && syst_list[isyst] != "") {break;}
                                else if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name && syst_list[isyst] != "") {break;} //Syst event weights only stored in nominal TTree
                                else if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefir") ) {continue;} //no prefire in 2018
                                else if((sample_list[isample].Contains("NPL") && syst_list[isyst] != "" && !syst_list[isyst].BeginsWith("FR")) || (!sample_list[isample].Contains("NPL") && syst_list[isyst].BeginsWith("FR"))) {continue;} //NPL <-> only fakes sytematics; all others <-> no fakes systematics
                                else if( (syst_list[isyst].BeginsWith("PDF") || syst_list[isyst].BeginsWith("ME") || syst_list[isyst].BeginsWith("alpha") || syst_list[isyst].BeginsWith("ISR") || syst_list[isyst].BeginsWith("FSR")) && !sample_list[isample].Contains("PrivMC") && sample_list[isample] != "tZq" && sample_list[isample] != "ttZ") {continue;} //Only considered for signal samples  //NB: no TH weights available in central tWZ V12
                                else if(syst_list[isyst].Contains("njets_tZq") && !sample_list[isample].Contains("PrivMC_tZq")) {continue;} //Only applies to LO tZq

        						// cout<<"-- sample "<<sample_list[isample]<<" / channel "<<channel_list[ichan]<<" / syst "<<syst_list[isyst]<<endl;

                                //-- Start from nominal weight (no syst) //Re-initialize for each syst
                                if(isPrivMC) {weight_tmp = w_SMpoint;}
        						else {weight_tmp = eventWeight*eventMCFactor;}

                                //No prefiring uncertainty for 2018 samples ; still write the histo (=nominal) to avoid troubles down the way... ?
                                // if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefir")) {*(v_double_systWeights[isyst]) = 1;}

                                // if(isyst>0) cout<<"syst_list[isyst] "<<syst_list[isyst]<<" : "<<*(v_double_systWeights[isyst])<<endl;
                                // cout<<"nominal : "<<weight_tmp<<endl;
                                if(syst_list[isyst] != "") //Syst weights were already divided by nominal weight
                                {
                                    if(syst_list[isyst].Contains("njets_tZq")) {weight_tmp*= Apply_nJets_SF(v_njets_SF_tZq, (int) *njets, iyear, syst_list[isyst]);}
                                    else if(abs(*(v_double_systWeights[isyst])) < 100) //Safe protection against huge, non-physical PDF/PSweights/... in PrivMC_tZq/ttZ samples (why?)
                                    {
                                        weight_tmp*= *(v_double_systWeights[isyst]);
                                    }
                                }
                                // cout<<"syst : "<<weight_tmp<<endl;

        						if(std::isnan(weight_tmp) || std::isinf(weight_tmp))
        						{
        							cout<<BOLD(FRED("* Found event with syst. weight = "<<weight_tmp<<" ; remove it..."))<<endl;
        							cout<<"(sample "<<sample_list[isample]<<" / channel "<<channel_list[ichan]<<" / year "<<v_lumiYears[iyear]<<" / syst "<<syst_list[isyst]<<")"<<endl;
        							continue;
        						}

                                //-- F i l l     h i s t o
                                if(isPrivMC && (systTree_list[itree] == "" || systTree_list[itree] == nominal_tree_name) && syst_list[isyst] == "" && (!this->make_fixedRegions_templates || total_var_list[ivar].Contains("ttZ4l")) && (include_EFTparam_inControlHist || !makeHisto_inputVars) ) //Only fill/store TH1EFTs for nominal private SMEFT samples in SRs //For all the rest, can use regular TH1s
                                {
                                    Fill_TH1EFT_UnderOverflow(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar], *total_var_pfloats[ivar], w_SMpoint, *eft_fit);

                                    // if(syst_list[isyst] == "") {Fill_TH1EFT_UnderOverflow(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar], total_var_floats[ivar], w_SMpoint, *eft_fit);} //Nominal --> need TH1EFT to store WCFit objects
                                    // else {Fill_TH1F_UnderOverflow(v3_histo_chan_syst_var[ichan][isyst][ivar], total_var_floats[ivar], weight_tmp);} //Weight systematics --> can use regular TH1F objects

                                    //-- Can count the signal nentries in each histo bin here (for xchecks)
                                    int bin = (int)(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->GetXaxis()->FindBin(*total_var_pfloats[ivar]))-1;
                                    if(bin<0) {bin=0;} else if(bin>=v4_nEntries_year_proc_var_bin[iyear][isample][ivar].size()) {bin = v4_nEntries_year_proc_var_bin[iyear][isample][ivar].size()-1;} //Deal with under/overflow
                                    v4_nEntries_year_proc_var_bin[iyear][isample][ivar][bin]++; //Count how many (signal) events end up in each bin
                                }
                                else //Use regular TH1 object
                                {
                                    Fill_TH1F_UnderOverflow(v3_histo_chan_syst_var[ichan][isyst][ivar], *total_var_pfloats[ivar], weight_tmp);
                                }
                            } //syst loop
    					} //var loop

    					if(channel_list[ichan] != "") {break;} //subcategories are orthogonal ; if already found, can break subcat. loop
    				} //subcat/chan loop
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
    				for(int isyst=0; isyst<syst_list.size(); isyst++)
    				{
                        //-- Protections : not all syst weights apply to all samples, etc.
                        if(sample_list[isample] == "DATA" && syst_list[isyst] != "") {break;}
                        else if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name && syst_list[isyst] != "") {break;} //Syst event weights only stored in nominal TTree
                        else if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefir") ) {continue;} //no prefire in 2018
                        else if((sample_list[isample].Contains("NPL") && syst_list[isyst] != "" && !syst_list[isyst].BeginsWith("FR")) || (!sample_list[isample].Contains("NPL") && syst_list[isyst].BeginsWith("FR"))) {continue;} //NPL <-> only fakes sytematics; all others <-> no fakes systematics
                        else if((syst_list[isyst].BeginsWith("PDF") || syst_list[isyst].BeginsWith("ME") || syst_list[isyst].BeginsWith("alpha") || syst_list[isyst].BeginsWith("ISR") || syst_list[isyst].BeginsWith("FSR")) && !sample_list[isample].Contains("PrivMC") && sample_list[isample] != "tZq" && sample_list[isample] != "ttZ") {continue;}  //NB: no TH weights available in central tWZ V12
                        else if(syst_list[isyst].Contains("njets_tZq") && !sample_list[isample].Contains("PrivMC_tZq")) {continue;} //Only applies to LO tZq

    					for(int ivar=0; ivar<total_var_list.size(); ivar++)
    					{
                            if(isPrivMC && (syst_list[isyst] != "" || (systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name)) && this->make_fixedRegions_templates && !total_var_list[ivar].Contains("ttZ4l") && store_emptyPrivMC_inCR) {continue;} //Don't write shifted PrivMC histos in CRs (only need nominals to avoid bug in combine)

                            // cout<<"channel "<<channel_list[ichan]<<" / systTree "<<systTree_list[isyst]<<" / syst "<<syst_list[isyst]<<" / var "<<total_var_list[ivar]<<endl;

    						TString output_histo_name;
                            output_histo_name = total_var_list[ivar];
                            if(channel_list[ichan] != "") {output_histo_name+= "_" + channel_list[ichan];}
                            if(region != "" && !makeHisto_inputVars && !use_predefined_EFT_strategy && !this->make_fixedRegions_templates) {output_histo_name+= "_" + region;}
                            output_histo_name+= "_" + v_lumiYears[iyear] + "__" + samplename;
							if(syst_list[isyst] != "" || (systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name)) {output_histo_name+= "__" + Get_Modified_SystName(syst_list[isyst]+systTree_list[itree], v_lumiYears[iyear], samplename);}
							else if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name) {output_histo_name+= "__" + systTree_list[itree];}

                            //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                            output_histo_name.ReplaceAll('-', 'm');

    						file_output->cd();

                            if(isPrivMC && (systTree_list[itree] == "" || systTree_list[itree] == nominal_tree_name) && syst_list[isyst] == "" && (!this->make_fixedRegions_templates || total_var_list[ivar].Contains("ttZ4l")) && (include_EFTparam_inControlHist || !makeHisto_inputVars)) //Only fill/store TH1EFTs for nominal private SMEFT samples in SRs //For all the rest, can use regular TH1s
                            {
                                if(sample_list[isample] != "NPL_MC" && sample_list[isample] != "DATA") //NB: NPL_MC *should* have negative bins (substraction); and allow empty data
                                {
                                    //-- Nominal: prevent negative bins, and error>content
                                    Avoid_Histogram_EmptyOrNegativeBins(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]);

                                    v_storeNominalBinContents_chan_var_bin[ichan][ivar].resize(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->GetNbinsX());
                                    v_storeNominalBinErrors_chan_var_bin[ichan][ivar].resize(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->GetNbinsX());
                                    for(int ibin=0; ibin<v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->GetNbinsX(); ibin++)
                                    {
                                        v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin] = v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1);
                                        v_storeNominalBinErrors_chan_var_bin[ichan][ivar][ibin] = v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->GetBinError(ibin+1);
                                        // cout<<"v_storeNominalBinContents_chan_var_bin["<<ichan<<"]["<<ivar<<"]["<<ibin<<"] "<<v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin]<<endl;
                                        // cout<<"v_storeNominalBinErrors_chan_var_bin["<<ichan<<"]["<<ivar<<"]["<<ibin<<"] "<<v_storeNominalBinErrors_chan_var_bin[ichan][ivar][ibin]<<endl;
                                    }
                                }

                                //-- Apply ad-hoc scale factor to private sample so that SM yield matches that of central sample
                                // if(sample_list[isample] == "PrivMC_ttZ") {v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Scale(SF_SMEFT_ttZ);}

                                //-- Normalize specific systematics (to nominal, before of after selection)
                                if(systTree_list[itree] == "" && syst_list[isyst] == "") {v_storeNominalIntegral_chan_var[ichan][ivar] = v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Integral();} //Store nominal norm.
                                else if(syst_list[isyst].BeginsWith("PDF") || syst_list[isyst].BeginsWith("ME") || syst_list[isyst].BeginsWith("alpha") || syst_list[isyst].BeginsWith("ISR") || syst_list[isyst].BeginsWith("FSR")) //Apply SF corresponding to (var/nom) *before selection* --> Remaining norm. effect is from acceptance
                                {
                                    // cout<<"Before: v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Integral() "<<v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Integral()<<endl;
                                    Scale_THSyst_toBeforeSelection(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar], v_SWE_beforeSel, syst_list[isyst]);
                                    // cout<<"After: v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Integral() "<<v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Integral()<<endl;
                                    // Scale_THSyst_toBeforeSelection(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar], h_sumWeights_beforeSel, syst_list[isyst]);

                                    // v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Scale(v_storeNominalIntegral_chan_var[ichan][ivar]/v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Integral());
                                }
                                else if(syst_list[isyst].BeginsWith("njets_tZq")) //Special case: rescale syst to nominal norm.
                                {
                                    v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Scale(v_storeNominalIntegral_chan_var[ichan][ivar]/v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Integral());
                                }

                                file_output->cd();
                                v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Write(output_histo_name);
                                nhistos++; //Count nof histos written to output file

                                // cout<<"Wrote TH1EFT histo : "<<output_histo_name<<" (Integral: "<<v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Integral()<<")"<<endl;

                                //-- BIN SPLITTING //For MVA-EFT: need to store each histogram bin separately so that they can be scaled independently in Combine
                                if(this->split_EFTtemplates_perBin && !makeHisto_inputVars && make_SMvsEFT_templates_plots && !this->make_fixedRegions_templates) {StoreEachHistoBinIndividually(file_output, v3_TH1EFT_chan_syst_var[ichan][isyst][ivar], output_histo_name, nhistos, !this->make_fixedRegions_templates);}

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
                            } //private SMEFT samples
                            else //Central SM samples, or systematics for private SMEFT samples --> Use regular TH1F objects
                            {
                                if(sample_list[isample] != "NPL_MC" && sample_list[isample] != "DATA") //NB: NPL_MC *should* have negative bins (substraction); and allow empty data
                                {
                                    //-- Nominal: prevent negative bins, and error>content
                                    if(syst_list[isyst] == "" && systTree_list[itree] == "")
                                    {
                                        Avoid_Histogram_EmptyOrNegativeBins(v3_histo_chan_syst_var[ichan][isyst][ivar]);

                                        v_storeNominalBinContents_chan_var_bin[ichan][ivar].resize(v3_histo_chan_syst_var[ichan][isyst][ivar]->GetNbinsX());
                                        v_storeNominalBinErrors_chan_var_bin[ichan][ivar].resize(v3_histo_chan_syst_var[ichan][isyst][ivar]->GetNbinsX());
                                        for(int ibin=0; ibin<v3_histo_chan_syst_var[ichan][isyst][ivar]->GetNbinsX(); ibin++)
                                        {
                                            v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin] = v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1);
                                            v_storeNominalBinErrors_chan_var_bin[ichan][ivar][ibin] = v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinError(ibin+1);
                                            // cout<<"v_storeNominalBinContents_chan_var_bin["<<ichan<<"]["<<ivar<<"]["<<ibin<<"] "<<v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin]<<endl;
                                            // cout<<"v_storeNominalBinErrors_chan_var_bin["<<ichan<<"]["<<ivar<<"]["<<ibin<<"] "<<v_storeNominalBinErrors_chan_var_bin[ichan][ivar][ibin]<<endl;
                                        }
                                    }
                                    //-- Others: if bin<0, or nominal~0, or (nom-var)<0, set bin to nominal; also adapt error
                                    else
                                    {
                                        for(int ibin=0; ibin<v3_histo_chan_syst_var[ichan][isyst][ivar]->GetNbinsX(); ibin++)
                                        {
                                            //-- bin<0, or nominal already set to ~0 --> set variation to nominal //Don't call 'Avoid_Histogram_EmptyOrNegativeBins' because we need to check case-by-case depending on the original bin contents
                                            if(v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1) <= 0 || v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin] < pow(10, -4))
                                            {
                                                v3_histo_chan_syst_var[ichan][isyst][ivar]->SetBinContent(ibin+1, v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin]);
                                                v3_histo_chan_syst_var[ichan][isyst][ivar]->SetBinError(ibin+1, pow(10,-6)); //Don't need MCstat error for shapeSysts
                                            }

                                            //-- Bin error > content --> set error=content
                                            // if(v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinError(ibin+1) > v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1) / 2.)
                                            if(v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinError(ibin+1) > v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1))
                                            {
                                                v3_histo_chan_syst_var[ichan][isyst][ivar]->SetBinError(ibin+1, v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1)-pow(10, -3));
                                                // v3_histo_chan_syst_var[ichan][isyst][ivar]->SetBinError(ibin+1, v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1)/2.);
                                            }

                                            //-- Careful ! Make sure not to remove any genuine shape effect here... (should be safe)
                                            //-- var/nom > truncation_factor or < (1/truncation_factor) --> indicates very low stat, and can cause problems when splitting per bin --> set var==nom //cf. https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/1086/1.html //Again verified that even when got up/nom>5, this led to unphysical postfit uncertainties in a bin !
                                            float truncation_factor = 2.;
                                            if(v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1) / v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin] > truncation_factor || v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1) / v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin] < 1/truncation_factor)
                                            // if(v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1) / v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin] > 10 || v3_histo_chan_syst_var[ichan][isyst][ivar]->GetBinContent(ibin+1) / v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin] < 0.1)
                                            {
                                                v3_histo_chan_syst_var[ichan][isyst][ivar]->SetBinContent(ibin+1, v_storeNominalBinContents_chan_var_bin[ichan][ivar][ibin]);
                                                v3_histo_chan_syst_var[ichan][isyst][ivar]->SetBinError(ibin+1, pow(10,-6)); //Don't need MCstat error for shapeSysts
                                            }
                                        }
                                    }
                                } //end protections for combine

                                //-- Normalize specific systematics (to nominal, before of after selection)
                                if(systTree_list[itree] == "" && syst_list[isyst] == "") {v_storeNominalIntegral_chan_var[ichan][ivar] = v3_histo_chan_syst_var[ichan][isyst][ivar]->Integral();} //Store nominal norm.
                                else if(syst_list[isyst].BeginsWith("PDF") || syst_list[isyst].BeginsWith("ME") || syst_list[isyst].BeginsWith("alpha") || syst_list[isyst].BeginsWith("ISR") || syst_list[isyst].BeginsWith("FSR")) //Apply SF corresponding to (var/nom) *before selection* --> Remaining norm. effect is from acceptance
                                {
                                    // cout<<"Before: v3_histo_chan_syst_var[ichan][isyst][ivar]->Integral() "<<v3_histo_chan_syst_var[ichan][isyst][ivar]->Integral()<<endl;
                                    Scale_THSyst_toBeforeSelection(v3_histo_chan_syst_var[ichan][isyst][ivar], v_SWE_beforeSel, syst_list[isyst]);
                                    // cout<<"After: v3_histo_chan_syst_var[ichan][isyst][ivar]->Integral() "<<v3_histo_chan_syst_var[ichan][isyst][ivar]->Integral()<<endl;
                                    // Scale_THSyst_toBeforeSelection(v3_histo_chan_syst_var[ichan][isyst][ivar], h_sumWeights_beforeSel, syst_list[isyst]);
                                    // v3_histo_chan_syst_var[ichan][isyst][ivar]->Scale(v_storeNominalIntegral_chan_var[ichan][ivar]/v3_histo_chan_syst_var[ichan][isyst][ivar]->Integral());
                                }
                                else if(syst_list[isyst].BeginsWith("njets_tZq")) //Special case: rescale syst to nominal norm.
                                {
                                    v3_histo_chan_syst_var[ichan][isyst][ivar]->Scale(v_storeNominalIntegral_chan_var[ichan][ivar]/v3_histo_chan_syst_var[ichan][isyst][ivar]->Integral());
                                }

                                file_output->cd();
                                v3_histo_chan_syst_var[ichan][isyst][ivar]->Write(output_histo_name);
                                nhistos++; //Count nof histos written to output file

                                // cout<<"Wrote histo : "<<output_histo_name<<" (Integral: "<<v3_histo_chan_syst_var[ichan][isyst][ivar]->Integral()<<")"<<endl;

                                //-- BIN SPLITTING //For MVA-EFT: need to store each histogram bin separately so that they can be scaled independently in Combine
                                if(this->split_EFTtemplates_perBin && !makeHisto_inputVars && make_SMvsEFT_templates_plots && !this->make_fixedRegions_templates) {StoreEachHistoBinIndividually(file_output, v3_histo_chan_syst_var[ichan][isyst][ivar], output_histo_name, nhistos, !this->make_fixedRegions_templates);}
                            } //central samples

                            //-- Delete TH1F/TH1EFT objects, if memory was allocated
    						if(v3_histo_chan_syst_var[ichan][isyst][ivar]) {delete v3_histo_chan_syst_var[ichan][isyst][ivar]; v3_histo_chan_syst_var[ichan][isyst][ivar] = NULL;}
                            if(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]) {delete v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]; v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = NULL;}
    					} //var loop
    				} //syst loop
    			} //chan loop

    			// cout<<"Done with "<<sample_list[isample]<<" sample"<<endl;

                if(v_wgts) {delete v_wgts; v_wgts = NULL;}
                if(v_ids) {delete v_ids; v_ids = NULL;}
                if(njets) {delete njets; njets = NULL;}
                if(channel) {delete channel; channel = NULL;}
                if(mTW) {delete mTW; mTW = NULL;} //Careful not to delete twice, can cause nasty segfaults
                if(eft_fit) {delete eft_fit; eft_fit = NULL;}

                file_input->Close(); file_input = NULL; //Also deletes associated TTrees, histos, etc.
            } //TTree loop
    	} //sample loop

        if(!makeHisto_inputVars)
        {
            if(classifier_name == "BDT")  {delete reader; reader = NULL;}
            else
            {
                if(clfy1 && use_specificMVA_eachYear) {delete clfy1; clfy1 = NULL; MVA_alreadyLoaded = false;} //Only delete if not going to be reused for another iteration
                //-- NB: do not delete clfy2 //Causes crash (because it uses the same TF session as clfy1, hence double-free... ?)
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
	cout<<FYEL("containing the "<<classifier_name<<" templates as histograms for : all samples / all channels")<<endl;
    cout<<DIM("("<<nhistos<<" histograms)")<<endl<<endl;
	file_output->Close(); file_output = NULL;
    if(make_njetstZqSF_file) {return;}

    //-- Can verify that the total nof processed entries computed from Count_Total_Nof_Entries() was effectively the nof processed entries
    // cout<<"total_nentries_toRead --> "<<total_nentries_toRead<<endl;
    // cout<<"tmp_compare --> "<<tmp_compare<<endl;

    //Restore potfile_outputentially modified variables
    classifier_name = restore_classifier_name;
    syst_list = restore_syst_list;
    systTree_list = restore_systTree_list;

    //NB: if get segfault when deleting a given array, may have allocated wrong size
    if(!makeHisto_inputVars || !noSysts_inputVars)
    {
        if(array_PU) {delete[] array_PU; array_PU = NULL;}
        if(array_prefiringWeight) {delete[] array_prefiringWeight; array_prefiringWeight = NULL;}
        if(array_Btag) {delete[] array_Btag; array_Btag = NULL;}
        if(array_jetPileupID) {delete[] array_jetPileupID; array_jetPileupID = NULL;}
        if(array_fakeFactor) {delete[] array_fakeFactor; array_fakeFactor = NULL;}
        if(array_ME) {delete[] array_ME; array_ME = NULL;}
        if(array_alphaS) {delete[] array_alphaS; array_alphaS = NULL;}
        if(array_PDFtotal) {delete[] array_PDFtotal; array_PDFtotal = NULL;}
        if(array_LepEffLoose_mu) {delete[] array_LepEffLoose_mu; array_LepEffLoose_mu = NULL;}
        if(array_partonShower) {delete[] array_partonShower; array_partonShower = NULL;}
        if(array_LepEffTight_mu) {delete[] array_LepEffTight_mu; array_LepEffTight_mu = NULL;}
        if(array_LepEffLoose_el) {delete[] array_LepEffLoose_el; array_LepEffLoose_el = NULL;}
        if(array_LepEffTight_el) {delete[] array_LepEffTight_el; array_LepEffTight_el = NULL;}
    }
    if(clfy1) {delete clfy1; clfy1 = NULL;}
    for(int ivar=0; ivar<var_list_pfloats.size(); ivar++) //Delete input vars
    {
        if(var_list_pfloats[ivar] && var_list[ivar] != "mTW" && var_list[ivar] != "njets" && var_list[ivar] != "channel") //Exception: mTW/njets/... are sepcial variables already deleted above //Ideally, use smart pointers (auto-delete shared mem, only once)
        {
            // cout<<"del var "<<var_list[ivar]<<" = "<<var_list_pfloats[ivar]<<endl;
            delete var_list_pfloats[ivar]; var_list_pfloats[ivar] = NULL;
        }
    }
    for(int ivar=0; ivar<var_list_pfloats_2.size(); ivar++)
    {
        if(var_list_pfloats_2[ivar] && v_varIndices_inMVA1[ivar]==-1 && var_list_NN2[ivar] != "mTW" && var_list_NN2[ivar] != "njets" && var_list_NN2[ivar] != "channel") //Be careful not to delete any previously-deleted variable ! //Ideally, use smart pointers (auto-delete shared mem, only once)
        {
            /*cout<<"del var "<<var_list_NN2[ivar]<<endl;*/
            delete var_list_pfloats_2[ivar]; var_list_pfloats_2[ivar] = NULL;
        }
    }

    //-- Not needed ? //Careful not to delete any variable twice
    for(int ivar=0; ivar<total_var_pfloats.size(); ivar++)
    {
        if(!total_var_list[ivar].Contains("mTW") && !total_var_list[ivar].Contains("njets") && !total_var_list[ivar].Contains("channel"))
        {
            /*cout<<"del var "<<total_var_list[ivar]<<endl;*/
            delete total_var_pfloats[ivar]; total_var_pfloats[ivar] = NULL;
        }
    }

    //-- Printout: count nof signal events in each bin //To understand/prevent low-stat issues
    if(!make_histos_forControlPlotsPaper && !make_njetstZqSF_file)
    {
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
            cout<<"year "<<v_lumiYears[iyear]<<endl;
            for(int isample=0; isample<sample_list.size(); isample++)
            {
                if(!sample_list[isample].Contains("PrivMC")) {continue;}
                cout<<"sample "<<sample_list[isample]<<endl;
                for(int ivar=0; ivar<total_var_list.size(); ivar++)
                {
                    cout<<"var "<<total_var_list[ivar]<<endl;
                    for(int ibin=0; ibin<v4_nEntries_year_proc_var_bin[iyear][isample][ivar].size(); ibin++)
                    {
                        cout<<"Bin "<<ibin+1<<" ---> Entries "<<v4_nEntries_year_proc_var_bin[iyear][isample][ivar][ibin]<<endl;
                    }
                }
            }
        }
    }

// #    # ###### #####   ####  ######
// ##  ## #      #    # #    # #
// # ## # #####  #    # #      #####
// #    # #      #####  #  ### #
// #    # #      #   #  #    # #
// #    # ###### #    #  ####  ######

    //-- For COMBINE fit, want to directly merge contributions from different processes into single histograms
    //-- For control histograms, only need to substract MC NPL from data-driven NPL
    if(!this->process_samples_byGroup) {MergeSplit_Templates(makeHisto_inputVars, output_file_name, total_var_list, template_name, region, true);} //If already considering group ntuples, no need to group anything more

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
void TopEFT_analysis::Draw_Templates(bool drawInputVars, TString channel, bool plot_onlyMaxNodeEvents, bool plot_onlyMVACutEvents, TString template_name, bool prefit, bool use_combine_file, TString EFTpoint, bool store_ymax_fixed, float* ymax_fixed1, float* ymax_fixed2)
{
//--- OPTIONS --------------------------------
//--------------------------------------------
    bool draw_errors = true; //true <-> superimpose error bands on plot/ratio plot

	bool draw_logarithm = false; //true <-> plot y-axis in log scale //May be overwriten for special cases below

    bool use_poisson_dataErrors = false; //true <-> use asymmetric poisson errors for data (recommended for low-stat templates), see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars //To be tested

    bool overlay_EFT_ratioPlot = false; //true <-> show (SM+EFT/SM) line(s) in ratio plot

    bool combineFile_fromHarvester = true; //true <-> use CombineHarvester conventions for combine file; else use FitDiagnostics conventions //NB: happened that data does not get drawn because of NaN in combine file (although GetBinContent loos normal...)

    bool superimpose_GENhisto = false; //true <-> superimpose corresponding GEN-level EFT histogram, for shape comparison... //Experimental

    //-- Hardcoded
    bool make_histos_forControlPlotsPaper = false; //true <-> plot the histograms produced when activating this option in Produce_Templates()

    bool superimpose_EFThist = true; //true <-> superimpose shape of EFT hists
        bool normalize_EFThist = false; //true <-> normalize EFT hists (arbitrary)
        bool superimpose_EFT_auto = true; //true <-> bypass vectors filled below, and automatically draw proc/WCs depending on region

        //-- HARDCODED BEST POSTFIT POINT
        TString bestfit_string = ""; //Hardcoded //String for rescaling the TH1EFT, corresponding to the bestfit point (for psotfit plots)
        // TString bestfit_string = "rwgt_ctz_-0.07_ctw_-0.01_cpq3_0.50_cpqm_-2.30_cpt_-15.69"; //Hardcoded //String for rescaling the TH1EFT, corresponding to the bestfit point (for psotfit plots)

        //-- Names of the private EFT samples to superimpose
        //NB: If empty, vectors are filled automatically depending on region, etc.
        vector<TString> v_EFT_samples;
        // v_EFT_samples.push_back("PrivMC_tZq");
        // v_EFT_samples.push_back("PrivMC_ttZ");
        // v_EFT_samples.push_back("PrivMC_tWZ");
        vector<TString> v_EFT_points; //Names of the EFT points at which to reweight the histos //Must follow naming convention used for private generation
        v_EFT_points.push_back("rwgt_ctz_0"); //<-> SM
        // v_EFT_points.push_back("rwgt_ctz_5");
        // v_EFT_points.push_back("rwgt_ctw_0");
        // v_EFT_points.push_back("rwgt_ctw_5");
        // v_EFT_points.push_back("rwgt_ctz_5");
        // v_EFT_points.push_back("rwgt_cpqm_10");
        // v_EFT_points.push_back("rwgt_cpt_5");
        // v_EFT_points.push_back("rwgt_cpt_-17");
        // v_EFT_points.push_back("rwgt_ctz_0_ctw_0_cpqm_0_cpq3_0_cpt_15");

        //-- AUTO:
        if(drawInputVars) {superimpose_EFT_auto = false;}
        if(superimpose_EFT_auto)
        {
            v_EFT_samples.clear();
            v_EFT_samples.push_back("PrivMC_tZq");
            v_EFT_samples.push_back("PrivMC_ttZ");

            v_EFT_points.clear();
            v_EFT_points.push_back("rwgt_ctz_0"); //SM
            v_EFT_points.push_back("rwgt_ctz_5"); //For NN-ctz

            v_EFT_points.push_back("rwgt_ctw_0"); //SM
            v_EFT_points.push_back("rwgt_ctw_5"); //For NN-ctw

            v_EFT_points.push_back("rwgt_cpq3_0"); //SM
            v_EFT_points.push_back("rwgt_cpq3_5"); //For NN-cpq3
        }

//--------------------------------------------
//--------------------------------------------

    if(template_name == "" && classifier_name != "BDT" && classifier_name != "NN") {cout<<BOLD(FRED("Error : classifier_name value ["<<classifier_name<<"] not supported !"))<<endl; return;}
    if(template_name == "") {template_name = classifier_name;}
    if(!drawInputVars && !make_SMvsEFT_templates_plots && categorization_strategy>0 && template_name!="BDT") {template_name = "NN_SM";} //Special case: if want to produce SM vs SM classifier plots, make sure we consider NN_SM templates (for BDT, must specify in main) //"NN"?
    if(this->make_fixedRegions_templates) {template_name = "";} //Irrelevant

    //-- For parameterized NN, can produce plots for *each EFT point* considered in the scan (adapt histo name accordingly, etc.)
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
		if(drawInputVars && !prefit) {cout<<"Error ! Can not draw postfit input vars yet !"<<endl; return;}
	}
    if(EFTpoint != "") {cout<<DIM("EFTpoint =  "<<EFTpoint<<"")<<endl;}

    //-- Obsolete (always only make templates/plots for private signal samples now)
    // if(make_SMvsEFT_templates_plots) {cout<<FBLU("[make_SMvsEFT_templates_plots==true] <-> will plot the private signal MC samples (not central)")<<endl;}
    // else {cout<<FRED("[make_SMvsEFT_templates_plots==false] <-> will plot the central signal MC samples (not private)")<<endl;}


//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

	//-- Can use 2 different types of combine files :
	//- the files containing the template histograms, produced with this code (-> only prefit plots)
	//- or, better, the file produced by Combine from the templates : contains the prefit distributions with total errors, and the postfit distribution
	//If want postfit plots, must use the Combine file. If want prefit plots, can use both of them (NB : errors will be different)

    // TString cat_tmp = (region=="") ? "SR" : region+"Cat";
    TString cat_tmp = region;

    bool use_predefined_EFT_strategy = false;
    if(!drawInputVars && this->categorization_strategy > 0) {use_predefined_EFT_strategy = true;}

    if(use_predefined_EFT_strategy && make_SMvsEFT_templates_plots && cat_tmp == "signal") {cat_tmp = "";} //Don't need this info in filename (default)
	if(!prefit) {use_combine_file = true;}
    if(!use_combine_file) {combineFile_fromHarvester = false;} //Necessary convention

    //-- Read input file (may be year-dependent)
    TFile* file_input = NULL;
    TString template_type = this->make_fixedRegions_templates? "otherRegions":template_name;
    TString inputFile_path = "";

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
	else //Templates
	{
        //-- Read the relevant NN info file, just to know if we are dealing with multiclass NN templates or not !
        if(template_name.Contains("NN") && !this->categorization_strategy && !use_combine_file)
        {
            TString NNinfo_input_path = Get_MVAFile_InputPath(template_name, signal_process, v_lumiYears[0], use_specificMVA_eachYear, make_SMvsEFT_templates_plots, true, this->categorization_strategy, this->scanOperators_paramNN, this->use_NN_cpq3_SRttZ);
            if(!Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN, v_NN_nodeLabels, NN_inputLayerName, NN_outputLayerName, NN_iMaxNode, NN_nNodes, minmax_bounds)) {return;} //Error: missing NN infos
        }

        //-- Can now fill the list of variables (templates names)
        Fill_Variables_List(total_var_list, use_predefined_EFT_strategy, template_name, this->region, this->scanOperators_paramNN, this->NN_nNodes, this->make_SMvsEFT_templates_plots, operator_scan1, operator_scan2, v_WCs_operator_scan1, v_WCs_operator_scan2, this->make_fixedRegions_templates, use_combine_file);

        //-- If reading a combine file with split bins, will need to build 'full' template from individual bins --> Must read info files for NN1 & NN2, so that we know about the ranges, etc. //But actually we can read this hardcoded info from Get_Template_Range()
        /*
        if(use_combine_file && split_EFTtemplates_perBin)
        {
            //-- Read info file for first NN
            TString MVA_input_path = Get_MVAFile_InputPath(total_var_list[0], "tZq", v_lumiYears[0], this->use_specificMVA_eachYear, make_SMvsEFT_templates_plots, false, this->categorization_strategy);
            TString NNinfo_input_path = Get_MVAFile_InputPath(total_var_list[0], signal_process, v_lumiYears[0], use_specificMVA_eachYear, make_SMvsEFT_templates_plots, true, this->categorization_strategy);
            var_list_NN.resize(0); NN_iMaxNode = -1; NN_strategy = ""; NN_inputLayerName = ""; NN_outputLayerName = ""; NN_nNodes = -1; minmax_bounds.clear();
            if(!Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN, v_NN_nodeLabels, NN_inputLayerName, NN_outputLayerName, NN_iMaxNode, NN_nNodes, minmax_bounds)) {return;} //Error: missing NN infos

            //-- Read info file for second NN //NB: it is assumed that the NN strategy (e.g. parametrized NN) is the same as for the first NN !
            MVA_input_path = Get_MVAFile_InputPath(total_var_list[1], "ttZ", v_lumiYears[0], this->use_specificMVA_eachYear, make_SMvsEFT_templates_plots, false, this->categorization_strategy);
            NNinfo_input_path = Get_MVAFile_InputPath(total_var_list[1], "ttZ", v_lumiYears[0], this->use_specificMVA_eachYear, make_SMvsEFT_templates_plots, true, this->categorization_strategy);
            var_list_NN2.resize(0); NN2_iMaxNode = -1; NN2_strategy = ""; NN2_inputLayerName = ""; NN2_outputLayerName = ""; NN2_nNodes = -1; minmax_bounds2.clear();
            if(Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN2, v_NN2_nodeLabels, NN2_inputLayerName, NN2_outputLayerName, NN2_iMaxNode, NN2_nNodes, minmax_bounds2)) {clfy2 = new TFModel(MVA_input_path.Data(), var_list_NN2.size(), NN2_inputLayerName.Data(), NN2_nNodes, NN2_outputLayerName.Data());} //Load neural network model
        }*/
	} //Templates
    if(make_histos_forControlPlotsPaper)
    {
        inputFile_path = "./outputs/ControlPlotsPaper.root";
        use_combine_file = false;

        //-- Hardcodeed list of variables
        total_var_list.clear();
        total_var_list.push_back("njets");
        total_var_list.push_back("nbjets");
        total_var_list.push_back("jPrimeAbsEta");
        total_var_list.push_back("jprime_Pt");
        total_var_list.push_back("metEt");
        total_var_list.push_back("mTW");
        total_var_list.push_back("lAsymmetry");
        total_var_list.push_back("recoZ_Pt");
        total_var_list.push_back("recoZ_dPhill");
        total_var_list.push_back("dR_tZ");
        total_var_list.push_back("maxDeepJet");
    }
    else if(!use_combine_file)
    {
        inputFile_path = Get_HistoFile_InputPath(!drawInputVars, template_type, cat_tmp, lumiName, use_combine_file, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, (EFTpoint!=""), combineFile_fromHarvester, prefit);
        if(inputFile_path == "") {cat_tmp = signal_process; inputFile_path = Get_HistoFile_InputPath(!drawInputVars, template_type, cat_tmp, lumiName, use_combine_file, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, (EFTpoint!=""), combineFile_fromHarvester, prefit);} //Retry with 'signal_process' as 'region' argument
    }
    if(inputFile_path == "") {cout<<"Get_HistoFile_InputPath --> file not found ! "<<endl; return;}
    else {file_input = TFile::Open(inputFile_path, "READ");}

    //Special case: for NN_5D, may want to superimpose all operators
    if(superimpose_EFT_auto && (total_var_list[0].Contains("NN_all") || total_var_list[0].Contains("NN_5D")))
    {
        v_EFT_points.clear();
        v_EFT_points.push_back("rwgt_ctz_0"); //SM

        if(!prefit && use_combine_file && bestfit_string != "") {v_EFT_points.push_back(bestfit_string);} //For postfit plot, should shown signal lines for SM and bestfit (bestfit point must have been hardcoded)
        else //For NN-5D, show signal lines for 3 relevant WCs
        {
            v_EFT_points.push_back("rwgt_ctz_5");
            v_EFT_points.push_back("rwgt_ctw_5");
            v_EFT_points.push_back("rwgt_cpq3_5");
        }
    }

    //-- TH1EFT not stored in combine output file; if need TH1EFT (to plot EFT prediction), need to read the template file given to combine for the fit
    TString inputFile_pathEFT = "";
    TFile* f_EFT = NULL;
    if(use_combine_file && superimpose_EFThist)
    {
        inputFile_pathEFT = Get_HistoFile_InputPath(!drawInputVars, template_type, cat_tmp, lumiName, false, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, (EFTpoint!=""), combineFile_fromHarvester, prefit);
        if(!Check_File_Existence(inputFile_pathEFT) ) {cout<<endl<<FRED("File "<<inputFile_pathEFT<<" not found! Can't superimpose EFT predictions !")<<endl; superimpose_EFThist = false;}
        else {f_EFT = TFile::Open(inputFile_pathEFT);}
    }

    //-- Define ranges of jet/bjets multiplicities -- for 'categ' templates only (modified in 'Get_Template_Range()')
    int nbjets_min = 1, nbjets_max=2, njets_min=2, njets_max=6;
    int nbins_tmp; float xmin_tmp, xmax_tmp;

    vector<float> v_yields_processes(sample_list.size()); //Can store and print yields per process, summed over years (e.g. to compare prefit/postfit)


// #       ####   ####  #####   ####
// #      #    # #    # #    # #
// #      #    # #    # #    #  ####
// #      #    # #    # #####       #
// #      #    # #    # #      #    #
// ######  ####   ####  #       ####

	for(int ivar=0; ivar<total_var_list.size(); ivar++)
	{
		if(drawInputVars) {cout<<endl<<FBLU("* Variable : "<<total_var_list[ivar]<<" ")<<endl<<endl;}

        if(this->make_SMvsEFT_templates_plots && !drawInputVars && !drawInputVars)
        {
            if(total_var_list[ivar].Contains("SRt") && !total_var_list[ivar].Contains("NN_SM") && !total_var_list[ivar].Contains("Zpt") && !total_var_list[ivar].Contains("countExp") && (!total_var_list[ivar].Contains("NN_cpq3_SRttZ") || use_NN_cpq3_SRttZ)) {draw_logarithm = true;} //Force log y-scale
            else {draw_logarithm = false;}
        }
        // draw_logarithm = false; //Force use of linear y-scale

        TH1F* h_tmp = NULL; //Tmp storing histo
        TH1F* hdata_tmp = NULL; //Tmp storing data histo
		TH1F* h_tzq = NULL; //Store tZq shape
		TH1F* h_ttz = NULL; //Store ttZ shape
		TH1F* h_sum_data = NULL; //Will store data histogram
		vector<TH1F*> v_MC_histo; //Will store all MC histograms (1 TH1F* per MC sample)
        TString data_histo_name = "";

		TGraphAsymmErrors* g_data = NULL; //If using Combine file, data are stored in TGAE
		TGraphAsymmErrors* g_tmp = NULL; //Tmp storing graph

		vector<TString> MC_samples_legend; //List the MC samples to mention in legend

		//-- Init error vectors
		double x, y, errory_low, errory_high;

		vector<double> v_eyl, v_eyh, v_exl, v_exh, v_x, v_y; //Contain the systematic errors (used to create the TGraphError)
		int nofbins=-1;

        vector<vector<TH1EFT*>> v2_th1eft(v_EFT_points.size()); //Store TH1EFT objects (inner: samples, outer: EFT points)
        vector<vector<TString>> v2_th1eft_labels(v_EFT_points.size()); //Corresponding names for legend (inner: samples, outer: EFT points)

        for(int ipt=0; ipt<v2_th1eft.size(); ipt++) {v2_th1eft[ipt].resize(v_EFT_samples.size()); v2_th1eft_labels[ipt].resize(v_EFT_samples.size());} //Init each element as a vector

        //-- Combine file: histos stored in subdirs -- define dir name
        TString dir_hist_prefix = "";

        //-- Only for FitDiagnostics output (not from CombineHarvester) ?
        if(!combineFile_fromHarvester)
        {
            if(prefit) {dir_hist_prefix = "shapes_prefit/";}
            else {dir_hist_prefix = "shapes_fit_s/";}
        }

        TString var_tmp = ""; //Dummy variable following expected naming convention --> check if exists
        //-- Combine file: if divided full histos into individual bins for fit, will now need to re-combine all the bins into full templates
        //-- First: check contents of the file to see whether templates were split or not
        bool combineIndividualBins=false; //Auto-set below
        int nIndivBins = 1; //Will update this number with the actual nof bins found in histos
        if(use_combine_file)
        {
            var_tmp = total_var_list[ivar];
            if(channel_list[0] != "") {var_tmp+= "_" + channel_list[0];}
            if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {var_tmp+= "_" + region;}
            // var_tmp+= "_" + v_lumiYears[0]; //Dummy year

            // cout<<dir_hist_prefix<<"bin"+Convert_Number_To_TString(nIndivBins)+"_"<<var_tmp<<endl;
            TString hpath_tmp = dir_hist_prefix; //change naming convention?
            if(combineFile_fromHarvester)
            {
                if(prefit) hpath_tmp+= "_prefit";
                else hpath_tmp+= "_postfit";
            }
            hpath_tmp+= "bin1_" + var_tmp;

            inputFile_path = Get_HistoFile_InputPath(!drawInputVars, "bin1_" + var_tmp, cat_tmp, lumiName, use_combine_file, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, (EFTpoint!=""), combineFile_fromHarvester, prefit, false);
            // if(inputFile_path == "") {cout<<"Get_HistoFile_InputPath --> file not found ! "<<endl; return;}
            if(inputFile_path != "") {combineIndividualBins = true;} //Seems like we have fitted individual bins, so in this chode we'll need to combine them all together again //Look for dummy bin (must be present if templates are split per bins)
            else
            {
                inputFile_path = Get_HistoFile_InputPath(!drawInputVars, var_tmp, cat_tmp, lumiName, use_combine_file, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, (EFTpoint!=""), combineFile_fromHarvester, prefit, false);
                if(inputFile_path == "") {cout<<"Get_HistoFile_InputPath --> file not found ! "<<endl; continue;}
            }

            file_input = TFile::Open(inputFile_path, "READ");

            if(combineIndividualBins) {Get_Template_Range(nIndivBins, xmin_tmp, xmax_tmp, total_var_list[ivar], this->make_SMvsEFT_templates_plots, this->categorization_strategy, plot_onlyMaxNodeEvents, nbjets_min, nbjets_max, njets_min, njets_max, minmax_bounds, use_NN_SRother, use_NN_cpq3_SRttZ);}
            // cout<<"nIndivBins "<<nIndivBins<<endl;
            // if(file_input->GetDirectory(hpath_tmp)) {combineIndividualBins = true;} //Seems like we have fitted individual bins, so in this chode we'll need to combine them all together again //Look for dummy bin (must be present if templates are split per bins)
            // if(combineIndividualBins) //Now, need to infer from the file how many individual bins will need to be combined
            // {
            //     while(file_input->GetDirectory(hpath_tmp))
            //     {
            //         nIndivBins++;
            //         hpath_tmp = dir_hist_prefix;
            //         if(combineFile_fromHarvester)
            //         {
            //             if(prefit) hpath_tmp+= "_prefit";
            //             else hpath_tmp+= "_postfit";
            //         }
            //         hpath_tmp+= "bin"+Convert_Number_To_TString(nIndivBins)+"_" + var_tmp;
            //     }
            //     nIndivBins--; //Last bin was not found --> update total nof bins
            //     cout<<"nIndivBins "<<nIndivBins<<endl;
            // }
        }

        float bin_width = -1; //Get bin width of histograms for current variable

        bool data_notEmpty = true; //Default=true <-> data found //NB: if blinded, will set all data bins to 0 manually below

        //-- All histos are for given lumiYears and sub-channels --> Need to sum them all for plots
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
        	//Need to rescale signal to fitted signal strength manually, and add its error in quadrature in each bin (verify)
        	double sigStrength = 0;
        	double sigStrength_Error = 0;
        	if(!combineFile_fromHarvester && !prefit)
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

    			//-- Combine file : histos stored in subdirs -- define dir name //Only used if dealing with 'full' templates
    			TString dir_hist = dir_hist_prefix + total_var_list[ivar];
                if(channel_list[ichan] != "") {dir_hist+= "_" + channel_list[ichan];}
                if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {dir_hist+= "_" + region;}
                dir_hist+= "_" + v_lumiYears[iyear];
                if(use_combine_file && combineFile_fromHarvester)
                {
                    if(prefit) {dir_hist+= "_prefit";}
                    else {dir_hist+= "_postfit";}
                }
                dir_hist+= "/";
                if(use_combine_file && !combineIndividualBins && !file_input->GetDirectory(dir_hist) ) {cout<<ITAL(DIM("ERROR: directory '"<<dir_hist<<"' : not found !"))<<endl; continue;}


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

                    // cout<<"sample_list[isample] "<<sample_list[isample]<<endl;
                    // cout<<"index_MC_sample "<<index_MC_sample<<endl;

    				//-- In Combine, some individual contributions are merged as "Rares"/"EWK", etc.
    				//-- If using Combine file, change the names of the samples we look for, and look only once for histogram of each "group"
    				TString samplename = sample_list[isample];
    				if(use_combine_file)
    				{
    					if(isample > 0 && sample_groups[isample] == sample_groups[isample-1]) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //if same group as previous sample, skip it
                        else if(make_SMvsEFT_templates_plots && (sample_groups[isample] == "tZq" || sample_groups[isample] == "ttZ" || sample_groups[isample] == "tWZ")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //SM vs EFT --> use private signal samples
                        else {samplename = sample_groups[isample];}
    				}

    				//-- Protections, special cases
    				if(sample_list[isample] == "DATA") {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}
                    // else if(!make_SMvsEFT_templates_plots && sample_list[isample].Contains("PrivMC")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //SM configuration --> only stack central samples (not private samples)
                    // else if(make_SMvsEFT_templates_plots && (sample_list[isample] == "tZq" || sample_list[isample] == "ttZ" || sample_list[isample] == "tWZ")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //EFT configuration --> only stack private samples (at SM point), not central samples
                    else if(sample_list[isample] == "tZq" || sample_list[isample] == "ttZ" || sample_list[isample] == "tWZ") {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //Default: only stack private samples (at SM point), not central samples
                    else if(sample_list[isample] == "NPL_DATA")  {samplename = "NPL";} //Instead of 'NPL_DATA' and 'NPL_MC', we only want to read the merged histo 'NPL'
                    else if(sample_list[isample] == "NPL_MC")  {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //NPL_MC gets substracted from NPL histograms and deleted --> Ignore this vector element //Remove ?

                    //-- Add sample name to list (used for legend) //NB: add even if histo was not found and skipped, because expect that it will be found for some other year/channel/... But if not found at all, legend will be wrong
                    if(iyear==0 && samplename != "DATA")
                    {
                        if(v_MC_histo.size() <=  index_MC_sample) {MC_samples_legend.push_back(samplename);}
                        // cout<<"v_MC_histo.size() "<<v_MC_histo.size()<<" / index_MC_sample "<<index_MC_sample<<" / MC_samples_legend.size() "<<MC_samples_legend.size()<<endl;
                    }
                    if(v_isSkippedSample[isample] == true) {continue;} //Skip this sample

    				// cout<<endl<<UNDL(FBLU("-- Sample : "<<sample_list[isample]<<" : "))<<endl;

    				h_tmp = NULL;
    				TString histo_name = "";
                    if(!drawInputVars && use_combine_file && combineIndividualBins) //Build 'full' template from individual bins
                    {
                        if(total_var_list[ivar].Contains("SRtZq")) {Get_Template_Range(nbins_tmp, xmin_tmp, xmax_tmp, total_var_list[ivar], this->make_SMvsEFT_templates_plots, this->categorization_strategy, plot_onlyMaxNodeEvents, nbjets_min, nbjets_max, njets_min, njets_max, minmax_bounds, use_NN_SRother, use_NN_cpq3_SRttZ);}
                        else {Get_Template_Range(nbins_tmp, xmin_tmp, xmax_tmp, total_var_list[ivar], this->make_SMvsEFT_templates_plots, this->categorization_strategy, plot_onlyMaxNodeEvents, nbjets_min, nbjets_max, njets_min, njets_max, minmax_bounds2, use_NN_SRother, use_NN_cpq3_SRttZ);}

                        if(template_name.Contains("NN") && make_SMvsEFT_templates_plots) {xmin_tmp = 0; xmax_tmp = nIndivBins;} //NN template range may have been auto-adapted based on the NN info (to avoid empty bins at boundaries), or it may have been 'discarded' voluntarily if we split the template into individual bins; since this info can't be propagated through Combine fit, we can't infer the initial range here --> Need to hardcode it case-by-case !

                        h_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
                        for(int ibin=1; ibin<nIndivBins+1; ibin++)
                        {
                            inputFile_path = Get_HistoFile_InputPath(!drawInputVars, "bin"+Convert_Number_To_TString(ibin)+"_" + var_tmp, cat_tmp, lumiName, use_combine_file, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, (EFTpoint!=""), combineFile_fromHarvester, prefit, false);
                            if(inputFile_path == "") {cout<<"Get_HistoFile_InputPath --> file not found ! "<<endl; continue;}
                            file_input = TFile::Open(inputFile_path, "READ");

                            TString dir_hist_tmp = dir_hist_prefix + "bin" + Convert_Number_To_TString(ibin) + "_" + total_var_list[ivar];
                            if(channel_list[ichan] != "") {dir_hist_tmp+= "_" + channel_list[ichan];}
                            if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {dir_hist_tmp+= "_" + region;}
                            dir_hist_tmp+= "_" + v_lumiYears[iyear];
                            if(combineFile_fromHarvester)
                            {
                                if(prefit) dir_hist_tmp+= "_prefit";
                                else dir_hist_tmp+= "_postfit";
                            }
                            dir_hist_tmp+= "/";

                            // cout<<"dir_hist_tmp/samplename "<<dir_hist_tmp<<samplename<<endl;
                            if(!file_input->GetDirectory(dir_hist_tmp) || !file_input->GetDirectory(dir_hist_tmp)->GetListOfKeys()->Contains(samplename) ) {cout<<FRED("Directory '"<<dir_hist_tmp<<"' or histogram '"<<dir_hist_tmp<<samplename<<"' not found ! Skip...")<<endl; continue;}

                            h_tmp->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist_tmp+samplename))->GetBinContent(1)); //Get content/error from individual bin
                            h_tmp->SetBinError(ibin, ((TH1F*) file_input->Get(dir_hist_tmp+samplename))->GetBinError(1));

                            // cout<<"ibin "<<ibin<<endl;
                            // cout<<"h_tmp->GetBinContent("<<ibin<<") "<<h_tmp->GetBinContent(ibin)<<endl;
                            // cout<<"h_tmp->GetBinError(1) "<<h_tmp->GetBinError(1)<<endl;

                            file_input->Close();
                        }
                    }
                    else //Get 'full' templates directly
                    {
                        if(use_combine_file) {histo_name = dir_hist + samplename;}
                        else
                        {
                            histo_name = total_var_list[ivar];
                            if(channel != "") {histo_name+= "_" + channel;}
                            if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {histo_name+= "_" + region;}
                            if(EFTpoint != "") {histo_name+= "_" + EFTpoint;}
                            histo_name+= "_" + v_lumiYears[iyear];
                            histo_name+= "__" + samplename;
                        }

                        // cout<<"histo_name "<<histo_name<<endl;

                        // Changed -- even if histo not found, still fill dummy vector element (may be just 1 year missing, still need to account for this sample in indices, etc.)
                        if((use_combine_file && !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(samplename)) || (!use_combine_file && !file_input->GetListOfKeys()->Contains(histo_name)) )
                        {
                            if(v_MC_histo.size() <=  index_MC_sample) {v_MC_histo.push_back(NULL);}
                            cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found ! Skip..."))<<endl;
                            continue;
                        }
                        //-- Obsolete
                        // if(use_combine_file && !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(samplename) ) {cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found ! Skip..."))<<endl; v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}
                        // else if(!use_combine_file && !file_input->GetListOfKeys()->Contains(histo_name) ) {cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found ! Skip..."))<<endl; v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}

                        h_tmp = (TH1F*) file_input->Get(histo_name);
                        // cout<<"h_tmp->Integral() = "<<h_tmp->Integral()<<endl;
    				}

                    if(use_combine_file && total_var_list[ivar].Contains("countExp")) {h_tmp->Rebin(h_tmp->GetNbinsX());} //Trick: incorrect nof bins in file


                    if(bin_width < 0) {bin_width = h_tmp->GetXaxis()->GetBinWidth(1);}
    				if(draw_errors)
    				{
    					//Initialize error vectors (only once at start)
    					if(nofbins == -1) //if not yet init, get histo parameters
    					{
    						nofbins = h_tmp->GetNbinsX();
    						for(int ibin=0; ibin<nofbins; ibin++)
    						{
    							v_eyl.push_back(0); v_eyh.push_back(0);
                                v_exl.push_back(bin_width / 2); v_exh.push_back(h_tmp->GetXaxis()->GetBinWidth(ibin+1) / 2);
                                // v_exl.push_back(h_tmp->GetXaxis()->GetBinWidth(ibin+1) / 2); v_exh.push_back(h_tmp->GetXaxis()->GetBinWidth(ibin+1) / 2);
    							v_x.push_back( (h_tmp->GetXaxis()->GetBinLowEdge(nofbins+1) - h_tmp->GetXaxis()->GetBinLowEdge(1) ) * ((ibin+1 - 0.5)/nofbins) + h_tmp->GetXaxis()->GetBinLowEdge(1));
    							v_y.push_back(0);
    						}
    					}

    					//-- Add errors in quadrature
    					for(int ibin=0; ibin<nofbins; ibin++) //Start at bin 1
    					{
                            // cout<<"h_tmp->GetBinContent(ibin+1) "<<h_tmp->GetBinContent(ibin+1)<<endl;
                            v_y[ibin]+= h_tmp->GetBinContent(ibin+1); //This vector is used to know where to draw the error zone on plot (= on top of stack)

                            //-- Protection //Prevent drawing of unphysical uncertainties from postfit combine file (should fix the source of the problem... often due to threshold effects near 0)
                            if(abs(h_tmp->GetBinError(ibin+1)) > h_tmp->GetBinContent(ibin+1)*50)
                            {
                                cout<<FRED("Warning: error>content*50, skip this error (bin "<<ibin+1<<" / "<<sample_list[isample]<<" / "<<v_lumiYears[iyear]<<")")<<endl;
                                continue;
                            }

    						//-- If using Combine output file (from MLF), bin error contains total error. Else if using template file directly, just MC stat. error
                            //NB: this is actually incorrect because we can't sum correlated errors simply in quadrature ! Works if summing only MCstat errors, else (e.g. for total error from combine) must do differently
                            v_eyl[ibin]+= pow(h_tmp->GetBinError(ibin+1), 2);
    						v_eyh[ibin]+= pow(h_tmp->GetBinError(ibin+1), 2);

                            //-- Debug printouts for a single bin //NB: ibin starts at 0 for first bin (ibin=6 -> look at histo bin 7)
    						// if(ibin != 0) {continue;} //cout only 1 bin
                            // cout<<"sample_list[isample] "<<sample_list[isample]<<endl;
                            // cout<<"v_lumiYears[iyear] "<<v_lumiYears[iyear]<<endl;
    						// cout<<"x = "<<v_x[ibin]<<endl;    cout<<", y = "<<v_y[ibin]<<endl;    cout<<", sqrt(eyl) = "<<sqrt(v_eyl[ibin])<<endl;    cout<<", sqrt(eyh) = "<<sqrt(v_eyh[ibin])<<endl; //cout<<", exl = "<<v_exl[ibin]<<endl;    cout<<", exh = "<<v_exh[ibin]<<endl;
    					} //loop on bins

    					//-- Draw all errors
    					//--------------------------------------------
    					if(!use_combine_file) //In Combine file, already accounted in binError
    					{
    						for(int itree=0; itree<systTree_list.size(); itree++)
    						{
                                if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name && (sample_list[isample] == "DATA" || sample_list[isample] == "DY" || sample_list[isample].Contains("NPL") || sample_list[isample].Contains("TTbar")) ) {continue;}

    							for(int isyst=0; isyst<syst_list.size(); isyst++)
    							{
    								//-- Protections : not all syst weights apply to all samples, etc.
    								if(syst_list[isyst] != "" && systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name) {break;} //JES,JER,... -> read first element only
                                    else if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefir") ) {continue;} //no prefire in 2018
                                    else if((syst_list[isyst].BeginsWith("PDF") || syst_list[isyst].BeginsWith("ME") || syst_list[isyst].BeginsWith("alpha") || syst_list[isyst].BeginsWith("ISR") || syst_list[isyst].BeginsWith("FSR")) && !sample_list[isample].Contains("PrivMC") && sample_list[isample] != "tZq" && sample_list[isample] != "ttZ") {continue;}  //NB: no TH weights available in central tWZ V12
                                    // else if(syst_list[isyst].BeginsWith("alpha") && sample_list[isample] == "tWZ") {continue;} //Was bugged -- can remove ?

    								// cout<<"sample "<<sample_list[isample]<<" / channel "<<channel_list[ichan]<<" / syst "<<syst_list[isyst]<<endl;

    								TString histo_name_syst = histo_name + "__";
                                    if(syst_list[isyst] != "" || (systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name)) {histo_name_syst+= Get_Modified_SystName(syst_list[isyst]+systTree_list[itree], v_lumiYears[iyear], sample_list[isample]);}
                                    else {histo_name_syst+= systTree_list[itree];}

    								if(!file_input->GetListOfKeys()->Contains(histo_name_syst)) {continue;} //No error messages if systematics histos not found

                                    TH1F* histo_syst = NULL; //Read the "systematics histograms"
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

                                        //-- Debug printouts
                                        if(overlay_EFT_ratioPlot)
                                        {
                                            // if(!total_var_list[ivar].Contains("4l")) {continue;}

                                            cout<<"//--------------------------------------------"<<endl;
                                            cout<<"Sample "<<sample_list[isample]<<" / Syst "<<syst_list[isyst]<< " / chan "<<channel_list[ichan]<< " / year "<<v_lumiYears[iyear]<<endl;
                                            cout<<"x = "<<v_x[ibin]<<endl;    cout<<", y = "<<v_y[ibin]<<endl;    cout<<", eyl = "<<v_eyl[ibin]<<endl;    cout<<", eyh = "<<v_eyh[ibin]<<endl; //cout<<", exl = "<<v_exl[ibin]<<endl;    cout<<", exh = "<<v_exh[ibin]<<endl;
                                            cout<<"(nominal value = "<<h_tmp->GetBinContent(ibin+1)<<" - shifted value = "<<histo_syst->GetBinContent(ibin+1)<<") = "<<h_tmp->GetBinContent(ibin+1)-histo_syst->GetBinContent(ibin+1)<<endl;
                                        }
    								}

    								delete histo_syst;
    							} //end syst loop
    						} //systTree_list loop
    					} //use combine file?
    				} //draw errors?

    				//-- Set histo style (use color vector filled in main) //NB: run for all sub-histos... for simplicity
                    //---------------------------------------------------
    				h_tmp->SetFillStyle(1001);
    				if(samplename == "Fakes") {h_tmp->SetFillStyle(3005);}
    		        else if(samplename == "QFlip" ) {h_tmp->SetFillStyle(3006);}

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
                    if(region == "twz" && samplename == "tWZ") {h_tmp->SetFillColor(kBlue);} //Different color
                    //---------------------------------------------------

                    //-- Fill vector of MC histograms
                    if(v_MC_histo.size() <=  index_MC_sample) {v_MC_histo.push_back((TH1F*) h_tmp->Clone());}
                    else if(!v_MC_histo[index_MC_sample] && h_tmp) {v_MC_histo[index_MC_sample] = (TH1F*) h_tmp->Clone();}
                    else {v_MC_histo[index_MC_sample]->Add((TH1F*) h_tmp->Clone());}
                    if(v_MC_histo[index_MC_sample]) {v_MC_histo[index_MC_sample]->SetDirectory(0);} //Dis-associate histo from TFile //https://root.cern.ch/root/htmldoc/guides/users-guide/ObjectOwnership.html

                    //-- Debug printouts
    				// cout<<"sample : "<<sample_list[isample]<<" / color = "<<color_list[isample]<<" fillstyle = "<<h_tmp->GetFillStyle()<<endl;
    				// cout<<"index_MC_sample "<<index_MC_sample<<endl;
    				// cout<<"v_MC_histo.size() "<<v_MC_histo.size()<<endl;
    				// cout<<"MC_samples_legend.size() "<<MC_samples_legend.size()<<endl<<endl;

                    // if(sample_list[isample].Contains("PrivMC_tZq")) {cout<<"h_tmp->Integral() "<<h_tmp->Integral()<<endl;}
                    v_yields_processes[isample]+= h_tmp->Integral(); //For printouts

    				delete h_tmp; h_tmp = NULL; //No crash ? (else only delete if new)
    			} //end sample loop


// #####    ##   #####   ##
// #    #  #  #    #    #  #
// #    # #    #   #   #    #
// #    # ######   #   ######
// #    # #    #   #   #    #
// #####  #    #   #   #    #

                //-- NB: implemented such that data histos should be read even if blinded; in the latter case, the data is then manually set to 0
                if(use_combine_file)
                {
                    TString dataname = "data"; //FitDiag convention ?
                    if(combineFile_fromHarvester) {dataname = "data_obs";} //CH convention

                    if(combineIndividualBins) //Build 'full' template from individual bins
                    {
                        if(combineFile_fromHarvester) //CH conventions : data = TH1F
                        {
                            hdata_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
                            if(use_poisson_dataErrors) {hdata_tmp->SetBinErrorOption(TH1::kPoisson);}

                            for(int ibin=1; ibin<nIndivBins+1; ibin++)
                            {
                                inputFile_path = Get_HistoFile_InputPath(!drawInputVars, "bin"+Convert_Number_To_TString(ibin)+"_" + var_tmp, cat_tmp, lumiName, use_combine_file, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, (EFTpoint!=""), combineFile_fromHarvester, prefit, false);
                                if(inputFile_path == "") {cout<<"Get_HistoFile_InputPath --> file not found ! "<<endl; continue;}
                                file_input = TFile::Open(inputFile_path, "READ");

                                // cout<<"ibin "<<ibin<<endl;
                                TString dir_hist_tmp = dir_hist_prefix + "bin" + Convert_Number_To_TString(ibin) + "_" + total_var_list[ivar];
                                if(channel_list[ichan] != "") {dir_hist_tmp+= "_" + channel_list[ichan];}
                                if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {dir_hist_tmp+= "_" + region;}
                                dir_hist_tmp+= "_" + v_lumiYears[iyear];
                                if(combineFile_fromHarvester)
                                {
                                    if(prefit) dir_hist_tmp+= "_prefit";
                                    else dir_hist_tmp+= "_postfit";
                                }
                                dir_hist_tmp+= "/";

                                // cout<<"dir_hist_tmp+dataname= "<<dir_hist_tmp+dataname<<endl;
                                if(!file_input->GetDirectory(dir_hist_tmp) || !file_input->GetDirectory(dir_hist_tmp)->GetListOfKeys()->Contains(dataname) ) {cout<<FRED("Directory '"<<dir_hist_tmp<<"' or histogram '"<<dir_hist_tmp+dataname<<"' not found ! Skip...")<<endl; continue;}

                                float content = ((TH1F*) file_input->Get(dir_hist_tmp+dataname))->GetBinContent(1), error = ((TH1F*) file_input->Get(dir_hist_tmp+dataname))->GetBinError(1);
                                // cout<<"content "<<content<<endl;
                                if(!isnan(content) && !isinf(content))
                                {
                                    hdata_tmp->SetBinContent(ibin, content); //Get content/error from individual bin
                                    hdata_tmp->SetBinError(ibin, error);
                                }
                                else //NaN data //Seems to happen in some cases for bins with no data
                                {
                                    cout<<FRED("ERROR: bin "<<ibin<<" is NaN/inf ! ("<<dir_hist_tmp+dataname<<")")<<endl;
                                }
                                // cout<<"hdata_tmp->GetBinContent(1) "<<hdata_tmp->GetBinContent(1)<<endl;
                                // cout<<"hdata_tmp->GetBinError(1) "<<hdata_tmp->GetBinError(1)<<endl;

                                file_input->Close();
                            }
                        }
                        else //FitDiag conventions : data = TGraph
                        {
                            double theErrorX_h[nIndivBins];
                            double theErrorY_h[nIndivBins];
                            double theErrorX_l[nIndivBins];
                            double theErrorY_l[nIndivBins];
                            double theY[nIndivBins];
                            double theX[nIndivBins];

                            for(int ibin=1; ibin<nIndivBins+1; ibin++)
                            {
                                // cout<<"ibin "<<ibin<<endl;
                                TString dir_hist_tmp = dir_hist_prefix + "bin" + Convert_Number_To_TString(ibin) + "_" + total_var_list[ivar];
                                if(channel_list[ichan] != "") {dir_hist_tmp+= "_" + channel_list[ichan];}
                                if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {dir_hist_tmp+= "_" + region;}
                                dir_hist_tmp+= "_" + v_lumiYears[iyear];
                                if(combineFile_fromHarvester) {dir_hist_tmp+= "_prefit";}
                                dir_hist_tmp+= "/";

                                // cout<<"dir_hist_tmp+dataname = "<<dir_hist_tmp+dataname<<endl;
                                if(!file_input->GetDirectory(dir_hist_tmp) || !file_input->GetDirectory(dir_hist_tmp)->GetListOfKeys()->Contains(dataname) ) {cout<<FRED("Directory '"<<dir_hist_tmp<<"' or histogram '"<<dir_hist_tmp+dataname<<"' not found ! Skip...")<<endl; continue;}
                                TGraphAsymmErrors* g_tmp2 = (TGraphAsymmErrors*) file_input->Get(dir_hist_tmp+dataname);
                                Double_t x, y; g_tmp2->GetPoint(0, x, y); //Read y coordinate (x is irrelevant)

                                theX[ibin-1] = ibin - bin_width/2.; theY[ibin-1] = y; //Set x/y coordinates for this point
                                theErrorY_h[ibin-1] = g_tmp2->GetErrorYhigh(0); theErrorY_l[ibin-1] = g_tmp2->GetErrorYlow(0); //Set y-errors
                                theErrorX_h[ibin-1] = 0; theErrorX_l[ibin-1] = 0; //Irrelevant
                            }
                            g_tmp = new TGraphAsymmErrors(nofbins,theX,theY,theErrorX_l,theErrorX_h,theErrorY_l,theErrorY_h); //Build full TGraph
                        }
                    }
                    else //Get 'full' templates directly
                    {
                        if(combineFile_fromHarvester) //CH conventions : data = TH1F
                        {
                            TString dir_hist_tmp = total_var_list[ivar];
                            if(channel_list[ichan] != "") {dir_hist_tmp+= "_" + channel_list[ichan];}
                            if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {dir_hist_tmp+= "_" + region;}
                            dir_hist_tmp+= "_" + v_lumiYears[iyear];
                            if(prefit) {dir_hist_tmp+= "_prefit";}
                            else {dir_hist_tmp+= "_postfit";}
                            dir_hist_tmp+= "/";

                            // cout<<"dir_hist_tmp+dataname "<<dir_hist_tmp+dataname<<endl;
                            hdata_tmp = (TH1F*) file_input->Get(dir_hist_tmp+dataname)->Clone();
                            if(total_var_list[ivar].Contains("countExp")) {hdata_tmp->Rebin(hdata_tmp->GetNbinsX());} //Trick: incorrect nof bins in file

                            hdata_tmp->SetDirectory(0); //Dis-associate from TFile //https://root.cern.ch/root/htmldoc/guides/users-guide/ObjectOwnership.html
                        }
                        else
                        {
                            data_histo_name = dir_hist + dataname;
                            // cout<<"data_histo_name "<<data_histo_name<<endl;
                            g_tmp = (TGraphAsymmErrors*) file_input->Get(data_histo_name); //stored as TGraph
                        }
                    }

                    if(!combineFile_fromHarvester)
                    {
                        //-- Remove X-axis error bars, not needed for plot
                        for(int ipt=0; ipt<g_tmp->GetN(); ipt++)
                        {
                            g_tmp->SetPointEXhigh(ipt, 0);
                            g_tmp->SetPointEXlow(ipt, 0);
                        }

                        if(!g_data) {g_data = (TGraphAsymmErrors*) g_tmp->Clone();}
                        else //Need to sum TGraphs content by hand
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
                        if(g_tmp) {delete g_tmp; g_tmp = NULL;} //No crash ? (else only delete if new)
                    }
                    else
                    {
                        if(h_sum_data == NULL) {h_sum_data = (TH1F*) hdata_tmp->Clone();}
                        else {h_sum_data->Add((TH1F*) hdata_tmp->Clone());}

                        h_sum_data->SetDirectory(0); //Dis-associate from TFile
                        delete hdata_tmp; hdata_tmp = NULL; //No crash ? (else only delete if new)
                    }
        		}
                else //If using template file made from this code
                // else if(!this->is_blind)
                {
                    data_histo_name = total_var_list[ivar];
                    if(channel != "") {data_histo_name+= "_" + channel;}
                    if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {data_histo_name+= "_" + region;}
                    if(EFTpoint != "") {data_histo_name+= "_" + EFTpoint;}
                    data_histo_name+= "_" + v_lumiYears[iyear] + "__data_obs";

        			if(!file_input->GetListOfKeys()->Contains(data_histo_name)) {cout<<data_histo_name<<" : not found"<<endl;}
        			else
        			{
        				h_tmp = (TH1F*) file_input->Get(data_histo_name);
        				if(h_sum_data == NULL) {h_sum_data = (TH1F*) h_tmp->Clone();}
        				else {h_sum_data->Add((TH1F*) h_tmp->Clone());} //not needed anymore (1 channel only)
                        h_sum_data->SetDirectory(0); //Dis-associate from TFile
        			}
        		}

        		if(use_combine_file && !combineFile_fromHarvester && !g_data && !this->is_blind) {cout<<endl<<BOLD(FRED("--- Empty data TGraph !"))<<endl<<endl; data_notEmpty = false;}
        		if((!use_combine_file || combineFile_fromHarvester) && !h_sum_data && !this->is_blind) {cout<<endl<<BOLD(FRED("--- Empty data histogram "<<data_histo_name<<" !"))<<endl<<endl; data_notEmpty = false;}

        		//-- Make sure there are no negative bins
                //-- Also blind the data manually if needed here
        		if(data_notEmpty)
        		{
        			if(use_combine_file && !combineFile_fromHarvester) //TGraph object
        			{
        				for(int ipt=0; ipt<g_data->GetN(); ipt++)
        				{
        					g_data->GetPoint(ipt, x, y);
                            // cout<<"x "<<x<<" / y "<<y<<endl;
                            if(this->is_blind || y<0) {g_data->SetPoint(ipt, x, 0); g_data->SetPointError(ipt,0,0,0,0);}
                        }
        			}
        			else //TH1F object
        			{
        				for(int ibin=1; ibin<h_sum_data->GetNbinsX()+1; ibin++)
        				{
                            if(h_sum_data && (this->is_blind || h_sum_data->GetBinContent(ibin) < 0)) {h_sum_data->SetBinContent(ibin, 0); h_sum_data->SetBinError(ibin, 0);}
                        }
        			}

                    //-- Obsolete for MC (or at least should skip NPL_MC, expected to be negative !)
        			// for(int k=0; k<v_MC_histo.size(); k++)
        			// {
        			// 	if(!v_MC_histo[k]) {continue;} //Fakes templates can be null
        			// 	for(int ibin=0; ibin<v_MC_histo[k]->GetNbinsX(); ibin++)
        			// 	{
        			// 		if(v_MC_histo[k]->GetBinContent(ibin) < 0) {v_MC_histo[k]->SetBinContent(ibin, 0);}
        			// 	}
        			// }
        		}
                if(use_poisson_dataErrors) {h_sum_data->SetBinErrorOption(TH1::kPoisson);}


// #####  #####  # #    #         ####    ##   #    # #####  #      ######
// #    # #    # # #    #        #       #  #  ##  ## #    # #      #
// #    # #    # # #    #         ####  #    # # ## # #    # #      #####
// #####  #####  # #    # ###         # ###### #    # #####  #      #
// #      #   #  #  #  #  ###    #    # #    # #    # #      #      #
// #      #    # #   ##   ###     ####  #    # #    # #      ###### ######

                //-- Protection: if want to plot private SMEFT samples, make sure they are included in the main sample list
                //-- NB: if using Combine file, treat private SMEFT samples like central samples (included in stack); can't rescale & superimpose on plot, because Combine does not store TH1EFT objects !
                bool PrivMC_sample_found = false;
                for(int isample=0; isample<sample_list.size(); isample++) {if(sample_list[isample].Contains("PrivMC")) {PrivMC_sample_found = true;}}
                if(!PrivMC_sample_found) {superimpose_EFThist = false;} //Special case: if considering combine file for SM vs SM, private signals are not considered
                // if(!PrivMC_sample_found || use_combine_file) {superimpose_EFThist = false;} //Special case: if considering combine file for SM vs SM, private signals are not considered

                if(superimpose_EFThist)
                {
                    // int icolor = 4; //Use different colors for each histo //Start from blue
                    int icolor = 0; //Index
                    vector<int> v_col_eft;
                    v_col_eft.push_back(kBlack); //SM
                    v_col_eft.push_back(kAzure+6);
                    v_col_eft.push_back(kMagenta+1);
                    v_col_eft.push_back(kOrange+6);
                    v_col_eft.push_back(kGreen+2);

                    for(int ipoint=0; ipoint<v_EFT_points.size(); ipoint++)
                    {
                        for(int isample=0; isample<v_EFT_samples.size(); isample++)
                        {
                            if(superimpose_EFT_auto)
                            {
                                if(!total_var_list[ivar].Contains("NN_all") && !total_var_list[ivar].Contains("NN_5D"))
                                {
                                    if(!total_var_list[ivar].Contains("ctz") && !total_var_list[ivar].Contains("ctw") && (!total_var_list[ivar].Contains("cpq3") || (total_var_list[ivar].Contains("SRttZ") && !this->use_NN_cpq3_SRttZ))) {continue;}
                                    if(total_var_list[ivar].Contains("cpq3") && !v_EFT_points[ipoint].Contains("cpq3")) {continue;}
                                    if(total_var_list[ivar].Contains("ctz") && !v_EFT_points[ipoint].Contains("ctz")) {continue;}
                                    if(total_var_list[ivar].Contains("ctw") && !v_EFT_points[ipoint].Contains("ctw")) {continue;}
                                }
                                // if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3")) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                                if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3") && !this->use_NN_cpq3_SRttZ) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                                if(total_var_list[ivar].Contains("SRtZq") && !v_EFT_samples[isample].Contains("tZq")) {continue;}
                                if(total_var_list[ivar].Contains("SRttZ") && !v_EFT_samples[isample].Contains("ttZ")) {continue;}
                            }

                            // cout<<"isample "<<isample<<" / ipoint "<<ipoint<<endl;

                            TH1EFT* th1eft_tmp = NULL;

                            /* -- Irrelevant, but save for now
                            if(use_combine_file && combineIndividualBins) //Build 'full' template from individual bins
                            {
                                   th1eft_tmp = new TH1EFT("", "", nIndivBins, -9, -8);
                                   for(int ibin=1; ibin<nIndivBins+1; ibin++)
                                   {
                                       // cout<<"ibin "<<ibin<<endl;
                                       TString dir_hist_tmp = dir_hist_prefix + "bin" + Convert_Number_To_TString(ibin) + "_" + total_var_list[ivar];
                                       if(channel_list[ichan] != "") {dir_hist_tmp+= "_" + channel_list[ichan];}
                                       if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {dir_hist_tmp+= "_" + region;}
                                       dir_hist_tmp+= "_" + v_lumiYears[iyear] + "/";

                                       // cout<<"dir_hist_tmp/samplename "<<dir_hist_tmp<<v_EFT_samples[isample]<<endl;
                                       if(!file_input->GetDirectory(dir_hist_tmp) || !file_input->GetDirectory(dir_hist_tmp)->GetListOfKeys()->Contains(v_EFT_samples[isample]) ) {cout<<FRED("Directory '"<<dir_hist_tmp<<"' or histogram '"<<dir_hist_tmp<<v_EFT_samples[isample]<<"' not found ! Skip...")<<endl; continue;}
                                       th1eft_tmp->SetBinContent(ibin, ((TH1EFT*) file_input->Get(dir_hist_tmp+v_EFT_samples[isample]))->GetBinContent(ibin));
                                       th1eft_tmp->SetBinError(ibin, ((TH1EFT*) file_input->Get(dir_hist_tmp+v_EFT_samples[isample]))->GetBinError(ibin));
                                       // th1eft_tmp->Set_WCFit_Bin(ibin, ((TH1EFT*) file_input->Get(dir_hist_tmp+v_EFT_samples[isample]))->GetBinFit(ibin));

                                       // cout<<"th1eft_tmp->GetBinContent(1) "<<th1eft_tmp->GetBinContent(1)<<endl;
                                       // cout<<"th1eft_tmp->GetBinError(1) "<<th1eft_tmp->GetBinError(1)<<endl;
                                   }
                            }*/

                            TString histo_name = total_var_list[ivar];
                            if(channel != "") {histo_name+= "_" + channel;}
                            if(region != "" && !drawInputVars && !use_predefined_EFT_strategy) {histo_name+= "_" + region;}
                            if(EFTpoint != "") {histo_name+= "_" + EFTpoint;}
                            histo_name+= "_" + v_lumiYears[iyear];
                            histo_name+= "__" + v_EFT_samples[isample];
                            // cout<<"histo_name "<<histo_name<<endl;

                            TH1EFT* th1eft_tmp2 = NULL; //Special case: if reading different TFiles for th1s/th1eft objects, need to use same x-axis (taken from combine file)
                            if(use_combine_file)
                            {
                                if(!f_EFT->GetListOfKeys()->Contains(histo_name) ) {cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found ! Skip..."))<<endl; v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}

                                // th1eft_tmp = (TH1EFT*) v_MC_histo[0]->Clone();

                                th1eft_tmp2 = new TH1EFT("", "", v_MC_histo[0]->GetNbinsX(), v_MC_histo[0]->GetXaxis()->GetXmin(), v_MC_histo[0]->GetXaxis()->GetXmax());

                                th1eft_tmp = (TH1EFT*) f_EFT->Get(histo_name);

                                // v_yields_processes[isample]+= th1eft_tmp->Integral(); //For printouts //Not needed here, only for combine file ?

                                //-- Rescale TH1EFT accordingly to current reweight //Pay attention to operator exact names !
                                WCPoint wcp = WCPoint((string) v_EFT_points[ipoint], 1.);
                                th1eft_tmp->Scale(wcp);

                                for(int ibin=1; ibin<th1eft_tmp2->GetNbinsX()+1; ibin++)
                                {
                                    // cout<<"ibin "<<ibin<<" / content "<<th1eft_tmp->GetBinContent(ibin)<<" / error "<<th1eft_tmp->GetBinError(ibin)<<endl;
                                    th1eft_tmp2->SetBinContent(ibin, th1eft_tmp->GetBinContent(ibin));
                                    th1eft_tmp2->SetBinError(ibin, th1eft_tmp->GetBinError(ibin));
                                }
                                th1eft_tmp = th1eft_tmp2;
                            }
                            else
                            {
                                if(!file_input->GetListOfKeys()->Contains(histo_name) ) {cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found ! Skip..."))<<endl; v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}
                                th1eft_tmp = (TH1EFT*) file_input->Get(histo_name);

                                //-- Rescale TH1EFT accordingly to current reweight //Pay attention to operator exact names !
                                WCPoint wcp = WCPoint((string) v_EFT_points[ipoint], 1.);

                                if((total_var_list[ivar] == "NN_5D_SRttZ" || total_var_list[ivar] == "NN_cpq3_SRttZ") && v_EFT_points[ipoint] == "rwgt_cpq3_5") {wcp = WCPoint((string) "rwgt_cpq3_10", 1.);} //Hardcoded: use larger value cpq3=10 for NN-5D and NN-cpq3 in SRttZ

                                // cout<<"th1eft_tmp->Integral() "<<th1eft_tmp->Integral()<<endl;
                                th1eft_tmp->Scale(wcp);
                                // cout<<"th1eft_tmp->Integral() "<<th1eft_tmp->Integral()<<endl;
                            }
                            // cout<<"histo_name "<<histo_name<<endl;
                            // cout<<"th1eft_tmp->Integral() "<<th1eft_tmp->Integral()<<endl;

                            if(iyear == 0 || v2_th1eft[ipoint].size() <= isample || v2_th1eft[ipoint][isample] == 0) //New point or sample --> Add new elements
                            {
                                // v2_th1eft[ipoint][isample] = (TH1EFT*) th1eft_tmp->Clone();
                                v2_th1eft[ipoint][isample] = new TH1EFT; //Create object
                                v2_th1eft[ipoint][isample]->CloneTH1EFT(th1eft_tmp); //Use class function to copy tmp TH1EFT contents

                                //-- Get/save corresponding label
                                vector<pair<TString,float>> v = Parse_EFTreweight_ID(v_EFT_points[ipoint]);
                                if((total_var_list[ivar] == "NN_5D_SRttZ" || total_var_list[ivar] == "NN_cpq3_SRttZ") && v_EFT_points[ipoint] == "rwgt_cpq3_5") {v = Parse_EFTreweight_ID("rwgt_cpq3_10");} //Hardcoded: use larger value cpq3=10 for NN-5D and NN-cpq3 in SRttZ
                                TString EFTpointlabel = "";
                                for(int i=0; i<v.size(); i++)
                                {
                                    if(v[i].second != 0)
                                    {
                                        if(EFTpointlabel != "") {EFTpointlabel+= ",";}
                                        EFTpointlabel+= Get_EFToperator_label((TString) v[i].first) + "=" + Convert_Number_To_TString(v[i].second);
                                    }
                                }
                                if(v_EFT_points[ipoint] == bestfit_string) {EFTpointlabel = "best fit";}
                                if(EFTpointlabel == "") {EFTpointlabel = "SM";} //Ex: 'ctz_3' -> 'SM'
                                // cout<<"EFTpointlabel "<<EFTpointlabel<<endl;

                                std::vector<std::string> words;
                                split_string((string) v_EFT_samples[isample], words, "_");
                                TString leg_name = "tZq";
                                if(words.at(1) == "ttZ") {leg_name = "t#bar{t}Z";}
                                else if(words.at(1) == "tWZ") {leg_name = "tWZ";}
                                leg_name+= "("+EFTpointlabel+")"; //NB: can't have space before bracket, else too wide labels
                                // TString leg_name = "#splitline{"+words.at(1)+"}{("+EFTpointlabel+")}"; //Process name + EFT point //Split over 2 lines
                                v2_th1eft_labels[ipoint][isample] = leg_name;
                            }
                            else  //Sum different years together
                            {
                                v2_th1eft[ipoint][isample]->Add(th1eft_tmp);
                            }

                            //-- Style
                            v2_th1eft[ipoint][isample]->SetLineColor(v_col_eft[icolor]);
                            v2_th1eft[ipoint][isample]->SetLineWidth(4);
                            icolor++;
                            // th1eft_tmp->SetLineStyle(9); //Dashed lines
                            // th1eft_tmp->SetLineColor(icolor);
                            // icolor++; if(icolor==5) {icolor++;} //Update color and skip yellow
                            // th1eft_tmp->SetLineWidth(5);

                            v2_th1eft[ipoint][isample]->SetDirectory(0); //Dis-associate from TFile
                            if(th1eft_tmp) {delete th1eft_tmp; th1eft_tmp = NULL;} //Needed ?
                        } //sample loop
                    } //EFT points loop
                } //Superimpose private SMEFT histograms

            } //channels loop
        } //years loop


// # #    # #####  ###### #    #
// # ##   # #    # #       #  #
// # # #  # #    # #####    ##
// # #  # # #    # #        ##
// # #   ## #    # #       #  #
// # #    # #####  ###### #    #

        //-- Get indices of particular samples, sum the others into 1 single histo (used for ratio subplot)
        TH1F* histo_total_MC = NULL; //Sum of all MC samples

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
    		if(!v_MC_histo[i]) {continue;} //Some templates may be null

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

        //-- If histo_total_MC is null, variable was not found --> skip it
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

            if(region=="twz" && MC_samples_legend[i] == "tWZ") {continue;} //Put tWZ on top in that region

			stack_MC->Add(v_MC_histo[i]);
			// cout<<"Stacking sample "<<MC_samples_legend[i]<<" / integral "<<v_MC_histo[i]->Integral()<<endl;
            // cout<<"stack bin 1 content = "<<((TH1*) stack_MC->GetStack()->Last())->GetBinContent(1)<<endl;
		}
        if(region=="twz") {stack_MC->Add(v_MC_histo[index_tWZ_sample]);} //Put tWZ on top in that region

        //-- Debug printout (THStack bins contents)
        // for(int ibin=1; ibin<last->GetNbinsX()+1; ibin++) {cout<<"stack bin "<<ibin<<" content = "<<((TH1*) stack_MC->GetStack()->Last())->GetBinContent(ibin)<<endl;}


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
                    if(superimpose_EFT_auto)
                    {
                        if(!total_var_list[ivar].Contains("NN_all") && !total_var_list[ivar].Contains("NN_5D"))
                        {
                            if(!total_var_list[ivar].Contains("ctz") && !total_var_list[ivar].Contains("ctw") && (!total_var_list[ivar].Contains("cpq3") || (total_var_list[ivar].Contains("SRttZ") && !this->use_NN_cpq3_SRttZ))) {continue;}
                            if(total_var_list[ivar].Contains("cpq3") && !v_EFT_points[ipoint].Contains("cpq3")) {continue;}
                            if(total_var_list[ivar].Contains("ctz") && !v_EFT_points[ipoint].Contains("ctz")) {continue;}
                            if(total_var_list[ivar].Contains("ctw") && !v_EFT_points[ipoint].Contains("ctw")) {continue;}
                        }
                        // if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3")) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                        if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3") && !this->use_NN_cpq3_SRttZ) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                        if(total_var_list[ivar].Contains("SRtZq") && !v_EFT_samples[isample].Contains("tZq")) {continue;}
                        if(total_var_list[ivar].Contains("SRttZ") && !v_EFT_samples[isample].Contains("ttZ")) {continue;}
                    }

                    if(normalize_EFThist) {v2_th1eft[ipoint][isample]->Scale(((TH1*) stack_MC->GetStack()->Last())->Integral()/(2.*v2_th1eft[ipoint][isample]->Integral()));} //Normalize to half-integral of stack (arbitrary)  //NB: access stack integral via summed object '((TH1*) stack_MC->GetStack()->Last())'

                    //-- Debugging: print each bin's content under each hypothesis
                    // for(int ibin=1; ibin<v2_th1eft[ipoint][isample]->GetNbinsX()+1; ibin++)
                    // {
                    //     cout<<"Sample "<<v_EFT_samples[isample]<<" / Point "<<v_EFT_points[ipoint]<<" / Bin "<<ibin<<" / Content "<<v2_th1eft[ipoint][isample]->GetBinContent(ibin)<<endl;
                    // }
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

        //NB: for best style, better to adjust legend size exactly to the entries (too much space --> non uniform)

        if(superimpose_EFThist)
        {
            if(superimpose_EFT_auto)
            {
                // if(total_var_list[ivar].Contains("NN_all")) {nSampleGroups+= 6;}
                if(total_var_list[ivar].Contains("NN_all") || total_var_list[0].Contains("NN_5D")) {nSampleGroups+= 4;}
                else {nSampleGroups+= 2;}
            }
            else {nSampleGroups+= v_EFT_points.size() * v_EFT_samples.size();}
        }

        int n_columns = ceil(nSampleGroups/2.) > 6 ? 6 : ceil(nSampleGroups/2.); //ceil = upper int
        float x_left = 0.94-n_columns*0.12; //Each column allocated same x-space //0.12 needed for most crowded plots
        if(superimpose_EFThist && (total_var_list[ivar].Contains("SRtZq") || total_var_list[ivar].Contains("SRttZ"))) //Only if really adding EFTpoints to legend
        {
            n_columns-= 1;
            x_left = 0.94-n_columns*0.1; //Need more space per entry
        }
        if(x_left < 0.4) {x_left = 0.4;} //Leave some space for region label

        TLegend* qw = NULL;
        float ylegend = 0.78; //Default
        if(use_paperStyle) {ylegend+= 0.01;}
        qw = new TLegend(x_left-0.04,ylegend,0.94,ylegend+0.15); //Default
        qw->SetTextSize(0.035);
        qw->SetNColumns(n_columns);
        qw->SetBorderSize(0);
        qw->SetFillStyle(0); //transparent
        qw->SetTextAlign(12); //align = 10*HorizontalAlign + VerticalAlign //Horiz: 1=left adjusted, 2=centered, 3=right adjusted //Vert: 1=bottom adjusted, 2=centered, 3=top adjusted
        // cout<<"x_left "<<x_left<<endl;
        // cout<<"ceil(nSampleGroups/2.) "<<ceil(nSampleGroups/2.)<<endl;
        /*
        if(superimpose_EFThist)
        {
            float ylegend = 0.72; //Default
            if(use_paperStyle) {ylegend+= 0.04;}
            qw = new TLegend(x_left-0.06,ylegend,0.94,ylegend+0.15); //Default
            qw->SetTextSize(0.03);
        }
        else
        {
            float ylegend = 0.78; //Default
            if(use_paperStyle) {ylegend+= 0.04;}
            qw = new TLegend(x_left,ylegend,0.94,ylegend+0.09); //y-space for 2 rows
            qw->SetTextSize(0.04);
        }*/

        //-- Dummy object, only used to display uncertainty band also in legend
        TH1F* h_uncert = new TH1F("h_uncert", "h_uncert", 1, 0, 1);
        h_uncert->SetFillStyle(3254); //3002 //3004
        h_uncert->SetFillColor(kBlack);
        h_uncert->SetLineWidth(0.);
        qw->AddEntry(h_uncert, "Uncert.", "F");
        // qw->AddEntry(h_uncert, "Uncertainty", "F");

		//--Data on top of legend
        if(!this->is_blind)
        {
            if(use_combine_file && !combineFile_fromHarvester && g_data != 0) {qw->AddEntry(g_data, "Data" , "ep");}
            else if((!use_combine_file || combineFile_fromHarvester) && h_sum_data != 0) {qw->AddEntry(h_sum_data, "Data" , "ep");}
            else {cout<<__LINE__<<BOLD(FRED(" : null data !"))<<endl;}
        }

		for(int i=0; i<v_MC_histo.size(); i++)
		{
            // cout<<"MC_samples_legend["<<i<<"] "<<MC_samples_legend[i]<<endl;
			if(!v_MC_histo[i]) {continue;} //Fakes templates can be null

            if(MC_samples_legend[i].Contains("tZq")) {qw->AddEntry(v_MC_histo[i], "tZq", "f");}
            else if(MC_samples_legend[i].EndsWith("ttZ") ) {qw->AddEntry(v_MC_histo[i], "t#bar{t}Z", "f");}
            else if(MC_samples_legend[i].EndsWith("tWZ") ) {qw->AddEntry(v_MC_histo[i], "tWZ", "f");}
            else if(MC_samples_legend[i] == "ttW" || MC_samples_legend[i] == "tX") {qw->AddEntry(v_MC_histo[i], "t(#bar{t})X", "f");}
            else if(MC_samples_legend[i] == "WZ") {qw->AddEntry(v_MC_histo[i], "WZ", "f");}
            else if(MC_samples_legend[i] == "WWZ" || MC_samples_legend[i] == "VVV") {qw->AddEntry(v_MC_histo[i], "VV(V)", "f");}
            else if(MC_samples_legend[i] == "TTGamma_Dilep" || MC_samples_legend[i] == "XG") {qw->AddEntry(v_MC_histo[i], "X#gamma", "f");}
            else if(MC_samples_legend[i] == "TTbar_DiLep" || MC_samples_legend[i] == "NPL" || MC_samples_legend[i] == "NPL_DATA") {qw->AddEntry(v_MC_histo[i], "NPL", "f");}
		}

        // qw->SetTextFont(43);
        // gStyle->SetTextFont(72);
        // qw->AddEntry(v_MC_histo[0], ((TString) "\\text{a}"), "L"); //test, to remove

        //-- Separate legend for EFT entries
        TLegend* leg_eft = NULL;
        ylegend = 0.78-0.10; //Default
        if(use_paperStyle) {ylegend+= 0.01;}
        leg_eft = new TLegend(x_left-0.04,ylegend,0.94,ylegend+0.15); //Default
        leg_eft->SetTextSize(0.03);
        leg_eft->SetBorderSize(0);
        leg_eft->SetFillStyle(0); //transparent
        leg_eft->SetTextAlign(12); //align = 10*HorizontalAlign + VerticalAlign //Horiz: 1=left adjusted, 2=centered, 3=right adjusted //Vert: 1=bottom adjusted, 2=centered, 3=top adjusted

        if(superimpose_EFThist)
        {
            for(int isample=0; isample<v_EFT_samples.size(); isample++)
            {
                for(int ipoint=0; ipoint<v_EFT_points.size(); ipoint++)
                {
                    if(superimpose_EFT_auto)
                    {
                        if(!total_var_list[ivar].Contains("NN_all") && !total_var_list[ivar].Contains("NN_5D"))
                        {
                            if(!total_var_list[ivar].Contains("ctz") && !total_var_list[ivar].Contains("ctw") && (!total_var_list[ivar].Contains("cpq3") || (total_var_list[ivar].Contains("SRttZ") && !this->use_NN_cpq3_SRttZ))) {continue;}
                            if(total_var_list[ivar].Contains("cpq3") && !v_EFT_points[ipoint].Contains("cpq3")) {continue;}
                            if(total_var_list[ivar].Contains("ctz") && !v_EFT_points[ipoint].Contains("ctz")) {continue;}
                            if(total_var_list[ivar].Contains("ctw") && !v_EFT_points[ipoint].Contains("ctw")) {continue;}
                        }
                        // if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3")) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                        if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3") && !this->use_NN_cpq3_SRttZ) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                        if(total_var_list[ivar].Contains("SRtZq") && !v_EFT_samples[isample].Contains("tZq")) {continue;}
                        if(total_var_list[ivar].Contains("SRttZ") && !v_EFT_samples[isample].Contains("ttZ")) {continue;}
                    }

                    leg_eft->AddEntry(v2_th1eft[ipoint][isample], v2_th1eft_labels[ipoint][isample], "L");
                    // qw->AddEntry(v2_th1eft[ipoint][isample], v2_th1eft_labels[ipoint][isample], "L");
                }
            }
        }

        if(leg_eft->GetNRows()>0) {leg_eft->SetNColumns(leg_eft->GetNRows());} //Display a row, not a column


// #####  #####    ##   #    #
// #    # #    #  #  #  #    #
// #    # #    # #    # #    #
// #    # #####  ###### # ## #
// #    # #   #  #    # ##  ##
// #####  #    # #    # #    #

		//-- Canvas definition
		Load_Canvas_Style(); //Default top/bottom/left/right margins: 0.07/0.13/0.16/0.03
		TCanvas* c1 = new TCanvas("c1", "c1", 1000, 800);

        //-- Override Load_Canvas_Style margins
        c1->SetTopMargin(0.06);
        c1->SetBottomMargin(0.30);
        c1->SetRightMargin(0.04);
        c1->SetLeftMargin(0.14);

		if(draw_logarithm) {c1->SetLogy();}

		//Draw stack
		stack_MC->Draw("hist");

		//Draw data
		if(data_notEmpty && !this->is_blind)
		{
			if(use_combine_file && !combineFile_fromHarvester)
			{
				g_data->SetMarkerStyle(20);
				g_data->Draw("e0psame");
			}
			else
			{
				h_sum_data->SetMarkerStyle(20);
				h_sum_data->SetMinimum(0.001); //0 fails in logscale
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
        leg_eft->Draw("same"); //Draw EFT legend

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
                    if(superimpose_EFT_auto)
                    {
                        if(!total_var_list[ivar].Contains("NN_all") && !total_var_list[ivar].Contains("NN_5D"))
                        {
                            if(!total_var_list[ivar].Contains("ctz") && !total_var_list[ivar].Contains("ctw") && (!total_var_list[ivar].Contains("cpq3") || (total_var_list[ivar].Contains("SRttZ") && !this->use_NN_cpq3_SRttZ))) {continue;}
                            if(total_var_list[ivar].Contains("cpq3") && !v_EFT_points[ipoint].Contains("cpq3")) {continue;}
                            if(total_var_list[ivar].Contains("ctz") && !v_EFT_points[ipoint].Contains("ctz")) {continue;}
                            if(total_var_list[ivar].Contains("ctw") && !v_EFT_points[ipoint].Contains("ctw")) {continue;}
                        }
                        // if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3")) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                        if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3") && !this->use_NN_cpq3_SRttZ) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                        if(total_var_list[ivar].Contains("SRtZq") && !v_EFT_samples[isample].Contains("tZq")) {continue;}
                        if(total_var_list[ivar].Contains("SRttZ") && !v_EFT_samples[isample].Contains("ttZ")) {continue;}
                    }

                    if(v2_th1eft[ipoint][isample]->GetMaximum() > ymax_EFT) {ymax_EFT = v2_th1eft[ipoint][isample]->GetMaximum();}
                }
            }
        }

        double ymax = 0;
        if(use_combine_file && !combineFile_fromHarvester)
		{
			Long64_t locmax = TMath::LocMax(g_data->GetN(), g_data->GetY()); //the corresponding x value can be obtained with double xmax = gr->GetX()[locmax];
            if(data_notEmpty) {ymax = g_data->GetY()[locmax];} //Data ymax
            if(ymax < stack_MC->GetMaximum()) {ymax = stack_MC->GetMaximum();} //MC ymax
		}
		else
		{
            if(data_notEmpty) {ymax = h_sum_data->GetMaximum();} //Data ymax
            if(ymax < stack_MC->GetMaximum()) {ymax = stack_MC->GetMaximum();} //MC ymax
            if(ymax < ymax_EFT) {ymax = ymax_EFT;} //EFT ymax
		}
        // ymax*= 1.4; //Previous default, for ratio pad 0.25
        ymax*= 1.6; //CHANGED -- larger ratio pad (0.3) -> must make more place for legend
        stack_MC->SetMaximum(ymax);
        // if(ymax > qw->GetY1()) {ymax = qw->GetY1();} //Avoid overlap with TLegend
        // cout<<"qw->GetY1() "<<qw->GetY1()<<endl;
        // cout<<"c1->GetUymax() "<<c1->GetUymax()<<endl;

		stack_MC->SetMinimum(0.0001); //Remove '0' label

		if(draw_logarithm) //Can't use 0
		{
            stack_MC->SetMinimum(1.2); //Default
            if(total_var_list[ivar].Contains("mTW")) {stack_MC->SetMinimum(2.5);}

			stack_MC->SetMaximum(stack_MC->GetMaximum()*20); //Must use higher threshold in log
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

        //-- CHANGED -- Read full error from combine here
        //NB: the only way to get the correct, total uncertainty in a given bin, is to make combine sum all errors by itself  (--total-shapes options) ==> Read corresponding bin error here, and assign it to corresponding histo bin
        if(use_combine_file)
        {
            for(int ibin=1; ibin<nIndivBins+1; ibin++)
            {

                inputFile_path = Get_HistoFile_InputPath(!drawInputVars, "bin"+Convert_Number_To_TString(ibin)+"_" + var_tmp, cat_tmp, lumiName, use_combine_file, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, (EFTpoint!=""), combineFile_fromHarvester, prefit, false);
                if(inputFile_path == "") {cout<<"Get_HistoFile_InputPath --> file not found ! "<<endl; continue;}
                file_input = TFile::Open(inputFile_path, "READ");
                TString dir_hist_tmp = prefit? "prefit/":"postfit/"; //Total dir
                // cout<<"dir_hist_tmp/TotalProcs "<<dir_hist_tmp+"TotalProcs"<<" / Total error = "<<((TH1F*) file_input->Get(dir_hist_tmp+"TotalProcs"))->GetBinError(1)<<endl;
                if(!file_input->GetDirectory(dir_hist_tmp) || !file_input->GetDirectory(dir_hist_tmp)->GetListOfKeys()->Contains("TotalProcs") ) {cout<<FRED("Directory '"<<dir_hist_tmp<<"' or histogram '"<<dir_hist_tmp<<"TotalProcs' not found ! Skip...")<<endl; continue;}
                v_eyl[ibin-1] = ((TH1F*) file_input->Get(dir_hist_tmp+"TotalProcs"))->GetBinError(1);
                v_eyh[ibin-1] = ((TH1F*) file_input->Get(dir_hist_tmp+"TotalProcs"))->GetBinError(1);
                file_input->Close();
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
		TGraphAsymmErrors* gr_error = NULL;

		gr_error = new TGraphAsymmErrors(nofbins,xx,yy,exl,exh,eyl,eyh);
        gr_error->SetFillStyle(3254); //3002 //3004
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
		TPad* pad_ratio = new TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 1.0);
		pad_ratio->SetFillColor(0);
		pad_ratio->SetFillStyle(0);
		pad_ratio->SetGridy(1);
		pad_ratio->Draw();
		pad_ratio->cd(0);

        //-- Override Load_Canvas_Style setting
        pad_ratio->SetTopMargin(0.70);
        pad_ratio->SetRightMargin(0.04);
        pad_ratio->SetLeftMargin(0.14);

		if(use_combine_file && !combineFile_fromHarvester && data_notEmpty) //Copy the content of the data graph into a TH1F (NB : symmetric errors... but anyway symmetric for data)
		{
			if(!v_MC_histo[0]) {cout<<__LINE__<<FRED("Error : v_MC_histo[0] is null ! Abort")<<endl; return;}

	        h_sum_data = (TH1F*) v_MC_histo[0]->Clone(); //Clone binning of the MC histos
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

		TH1F* histo_ratio_data = NULL;
		if(data_notEmpty)
		{
			histo_ratio_data = (TH1F*) h_sum_data->Clone();
            if(use_poisson_dataErrors) {histo_ratio_data->SetBinErrorOption(TH1::kPoisson);}

			if(!show_pulls_ratio) //Compute ratios (with error bars)
			{
				//To get correct error bars in ratio plot, must only account for errors from data, not MC ! (MC error shown as separate band)
				for(int ibin=1; ibin<histo_total_MC->GetNbinsX()+1; ibin++)
				{
					histo_total_MC->SetBinError(ibin, 0.);
				}

				histo_ratio_data->Divide(histo_total_MC);

                //Debug printouts
                // int ibin = 7;
    			// cout<<"h_sum_data->GetBinContent("<<ibin<<") "<<h_sum_data->GetBinContent(ibin)<<endl;
    			// cout<<"h_sum_data->GetBinError("<<ibin<<") "<<h_sum_data->GetBinError(ibin)<<endl;
    			// cout<<"histo_total_MC->GetBinContent("<<ibin<<") "<<histo_total_MC->GetBinContent(ibin)<<endl;
    			// cout<<"histo_total_MC->GetBinError("<<ibin<<") "<<histo_total_MC->GetBinError(ibin)<<endl;
                // cout<<"histo_ratio_data->GetBinError("<<ibin<<") "<<histo_ratio_data->GetBinError(ibin)<<endl;
                // cout<<"histo_ratio_data->GetBinError("<<ibin<<") "<<histo_ratio_data->GetBinError(ibin)<<endl;
			}
		 	else //-- Compute pulls (no error bars)
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

			//-- Debug printouts
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
		else
        {
            histo_ratio_data->GetYaxis()->SetTitle("Data/MC");
            // histo_ratio_data->GetYaxis()->SetTitle("Data/Prediction");
        }
		histo_ratio_data->GetYaxis()->SetTickLength(0.);
        histo_ratio_data->GetYaxis()->SetTitleOffset(1.15);
        // histo_ratio_data->GetYaxis()->SetTitleOffset(1.15);
        histo_ratio_data->GetXaxis()->SetTitleOffset(1.05);
        histo_ratio_data->GetYaxis()->SetLabelSize(0.04);
        histo_ratio_data->GetXaxis()->SetLabelSize(0.045);
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
        histo_ratio_data->GetXaxis()->SetTitleSize(0.05);

        //-- If a point is outside the y-range of the ratio pad defined by SetMaximum/SetMinimum(), it disappears with its error
        //-- Trick: fill 2 histos with points either above/below y-range, to plot some markers indicating missing points (cleaner)
        //NB: only for ratio plot, not pulls
        float ratiopadmin = 0.4, ratiopadmax = 1.6; //Define ymin/ymax for ratio plot

        //-- Zoom in ratio plot for CR plots (better agreement)
        if(total_var_list[ivar].Contains("CR") || total_var_list[ivar].Contains("SRttZ4l")) {ratiopadmin = 0.8; ratiopadmax = 1.2;}

        TH1F* h_pointsAboveY = (TH1F*) histo_ratio_data->Clone();
        h_pointsAboveY->SetMarkerStyle(26); //Open triangle pointing up
        h_pointsAboveY->SetMarkerSize(1.5);
        TH1F* h_pointsBelowY = (TH1F*) histo_ratio_data->Clone();
        h_pointsBelowY->SetMarkerStyle(32); //Open triangle pointing down
        h_pointsBelowY->SetMarkerSize(1.5);
        if(show_pulls_ratio)
		{
			histo_ratio_data->SetMinimum(-2.99);
			histo_ratio_data->SetMaximum(2.99);
		}
		else
		{
            //-- Default
            histo_ratio_data->SetMinimum(ratiopadmin); //NB: removes error bars if data point is below ymin...?
            histo_ratio_data->SetMaximum(ratiopadmax);

            //-- Fill histos with points outside yrange
            for(int ibin=1; ibin<histo_ratio_data->GetNbinsX()+1; ibin++)
            {
                //-- Debug printouts
                // cout<<"histo_ratio_data->GetBinContent("<<ibin<<") "<<histo_ratio_data->GetBinContent(ibin)<<endl;
                // cout<<"histo_ratio_data->GetBinError("<<ibin<<") "<<histo_ratio_data->GetBinError(ibin)<<endl;

                //-- Default: make point invisible
                h_pointsAboveY->SetBinContent(ibin, -999);
                h_pointsBelowY->SetBinContent(ibin, -999);

                if(histo_ratio_data->GetBinContent(ibin) > ratiopadmax && h_sum_data->GetBinContent(ibin) >= 1)
                {
                    //Adjust error
                    float initial_y = histo_ratio_data->GetBinContent(ibin);
                    float initial_err = histo_ratio_data->GetBinError(ibin);
                    float new_err = initial_err - (initial_y-ratiopadmax);
                    if(new_err<0) {new_err=0.;}

                    //-- Debug printouts
                    // cout<<"ABOVE, bin "<<ibin<<endl;
                    // cout<<"initial_y "<<initial_y<<endl;
                    // cout<<"initial_err "<<initial_err<<endl;
                    // cout<<"new_err "<<new_err<<endl;

                    h_pointsAboveY->SetBinContent(ibin, ratiopadmax-0.05); //Add some padding
                    h_pointsAboveY->SetBinError(ibin, new_err);
                }
                else if(histo_ratio_data->GetBinContent(ibin) < ratiopadmin && h_sum_data->GetBinContent(ibin) >= 1)
                {
                    //Adjust error
                    float initial_y = histo_ratio_data->GetBinContent(ibin);
                    float initial_err = histo_ratio_data->GetBinError(ibin);
                    float new_err = initial_err - (ratiopadmin-initial_y);
                    if(new_err<0) {new_err=0.;}

                    //-- Debug printouts
                    // cout<<"BELOW, bin "<<ibin<<endl;
                    // cout<<"histo_total_MC->GetBinContent(ibin) "<<histo_total_MC->GetBinContent(ibin)<<endl;
                    // cout<<"initial_y "<<initial_y<<endl;
                    // cout<<"initial_err "<<initial_err<<endl;
                    // cout<<"new_err "<<new_err<<endl;

                    h_pointsBelowY->SetBinContent(ibin, ratiopadmin+(ratiopadmin/10.)); //Add some padding
                    h_pointsBelowY->SetBinError(ibin, new_err);
                }
            }
		}

        //-- SET X_AXIS TITLES
		if(drawInputVars) {histo_ratio_data->GetXaxis()->SetTitle(Get_Variable_Name(total_var_list[ivar]));}
		else
		{
            //-- Hardcoded x-titles
            TString varname = total_var_list[ivar];
            if(use_NN_SRother && varname.Contains("SRother")) {varname = "NN_SM_SRother";} //Special case
            TString xtitle = Get_Template_XaxisTitle(varname, this->use_paperStyle, this->use_NN_cpq3_SRttZ);
            histo_ratio_data->GetXaxis()->SetTitle(xtitle);

			if(total_var_list[ivar].Contains("categ")) //Vertical text X labels (categories names)
			{
                histo_ratio_data->GetXaxis()->SetTitle("");
                histo_ratio_data->GetXaxis()->SetLabelSize(0.06); //Increase x-label size
                histo_ratio_data->GetXaxis()->SetLabelOffset(0.02); //Add some x-label offset
                histo_ratio_data->LabelsOption("v", "X"); //Make X labels vertical

                //Hard-coded
                {
                    const char *labels[10]  = {"1bj,2j","1bj,3j","1bj,4j","1bj,5j","1bj,6j","2bj,2j","2bj,3j","2bj,4j","2bj,5j","2bj,6j"};
                    for(int i=1;i<=10;i++) {histo_ratio_data->GetXaxis()->SetBinLabel(i,labels[i-1]);}
                }
			}
            else if(total_var_list[ivar].Contains("countExp")) {histo_ratio_data->GetXaxis()->SetLabelSize(0.);} //Counting exp. <-> single bin <-> x labels useless
            else if(use_combine_file) //Need to hard-code x-axis labels
            {
                int nbins=1; float xmin=0, xmax=0;
                int nbjets_min = 1, nbjets_max=2, njets_min=2, njets_max=6;
                vector<float> dummy;
                Get_Template_Range(nbins, xmin, xmax, total_var_list[ivar], this->make_SMvsEFT_templates_plots, this->categorization_strategy, plot_onlyMaxNodeEvents, nbjets_min, nbjets_max, njets_min, njets_max, dummy, use_NN_SRother, use_NN_cpq3_SRttZ);

                //Arbitrary bins
                if(total_var_list[ivar].Contains("NN") || total_var_list[ivar].Contains("BDT"))
                {
                    // histo_ratio_data->GetYaxis()->SetMoreLogLabels();
                    histo_ratio_data->GetXaxis()->SetLabelSize(0.05); //Increase x-label size
                    histo_ratio_data->GetXaxis()->SetLabelOffset(0.02); //Add some x-label offset

                    for(int i=1;i<=nbins;i++)
                    {
                        // TString label = Convert_Number_To_TString(xmin + (i-1)*((xmax-xmin)/nbins));
                        TString label = "Bin " + Convert_Number_To_TString(i);
                        histo_ratio_data->GetXaxis()->SetBinLabel(i, label);
                        // histo_ratio_data->LabelsOption("v", "X"); //X labels vertical

                        // custom_axis = new TGaxis(0., 0.,  1., histo_ratio_data->GetMaximum(), xmin, xmax, 503, "-S");
                        // custom_axis->Draw("same");
                    }
                }
                else if(total_var_list[ivar].Contains("mTW"))
                {
                    histo_ratio_data->GetXaxis()->SetNdivisions(15); //Must set Ndivisions=nbins

                    vector<TString> v_labels { "0"," "," "," "," ","50"," "," "," "," ","100"," "," "," "," ","150"}; //NB: empty entry <-> " " (space)
                    // vector<TString> v_labels { "0","","20","","40","","60","","80","","100","","120","0","140", "" };

                    for(int ibin=0;ibin<nbins+1;ibin++)
                    {
                        // cout<<"ibin "<<ibin<<" / v_labels[i] "<<v_labels[ibin]<<endl;
                        // TString label = Convert_Number_To_TString(xmin + (i-1)*((xmax-xmin)/nbins));
                        // ChangeLabel (Int_t labNum=0, Double_t labAngle=-1., Double_t labSize=-1., Int_t labAlign=-1, Int_t labColor=-1, Int_t labFont=-1, TString labText="")
                        histo_ratio_data->GetXaxis()->ChangeLabel(ibin+1, -1, -1, -1, -1, -1, v_labels[ibin]);
                    }
                }
            }
		}

		pad_ratio->cd(0);
		if(show_pulls_ratio) {histo_ratio_data->Draw("HIST P");} //Draw ratio points
        else
        {
            histo_ratio_data->Draw("E1 X0 P"); //Draw ratio points ; E1 : perpendicular lines at end ; X0 : suppress x errors

            h_pointsAboveY->Draw("E1 X0 P same");
            h_pointsBelowY->Draw("E1 X0 P same");
        }

        //-- Compute a chi2-style quantity
        // float total_chi2 = 0;
        // for(int ibin=0; ibin<h_sum_data->GetNbinsX(); ibin++)
        // {
        //     float chi2 = abs(h_sum_data->GetBinContent(ibin+1) - histo_total_MC->GetBinContent(ibin+1)) / (h_sum_data->GetBinError(ibin+1)+histo_total_MC->GetBinError(ibin+1));
        //     total_chi2+= chi2;
        //     cout<<"-- ibin "<<ibin<<"--> chi2 "<<chi2<<endl;
        // }
        // cout<<"--> TOTAL chi2 "<<total_chi2/h_sum_data->GetNbinsX()<<endl;



// ###### #####  #####   ####  #####   ####     #####    ##   ##### #  ####
// #      #    # #    # #    # #    # #         #    #  #  #    #   # #    #
// #####  #    # #    # #    # #    #  ####     #    # #    #   #   # #    #
// #      #####  #####  #    # #####       #    #####  ######   #   # #    #
// #      #   #  #   #  #    # #   #  #    #    #   #  #    #   #   # #    #
// ###### #    # #    #  ####  #    #  ####     #    # #    #   #   #  ####

		TGraphAsymmErrors* gr_ratio_error = NULL;
		if(draw_errors)
		{
			//Copy previous TGraphAsymmErrors, then modify it -> error TGraph for ratio plot
			TGraphAsymmErrors *thegraph_tmp = NULL;
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
            gr_ratio_error->SetFillStyle(3254); //3002 //3004
            gr_ratio_error->SetFillColor(kBlack); //kBlue+2 //kCyan

			pad_ratio->cd(0);
			if(!show_pulls_ratio) {gr_ratio_error->Draw("e2 same");} //Draw error bands in ratio plot

            //-- Add sub-legend here ? (stat or stat+syst)
		} //draw errors

        //-- Testing EFT scenario in ratio plot
        vector<TH1F*> v_hEFT_ratioPlot;
        if(false)
        {
            int icolor = 0; //Index
            vector<int> v_col_eft;
            v_col_eft.push_back(kBlack); //SM
            v_col_eft.push_back(kAzure+6);
            v_col_eft.push_back(kMagenta+1);
            v_col_eft.push_back(kOrange+6);
            v_col_eft.push_back(kGreen+2);

            for(int ipoint=0; ipoint<v_EFT_points.size(); ipoint++)
            {
                for(int isample=0; isample<v_EFT_samples.size(); isample++)
                {
                    if(superimpose_EFT_auto)
                    {
                        if(!total_var_list[ivar].Contains("NN_all") && !total_var_list[ivar].Contains("NN_5D"))
                        {
                            if(!total_var_list[ivar].Contains("ctz") && !total_var_list[ivar].Contains("ctw") && (!total_var_list[ivar].Contains("cpq3") || (total_var_list[ivar].Contains("SRttZ") && !this->use_NN_cpq3_SRttZ))) {continue;}
                            if(total_var_list[ivar].Contains("cpq3") && !v_EFT_points[ipoint].Contains("cpq3")) {continue;}
                            if(total_var_list[ivar].Contains("ctz") && !v_EFT_points[ipoint].Contains("ctz")) {continue;}
                            if(total_var_list[ivar].Contains("ctw") && !v_EFT_points[ipoint].Contains("ctw")) {continue;}
                        }
                        // if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3")) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                        if(total_var_list[ivar].Contains("NN_5D") && total_var_list[ivar].Contains("SRttZ") && v_EFT_points[ipoint].Contains("rwgt_cpq3") && !this->use_NN_cpq3_SRttZ) {continue;} //NB: specify full substring 'rwgt_cpq3...' because only want to ignore cases where cpq3 is the only operator in the string (not e.g. strings containing the 5 WCs)
                        if(total_var_list[ivar].Contains("SRtZq") && !v_EFT_samples[isample].Contains("tZq")) {continue;}
                        if(total_var_list[ivar].Contains("SRttZ") && !v_EFT_samples[isample].Contains("ttZ")) {continue;}
                    }

                    v_hEFT_ratioPlot.push_back((TH1F*) histo_total_MC->Clone());

                    v_hEFT_ratioPlot[v_hEFT_ratioPlot.size()-1]->SetLineColor(v_col_eft[icolor]);
                    // v_hEFT_ratioPlot[v_hEFT_ratioPlot.size()-1]->SetLineStyle(9); //Dashed lines
                    v_hEFT_ratioPlot[v_hEFT_ratioPlot.size()-1]->SetLineWidth(3);
                    v_hEFT_ratioPlot[v_hEFT_ratioPlot.size()-1]->SetFillStyle(0);
                    icolor++;

                    for(int ibin=1; ibin<v_hEFT_ratioPlot[v_hEFT_ratioPlot.size()-1]->GetNbinsX()+1; ibin++)
                    {
                        float content = v_hEFT_ratioPlot[v_hEFT_ratioPlot.size()-1]->GetBinContent(ibin);
                        if(normalize_EFThist) {v2_th1eft[ipoint][isample]->Scale(((TH1*) stack_MC->GetStack()->Last())->Integral()/(2.*v2_th1eft[ipoint][isample]->Integral()));} //Normalize to half-integral of stack (arbitrary)  //NB: access stack integral via summed object '((TH1*) stack_MC->GetStack()->Last())'
                        v2_th1eft[ipoint][isample]->Scale("SM");
                        content-= v2_th1eft[ipoint][isample]->GetBinContent(ibin);
                        WCPoint wcp = WCPoint((string) v_EFT_points[ipoint], 1.);
                        v2_th1eft[ipoint][isample]->Scale(wcp);
                        content+= v2_th1eft[ipoint][isample]->GetBinContent(ibin);
                        v_hEFT_ratioPlot[v_hEFT_ratioPlot.size()-1]->SetBinContent(ibin, content);
                        // cout<<"bin "<<ibin<<" / content "<<content<<endl;
                    }
                    v_hEFT_ratioPlot[v_hEFT_ratioPlot.size()-1]->Divide(histo_total_MC);
                    v_hEFT_ratioPlot[v_hEFT_ratioPlot.size()-1]->Draw("hist same");
                } //EFT samples
            } //EFT points
        } //Draw EFT in ratio plot


//  ####   ####   ####  #    # ###### ##### #  ####   ####
// #    # #    # #      ##  ## #        #   # #    # #
// #      #    #  ####  # ## # #####    #   # #       ####
// #      #    #      # #    # #        #   # #           #
// #    # #    # #    # #    # #        #   # #    # #    #
//  ####   ####   ####  #    # ######   #   #  ####   ####

		//-- Draw ratio y-lines manually
		TH1F *h_line1 = NULL;
		TH1F *h_line2 = NULL;
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

        // TString Y_label = ""; //cf. CMS guidelines
		double xmax_stack = stack_MC->GetXaxis()->GetXmax();
		double xmin_stack = stack_MC->GetXaxis()->GetXmin();
        TString Y_label = "Events / bin"; //Default
		if(!use_combine_file) //Compute bin width
		{
			double xmax = histo_total_MC->GetXaxis()->GetXmax();
			double xmin = histo_total_MC->GetXaxis()->GetXmin();
			Y_label = "Events / " + Convert_Number_To_TString( (xmax - xmin) / histo_total_MC->GetNbinsX(), 3); //Automatically get the Y label depending on binning
            Y_label+= Get_Unit_Variable(total_var_list[ivar]);
        }
        else //combine ruins x-axis info, so can't rely on binning of the TH1F; call same functions as used to produce the histos to get the same binning
        {
            int nbins = 0; float xmin = 0, xmax = 0;
            if(drawInputVars) {Get_Variable_Range(total_var_list[ivar], nbins, xmin, xmax);}
            else {int nbjets_min=0,nbjets_max=0,njets_min=0,njets_max=0; vector<float> minmax_bounds; Get_Template_Range(nbins, xmin, xmax, total_var_list[ivar], this->make_SMvsEFT_templates_plots, this->categorization_strategy, plot_onlyMaxNodeEvents, nbjets_min, nbjets_max, njets_min, njets_max, minmax_bounds, this->use_NN_SRother, this->use_NN_cpq3_SRttZ);} //NB: for NNs the binning is anyway arbitrary (would need to read NN_settings file, etc.)

            if(xmax != 0)
            {
                Y_label = "Events / " + Convert_Number_To_TString( (xmax - xmin) / nbins, 3); //Automatically get the Y label depending on binning
                Y_label+= Get_Unit_Variable(total_var_list[ivar]);
            }
        }
        if(total_var_list[ivar].Contains("countExp")) {Y_label = "Events";} //No unit for simple counting experiment

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
            stack_MC->GetYaxis()->SetTitleOffset(1.15);
            // stack_MC->GetYaxis()->SetTitleOffset(1.28);
			stack_MC->GetYaxis()->SetTitle(Y_label);

            stack_MC->GetXaxis()->SetTickLength(0.);
		}

    	//----------------
    	// CAPTIONS //
    	//----------------
    	// -- using https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines
        // -- About fonts: https://root.cern.ch/doc/master/classTAttText.html#T5

		float l = c1->GetLeftMargin();
		float t = c1->GetTopMargin();

		TString cmsText = "CMS";
		TLatex latex;
		latex.SetNDC();
		latex.SetTextColor(kBlack);
        // latex.SetTextFont(61);
        latex.SetTextFont(62); //Changed
		latex.SetTextAlign(11);
		latex.SetTextSize(0.06);
        if(use_paperStyle) {latex.DrawLatex(l + 0.04, 0.87, cmsText);} //CMS guideline: within frame
		else {latex.DrawLatex(l + 0.01, 0.95, cmsText);} //Default: outside frame

		TString extraText = "Preliminary";
		latex.SetTextFont(52);
		latex.SetTextSize(0.05);
		if(!use_paperStyle) {latex.DrawLatex(l + 0.12, 0.91+0.04, extraText);}

		float lumi = lumiValue;
		TString lumi_ts = Convert_Number_To_TString(lumi);
		lumi_ts += " fb^{-1} (13 TeV)";
		latex.SetTextFont(42);
		latex.SetTextAlign(31);
		latex.SetTextSize(0.04);
        if(use_paperStyle) {latex.DrawLatex(0.96, 0.95,lumi_ts);}
        else {latex.DrawLatex(0.96, 0.95,lumi_ts);}

		//------------------
		//-- channel info
		TLatex text2 ;
		text2.SetNDC();
		text2.SetTextAlign(13);
		text2.SetTextSize(0.045);
		text2.SetTextFont(42);

        TString info_data = Get_Region_Label(region, total_var_list[ivar]);
        if(info_data != "")
        {
            if(use_paperStyle) {text2.DrawLatex(l + 0.04,0.83,info_data);}
            else {text2.DrawLatex(0.20,0.90,info_data);} //Default: outside frame
        }


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
            if(make_histos_forControlPlotsPaper) {outdir = "plots/tmp/";}
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
            if(cat_tmp != "")
            {
                outdir+= cat_tmp + "/";
                mkdir(outdir.Data(), 0777);
            }
            outdir+= lumiName;
            mkdir(outdir.Data(), 0777);
            if(prefit) {outdir+= "/prefit";}
            else {outdir+= "/postfit";}
            mkdir(outdir.Data(), 0777);
            if(this->categorization_strategy>0 && make_SMvsEFT_templates_plots) {outdir+= "/strategy" + Convert_Number_To_TString(this->categorization_strategy); if(EFTpoint!="") {outdir+= "param";}}
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
            output_plot_name = outdir + total_var_list[ivar] + (EFTpoint==""? "":"_"+EFTpoint) + "_template";
            // output_plot_name = outdir + total_var_list[ivar] + (EFTpoint==""? "":"_"+EFTpoint) + "_template_" + signal_process;
		}
		if(channel != "") {output_plot_name+= "_" + channel;}
        if(!drawInputVars && categorization_strategy == 2 && (make_SMvsEFT_templates_plots || (!make_SMvsEFT_templates_plots && plot_onlyMaxNodeEvents)) ) {output_plot_name+= "_maxNode";} //Cases for which we want to cut on a multiclass MVA-SM
        else if(!drawInputVars && categorization_strategy == 1 && (make_SMvsEFT_templates_plots || (!make_SMvsEFT_templates_plots && plot_onlyMVACutEvents)) ) {output_plot_name+= "_MVAcut";}
		output_plot_name+= this->filename_suffix;
		if(draw_logarithm) {output_plot_name+= "_log";}
        if(this->is_blind) {output_plot_name+= "_blind";}
        if(use_paperStyle) {output_plot_name+= ".eps";} //Use .eps to have correct \ell labels
		else {output_plot_name+= this->plot_extension;}

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
        for(int i=0; i<v_hEFT_ratioPlot.size(); i++) {delete v_hEFT_ratioPlot[i];}
		if(use_combine_file && !combineFile_fromHarvester) {delete g_data; g_data = NULL;}
        if(h_uncert) {delete h_uncert; h_uncert = NULL;}
	} //Var loop

	file_input->Close();
    if(f_EFT) {f_EFT->Close();}

    //-- Activate this function to printout yields in each subregion //NB: only works with my own template files (not those from combine)
    vector<TString> vregions;
    vregions.push_back("SRtZq");
    vregions.push_back("SRttZ");
    vregions.push_back("SRother");
    // vregions.push_back("WZCR");
    // vregions.push_back("ZZCR");
    // Print_Yields_fromHistograms(inputFile_path, template_name, v_lumiYears, vregions, sample_list);

    //-- Debug printouts
    // for(int isample=0; isample<sample_list.size(); isample++)
    // {
    //     cout<<"sample "<<sample_list[isample]<<" --> "<<v_yields_processes[isample]<<endl;
    // }

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
void TopEFT_analysis::Compare_TemplateShapes_Processes(TString template_name, TString channel)
{
//--- OPTIONS --------------------------------
//--------------------------------------------
	bool drawInputVars = true;

	bool normalize = false; //Ratio plot incorrect if false, why ?

    TString type = "ttZ"; //'' / 'tZq' / 'ttZ' / 'tWZ' --> Compare corresponding private/central samples

    vector<TString> v_years; float lumi = 0.; TString luminame_tmp = ""; //Select 1 or multiple years
    v_years.push_back("2016"); lumi+= 36.33; luminame_tmp = "2016";
    v_years.push_back("2017"); lumi+= 41.53; luminame_tmp = luminame_tmp=="2016"? "201617":"2017";
    v_years.push_back("2018"); lumi+= 59.74; luminame_tmp = (luminame_tmp=="2016"? "201618":(luminame_tmp=="201617"? "Run2":(luminame_tmp=="2017"? "201718":"2018")) );

    vector<TString> total_var_list;
    //* Default input vars
    for(int i=0; i<var_list.size(); i++) {total_var_list.push_back(var_list[i]);}
    for(int i=0; i<v_add_var_names.size(); i++) {total_var_list.push_back(v_add_var_names[i]);}

    //* Additional input vars
    /*
    total_var_list.push_back("mTW");
    total_var_list.push_back("mHT");
    total_var_list.push_back("Mass_3l");
    total_var_list.push_back("maxEtaJet");
    total_var_list.push_back("jPrimeAbsEta");
    total_var_list.push_back("channel");
    total_var_list.push_back("njets");
    total_var_list.push_back("nbjets");
    total_var_list.push_back("metEt");*/


    //-- Hardcode samples here... or could filter the main sample list
	vector<TString> v_samples; vector<TString> v_groups; vector<int> v_colors;
    // v_samples.push_back("tZq"); v_groups.push_back("tZq (Central)"); v_colors.push_back(kRed);
    // v_samples.push_back("PrivMC_tZq"); v_groups.push_back("tZq (Private)"); v_colors.push_back(kBlue);
    // v_samples.push_back("PrivMC_tZq_v3"); v_groups.push_back("tZq (Private v3)"); v_colors.push_back(kMagenta);
    // v_samples.push_back("ttZ"); v_groups.push_back("ttZ (Central)"); v_colors.push_back(kRed);
    // v_samples.push_back("PrivMC_ttZ"); v_groups.push_back("ttZ (Private)"); v_colors.push_back(kBlue);

    // v_samples.push_back("PrivMC_tZq"); v_groups.push_back("PrivMC_tZq"); v_colors.push_back(kRed);
    // v_samples.push_back("PrivMC_ttZ"); v_groups.push_back("PrivMC_ttZ"); v_colors.push_back(kBlue);
    // v_samples.push_back("PrivMC_tWZ"); v_groups.push_back("PrivMC_tWZ"); v_colors.push_back(kMagenta);
    // v_samples.push_back("tX"); v_groups.push_back("tX"); v_colors.push_back(kRed);
    // v_samples.push_back("VVV"); v_groups.push_back("VVV"); v_colors.push_back(kBlue);
    // v_samples.push_back("XG"); v_groups.push_back("XG"); v_colors.push_back(kMagenta);

    if(type == "tZq")
    {
        v_samples.resize(0); v_groups.resize(0); v_colors.resize(0);
        v_samples.push_back("tZq"); v_groups.push_back("tZq (Central)"); v_colors.push_back(kBlack);
        // v_samples.push_back("PrivMC_tZq_v3"); v_groups.push_back("tZq (Private v3)"); v_colors.push_back(kBlue);
        v_samples.push_back("PrivMC_tZq"); v_groups.push_back("tZq (Private)"); v_colors.push_back(kRed);
        // v_samples.puESUp");
    // v_syst.push_back("Jsh_back("PrivMC_tZq_TOP19001"); v_groups.push_back("tZq (TOP19001)"); v_colors.push_back(kOrange);
    }
    else if(type == "ttZ")
    {
        v_samples.resize(0); v_groups.resize(0); v_colors.resize(0);
        v_samples.push_back("ttZ"); v_groups.push_back("ttZ (Central)"); v_colors.push_back(kBlack);
        v_samples.push_back("PrivMC_ttZ"); v_groups.push_back("ttZ (Private)"); v_colors.push_back(kRed);
        // v_samples.push_back("PrivMC_ttZ_TOP19001"); v_groups.push_back("ttZ (TOP19001)"); v_colors.push_back(kOrange);
    }
    else if(type == "tWZ")
    {
        v_samples.resize(0); v_groups.resize(0); v_colors.resize(0);
        v_samples.push_back("tWZ"); v_groups.push_back("tWZ (Central)"); v_colors.push_back(kBlack);
        v_samples.push_back("PrivMC_tWZ"); v_groups.push_back("tWZ (Private)"); v_colors.push_back(kRed);
    }

    vector<TString> v_syst;
    v_syst.push_back("");
    // v_syst.push_back("JESUp");
    // v_syst.push_back("JESDown");
    // v_syst.push_back("JER2017Up");
    // v_syst.push_back("JER2017Down");
    // v_syst.push_back("MET2017Up");
    // v_syst.push_back("MET2017Down");

//--------------------------------------------
//--------------------------------------------

    bool use_predefined_EFT_strategy = false;
    if(!drawInputVars && categorization_strategy > 0 && !this->make_fixedRegions_templates) {use_predefined_EFT_strategy = true;}

    if(!drawInputVars)
    {
        total_var_list.clear();
        // if(template_name == "BDT") {total_var_list.push_back(classifier_name);}
        // else
        // {
        //     if(NN_nNodes==1) {total_var_list.push_back(template_name);}
        //     else
        //     {
        //         for(int inode=0; inode<NN_nNodes; inode++) {total_var_list.push_back(template_name + Convert_Number_To_TString(inode));}
        //     }
        // }

        Fill_Variables_List(total_var_list, use_predefined_EFT_strategy, template_name, this->region, this->scanOperators_paramNN, this->NN_nNodes, this->make_SMvsEFT_templates_plots, operator_scan1, operator_scan2, v_WCs_operator_scan1, v_WCs_operator_scan2, this->make_fixedRegions_templates);
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
    if(cat_tmp=="") {cat_tmp = "signal";} //Default -> SR3l

    //-- Read input file (may be year-dependent)
    TString input_name;
    TString template_type = this->make_fixedRegions_templates? "otherRegions":template_name;
    input_name = Get_HistoFile_InputPath(!drawInputVars, template_type, cat_tmp, luminame_tmp, false, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, false);
    if(input_name == "") {cat_tmp = "signal"; input_name = Get_HistoFile_InputPath(!drawInputVars, template_type, cat_tmp, "Run2", false, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, false);} //Retry with 'signal_process' as 'region' argument
    // if(input_name == "") {cat_tmp = signal_process; input_name = Get_HistoFile_InputPath(!drawInputVars, template_type, cat_tmp, "Run2", false, this->filename_suffix, make_SMvsEFT_templates_plots, categorization_strategy, this->make_fixedRegions_templates, false);} //Retry with 'signal_process' as 'region' argument
    if(input_name == "") {return;}
    TFile* file_input = TFile::Open(input_name, "READ");
    cout<<DIM("Opening file: "<<input_name<<"")<<endl;

	//To retrieve/store distributions
	TH1F* h_tmp = NULL; //Tmp storing histo

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

    			TString samplename = v_samples[isample];

                for(int isyst=0; isyst<v_syst.size(); isyst++)
                {
    				// cout<<"syst "<<v_syst[isyst]<<endl;

                    v3_histos_var_sample_syst[ivar][isample][isyst] = NULL; //Init

                    for(int iyear=0; iyear<v_years.size(); iyear++)
                    {
                        TString histo_name = total_var_list[ivar];
                        if(channel_list[ichan] != "") {histo_name+= "_" + channel_list[ichan];}
                        histo_name+= "_" + v_years[iyear];
                        histo_name+= "__" + samplename;
                        if(v_syst[isyst] != "") {histo_name+= "__" + v_syst[isyst];}

            			if(!file_input->GetListOfKeys()->Contains(histo_name)) {cout<<ITAL("Histogram '"<<histo_name<<"' : not found ! Skip...")<<endl; continue;}

            			h_tmp = (TH1F*) file_input->Get(histo_name);
        				h_tmp->SetDirectory(0); //Dis-associate from TFile
            			// cout<<"histo_name "<<histo_name<<" / h_tmp->Integral() = "<<h_tmp->Integral()<<endl;


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
                        // if(samplename.Contains("PrivMC")) {h_tmp->SetLineStyle(2);}

        				//HARDCODED
        				// if(v_syst[isyst] == "JESUp") {h_tmp->SetLineColor(kRed);}
        				// else if(v_syst[isyst] == "JESDown") {h_tmp->SetLineColor(kBlue);}
                        // h_tmp->SetLineColor(v_colors[isample]+isyst);
                        if(v_syst[isyst] != "") {h_tmp->SetLineColor(2+isyst);}
        				// cout<<"v_colors[isample] "<<v_colors[isample]<<endl;

            			h_tmp->SetLineWidth(3);

                        if(v_syst[isyst] != "") {h_tmp->SetLineStyle(2);}

            			if(normalize) {h_tmp->Scale(1./h_tmp->Integral() );}

                        if(v3_histos_var_sample_syst[ivar][isample][isyst]==NULL) {v3_histos_var_sample_syst[ivar][isample][isyst] = (TH1F*) h_tmp->Clone();}
                        else {v3_histos_var_sample_syst[ivar][isample][isyst]->TH1::Add((TH1F*) h_tmp);} //Sum years
        				// cout<<"v3_histos_var_sample_syst[ivar]["<<isample<<"]["<<isyst<<"]->Integral() "<<v3_histos_var_sample_syst[ivar][isample][isyst]->Integral()<<endl;

            			delete h_tmp; h_tmp = NULL;
                    } //end year loop
                } //end syst loop
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
        if(v3_histos_var_sample_syst[ivar][0][0] == NULL) {continue;} //In case this variable was not found

        //Canvas definition
        Load_Canvas_Style();
        gStyle->SetOptTitle(1);
        TCanvas* c = new TCanvas("", "", 1000, 800);
        c->SetTopMargin(0.1); //Overrides Load_Canvas_Style setting
        c->SetBottomMargin(0.25); //Overrides Load_Canvas_Style setting
        c->cd();

        // c->SetLogy();

        TLegend* qw = new TLegend(0.75,.70,1.,1.);

        //-- Find y-axis max threshold
        float ymax = 0.;
        for(int isample=0; isample<v3_histos_var_sample_syst[ivar].size(); isample++)
    	{
    		for(int isyst=0; isyst<v_syst.size(); isyst++)
    		{
                // cout<<"var "<<total_var_list[ivar]<<" / isample "<<isample<<" / isyst "<<isyst<<endl;
                if(v3_histos_var_sample_syst[ivar][isample][isyst]->GetMaximum() > ymax) {ymax = v3_histos_var_sample_syst[ivar][isample][isyst]->GetMaximum();}
            }
        }

        //-- Draw
    	for(int isample=0; isample<v3_histos_var_sample_syst[ivar].size(); isample++)
    	{
    		for(int isyst=0; isyst<v_syst.size(); isyst++)
    		{
    			if(v_samples[isample].Contains("Fake") && !v_syst[isyst].Contains("Clos") && !v_syst[isyst].Contains("FR") && v_syst[isyst] != "") {continue;}

    			if(v3_histos_var_sample_syst[ivar][isample][isyst] == 0) {cout<<"Null histo ! Skip"<<endl; continue;}

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
                else {v3_histos_var_sample_syst[ivar][isample][isyst]->GetYaxis()->SetTitle("Events / bin");}

                //-- CHANGED
                v3_histos_var_sample_syst[ivar][isample][isyst]->SetMaximum(ymax*1.4);
                // v3_histos_var_sample_syst[ivar][isample][isyst]->SetMaximum(v3_histos_var_sample_syst[ivar][isample][isyst]->GetMaximum()*1.4);
                // if(normalize) {v3_histos_var_sample_syst[ivar][isample][isyst]->SetMaximum(0.5);}
                if(normalize) {v3_histos_var_sample_syst[ivar][isample][isyst]->SetMaximum(v3_histos_var_sample_syst[ivar][isample][isyst]->GetMaximum()*1.2);}
                v3_histos_var_sample_syst[ivar][isample][isyst]->SetMinimum(0.001);

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
                    else if(v_groups[isample] == "TTGamma_Dilep") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "X#gamma", "L");}
                    else if(v_groups[isample] == "DY" ) {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "V+jets", "L");}
                    else if(v_groups[isample] == "TTbar_DiLep" || v_groups[isample] == "NPL") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "t#bar{t}", "L");}
                    else {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], v_groups[isample], "L");}
    			}

    			//-- Syst (can hardcode names)
                if(v_syst[isyst] != "") qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], v_syst[isyst], "L");
    			// if(v_syst[isyst] == "JESUp") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "JES Up", "L");}
    			// if(v_syst[isyst] == "JESDown") {qw->AddEntry(v3_histos_var_sample_syst[ivar][isample][isyst], "JES Down", "L");}

                //-- Debug printouts
                // for(int ibin=1; ibin<v3_histos_var_sample_syst[ivar][isample][isyst]->GetNbinsX()+1; ibin++) {cout<<"var="<<total_var_list[ivar]<<"/sample="<<v_samples[isample]<<"/syst="<<v_syst[isyst]<<"/ BIN("<<ibin<<") "<<v3_histos_var_sample_syst[ivar][isample][isyst]->GetBinContent(ibin)<<endl;}
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

        vector<TH1F*> v_histos_ratio; //For each sample *but central SM sample*, create 1 ratio histogram
        TH1F* histo_ratio_denominator = NULL;
        TPad* pad_ratio = NULL;
        if(v3_histos_var_sample_syst[ivar].size()>=2 || v_syst.size()>1)
        {
    		//-- create subpad to plot ratio
    		pad_ratio = new TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 1.0);
    		pad_ratio->SetTopMargin(0.75);
    		pad_ratio->SetFillColor(0);
    		pad_ratio->SetFillStyle(0);
    		pad_ratio->SetGridy(1);
    		pad_ratio->Draw();
    		pad_ratio->cd(0);

            histo_ratio_denominator = (TH1F*) v3_histos_var_sample_syst[ivar][0][0]->Clone();
            // TH1F* histo_ratio = (TH1F*) v3_histos_var_sample_syst[ivar][1][0]->Clone();
            // histo_ratio->Divide(histo_ratio_denominator);

            if(v_samples.size()>=2)
            {
                for(int ihisto=0; ihisto<v_samples.size(); ihisto++)
                {
                    if(!ihisto) {continue;} //No SM/SM ratio histo
                    v_histos_ratio.push_back((TH1F*) v3_histos_var_sample_syst[ivar][ihisto][0]->Clone()); //EFT sample
                    v_histos_ratio[ihisto-1]->Divide(histo_ratio_denominator); //Divide by central SM sample
                }
            }
            else if(v_syst.size()>1)
            {
                for(int isyst=0; isyst<v_syst.size(); isyst++)
                {
                    if(!isyst) {continue;} //No SM/SM ratio histo
                    v_histos_ratio.push_back((TH1F*) v3_histos_var_sample_syst[ivar][0][isyst]->Clone()); //EFT sample
                    v_histos_ratio[isyst-1]->Divide(histo_ratio_denominator); //Divide by central SM sample

                    //-- Debug printouts
                    // for(int ibin=1; ibin<v_histos_ratio[isyst-1]->GetNbinsX()+1; ibin++) {cout<<"var="<<total_var_list[ivar]<<"/syst="<<v_syst[isyst]<<"/ BIN("<<ibin<<") "<<v_histos_ratio[isyst-1]->GetBinContent(ibin)<<endl;}
                }
            }

    		v_histos_ratio[0]->GetYaxis()->SetTitle("Ratio");
    		v_histos_ratio[0]->GetYaxis()->SetTickLength(0.);
            v_histos_ratio[0]->GetYaxis()->SetTitleOffset(1.15);
            v_histos_ratio[0]->GetYaxis()->SetLabelSize(0.04);
            v_histos_ratio[0]->GetXaxis()->SetLabelSize(0.04);
    		v_histos_ratio[0]->GetXaxis()->SetLabelFont(42);
    		v_histos_ratio[0]->GetYaxis()->SetLabelFont(42);
    		v_histos_ratio[0]->GetXaxis()->SetTitleFont(42);
    		v_histos_ratio[0]->GetYaxis()->SetTitleFont(42);
            v_histos_ratio[0]->GetYaxis()->SetNdivisions(507); //grid drawn on primary tick marks only
    		v_histos_ratio[0]->GetXaxis()->SetNdivisions(505);
    		v_histos_ratio[0]->GetYaxis()->SetTitleSize(0.06);
    		v_histos_ratio[0]->GetXaxis()->SetTickLength(0.04);

    		//NB : when using SetMaximum(), points above threshold are simply not drawn
    		//So for ratio plot, even if point is above max but error compatible with 1, point/error bar not represented!
            for(int iratio=0; iratio<v_histos_ratio.size(); iratio++)
            {
                v_histos_ratio[iratio]->SetMinimum(0.7);
                v_histos_ratio[iratio]->SetMaximum(1.3);
                if(v_samples.size()==1 && v_syst.size()>1) {v_histos_ratio[iratio]->SetMinimum(0.90); v_histos_ratio[iratio]->SetMaximum(1.10);} //zoom
            }

    		if(drawInputVars) {v_histos_ratio[0]->GetXaxis()->SetTitle(Get_Variable_Name(total_var_list[ivar]));}
    		else
    		{
                v_histos_ratio[0]->GetXaxis()->SetTitle(Get_Template_XaxisTitle(total_var_list[ivar], this->use_paperStyle, this->use_NN_cpq3_SRttZ)); //Default

    			if(template_name == "categ") //Vertical text X labels (categories names)
    			{
                    //Hard-coded -- should automate labels
                    {
                        const char *labels[10]  = {"1bj,2j","1bj,3j","1bj,4j","1bj,5j","1bj,6j","2bj,2j","2bj,3j","2bj,4j","2bj,5j","2bj,6j"};
                        for(int i=1;i<=10;i++) {v_histos_ratio[0]->GetXaxis()->SetBinLabel(i,labels[i-1]);}
                    }

                    v_histos_ratio[0]->GetXaxis()->SetTitle("");
                    // histo_ratio->GetXaxis()->SetTitle("Categ.");
    				v_histos_ratio[0]->GetXaxis()->SetLabelSize(0.06);
    				v_histos_ratio[0]->GetXaxis()->SetLabelOffset(0.02);
                    v_histos_ratio[0]->LabelsOption("v", "X"); //X labels vertical
    			}
    		}

            //-- Draw SMEFT/SM histos
    		pad_ratio->cd(0);
            if(v_samples.size()>=2)
            {
                for(int isample=0; isample<v_histos_ratio.size(); isample++)
                {
                    v_histos_ratio[isample]->Draw("hist E same");
                }
            }
            else if(v_syst.size()>1)
            {
                for(int isyst=0; isyst<v_histos_ratio.size(); isyst++)
                {
                    v_histos_ratio[isyst]->Draw("hist same");

                    //-- Debug printouts
                    // cout<<"v_histos_ratio[isyst]->GetMaximum() "<<v_histos_ratio[isyst]->GetMaximum()<<endl;
                    // v_histos_ratio[isyst]->SetMinimum(0.7);
                    // v_histos_ratio[isyst]->SetMaximum(1.3);
                    // for(int ibin=1; ibin<v_histos_ratio[isyst]->GetNbinsX()+1; ibin++) {cout<<"var="<<total_var_list[ivar]<<" / syst="<<v_syst[isyst+1]<<" / BIN("<<ibin<<") "<<v_histos_ratio[isyst]->GetBinContent(ibin)<<endl;}
                }
            }

            //Draw SM/SM markers (centered at 1)
            histo_ratio_denominator->Divide(histo_ratio_denominator);
            histo_ratio_denominator->Draw("E same");
            // for(int ibin=1; ibin<histo_ratio_denominator->GetNbinsX()+1; ibin++)
            // {
            //     cout<<"bin "<<ibin<<" / "<<histo_ratio_denominator->GetBinContent(ibin)<<" / "<<histo_ratio_denominator->GetBinError(ibin)<<" / "<<histo_ratio_denominator->GetBinError(ibin)/histo_ratio_denominator->GetBinContent(ibin)<<endl;
            // }
        }


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
    	latex.DrawLatex(l + 0.01, 0.91, cmsText);

    	TString extraText = "Preliminary";
    	latex.SetTextFont(52);
    	latex.SetTextSize(0.05);
    	if(!use_paperStyle) {latex.DrawLatex(l + 0.12, 0.91, extraText);}

		TString lumi_ts = Convert_Number_To_TString(lumi);
		lumi_ts += " fb^{-1} (13 TeV)";
		latex.SetTextFont(42);
		latex.SetTextAlign(31);
		latex.SetTextSize(0.04);
        latex.DrawLatex(0.72, 0.91,lumi_ts);

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
        outdir+= luminame_tmp + "/";
        mkdir(outdir.Data(), 0777);
        if(cat_tmp != "")
        {
            outdir+= cat_tmp + "/";
            mkdir(outdir.Data(), 0777);
        }

    	//Output
    	TString output_plot_name = outdir + total_var_list[ivar] +"_templatesShapes";
    	if(channel != "") {output_plot_name+= "_" + channel;}
        output_plot_name+= this->filename_suffix;
        if(use_paperStyle) {output_plot_name+= ".eps";} //Use .eps to have correct \ell labels
    	else {output_plot_name+= this->plot_extension;}

    	c->SaveAs(output_plot_name);

        delete c; c = NULL;
        delete qw; qw = NULL;
        // delete histo_ratio; histo_ratio = NULL;
        delete histo_ratio_denominator; histo_ratio_denominator = NULL;

        if(v_samples.size()>=2 || v_syst.size()>1)
        {
            for(int isample=0; isample<v_histos_ratio.size(); isample++)
            {
                delete v_histos_ratio[isample]; v_histos_ratio[isample] = NULL;
            }
        }
    } //var loop

    for(int ivar=0; ivar<v3_histos_var_sample_syst.size(); ivar++)
    {
        if(v3_histos_var_sample_syst[ivar][0][0] == NULL) {continue;} //In case this variable was not found
    	for(int isample=0; isample<v3_histos_var_sample_syst[ivar].size(); isample++)
    	{
    		for(int isyst=0; isyst<v_syst.size(); isyst++) {delete v3_histos_var_sample_syst[ivar][isample][isyst];}
        }
    }

    file_input->Close();

	return;
} //Compare_TemplateShapes_Processes()

















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
    else if(systname.BeginsWith("prefir"))
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
    else if(systname.BeginsWith("jetPUID"))
    {
        address_memberArray = array_jetPileupID;
        array_name = "varWeightJetPileupID";
        if(systname.EndsWith("EffDown")) {index = 1;}
        else if(systname.EndsWith("EffUp")) {index = 0;}
        else if(systname.EndsWith("MTDown")) {index = 3;}
        else if(systname.EndsWith("MTUp")) {index = 2;}
    }
    else if(systname.BeginsWith("FR"))
    {
        address_memberArray = array_fakeFactor;
        array_name = "varWeightFakeFactor";
        if(systname.Contains("_pt"))
        {
            if(systname.EndsWith("Down")) {index = 2;}
            else if(systname.EndsWith("Up")) {index = 3;}
        }
        else if(systname.Contains("_be"))
        {
            if(systname.EndsWith("Down")) {index = 4;}
            else if(systname.EndsWith("Up")) {index = 5;}
        }
        else //FR norm. variations //For David's FR, this is the only set of variations (first 2 array elements)
        {
            if(systname.EndsWith("Down")) {index = 0;}
            else if(systname.EndsWith("Up")) {index = 1;}
        }

        if(systname.Contains("m_")) {index+= 6;} //Convention: first 6 elements correspond to FR ele, last 6 to FR mu
    }
    else if(systname.BeginsWith("PDF"))
    {
        address_memberArray = array_PDFtotal;
        array_name = "varWeightPdfTotal";
        if(systname.EndsWith("Down")) {index = 1;}
        else if(systname.EndsWith("Up")) {index = 0;}
    }
    else if(systname.BeginsWith("ME"))
    {
        //-- NB: there are 6 physical variations stored; should plot them, and use only that with the largest effect (enveloppe)
        // var_weight_me[0] = Event.weightMEScaleUp;
        // var_weight_me[1] = Event.weightMEScaleDown;
        // var_weight_me[2] = Event.weightMEFacScaleUp;
        // var_weight_me[3] = Event.weightMEFacScaleDown;
        // var_weight_me[4] = Event.weightMERenScaleUp;
        // var_weight_me[5] = Event.weightMERenScaleDown;

        address_memberArray = array_ME;
        array_name = "varWeightMe";
        if(systname.EndsWith("Down")) {index = 1;}
        else if(systname.EndsWith("Up")) {index = 0;}
    }
    else if(systname.BeginsWith("alphas"))
    {
        address_memberArray = array_alphaS;
        array_name = "varWeightAlphas";
        if(systname.EndsWith("Down")) {index = 1;}
        else if(systname.EndsWith("Up")) {index = 0;}
    }
    else if(systname.BeginsWith("ISR"))
    {
        address_memberArray = array_partonShower;
        array_name = "varWeightPs";
        if(systname.EndsWith("Down")) {index = 1;}
        else if(systname.EndsWith("Up")) {index = 0;}
    }
    else if(systname.BeginsWith("FSR"))
    {
        address_memberArray = array_partonShower;
        array_name = "varWeightPs";
        if(systname.EndsWith("Down")) {index = 3;}
        else if(systname.EndsWith("Up")) {index = 2;}
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
        if(systname.Contains("Loose")) {address_memberArray = array_LepEffLoose_el; array_name = "varWeightElectronLoose";}
        else {address_memberArray = array_LepEffTight_el; array_name = "varWeightElectronTight";}
        if(systname.EndsWith("Down")) {index = 0;}
        else if(systname.EndsWith("Up")) {index = 1;}
    }

    else if(systname.Contains("njets_tZq")) {return;} //Handled differently (using class member SF variable)

    else {cout<<FRED("ERROR ! Systematic '"<<systname<<"' not included in function SetBranchAddress_SystVariation() from Helper.cxx ! Can *not* compute it !")<<endl; return;}

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
void TopEFT_analysis::MergeSplit_Templates(bool makeHisto_inputVars, TString filename, vector<TString> total_var_list, TString template_name, TString category, bool force_normTemplate_positive)
{
    bool store_countExp_SMvsEFT = true; //true <-> store full histograms as single bins for later comparison with counting experiment in Combine

//--------------------------------------------
//-- Automated

    if(makeHisto_inputVars || this->scanOperators_paramNN) {store_countExp_SMvsEFT = false;} //Un-necessary

//--------------------------------------------

    cout<<endl<<FYEL("==> Merging/splitting histograms in TFile : ")<<filename<<endl;

	if(!Check_File_Existence(filename) ) {cout<<endl<<FRED("File "<<filename<<" not found! Abort merging procedure !")<<endl; return;}
	TFile* f = TFile::Open(filename, "UPDATE");

    // int n_total_histos = Count_nofHistos_inTFile(f); //Too slow for huge files
    int counter = 0; //Count nof read histos

	//NB :here the order of loops is important because we sum histograms recursively ! The 'sample_list' loop *must be the most nested one* !
    for(int ivar=0; ivar<total_var_list.size(); ivar++) //There may be more than 1 template (e.g. several NN output nodes)
    {
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
            // cout<<"v_lumiYears[iyear] "<<v_lumiYears[iyear]<<endl;

            for(int ichan=0; ichan<channel_list.size(); ichan++)
        	{
                // cout<<"channel_list[ichan] "<<channel_list[ichan]<<endl;

        		for(int itree=0; itree<systTree_list.size(); itree++)
        		{
        			if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name && channel_list.size() > 1 && channel_list[ichan] == "") {continue;}

        			for(int isyst=0; isyst<syst_list.size(); isyst++)
        			{
                        // cout<<"syst_list[isyst] "<<syst_list[isyst]<<endl;

        				// if(((channel_list.size() > 1 && channel_list[ichan] == "") || (systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name)) && syst_list[isyst] != "") {continue;}
        				if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name && syst_list[isyst] != "") {break;}
                        else if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefir") ) {continue;} //no prefire in 2018

        				TH1F* h_merging = NULL; //Tmp merged histogram

                        //-- For EFT strategies, need to store all histogram bins separately
                        int n_singleBins = 0; //Default: will only merge: a) entire histograms, b) histograms stored as single bin //Gets updated below depending on strategy
                        for(int ibin=-1; ibin<n_singleBins+1; ibin++) //Convention (depending on NN_strategy): bin==-1 <-> merge full histogram; bin==0 <-> merge histo stored as single bin (counting exp.); bin>0 <-> merge corresponding template bin (individually)
                        {
                            // cout<<"ibin "<<ibin<<" ("<<"n_singleBins "<<n_singleBins<<")"<<endl;

                            if(!store_countExp_SMvsEFT && ibin==0) {continue;} //Choose not to store histograms as single bin (countExp, for comparisons) for speed
                            else if(!make_SMvsEFT_templates_plots && ibin == 0) {break;} //SM vs SM templates: don't need per-bin nor single-bin (countExp) histograms
                            else if(ibin==0 && (this->make_fixedRegions_templates || total_var_list[ivar].Contains("countExp"))) {break;} //countExp: no need to split per bin !
                            else if(ibin==0 && (n_singleBins == 1 || total_var_list[ivar].Contains("countExp")) ) {break;} //Idem if the full template has only 1 bin
                            else if(makeHisto_inputVars && ibin==0) {break;} //Input features: don't need to split per bin (?)

            				for(int isample=0; isample<sample_list.size(); isample++)
            				{
                            	//-- Protections : not all syst weights apply to all samples, etc.
                                if(sample_list[isample] == "DATA" && ((systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name) || syst_list[isyst] != "")) {continue;} //nominal data only
                                else if(makeHisto_inputVars && sample_groups[isample] != "NPL") {continue;} //For control plots, only need to substract prompt NPL from data-driven NPL
                                else if((sample_list[isample].Contains("NPL") && syst_list[isyst] != "" && !syst_list[isyst].BeginsWith("FR")) || (!sample_list[isample].Contains("NPL") && syst_list[isyst].BeginsWith("FR"))) {continue;} //NPL <-> only fakes sytematics; all others <-> no fakes systematics
                                else if((syst_list[isyst].BeginsWith("PDF") || syst_list[isyst].BeginsWith("ME") || syst_list[isyst].BeginsWith("alpha") || syst_list[isyst].BeginsWith("ISR") || syst_list[isyst].BeginsWith("FSR")) && !sample_list[isample].Contains("PrivMC") && sample_list[isample] != "tZq" && sample_list[isample] != "ttZ") {continue;}  //NB: no TH weights available in central tWZ V12
                                else if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name && (sample_list[isample] == "DY" || sample_list[isample].Contains("TTbar") || sample_list[isample].Contains("NPL")) ) {continue;}
                                else if(syst_list[isyst].Contains("njets_tZq") && !sample_list[isample].Contains("PrivMC_tZq")) {continue;} //Only applies to LO tZq
                                else if(sample_list[isample].Contains("PrivMC") && (syst_list[isyst] != "" || (systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name)) && this->make_fixedRegions_templates && !total_var_list[ivar].Contains("ttZ4l")) {continue;} //Don't write shifted PrivMC histos in CRs (only need nominals to avoid bug in combine)

            					//-- Check if this sample needs to be merged, i.e. if the samples before/after belong to the same "group of samples"
            					bool merge_this_sample = false;
            					if(!isample && sample_groups.size() > 1 && sample_groups[isample+1] == sample_groups[isample]) {merge_this_sample = true;}
            					else if(isample > 0 && isample == sample_list.size()-1 && sample_groups[isample-1] == sample_groups[isample]) {merge_this_sample = true;}
            					else if(isample > 0 && isample < sample_list.size()-1 && (sample_groups[isample+1] == sample_groups[isample] || sample_groups[isample-1] == sample_groups[isample])) {merge_this_sample = true;}

                                // cout<<endl<<"Year "<<v_lumiYears[iyear]<< " / tree "<<systTree_list[itree]<<" / chan "<<channel_list[ichan]<<" / syst "<<syst_list[isyst]<<" / sample "<<sample_list[isample]<<endl;
            					// cout<<"merge_this_sample "<<merge_this_sample<<endl;
            					if(!merge_this_sample) {continue;} //Only care about samples to merge : others are already stored in file

                                TString histoname = ""; //Default: merge full histograms
                                if(ibin==0) {histoname+= "countExp_";} //Merge full histos stored as single bins (counting exp.)
                                else if(n_singleBins > 1 && ibin>0) {histoname+= (TString) "bin"+Form("%d",ibin)+"_";} //Merge bin per bin
                                histoname+= total_var_list[ivar];
            					if(channel_list[ichan] != "") {histoname+= "_" + channel_list[ichan];}
                                if(region != "" && !makeHisto_inputVars && categorization_strategy==0) {histoname+= "_" + region;}
                                histoname+= "_" + v_lumiYears[iyear] + "__" + (sample_list[isample] == "DATA"? "data_obs":sample_list[isample]);
            					if(syst_list[isyst] != "" || (systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name)) {histoname+= "__" + Get_Modified_SystName(syst_list[isyst]+systTree_list[itree], v_lumiYears[iyear], sample_list[isample]);}
            					else if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name) {histoname+= "__" + systTree_list[itree];}

                                //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                                histoname.ReplaceAll('-', 'm');

                                // cout<<"histoname = "<<histoname<<endl;

                                if(!f->GetListOfKeys()->Contains(histoname) )
            					{
            						if(systTree_list[itree] == "" && syst_list[isyst] == "") {cout<<endl<<DIM("Histo "<<histoname<<" not found in file "<<filename<<" !")<<endl;}
                                    continue;
                                }

                                //NB -- very slow for large files !
            					TH1F* h_tmp = (TH1F*) f->Get(histoname); //Get individual histograms

                                counter++; //Increment counter of read histos
                                if(counter % 1000 == 0) {cout<<DIM("Read "<<counter<<" histograms...")<<endl;}
                                // if(counter % 1000 == 0) {cout<<DIM("Read "<<counter<<"/"<<n_total_histos<<" histograms...")<<endl;}

                                if(make_SMvsEFT_templates_plots) //For EFT templates, will also consider inidividual-bin histos
                                {
                                    if(n_singleBins==0 && ibin==-1) //Only read the binning once, from full histograms
                                    {
                                        n_singleBins = h_tmp->GetNbinsX(); //Update the limit for the for-loop *within the loop* (--> first read full histo to infer the correct binning)
                                    }
                                }

                                // cout<<"h_merging "<<h_merging<<endl;
            					if(h_tmp)
            					{
            						if(!h_merging) {h_merging = (TH1F*) h_tmp->Clone();}
            						else {h_merging->Add(h_tmp);}
                                    // cout<<"h_tmp->Integral() = "<<h_tmp->Integral()<<endl;
            					}
            					// else {cout<<DIM("h_tmp is null ! ("<<histoname<<")")<<endl;}
            					// cout<<"h_merging->Integral() = "<<h_merging->Integral()<<endl;

            					delete h_tmp; h_tmp = NULL;
            					if(!h_merging) {cout<<"Syst "<<syst_list[isyst]<<systTree_list[itree]<<" / chan "<<channel_list[ichan]<<" / sample "<<sample_list[isample]<<endl; cout<<"h_merging is null ! Latest histo read "<<histoname<<"... You should check this !"<<endl; continue;}

            					//Check if next sample will be merged with this one, or else if must write the histogram
            					if(isample < sample_list.size()-1 && sample_groups[isample+1] == sample_groups[isample]) {continue;}
            					else
            					{
                                    TString histoname_new = "";
                                    if(ibin==0) {histoname_new+= "countExp_";} //Merge full histos stored as single bins (counting exp.)
                                    else if(n_singleBins > 1 && ibin>0) {histoname_new+= (TString) "bin"+Form("%d",ibin)+"_";} //Merge bin per bin
                                    histoname_new+= total_var_list[ivar];
                                    // if(category != "") {histoname_new+= "_" + category;}
                                    if(channel_list[ichan] != "") {histoname_new+="_"  + channel_list[ichan];}
                                    if(region != "" && !makeHisto_inputVars && categorization_strategy==0) {histoname_new+= "_" + region;}
                                    histoname_new+= "_" + v_lumiYears[iyear] + "__" + sample_groups[isample];
                                    if(syst_list[isyst] != "" || (systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name)) {histoname_new+= "__" + Get_Modified_SystName(syst_list[isyst]+systTree_list[itree], v_lumiYears[iyear], sample_groups[isample]);}
                                    else if(systTree_list[itree] != "" && systTree_list[itree] != nominal_tree_name) {histoname_new+= "__" + systTree_list[itree];}

                                    //-- Protection: replace '-' (hyphen) with 'm' character (hyphen in histo name causes errors at reading)
                                    histoname_new.ReplaceAll('-', 'm');

            						if(force_normTemplate_positive && sample_list[isample] != "DATA")
            						{
                                        Avoid_Histogram_EmptyOrNegativeBins(h_merging);
            						}
            						// cout<<"h_merging->Integral() = "<<h_merging->Integral()<<endl;

            						f->cd();
            						h_merging->Write(histoname_new, TObject::kOverwrite);
            						// cout<<"-- Writing merged histo "<<histoname_new<<" with integral "<<h_merging->Integral()<<endl;

            						delete h_merging; h_merging = NULL;

                                    //-- Special case: for control histograms/plots, want to substract NPL_MC from (data-driven) NPL --> then overwrite "NPL" and delete "NPL_MC" (to avoid ambiguities)
                                    //-- Necessary ? Slow ? Could also delete NPL_DATA... ?
                                    // if(sample_list[isample] == "NPL_MC")
                                    // {
                                    //     f->Delete(histoname+";1"); //Delete (first cycle of) histogram
                                    //     cout<<DIM("Merged and deleted histogram "<<histoname+";1"<<"")<<endl;
                                    // }

                                } //write histo
                            } //sample loop
                        } //bin loop

        			} //syst loop
        		} //tree loop
        	} //channel loop
        } //years loop
    } //vars loop

    // cout<<DIM("(Total: "<<n_total_histos<<" histograms)")<<endl;

	f->Close();

    // cout<<endl<<FYEL("Updated file: " )<<filename<<endl<<endl<<endl;
    cout<<FYEL("... Done")<<" (file: "<<filename<<")"<<endl<<endl<<endl;

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
 * NB: partially rewrote the NN implementation (avoid un-necessary copies) for speed up //NB: not much faster
 */
bool TopEFT_analysis::Get_VectorAllEvents_passMVACut(vector<int>& v_maxNode, vector<vector<float>>& v_values, TString signal, TString classifier_name, TTree* tree, TString year, float cut_value, bool keep_aboveCut, bool use_specificMVA_eachYear, int categorization_strategy, bool MVA_EFT, int nentries_max, TString event_cat, bool also_applyCut_onMaxNodeValue, bool isFake)
// bool TopEFT_analysis::Get_VectorAllEvents_passMVACut(vector<int>& v_maxNode, vector<vector<float>>& v_values, TString signal, TString classifier_name, TString tree_name, TString input_file_path, TString year, float cut_value, bool keep_aboveCut, bool use_specificMVA_eachYear, int categorization_strategy, bool MVA_EFT, int nentries_max, TString event_cat, bool also_applyCut_onMaxNodeValue, bool isFake)
{
    cout<<endl<<FYEL("=== Filling vector with specific MVA information (pass/fail cut or max. node)... ")<<endl;
    // cout<<DIM("(File : "<<input_file_path<<")")<<endl;

    if(classifier_name == "") {classifier_name = "NN";}
    if(classifier_name != "BDT" && classifier_name != "NN") {cout<<BOLD(FRED("ERROR: wrong [classifier_name] option !"))<<endl; return false;}

    int nevents_passingCut = 0;
    TString BDT_method_name = "BDT";
    vector<TString> var_list_tmp;
    vector<float> var_floats_tmp;

    TString MVA_input_path = Get_MVAFile_InputPath(classifier_name, signal, year, use_specificMVA_eachYear, MVA_EFT, false, categorization_strategy, this->scanOperators_paramNN, this->use_NN_cpq3_SRttZ);
    if(MVA_input_path == "") {cout<<"MVA input file not found ! "<<endl; return false;} //MVA file not found

    //Use local variables to avoid conflict with class members (different MVAs)
    vector<TString> var_list_NN; TString NN_inputLayerName = ""; TString NN_outputLayerName = ""; int NN_iMaxNode = -1; int NN_nNodes = -1;

    TMVA::Reader* reader_tmp = NULL; //TMVA BDT reader
    TFModel* clfy_tmp = NULL; //TF NN reader

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
        TString NNinfo_input_path = Get_MVAFile_InputPath(classifier_name, signal, year, use_specificMVA_eachYear, MVA_EFT, true, categorization_strategy, this->scanOperators_paramNN, this->use_NN_cpq3_SRttZ);
        if(NNinfo_input_path == "") {cout<<"MVA info file not found ! "<<endl; return false;} //MVA file not found

        //-- Load neural network model
        if(Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN, v_NN_nodeLabels, NN_inputLayerName, NN_outputLayerName, NN_iMaxNode, NN_nNodes, minmax_bounds))
        {
            // cout<<"clfy_tmp "<<clfy_tmp<<endl;
            // cout<<"MVA_input_path "<<MVA_input_path<<endl;
            // cout<<"var_list_NN.size() "<<var_list_NN.size()<<endl;
            // cout<<"NN_inputLayerName "<<NN_inputLayerName<<endl;
            // cout<<"NN_outputLayerName "<<NN_outputLayerName<<endl;
            // cout<<"NN_nNodes "<<NN_nNodes<<endl;
            clfy_tmp = new TFModel(MVA_input_path.Data(), var_list_NN.size(), NN_inputLayerName.Data(), NN_nNodes, NN_outputLayerName.Data());
        }
        else {cout<<"Missing NN information ! "<<endl; return false;} //Error: missing NN infos

        var_list_tmp = var_list_NN;
        var_floats_tmp.resize(var_list_tmp.size());
    }

    //-- Create an input tensor
    long long int n_inputs = var_list_NN.size() > 0? var_list_NN.size():1;
    tensorflow::Tensor input(tensorflow::DT_FLOAT, { 1, n_inputs }); // single batch of dimension 10
    std::vector<tensorflow::Tensor> outputs; //Store outputs

    /* //CHANGED
	if(!Check_File_Existence(input_file_path) ) {cout<<BOLD(FRED("ERROR: "<<input_file_path<<" not found!"))<<endl; return false;}
	TFile* file_input = TFile::Open(input_file_path, "READ");

	TTree* tree = (TTree*) file_input->Get(tree_name);
    if(!tree) {cout<<BOLD(FRED("ERROR :"))<<" file "<<input_file_path<<", tree "<<tree_name<<" is NULL !"<<endl; return false;}
    if(!tree->GetEntries()) {cout<<"File "<<input_file_path<<" / tree "<<tree_name<<" is EMPTY !"<<endl; return false;}
    */

    //-- Set addresses of input features
    tree->SetBranchStatus("*", 0); //Disable all branches by default, speeds up considerably
	for(int ivar=0; ivar<var_list_tmp.size(); ivar++)
    {
        if(var_list_tmp[ivar] == "ctz" || var_list_tmp[ivar] == "ctw" || var_list_tmp[ivar] == "cpq3" || var_list_tmp[ivar] == "cpqm" || var_list_tmp[ivar] == "cpt") {continue;} //WC input values are arbitrary, there is no address to set !
        tree->SetBranchStatus(var_list_tmp[ivar], 1); //Activate only necessary branches
        tree->SetBranchAddress(var_list_tmp[ivar], &var_floats_tmp[ivar]); //FIXCMSSW
        // tree->SetBranchAddress(var_list_tmp[ivar], &input.matrix<float>()(0, ivar)); //Fill tensor directly for speed up //FIXLOCAL
    }

    //-- May cut on an 'event category flag' whose name is given as argument (<-> no need to evaluate MVA for events which do not enter the region of interest)
    Char_t is_goodCategory = true;
    TString cat_name = Get_Category_Boolean_Name(event_cat, isFake);
    if(cat_name != "") {tree->SetBranchStatus(cat_name, 1); tree->SetBranchAddress(cat_name, &is_goodCategory);}

	int nentries = tree->GetEntries();
    if(nentries_max > 0 && nentries > nentries_max) {nentries = nentries_max;}
    v_maxNode.clear(); v_maxNode.resize(nentries); std::fill(v_maxNode.begin(),v_maxNode.end(),0); //Fill with '0' for all tree entries
    v_values.clear(); v_values.resize(nentries);
    for(int jentry=0; jentry<v_values.size(); jentry++)
    {
        v_values[jentry].resize(NN_nNodes);
    }

    std::vector<float> clfy_outputs;

	for(int ientry=0; ientry<nentries; ientry++)
	{
		if(ientry && ientry%50000==0) {cout<<DIM(" --- "<<ientry<<" / "<<nentries<<"")<<endl;}
		tree->GetEntry(ientry);

        if(!is_goodCategory) {continue;}

        // for(int ivar=0; ivar<var_list_tmp.size(); ivar++) {cout<<ivar<<" / "<<var_list_tmp[ivar]<<" / "<<var_floats_tmp[ivar]<<endl;} //Debug printouts

        float mva_output = 0.;
        if(classifier_name == "BDT") {mva_output = reader_tmp->EvaluateMVA(BDT_method_name);}
        else //NN
        {
            //-- NB: if get segfault here of type 'Incompatible shapes: [1,105] vs. [35]' <-> means that the NN info file and actual .pb model are incompatible (need to retrain NN)
            clfy_outputs = clfy_tmp->evaluate(var_floats_tmp); //Evaluate output node(s) value(s) //Slow... ! //FIXCMSSW
            // clfy_tmp->evaluate_fast(input, outputs); //Evaluate output node(s) value(s) //CHANGED -- overloaded function avoids un-necessary copies //FIXLOCAL

            NN_iMaxNode = -1;
            for(int inode=0; inode<NN_nNodes; inode++)
            {
                v_values[ientry][inode] = clfy_outputs[inode];
                if(clfy_outputs[inode] > mva_output) {mva_output = clfy_outputs[inode]; NN_iMaxNode = inode;} //FIXCMSSW
                // cout<<"clfy_outputs[inode] "<<clfy_outputs[inode]<<" / mva_output "<<mva_output<<" / NN_iMaxNode "<<NN_iMaxNode<<endl;

                // if(outputs[0].matrix<float>()(0,inode) > mva_output) {mva_output = outputs[0].matrix<float>()(0,inode); NN_iMaxNode = inode;} //FIXLOCAL
                // cout<<"outputs[0].matrix<float>()(0,inode) "<<outputs[0].matrix<float>()(0,inode)<<" / mva_output "<<mva_output<<" / NN_iMaxNode "<<NN_iMaxNode<<endl;
            }

            //Multiclass --> Store max. node information
            if(NN_nNodes > 1)
            {
                // if(!also_applyCut_onMaxNodeValue || ((keep_aboveCut && mva_output >= cut_value) || (!keep_aboveCut && mva_output < cut_value))) {v_maxNode[ientry] = NN_iMaxNode;} //Don't use this function for now, for speed up
                v_maxNode[ientry] = NN_iMaxNode;
                continue;
            }
        }

        //Binary classifier --> Determine/store whether the event passes the MVA cut or not
        bool pass_cut = false;
        if((keep_aboveCut && mva_output >= cut_value) || (!keep_aboveCut && mva_output < cut_value)) {pass_cut = true; nevents_passingCut++;}
        v_maxNode[ientry] = pass_cut;

        if(isnan(mva_output)) {cout<<BOLD(FRED("ERROR: MVA = NaN ! Problem with NN evaluation ? Check it please !"))<<endl;} //Protection
	} //loop on entries

    if(v_maxNode.size() != nentries) {cout<<BOLD(FRED("Wrong number of entries in MVA cut vector ! Check it please !"))<<endl; return false;}
	// file_input->Close();

	// cout<<FMAG("---- Vector containing BDTfakeSR cut results is filled !")<<endl;
	if(NN_nNodes <= 1) {cout<<DIM("(Number of events passing the MVA cut : "<<nevents_passingCut<<")")<<endl;} //Meaningless for multiclass NN

    if(reader_tmp) {delete reader_tmp; reader_tmp = NULL;}
    if(clfy_tmp) {delete clfy_tmp; clfy_tmp = NULL;}

    cout<<FYEL("... Done ! ===")<<endl<<endl;

	return true;
}




















//--------------------------------------------
// ########     ###    ########  ######## ########     ########  ##        #######  ########  ######
// ##     ##   ## ##   ##     ## ##       ##     ##    ##     ## ##       ##     ##    ##    ##    ##
// ##     ##  ##   ##  ##     ## ##       ##     ##    ##     ## ##       ##     ##    ##    ##
// ########  ##     ## ########  ######   ########     ########  ##       ##     ##    ##     ######
// ##        ######### ##        ##       ##   ##      ##        ##       ##     ##    ##          ##
// ##        ##     ## ##        ##       ##    ##     ##        ##       ##     ##    ##    ##    ##
// ##        ##     ## ##        ######## ##     ##    ##        ########  #######     ##     ######
//--------------------------------------------


/**
 * Hard-coded function to produce figures for the paper.
 * Plot showing the postfit distributions in the regions common to all fits (CRs, SRttZ4l, SRother).
 */
void TopEFT_analysis::Make_PaperPlot_CommonRegions()
{
//--------------------------------------------

    cout<<endl<<BYEL("                          ")<<endl<<endl;
	cout<<FYEL("--- Producing paper figure [COMMON REGIONS] ---")<<endl;
    cout<<endl<<BYEL("                          ")<<endl<<endl;


//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

    //--------------------------
    // DIVIDE CANVAS IN 4 PADS
    //--------------------------

    //-- Canvas definition
    Load_Canvas_Style(); //Default top/bottom/left/right margins: 0.07/0.13/0.16/0.03
    TCanvas* c1 = new TCanvas("c1","c1", 1400, 600);
    c1->Divide(2, 1, 1E-11, 1E-11); //(x,y)
    // c1->Divide(2, 1, 0.01, 0.00001, 0); //See: https://root.cern/doc/master/classTPad.html#TPad:Divide

    //--------------------------
    // DEFINE 3 VARIABLES TO PLOT
    //--------------------------

    int nbins_tmp, nIndivBins; float xmin_tmp, xmax_tmp;

    int nvar = 2;
    vector<TString> total_var_list(nvar);
    total_var_list[0] = "countExp";
    total_var_list[1] = "mTW";

    //--------------------------
    // CREATE VECTORS OF OBJECTS
    //--------------------------

    TH1F *h_tmp = NULL;

    vector<TH1F*> v_hdata(nvar);
    vector<THStack*> v_stack(nvar);
    vector<TH1F*> v_histo_total_MC(nvar);
    vector<vector<TH1F*> > v_vector_MC_histo(nvar); //Store separately the histos for each MC sample --> stack them after loops
    vector<TH1F*> v_histo_ratio_data(nvar);
    vector<TH1F*> v_hlines1(nvar), v_hlines2(nvar); //Draw TLinesin ratio plots

    vector<TGraphAsymmErrors*> v_gr_error(nvar);
    vector<TGraphAsymmErrors*> v_gr_ratio_error(nvar);

    vector<TPad*> v_tpad_ratio(nvar);

    vector<TString> MC_samples_legend; //List the MC samples to mention in legend

    TLine* tline1; TLine* tline2;

    TString path_dir_shapes_binContent = "./outputs/dir_shapes_tmp/binContent/"; //Path of dir. containing shapes (with correct bin contents <-> did not freeze any parameter)
    TString path_dir_shapes_binError = "./outputs/dir_shapes_tmp/binError/"; //Path of dir. containing shapes (with correct bin errors <-> froze WC and problematic split JEC)


// #       ####   ####  #####   ####
// #      #    # #    # #    # #
// #      #    # #    # #    #  ####
// #      #    # #    # #####       #
// #      #    # #    # #      #    #
// ######  ####   ####  #       ####

//--------------------------------------------
	for(int ivar=0; ivar<nvar; ivar++)
	{
        cout<<endl<<"== VARIABLE: "<<total_var_list[ivar]<<endl;

    	TString primitive_name = "c1_" + Convert_Number_To_TString(ivar+1);
        // c1->GetListOfPrimitives()->Print();
        TPad* pad = (TPad*) c1->GetPrimitive(primitive_name);
        TString inputfilename = "";
    	if(total_var_list[ivar].Contains("countExp"))
    	{
            // inputfilename = "./outputs/shapes_postfit_otherRegions.root";
    		pad->SetPad(0, 0., 0.50, 1); //xlow, ylow, xup, yup
            pad->SetLogy();
    	}
    	else if(total_var_list[ivar].Contains("mTW"))
    	{
            // inputfilename = "./outputs/shapes_postfit_NN_5D.root";
    		pad->SetPad(0.50, 0., 1., 1); //xlow, ylow, xup, yup
            pad->SetRightMargin(0.04);
    	}
        pad->SetBottomMargin(0.30);

        // cout<<"-- Open "<<inputfilename<<endl;
        // TFile* file_input = TFile::Open(inputfilename);
        TFile* file_input = NULL;

        TH1F* h_tmp = NULL; //Tmp storing histo
        TH1F* hdata_tmp = NULL; //Tmp storing data histo
		// TH1F* h_sum_data = NULL; //Will store data histogram
		// vector<TH1F*> v_MC_histo; //Will store all MC histograms (1 TH1F* per MC sample)
        // TString data_histo_name = "";
		// TGraphAsymmErrors* g_data = NULL; //If using Combine file, data are stored in TGAE
		// TGraphAsymmErrors* g_tmp = NULL; //Tmp storing graph

		//-- Init error vectors
		double x, y, errory_low, errory_high;

		vector<double> v_eyl, v_eyh, v_exl, v_exh, v_x, v_y; //Contain the systematic errors (used to create the TGraphError)

        float bin_width = -1; //Get bin width of histograms for current variable

        //-- All histos are for given lumiYears and sub-channels --> Need to sum them all for plots
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {


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

                // cout<<"sample_list[isample] "<<sample_list[isample]<<endl;
                // cout<<"index_MC_sample "<<index_MC_sample<<endl;

				//-- In Combine, some individual contributions are merged as "Rares"/"EWK", etc.
				//-- If using Combine file, change the names of the samples we look for, and look only once for histogram of each "group"
				TString samplename = sample_list[isample];
				if(isample > 0 && sample_groups[isample] == sample_groups[isample-1]) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //if same group as previous sample, skip it
                else if(make_SMvsEFT_templates_plots && (sample_groups[isample] == "tZq" || sample_groups[isample] == "ttZ" || sample_groups[isample] == "tWZ")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //SM vs EFT --> use private signal samples
                else {samplename = sample_groups[isample];}

				//-- Protections, special cases
				if(sample_list[isample] == "DATA") {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}
                else if(!make_SMvsEFT_templates_plots && sample_list[isample].Contains("PrivMC")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //SM configuration --> only stack central samples (not private samples)
                else if(make_SMvsEFT_templates_plots && (sample_list[isample] == "tZq" || sample_list[isample] == "ttZ" || sample_list[isample] == "tWZ")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //EFT configuration --> only stack private samples (at SM point), not central samples
                else if(sample_list[isample] == "NPL_DATA")  {samplename = "NPL";} //Instead of 'NPL_DATA' and 'NPL_MC', we only want to read the merged histo 'NPL'
                else if(sample_list[isample] == "NPL_MC")  {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //NPL_MC gets substracted from NPL histograms and deleted --> Ignore this vector element //Remove ?

                //-- Add sample name to list (used for legend) //NB: add even if histo was not found and skipped, because expect that it will be found for some other year/channel/... But if not found at all, legend will be wrong
                if(iyear==0 && ivar==0 && samplename != "DATA")
                {
                    if(v_vector_MC_histo[ivar].size() <=  index_MC_sample) {MC_samples_legend.push_back(samplename);}
                }
                if(v_isSkippedSample[isample] == true) {continue;} //Skip this sample

				// cout<<endl<<UNDL(FBLU("-- Sample : "<<sample_list[isample]<<" : "))<<endl;

				h_tmp = NULL;
                if(total_var_list[ivar].Contains("mTW"))
                {
                    xmin_tmp = 0; xmax_tmp = 150; nIndivBins = 15;
                    h_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
                    v_eyl.resize(nIndivBins); v_eyh.resize(nIndivBins); v_exl.resize(nIndivBins); v_exh.resize(nIndivBins); v_y.resize(nIndivBins); v_x.resize(nIndivBins);
                    std::fill(v_y.begin(), v_y.end(), -1); //Init errors positions to <0 (invisible)

                    for(int ibin=1; ibin<nIndivBins+1; ibin++)
                    {
                        //-- Read bin content
                        //--------------------------------------------
                        inputfilename = path_dir_shapes_binContent + "shapes_postfit_datacard_bin"+Convert_Number_To_TString(ibin)+"_mTW_SRother.root";
                        if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                        file_input = TFile::Open(inputfilename, "READ");
                        // cout<<"inputfilename "<<inputfilename<<endl;
                        TString dir_hist = "bin" + Convert_Number_To_TString(ibin) + "_" + total_var_list[ivar] + "_SRother_" + v_lumiYears[iyear] + "_postfit/";
                        // cout<<"dir_hist/samplename "<<dir_hist<<samplename<<endl;
                        if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(samplename) ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<samplename<<"' not found ! Skip...")<<endl; continue;}
                        h_tmp->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+samplename))->GetBinContent(1)); //Get content/error from individual bin
                        // h_tmp->SetBinError(ibin, ((TH1F*) file_input->Get(dir_hist+samplename))->GetBinError(1));
                        file_input->Close();
                        //--------------------------------------------

                        //-- Read bin error
                        //--------------------------------------------
                        dir_hist = "postfit/"; //Reminder: read per-year histograms to get individual contributions from processes, *but* read total-sum histogram to get full postfit error
                        inputfilename = path_dir_shapes_binError + "shapes_postfit_datacard_bin"+Convert_Number_To_TString(ibin)+"_mTW_SRother.root";
                        if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                        file_input = TFile::Open(inputfilename, "READ");
                        if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains("TotalProcs") ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<"TotalProcs' not found ! Skip...")<<endl; continue;}
                        float bincontent = ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinContent(1); //Now read total yield (for error vector)
                        float binerror = ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinError(1);
                        file_input->Close();
                        //--------------------------------------------

                        //-- Errors
                        if(v_y[ibin-1] < 0) //Need to fill total error only once (not for each process)
                        {
                            float bin_width = (xmax_tmp-xmin_tmp) / nIndivBins;
                            v_x[ibin-1] = xmin_tmp + ((ibin-1)*bin_width) + (bin_width/2.);
                            v_y[ibin-1] = bincontent;
                            v_eyl[ibin-1] = binerror;
                            v_eyh[ibin-1] = binerror;
                            v_exl[ibin-1] = bin_width / 2; v_exh[ibin-1] = bin_width / 2;
                            // cout<<"bin "<<ibin<<" / "<<v_y[ibin-1]<<" / "<<v_eyl[ibin-1]<<endl;
                        }
                    } //nbins
                }
                else if(total_var_list[ivar].Contains("countExp"))
                {
                    xmin_tmp = 0; xmax_tmp = 3; nIndivBins = 3;
                    h_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
                    v_eyl.resize(nIndivBins); v_eyh.resize(nIndivBins); v_exl.resize(nIndivBins); v_exh.resize(nIndivBins); v_y.resize(nIndivBins); v_x.resize(nIndivBins);
                    std::fill(v_y.begin(), v_y.end(), -1); //Init errors positions to <0 (invisible)

                    for(int ibin=1; ibin<nIndivBins+1; ibin++)
                    {
                        TString cat = "";
                        if(ibin==1) {cat = "CRWZ";}
                        else if(ibin==2) {cat = "CRZZ";}
                        else if(ibin==3) {cat = "SRttZ4l";}

                        //-- Read bin content
                        //--------------------------------------------
                        inputfilename = path_dir_shapes_binContent + "shapes_postfit_datacard_countExp_"+cat+".root";
                        if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                        file_input = TFile::Open(inputfilename, "READ");
                        // cout<<"inputfilename "<<inputfilename<<endl;
                        TString dir_hist = total_var_list[ivar] + "_" + cat + "_" + v_lumiYears[iyear] + "_postfit/";
                        // cout<<"dir_hist/samplename "<<dir_hist<<samplename<<endl;
                        if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(samplename) ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<samplename<<"' not found ! Skip...")<<endl; continue;}
                        h_tmp->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+samplename))->GetBinContent(1)); //Get content/error from individual bin
                        // h_tmp->SetBinError(ibin, ((TH1F*) file_input->Get(dir_hist+samplename))->GetBinError(1));
                        file_input->Close();
                        //--------------------------------------------

                        //-- Read bin error
                        //--------------------------------------------
                        dir_hist = "postfit/"; //Reminder: read per-year histograms to get individual contributions from processes, *but* read total-sum histogram to get full postfit error
                        inputfilename = path_dir_shapes_binError + "shapes_postfit_datacard_countExp_"+cat+".root";
                        if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                        file_input = TFile::Open(inputfilename, "READ");
                        if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains("TotalProcs") ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<"TotalProcs' not found ! Skip...")<<endl; continue;}
                        float bincontent = ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinContent(1); //Now read total yield (for error vector)
                        float binerror = ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinError(1);
                        file_input->Close();
                        //--------------------------------------------

                        //-- Errors
                        if(v_y[ibin-1] < 0) //Need to fill total error only once (not for each process)
                        {
                            float bin_width = (xmax_tmp-xmin_tmp) / nIndivBins;
                            v_x[ibin-1] = xmin_tmp + ((ibin-1)*bin_width) + (bin_width/2.);
                            v_y[ibin-1] = bincontent;
                            v_eyl[ibin-1] = binerror;
                            v_eyh[ibin-1] = binerror;
                            v_exl[ibin-1] = bin_width / 2; v_exh[ibin-1] = bin_width / 2;
                            // cout<<"bin "<<ibin<<" / "<<v_y[ibin-1]<<endl;
                        }
                    } //nbins
                }

				//-- Set histo style (use color vector filled in main) //NB: run for all sub-histos... for simplicity
                //---------------------------------------------------
				h_tmp->SetFillStyle(1001);
				if(samplename == "Fakes") {h_tmp->SetFillStyle(3005);}
		        else if(samplename == "QFlip" ) {h_tmp->SetFillStyle(3006);}

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
                //---------------------------------------------------

                //-- Fill vector of MC histograms
                if(v_vector_MC_histo[ivar].size() <=  index_MC_sample) {v_vector_MC_histo[ivar].push_back((TH1F*) h_tmp->Clone());}
                else if(!v_vector_MC_histo[ivar][index_MC_sample] && h_tmp) {v_vector_MC_histo[ivar][index_MC_sample] = (TH1F*) h_tmp->Clone();}
                else {v_vector_MC_histo[ivar][index_MC_sample]->Add((TH1F*) h_tmp->Clone());}
                if(v_vector_MC_histo[ivar][index_MC_sample]) {v_vector_MC_histo[ivar][index_MC_sample]->SetDirectory(0);} //Dis-associate histo from TFile //https://root.cern.ch/root/htmldoc/guides/users-guide/ObjectOwnership.html

				delete h_tmp; h_tmp = NULL; //No crash ? (else only delete if new)
			} //end sample loop


// #####    ##   #####   ##
// #    #  #  #    #    #  #
// #    # #    #   #   #    #
// #    # ######   #   ######
// #    # #    #   #   #    #
// #####  #    #   #   #    #

            TString dataname = "data_obs";
            if(total_var_list[ivar].Contains("mTW"))
            {
                xmin_tmp = 0; xmax_tmp = 150; nIndivBins = 15;
                hdata_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
                for(int ibin=1; ibin<nIndivBins+1; ibin++)
                {
                    inputfilename = path_dir_shapes_binContent + "shapes_postfit_datacard_bin"+Convert_Number_To_TString(ibin)+"_mTW_SRother.root";
                    if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                    file_input = TFile::Open(inputfilename, "READ");
                    // cout<<"inputfilename "<<inputfilename<<endl;

                    TString dir_hist = "bin" + Convert_Number_To_TString(ibin) + "_" + total_var_list[ivar] + "_SRother_" + v_lumiYears[iyear] + "_postfit/";;
                    // cout<<"dir_hist/dataname "<<dir_hist<<dataname<<endl;
                    if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(dataname) ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<dataname<<"' not found ! Skip...")<<endl; continue;}

                    if(isnan(((TH1F*) file_input->Get(dir_hist+dataname))->GetBinContent(1))) {cout<<FRED("ERROR: NaN data (dir_hist="<<dir_hist<<") !")<<endl; continue;} //Can happen when input data is 0 (?)

                    hdata_tmp->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+dataname))->GetBinContent(1)); //Get content/error from individual bin
                    hdata_tmp->SetBinError(ibin, ((TH1F*) file_input->Get(dir_hist+dataname))->GetBinError(1));

                    file_input->Close();
                } //nbins
            }
            else if(total_var_list[ivar].Contains("countExp"))
            {
                xmin_tmp = 0; xmax_tmp = 3; nIndivBins = 3;
                hdata_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
                for(int ibin=1; ibin<nIndivBins+1; ibin++)
                {
                    TString cat = "";
                    if(ibin==1) {cat = "CRWZ";}
                    else if(ibin==2) {cat = "CRZZ";}
                    else if(ibin==3) {cat = "SRttZ4l";}
                    inputfilename = path_dir_shapes_binContent + "shapes_postfit_datacard_countExp_"+cat+".root";
                    if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                    file_input = TFile::Open(inputfilename, "READ");
                    // cout<<"inputfilename "<<inputfilename<<endl;

                    TString dir_hist = total_var_list[ivar] + "_" + cat + "_" + v_lumiYears[iyear] + "_postfit/";
                    // cout<<"dir_hist/dataname "<<dir_hist<<dataname<<endl;
                    if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(dataname) ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<dataname<<"' not found ! Skip...")<<endl; continue;}

                    hdata_tmp->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+dataname))->GetBinContent(1)); //Get content/error from individual bin
                    hdata_tmp->SetBinError(ibin, ((TH1F*) file_input->Get(dir_hist+dataname))->GetBinError(1));

                    file_input->Close();
                } //nbins
            }

            if(v_hdata[ivar] == NULL) {v_hdata[ivar] = (TH1F*) hdata_tmp->Clone();}
            else {v_hdata[ivar]->Add((TH1F*) hdata_tmp->Clone());}
            v_hdata[ivar]->SetMarkerStyle(20);

            v_hdata[ivar]->SetDirectory(0); //Dis-associate from TFile
            delete hdata_tmp; hdata_tmp = NULL; //No crash ? (else only delete if new)
        } //years loop

        //-- Protection against nan/inf data points
        for(int ibin=1; ibin<v_hdata[ivar]->GetNbinsX()+1; ibin++)
        {
            // cout<<"histo_ratio_data["<<ibin<<"] = "<<histo_ratio_data->GetBinContent(ibin)<<endl;
            if(std::isnan(v_hdata[ivar]->GetBinContent(ibin)) || std::isinf(v_hdata[ivar]->GetBinContent(ibin))) {cout<<FRED("ERROR: v_hdata["<<ivar<<"]->GetBinContent("<<ibin<<")="<<v_hdata[ivar]->GetBinContent(ibin)<<" ! May cause plotting bugs !")<<endl; v_hdata[ivar]->SetBinContent(ibin, 0.);}
        }


// ##### #    #  ####  #####   ##    ####  #    #
//   #   #    # #        #    #  #  #    # #   #
//   #   ######  ####    #   #    # #      ####
//   #   #    #      #   #   ###### #      #  #
//   #   #    # #    #   #   #    # #    # #   #
//   #   #    #  ####    #   #    #  ####  #    #

    	//-- Add legend entries -- iterate backwards, so that last histo stacked is on top of legend
        v_stack[ivar] = new THStack;

		for(int i=v_vector_MC_histo[ivar].size()-1; i>=0; i--)
		{
			if(!v_vector_MC_histo[ivar][i]) {continue;} //Some templates may be null
			v_stack[ivar]->Add(v_vector_MC_histo[ivar][i]);

            if(v_histo_total_MC[ivar] == NULL) {v_histo_total_MC[ivar] = (TH1F*) v_vector_MC_histo[ivar][i]->Clone();}
            else {v_histo_total_MC[ivar]->Add(v_vector_MC_histo[ivar][i]);}

			// cout<<"Stacking sample "<<MC_samples_legend[i]<<" / integral "<<v_vector_MC_histo[ivar][i]->Integral()<<endl;
            // cout<<"stack bin 1 content = "<<((TH1*) v_stack[ivar]->GetStack()->Last())->GetBinContent(1)<<endl;
		}


// ###### #####  #####   ####  #####   ####      ####  #####   ##    ####  #    #
// #      #    # #    # #    # #    # #         #        #    #  #  #    # #   #
// #####  #    # #    # #    # #    #  ####      ####    #   #    # #      ####
// #      #####  #####  #    # #####       #         #   #   ###### #      #  #
// #      #   #  #   #  #    # #   #  #    #    #    #   #   #    # #    # #   #
// ###### #    # #    #  ####  #    #  ####      ####    #   #    #  ####  #    #

        //-- Use pointers to vectors : need to give the adress of first element (all other elements can then be accessed iteratively)
        double* eyl = &v_eyl[0];
        double* eyh = &v_eyh[0];
        double* exl = &v_exl[0];
        double* exh = &v_exh[0];
        double* xx = &v_x[0];
        double* yy = &v_y[0];
        // cout<<"v_eyl[0] "<<v_eyl[0]<<endl;
        // cout<<"v_eyh[0] "<<v_eyh[0]<<endl;
        // cout<<"v_exl[0] "<<v_exl[0]<<endl;
        // cout<<"v_exh[0] "<<v_exh[0]<<endl;
        // cout<<"v_x[0] "<<v_x[0]<<endl;
        // cout<<"v_y[0] "<<v_y[0]<<endl;

        v_gr_error[ivar] = new TGraphAsymmErrors(nIndivBins,xx,yy,exl,exh,eyl,eyh);
        v_gr_error[ivar]->SetFillStyle(3254); //3002 //3004
        v_gr_error[ivar]->SetFillColor(kBlack);


// #####  #####    ##   #    #
// #    # #    #  #  #  #    #
// #    # #    # #    # #    #
// #    # #####  ###### # ## #
// #    # #   #  #    # ##  ##
// #####  #    # #    # #    #

        c1->cd(ivar+1);

        //Draw stack
        v_stack[ivar]->Draw("hist");

        v_hdata[ivar]->Draw("e0p same");

        v_gr_error[ivar]->Draw("e2 same"); //Superimposes the uncertainties on stack


// #   # #    #   ##   #    #
//  # #  ##  ##  #  #   #  #
//   #   # ## # #    #   ##
//   #   #    # ######   ##
//   #   #    # #    #  #  #
//   #   #    # #    # #    #

        //-- Set minimum
        v_stack[ivar]->SetMinimum(0.0001); //Remove '0' label
        if(total_var_list[ivar].Contains("countExp")) {v_stack[ivar]->SetMinimum(10.);}

        //-- Set Yaxis maximum
        double ymax = 0;
        ymax = v_hdata[ivar]->GetMaximum(); //Data ymax
        if(ymax < v_stack[ivar]->GetMaximum()) {ymax = v_stack[ivar]->GetMaximum();} //MC ymax
        ymax*= 1.4;
        v_stack[ivar]->SetMaximum(ymax);
		if(total_var_list[ivar].Contains("countExp")) //Can't use 0
		{

			v_stack[ivar]->SetMaximum(v_stack[ivar]->GetMaximum()*20); //Must use higher threshold in log
		}
        // else if(total_var_list[ivar].Contains("mTW")) {v_stack[ivar]->SetMinimum(2.5);}
        c1->Modified();


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
        v_tpad_ratio[ivar] = new TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 1.0);
        v_tpad_ratio[ivar]->SetTopMargin(0.70);
        if(total_var_list[ivar].Contains("mTW")) {v_tpad_ratio[ivar]->SetRightMargin(0.04);}
        v_tpad_ratio[ivar]->SetFillColor(0);
    	v_tpad_ratio[ivar]->SetFillStyle(0);
    	v_tpad_ratio[ivar]->SetGridy(1);
    	v_tpad_ratio[ivar]->Draw();
    	v_tpad_ratio[ivar]->cd(0);

        v_histo_ratio_data[ivar] = (TH1F*) v_hdata[ivar]->Clone();

    	if(!show_pulls_ratio) //Compute ratios (with error bars)
    	{
    		//To get correct error bars in ratio plot, must only account for errors from data, not MC ! (MC error shown as separate band)
    		for(int ibin=1; ibin<v_histo_total_MC[ivar]->GetNbinsX()+1; ibin++)
    		{
    			v_histo_total_MC[ivar]->SetBinError(ibin, 0.);
    		}

    		v_histo_ratio_data[ivar]->Divide(v_histo_total_MC[ivar]);
    	}
     	else //-- Compute pulls (no error bars)
    	{
    		for(int ibin=1; ibin<v_histo_ratio_data[ivar]->GetNbinsX()+1; ibin++)
    		{
    			//Add error on signal strength (since we rescale signal manually)
    			// double bin_error_mu = v_vector_MC_histo[ivar].at(index_tZq_sample)->GetBinError(ibin) * sig_strength_err;
    			// cout<<"bin_error_mu = "<<bin_error_mu<<endl;

    			double bin_error_mu = 0; //No sig strength uncert. for prefit ! //-- postfit -> ?

    			//Quadratic sum of systs, stat error, and sig strength error
    			double bin_error = pow(pow(v_histo_total_MC[ivar]->GetBinError(ibin), 2) + pow(v_histo_ratio_data[ivar]->GetBinError(ibin), 2) + pow(bin_error_mu, 2), 0.5);

    			if(!v_histo_total_MC[ivar]->GetBinError(ibin)) {v_histo_ratio_data[ivar]->SetBinContent(ibin,-99);} //Don't draw null markers
    			else{v_histo_ratio_data[ivar]->SetBinContent(ibin, (v_histo_ratio_data[ivar]->GetBinContent(ibin) - v_histo_total_MC[ivar]->GetBinContent(ibin)) / bin_error );}
    		}
    	}

    	if(show_pulls_ratio) {v_histo_ratio_data[ivar]->GetYaxis()->SetTitle("Pulls");}
    	// else {v_histo_ratio_data[ivar]->GetYaxis()->SetTitle("Data/MC");}
        else {v_histo_ratio_data[ivar]->GetYaxis()->SetTitle("#frac{Data}{Pred.}");}
        // v_histo_ratio_data[ivar]->GetYaxis()->SetTickLength(0.);
        v_histo_ratio_data[ivar]->GetYaxis()->SetTitleOffset(1.15);
        v_histo_ratio_data[ivar]->GetXaxis()->SetTitleOffset(1.05);
        if(total_var_list[ivar].Contains("countExp")) {v_histo_ratio_data[ivar]->GetXaxis()->SetTitleOffset(1.0);}
        // v_histo_ratio_data[ivar]->GetXaxis()->SetTitleOffset(1.1);
        v_histo_ratio_data[ivar]->GetYaxis()->SetLabelSize(0.045);
        v_histo_ratio_data[ivar]->GetXaxis()->SetLabelSize(0.045);
    	v_histo_ratio_data[ivar]->GetXaxis()->SetLabelFont(42);
    	v_histo_ratio_data[ivar]->GetYaxis()->SetLabelFont(42);
    	v_histo_ratio_data[ivar]->GetXaxis()->SetTitleFont(42);
    	v_histo_ratio_data[ivar]->GetYaxis()->SetTitleFont(42);
        v_histo_ratio_data[ivar]->GetYaxis()->SetNdivisions(503); //grid draw on primary tick marks only
    	v_histo_ratio_data[ivar]->GetXaxis()->SetNdivisions(505);
        v_histo_ratio_data[ivar]->GetYaxis()->SetTitleSize(0.06);
    	v_histo_ratio_data[ivar]->GetXaxis()->SetTickLength(0.04);
    	v_histo_ratio_data[ivar]->SetMarkerStyle(20);
    	v_histo_ratio_data[ivar]->SetMarkerSize(1.2); //changed from 1.4
        v_histo_ratio_data[ivar]->GetXaxis()->SetTitleSize(0.06);

        v_histo_ratio_data[ivar]->GetYaxis()->SetTickLength(0.15);
        v_histo_ratio_data[ivar]->GetYaxis()->SetNdivisions(303);

        if(total_var_list[ivar].Contains("countExp")) {v_histo_ratio_data[ivar]->GetXaxis()->SetLabelSize(0.);}

        //-- If a point is outside the y-range of the ratio pad defined by SetMaximum/SetMinimum(), it disappears with its error
        //-- Trick: fill 2 histos with points either above/below y-range, to plot some markers indicating missing points (cleaner)
        //NB: only for ratio plot, not pulls
        float ratiopadmin = 0.4, ratiopadmax = 1.6; //Define ymin/ymax for ratio plot

        //-- Zoom in ratio plot for CR plots (better agreement)
        if(total_var_list[ivar].Contains("countExp")) {ratiopadmin = 0.8; ratiopadmax = 1.2;}

        TH1F* h_pointsAboveY = (TH1F*) v_histo_ratio_data[ivar]->Clone();
        h_pointsAboveY->SetMarkerStyle(26); //Open triangle pointing up
        h_pointsAboveY->SetMarkerSize(1.5);
        TH1F* h_pointsBelowY = (TH1F*) v_histo_ratio_data[ivar]->Clone();
        h_pointsBelowY->SetMarkerStyle(32); //Open triangle pointing down
        h_pointsBelowY->SetMarkerSize(1.5);
        if(show_pulls_ratio)
    	{
    		v_histo_ratio_data[ivar]->SetMinimum(-2.99);
    		v_histo_ratio_data[ivar]->SetMaximum(2.99);
    	}
    	else
    	{
            //-- Default
            v_histo_ratio_data[ivar]->SetMinimum(ratiopadmin); //NB: removes error bars if data point is below ymin...?
            v_histo_ratio_data[ivar]->SetMaximum(ratiopadmax);

            //-- Fill histos with points outside yrange
            for(int ibin=1; ibin<v_histo_ratio_data[ivar]->GetNbinsX()+1; ibin++)
            {
                //-- Default: make point invisible
                h_pointsAboveY->SetBinContent(ibin, -999);
                h_pointsBelowY->SetBinContent(ibin, -999);

                if(v_histo_ratio_data[ivar]->GetBinContent(ibin) > ratiopadmax && v_hdata[ivar]->GetBinContent(ibin) >= 1)
                {
                    //Adjust error
                    float initial_y = v_histo_ratio_data[ivar]->GetBinContent(ibin);
                    float initial_err = v_histo_ratio_data[ivar]->GetBinError(ibin);
                    float new_err = initial_err - (initial_y-ratiopadmax);
                    if(new_err<0) {new_err=0.;}

                    h_pointsAboveY->SetBinContent(ibin, ratiopadmax-0.05); //Add some padding
                    h_pointsAboveY->SetBinError(ibin, new_err);
                }
                else if(v_histo_ratio_data[ivar]->GetBinContent(ibin) < ratiopadmin && v_hdata[ivar]->GetBinContent(ibin) >= 1)
                {
                    //Adjust error
                    float initial_y = v_histo_ratio_data[ivar]->GetBinContent(ibin);
                    float initial_err = v_histo_ratio_data[ivar]->GetBinError(ibin);
                    float new_err = initial_err - (ratiopadmin-initial_y);
                    if(new_err<0) {new_err=0.;}

                    h_pointsBelowY->SetBinContent(ibin, ratiopadmin+(ratiopadmin/10.)); //Add some padding
                    h_pointsBelowY->SetBinError(ibin, new_err);
                }
            }
    	}

        //-- SET X_AXIS TITLES
        if(total_var_list[ivar].Contains("mTW")) {v_histo_ratio_data[ivar]->GetXaxis()->SetTitle("m_{T}^{W}");}

        if(total_var_list[ivar].Contains("countExp")) //Vertical text X labels
		{
            v_histo_ratio_data[ivar]->GetXaxis()->SetTitle("Counting experiments");
            v_histo_ratio_data[ivar]->GetXaxis()->SetLabelSize(0.07); //Increase x-label size
            v_histo_ratio_data[ivar]->GetXaxis()->SetLabelOffset(0.025); //Add some x-label offset
            const char *labels[3]  = {"CR WZ", "CR ZZ", "\\text{SR-t}\\bar{\\text{t}}\\text{Z-4}\\ell"};
            // const char *labels[3]  = {"CR WZ", "CR ZZ", ""};
            // const char *labels[10]  = {Get_Region_Label("wz", ""), Get_Region_Label("zz", ""), Get_Region_Label("ttz4l", "")};
            for(int i=1;i<=3;i++) {v_histo_ratio_data[ivar]->GetXaxis()->SetBinLabel(i,labels[i-1]);}
            // v_histo_ratio_data[ivar]->LabelsOption("v", "X"); //Make X labels vertical
		}

    	if(show_pulls_ratio) {v_histo_ratio_data[ivar]->Draw("HIST P");} //Draw ratio points
        else
        {
            v_histo_ratio_data[ivar]->Draw("E1 X0 P"); //Draw ratio points ; E1 : perpendicular lines at end ; X0 : suppress x errors

            h_pointsAboveY->Draw("E1 X0 P same");
            h_pointsBelowY->Draw("E1 X0 P same");
        }


// ###### #####  #####   ####  #####   ####     #####    ##   ##### #  ####
// #      #    # #    # #    # #    # #         #    #  #  #    #   # #    #
// #####  #    # #    # #    # #    #  ####     #    # #    #   #   # #    #
// #      #####  #####  #    # #####       #    #####  ######   #   # #    #
// #      #   #  #   #  #    # #   #  #    #    #   #  #    #   #   # #    #
// ###### #    # #    #  ####  #    #  ####     #    # #    #   #   #  ####


		//Copy previous TGraphAsymmErrors, then modify it -> error TGraph for ratio plot
		TGraphAsymmErrors *thegraph_tmp = NULL;
		double *theErrorX_h;
		double *theErrorY_h;
		double *theErrorX_l;
		double *theErrorY_l;
		double *theY;
		double *theX;

		thegraph_tmp = (TGraphAsymmErrors*) v_gr_error[ivar]->Clone();
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

		v_gr_ratio_error[ivar] = new TGraphAsymmErrors(thegraph_tmp->GetN(), theX , theY ,  theErrorX_l, theErrorX_h, theErrorY_l, theErrorY_h);
        v_gr_ratio_error[ivar]->SetFillStyle(3254); //3002 //3004
        v_gr_ratio_error[ivar]->SetFillColor(kBlack); //kBlue+2 //kCyan

		if(!show_pulls_ratio) {v_gr_ratio_error[ivar]->Draw("e2 same");} //Draw error bands in ratio plot


//  ####   ####   ####  #    # ###### ##### #  ####   ####
// #    # #    # #      ##  ## #        #   # #    # #
// #      #    #  ####  # ## # #####    #   # #       ####
// #      #    #      # #    # #        #   # #           #
// #    # #    # #    # #    # #        #   # #    # #    #
//  ####   ####   ####  #    # ######   #   #  ####   ####

    	//-- Draw ratio y-lines manually
    	v_hlines1[ivar] = new TH1F("","",this->nbins, v_hdata[ivar]->GetXaxis()->GetXmin(), v_hdata[ivar]->GetXaxis()->GetXmax());
    	v_hlines2[ivar] = new TH1F("","",this->nbins, v_hdata[ivar]->GetXaxis()->GetXmin(), v_hdata[ivar]->GetXaxis()->GetXmax());
    	for(int ibin=1; ibin<this->nbins +1; ibin++)
    	{
    		if(show_pulls_ratio)
    		{
    			v_hlines1[ivar]->SetBinContent(ibin, -1);
    			v_hlines2[ivar]->SetBinContent(ibin, 1);
    		}
    		else
    		{
                if(total_var_list[ivar].Contains("countExp"))
                {
                    v_hlines1[ivar]->SetBinContent(ibin, 0.90);
                    v_hlines2[ivar]->SetBinContent(ibin, 1.10);
                }
                else
                {
                    v_hlines1[ivar]->SetBinContent(ibin, 0.75);
                    v_hlines2[ivar]->SetBinContent(ibin, 1.25);
                }
    		}
    	}
    	v_hlines1[ivar]->SetLineStyle(6); v_hlines2[ivar]->SetLineStyle(6);
    	// v_hlines1[ivar]->Draw("hist same"); v_hlines2[ivar]->Draw("hist same"); //Removed extra lines

        TString Y_label = "Events"; //Default
        if(total_var_list[ivar].Contains("mTW")) {Y_label = "Events / " + Convert_Number_To_TString( (v_stack[ivar]->GetXaxis()->GetXmax() - v_stack[ivar]->GetXaxis()->GetXmin()) / v_histo_total_MC[ivar]->GetNbinsX(), 2) + " GeV";}

    	if(v_stack[ivar]!= 0)
    	{
    		v_stack[ivar]->GetXaxis()->SetLabelFont(42);
    		v_stack[ivar]->GetYaxis()->SetLabelFont(42);
    		v_stack[ivar]->GetYaxis()->SetTitleFont(42);
    		v_stack[ivar]->GetYaxis()->SetTitleSize(0.06);
            v_stack[ivar]->GetYaxis()->SetTickLength(0.04);
    		v_stack[ivar]->GetXaxis()->SetLabelSize(0.0);
    		v_stack[ivar]->GetYaxis()->SetLabelSize(0.048);
    		v_stack[ivar]->GetXaxis()->SetNdivisions(505);
    		v_stack[ivar]->GetYaxis()->SetNdivisions(506);
            v_stack[ivar]->GetYaxis()->SetTitleOffset(1.15);
            // v_stack[ivar]->GetYaxis()->SetTitleOffset(1.28);
    		v_stack[ivar]->GetYaxis()->SetTitle(Y_label);
            v_stack[ivar]->GetXaxis()->SetTickLength(0.);
    	}

    	//----------------
    	// CAPTIONS //
    	//----------------
    	// -- using https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines
        // -- About fonts: https://root.cern.ch/doc/master/classTAttText.html#T5

    	float l = pad->GetLeftMargin();
    	float t = pad->GetTopMargin();

    	TString cmsText = "CMS";
    	TLatex latex;
    	latex.SetNDC();
    	latex.SetTextColor(kBlack);
        latex.SetTextFont(62); //Changed
    	latex.SetTextAlign(11);
    	latex.SetTextSize(0.06);
        // latex.DrawLatex(l + 0.04, 0.87, cmsText);
        if(use_paperStyle) {latex.DrawLatex(l + 0.04, 0.87, cmsText);} //CMS guideline: within frame
        else {latex.DrawLatex(l + 0.01, 0.94, cmsText);} //Default: outside frame

    	float lumi = lumiValue;
    	TString lumi_ts = Convert_Number_To_TString(lumi);
    	lumi_ts += " fb^{-1} (13 TeV)";
    	latex.SetTextFont(42);
    	latex.SetTextAlign(31);
    	latex.SetTextSize(0.04);
        latex.DrawLatex(0.96, 0.94,lumi_ts);

        TString extraText = "Preliminary";
        if(!use_paperStyle) //Default is without (for paper)
        {
            latex.SetTextFont(52);
            latex.SetTextSize(0.05);
            latex.DrawLatex(l + 0.32, 0.94, extraText);
        }

        if(total_var_list[ivar].Contains("countExp"))
        {
            float y1 = 0.8; //Based on yvalues of ratio pad //1.2
            float y2 = 2.3;
            tline1 = new TLine(1.,y1,1.,y2);
            tline1->SetLineColor(kBlack);
            tline1->SetLineWidth(2);
            // tline1->SetLineStyle(2);
            tline1->Draw("same");

            tline2 = new TLine(2.,y1,2.,y2);
            tline2->SetLineColor(kBlack);
            tline2->SetLineWidth(2);
            // tline2->SetLineStyle(2);
            tline2->Draw("same");
        }
        else
        {
            TString info_data = Get_Region_Label(region, "mTW_SRother");
            TLatex text2;
            text2.SetNDC();
            text2.SetTextAlign(13);
            text2.SetTextSize(0.04);
            text2.SetTextFont(42);
            float left = l+0.04;
            if(use_paperStyle) {text2.DrawLatex(left,0.85,info_data);} //0.83
            else {text2.DrawLatex(left,0.90,info_data);} //Default: outside frame //0.90
        }

        pad = NULL;
    } //Var loop
//--------------------------------------------


// ####### #
//    #    #       ######  ####  ###### #    # #####
//    #    #       #      #    # #      ##   # #    #
//    #    #       #####  #      #####  # #  # #    #
//    #    #       #      #  ### #      #  # # #    #
//    #    #       #      #    # #      #   ## #    #
//    #    ####### ######  ####  ###### #    # #####

    int ivar = 0; //Read first element by default

    int n_columns = ceil(nSampleGroups/2.) > 6 ? 6 : ceil(nSampleGroups/2.); //ceil = upper int
    float x_left = 0.94-n_columns*0.12; //Each column allocated same x-space //0.12 needed for most crowded plots
    if(x_left < 0.4) {x_left = 0.4;} //Leave some space for region label

    TLegend* qw = new TLegend(x_left-0.05,0.80,0.94,0.92); //Default /
    qw->SetTextSize(0.04);
    qw->SetNColumns(n_columns);
    qw->SetBorderSize(0);
    qw->SetFillStyle(0); //transparent
    qw->SetTextAlign(12); //align = 10*HorizontalAlign + VerticalAlign //Horiz: 1=left adjusted, 2=centered, 3=right adjusted //Vert: 1=bottom adjusted, 2=centered, 3=top adjusted
    // cout<<"x_left "<<x_left<<endl;
    // cout<<"ceil(nSampleGroups/2.) "<<ceil(nSampleGroups/2.)<<endl;

    //-- Dummy object, only used to display uncertainty band also in legend
    TH1F* h_uncert = new TH1F("h_uncert", "h_uncert", 1, 0, 1);
    h_uncert->SetFillStyle(3254); //3002 //3004
    h_uncert->SetFillColor(kBlack);
    h_uncert->SetLineWidth(0.);
    qw->AddEntry(h_uncert, "Unc.", "F");

	//-- Data on top of legend
    qw->AddEntry(v_hdata[0], "Data" , "ep");

    TLegend* qw_reduced = (TLegend*) qw->Clone();

	for(int i=0; i<v_vector_MC_histo[ivar].size(); i++)
	{
        if(MC_samples_legend[i].Contains("tZq")) {qw->AddEntry(v_vector_MC_histo[ivar][i], "tZq", "f");}
        else if(MC_samples_legend[i].EndsWith("ttZ") ) {qw->AddEntry(v_vector_MC_histo[ivar][i], "t#bar{t}Z", "f");}
        else if(MC_samples_legend[i].EndsWith("tWZ") ) {qw->AddEntry(v_vector_MC_histo[ivar][i], "tWZ", "f");}
        else if(MC_samples_legend[i] == "ttW" || MC_samples_legend[i] == "tX") {qw->AddEntry(v_vector_MC_histo[ivar][i], "t(#bar{t})X", "f");}
        else if(MC_samples_legend[i] == "WZ") {qw->AddEntry(v_vector_MC_histo[ivar][i], "WZ", "f");}
        else if(MC_samples_legend[i] == "WWZ" || MC_samples_legend[i] == "VVV") {qw->AddEntry(v_vector_MC_histo[ivar][i], "VV(V)", "f");}
        else if(MC_samples_legend[i] == "TTGamma_Dilep" || MC_samples_legend[i] == "XG") {qw->AddEntry(v_vector_MC_histo[ivar][i], "X#gamma", "f");}
        else if(MC_samples_legend[i] == "TTbar_DiLep" || MC_samples_legend[i] == "NPL" || MC_samples_legend[i] == "NPL_DATA") {qw->AddEntry(v_vector_MC_histo[ivar][i], "NPL", "f");}
	}

    for(int i=0; i<v_vector_MC_histo[ivar].size(); i++)
	{
        if(MC_samples_legend[i] == "WZ") {qw_reduced->AddEntry(v_vector_MC_histo[ivar][i], "WZ", "f");}
        else if(MC_samples_legend[i] == "VVV") {qw_reduced->AddEntry(v_vector_MC_histo[ivar][i], "VV(V)", "f");}
        else if(MC_samples_legend[i] == "tX") {qw_reduced->AddEntry(v_vector_MC_histo[ivar][i], "t(#bar{t})X", "f");}
        else if(MC_samples_legend[i] == "TTGamma_Dilep" || MC_samples_legend[i] == "XG") {qw_reduced->AddEntry(v_vector_MC_histo[ivar][i], "X#gamma", "f");}
        else if(MC_samples_legend[i].Contains("NPL")) {qw_reduced->AddEntry(v_vector_MC_histo[ivar][i], "NPL", "f");}
        else if(MC_samples_legend[i].EndsWith("ttZ") ) {qw_reduced->AddEntry(v_vector_MC_histo[ivar][i], "t#bar{t}Z", "f");}
        else if(MC_samples_legend[i].EndsWith("tWZ") ) {qw_reduced->AddEntry(v_vector_MC_histo[ivar][i], "tWZ", "f");}
	}

    for(int ivar=0; ivar<nvar; ivar++)
    {
        c1->cd(ivar+1);
        if(ivar==0) {qw_reduced->Draw("same");}
        else if(ivar==1) {qw->Draw("same");} //Draw legend
    }

    // TLatex text2; //Different offset for SRttz4l label ?
    // text2.DrawTextNDC(-0.4,0.03, "\\text{SR-t}\\bar{\\text{t}}\\text{Z-4}\\ell");


// #    # #####  # ##### ######     ####  #    # ##### #####  #    # #####
// #    # #    # #   #   #         #    # #    #   #   #    # #    #   #
// #    # #    # #   #   #####     #    # #    #   #   #    # #    #   #
// # ## # #####  #   #   #         #    # #    #   #   #####  #    #   #
// ##  ## #   #  #   #   #         #    # #    #   #   #      #    #   #
// #    # #    # #   #   ######     ####   ####    #   #       ####    #

    TString outdir = "plots/paperPlots/";
    mkdir(outdir.Data(), 0777);
    TString output_plot_name = outdir + "PostfitTemplates_commonRegions";
    if(!use_paperStyle) {output_plot_name+= "_prelim";}

    c1->SaveAs(output_plot_name + ".png");
    c1->SaveAs(output_plot_name + ".eps");
    c1->SaveAs(output_plot_name + ".pdf");

    delete c1; c1 = NULL;
    delete qw; qw = NULL;
    // delete qw_reduced; qw_reduced = NULL;
    if(h_uncert) {delete h_uncert; h_uncert = NULL;}
    if(tline1) {delete tline1; tline1 = NULL;}
    if(tline2) {delete tline2; tline2 = NULL;}

    for(int ivar=0; ivar<nvar; ivar++)
    {
        if(v_gr_error[ivar]) {delete v_gr_error[ivar]; v_gr_error[ivar] = NULL;}
        if(v_stack[ivar]) {delete v_stack[ivar]; v_stack[ivar] = NULL;}
        if(v_gr_ratio_error[ivar]) {delete v_gr_ratio_error[ivar]; v_gr_ratio_error[ivar] = NULL;}
        if(v_hlines1[ivar]) {delete v_hlines1[ivar]; v_hlines1[ivar] = NULL;}
        if(v_hlines2[ivar]) {delete v_hlines2[ivar]; v_hlines2[ivar] = NULL;}
        // delete v_tpad_ratio[ivar]; v_tpad_ratio[ivar] = NULL; //Deleteed with c1
    }

    return;
}











/**
 * Hard-coded function to produce figures for the paper.
 * Plot showing the postfit distributions in the regions common to all fits (CRs, SRttZ4l, SRother).
 * NB: for NN-SM distributions, not much sense to show postfit (which fit ? what EFT ratio to show ? etc.) --> Show prefit instead, without EFT ratio, and all 3 nodes. Hence, get treated separately throughout this function.
 */
void TopEFT_analysis::Make_PaperPlot_SignalRegions(TString template_name)
{
//--------------------------------------------

    cout<<endl<<BYEL("                          ")<<endl<<endl;
	cout<<FYEL("--- Producing paper figure [SIGNAL REGIONS] / "<<template_name<<" ---")<<endl;
    cout<<endl<<BYEL("                          ")<<endl<<endl;


//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

    //--------------------------
    // DEFINE VARIABLES TO PLOT
    //--------------------------

    vector<TString> total_var_list;
    total_var_list.push_back(template_name + "_SRttZ");
    total_var_list.push_back(template_name + "_SRtZq");
    if(template_name.Contains("SM")) {total_var_list.push_back(template_name + "_SRother");} //Special case: for NN-SM, plot all 3 output nodes (prefit)

    int nvar = total_var_list.size();

    //--------------------------
    // DIVIDE CANVAS IN 4 PADS
    //--------------------------

    //-- Canvas definition
    Load_Canvas_Style(); //Default top/bottom/left/right margins: 0.07/0.13/0.16/0.03
    TCanvas* c1 = NULL;
    if(template_name.Contains("SM"))
    {
        c1 = new TCanvas("c1","c1", 1400, 600);
    }
    else
    {
        c1 = new TCanvas("c1","c1", 1400, 950); //Optimized to use full page with 2 figures
        // c1 = new TCanvas("c1","c1", 1400, 800); //Default
    }
    c1->Divide(nvar, 1, 1E-11, 1E-11); //(x,y)

    //--------------------------
    // CREATE VECTORS OF OBJECTS
    //--------------------------

    vector<TH1F*> v_hdata(nvar);
    vector<THStack*> v_stack(nvar);
    vector<TH1F*> v_histo_total_MC(nvar);
    vector<vector<TH1F*> > v_vector_MC_histo(nvar); //Store separately the histos for each MC sample --> stack them after loops
    vector<TH1F*> v_histo_ratio_data(nvar);
    vector<TH1F*> v_hlines1(nvar), v_hlines2(nvar); //Draw TLinesin ratio plots

    vector<TGraphAsymmErrors*> v_gr_error(nvar);
    vector<TGraphAsymmErrors*> v_gr_ratio_error(nvar);

    vector<TPad*> v_tpad_ratio(nvar);

    //-- Ratio EFT/SM
    vector<vector<TH1F*>> v2_histo_ratio_eft(nvar); //1 vector per variable; 1 vector element per scenario //Store EFT points
    for(int i=0; i<v2_histo_ratio_eft.size(); i++) {v2_histo_ratio_eft[i].resize(2);}
    vector<TH1F*> v_histo_ratio_sm(nvar); //1 vector per variable; 1 vector element per scenario //Store SM point

    //-- Can also plot SM+EFT/SM ratio when considering total sum of MC, not only signal
    vector<vector<TH1F*>> v2_histo_ratio_eft_totalProcs_sm(nvar); //1 vector per variable; 1 vector element per scenario //Store EFT points
    vector<vector<TH1F*>> v2_histo_ratio_eft_totalProcs_smeft(nvar); //1 vector per variable; 1 vector element per scenario //Store EFT points
    vector<vector<vector<float>>> v3_EFTsig_binContents(nvar); //for each var and each EFT hypothesis, store 2 float per histo bin
    for(int i=0; i<nvar; i++)
    {
        v2_histo_ratio_eft_totalProcs_sm[i].resize(2);
        v2_histo_ratio_eft_totalProcs_smeft[i].resize(2);
        v3_EFTsig_binContents[i].resize(2);
    }

    vector<TPad*> v_tpad_ratio_eft(nvar);

    vector<TLegend*> v_tlegend_ratio_eft(nvar);
    vector<vector<TH1F*>> v_tlegend_dummy(nvar); //Dummy object to fill legend
    for(int i=0; i<v_tlegend_dummy.size(); i++) {v_tlegend_dummy[i].resize(2);}

    vector<TLatex*> v_tlatex_legend_eft(nvar);

    vector<TString> MC_samples_legend; //List the MC samples to mention in legend

    int nbjets_min =-1, nbjets_max=-1, njets_min=-1, njets_max=-1;
    int nIndivBins; float xmin_tmp, xmax_tmp;

    TString fit_type = "postfit"; //Default
    if(template_name.Contains("SM")) {fit_type = "prefit";} //Do prefit SM plots !

    TString path_dir_shapes_binContent = "./outputs/dir_shapes_tmp/binContent/"; //Path of dir. containing shapes (with correct bin contents <-> did not freeze any parameter)
    TString path_dir_shapes_binError = "./outputs/dir_shapes_tmp/binError/"; //Path of dir. containing shapes (with correct bin errors <-> froze WC and problematic split JEC)

    if(template_name.Contains("SM")) //For prefit plots, get bin content/error from same file in common folder
    {
        path_dir_shapes_binContent = "./outputs/dir_shapes_tmp/";
        path_dir_shapes_binError = path_dir_shapes_binContent;
    }


// #       ####   ####  #####   ####
// #      #    # #    # #    # #
// #      #    # #    # #    #  ####
// #      #    # #    # #####       #
// #      #    # #    # #      #    #
// ######  ####   ####  #       ####

//--------------------------------------------
	for(int ivar=0; ivar<nvar; ivar++)
	{
        cout<<endl<<FBLU("== VARIABLE: "<<total_var_list[ivar]<<"")<<endl;

    	TString primitive_name = "c1_" + Convert_Number_To_TString(ivar+1);
        TPad* pad = (TPad*) c1->GetPrimitive(primitive_name);
        TString inputfilename = "";
        float rightmargin = -1; //Default (use default right margin)
        // if(total_var_list[ivar].Contains("5D") || total_var_list[ivar].Contains("cpq3_SRttZ")) {rightmargin = 0.04;} //More bins -> need more space
        // if(total_var_list[ivar].Contains("SM")) {rightmargin = 0.05;} //Need more space

        // if(total_var_list[ivar].Contains("SRtZq")) //Left
        if(ivar==0) //Left
    	{
    		pad->SetPad(0, 0., 0.50, 1); //xlow, ylow, xup, yup
            if(template_name.Contains("SM")) {pad->SetPad(0, 0., 0.33, 1);}
    	}
    	else if(ivar==1) //Right
    	{
    		pad->SetPad(0.50, 0., 1., 1); //xlow, ylow, xup, yup
            if(template_name.Contains("SM")) {pad->SetPad(0.33, 0., 0.66, 1);}
        }
        else if(ivar==1)
        {
            if(template_name.Contains("SM")) {pad->SetPad(0.66, 0., 1., 1);}
        }
        if(rightmargin>0) {pad->SetRightMargin(rightmargin);}
        float top_panel_size = 0.6;
        if(template_name.Contains("SM")) {pad->SetBottomMargin(0.30);}
        else {pad->SetBottomMargin(1-top_panel_size);}

        bool isLogY = true; //Default
        if(total_var_list[ivar].Contains("SM") || (total_var_list[ivar].Contains("cpq3") && total_var_list[ivar].Contains("ttZ") && !this->use_NN_cpq3_SRttZ)) {isLogY = false;}
        if(isLogY) {pad->SetLogy();}

        TFile* file_input = NULL;
        TH1F* h_tmp = NULL; //Tmp storing histo
        TH1F* hdata_tmp = NULL; //Tmp storing data histo

		//-- Init error vectors
		double x, y, errory_low, errory_high;

		vector<double> v_eyl, v_eyh, v_exl, v_exh, v_x, v_y; //Contain the systematic errors (used to create the TGraphError)

        float bin_width = -1; //Get bin width of histograms for current variable

        string EFTpoint_name1 = "", EFTpoint_name2 = "";

        //-- All histos are for given lumiYears and sub-channels --> Need to sum them all for plots
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {


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

                // cout<<"sample_list[isample] "<<sample_list[isample]<<endl;
                // cout<<"index_MC_sample "<<index_MC_sample<<endl;

				//-- In Combine, some individual contributions are merged as "Rares"/"EWK", etc.
				//-- If using Combine file, change the names of the samples we look for, and look only once for histogram of each "group"
				TString samplename = sample_list[isample];
				if(isample > 0 && sample_groups[isample] == sample_groups[isample-1]) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //if same group as previous sample, skip it
                else if(make_SMvsEFT_templates_plots && (sample_groups[isample] == "tZq" || sample_groups[isample] == "ttZ" || sample_groups[isample] == "tWZ")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //SM vs EFT --> use private signal samples
                else {samplename = sample_groups[isample];}

				//-- Protections, special cases
				if(sample_list[isample] == "DATA") {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}
                else if(!make_SMvsEFT_templates_plots && sample_list[isample].Contains("PrivMC")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //SM configuration --> only stack central samples (not private samples)
                else if(make_SMvsEFT_templates_plots && (sample_list[isample] == "tZq" || sample_list[isample] == "ttZ" || sample_list[isample] == "tWZ")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //EFT configuration --> only stack private samples (at SM point), not central samples
                else if(sample_list[isample] == "NPL_DATA")  {samplename = "NPL";} //Instead of 'NPL_DATA' and 'NPL_MC', we only want to read the merged histo 'NPL'
                else if(sample_list[isample] == "NPL_MC")  {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //NPL_MC gets substracted from NPL histograms and deleted --> Ignore this vector element //Remove ?

                //-- Add sample name to list (used for legend) //NB: add even if histo was not found and skipped, because expect that it will be found for some other year/channel/... But if not found at all, legend will be wrong
                if(iyear==0 && ivar==0 && samplename != "DATA")
                {
                    if(v_vector_MC_histo[ivar].size() <=  index_MC_sample) {MC_samples_legend.push_back(samplename);}
                }
                if(v_isSkippedSample[isample] == true) {continue;} //Skip this sample

				// cout<<endl<<UNDL(FBLU("-- Sample : "<<sample_list[isample]<<" : "))<<endl;

				h_tmp = NULL;
                Get_Template_Range(nIndivBins, xmin_tmp, xmax_tmp, total_var_list[ivar], true, 2, true, nbjets_min, nbjets_max, njets_min, njets_max, minmax_bounds, this->use_NN_SRother, this->use_NN_cpq3_SRttZ);
                h_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
                v_eyl.resize(nIndivBins); v_eyh.resize(nIndivBins); v_exl.resize(nIndivBins); v_exh.resize(nIndivBins); v_y.resize(nIndivBins); v_x.resize(nIndivBins);
                std::fill(v_y.begin(), v_y.end(), -1); //Init errors positions to <0 (invisible)

                for(int ibin=1; ibin<nIndivBins+1; ibin++)
                {
                    int bin_to_read = 1; //Default (for split-per-bin postfit plots, always need to read bin1); but for prefit NN_SM templates, need to read current ibin (because we are reading full templates directly, not 1 bin at a time)
                    if(template_name.Contains("SM")) {bin_to_read=ibin;}

                    //-- Read bin content
                    //--------------------------------------------
                    inputfilename = path_dir_shapes_binContent + "shapes_"+fit_type+"_datacard_bin"+Convert_Number_To_TString(ibin)+"_"+total_var_list[ivar]+".root";
                    if(template_name.Contains("SM")) {inputfilename = path_dir_shapes_binContent + "shapes_"+fit_type+"_COMBINED_Datacard_TemplateFit_"+total_var_list[ivar]+"_Run2.root";} //Special case: for prefit NN-SM plots, consider full templates directly
                    if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                    file_input = TFile::Open(inputfilename, "READ");
                    // cout<<"inputfilename "<<inputfilename<<endl;
                    TString dir_hist = "bin" + Convert_Number_To_TString(ibin) + "_" + total_var_list[ivar] + "_" + v_lumiYears[iyear] + "_"+fit_type+"/";
                    if(template_name.Contains("SM")) {dir_hist = total_var_list[ivar] + "_" + v_lumiYears[iyear] + "_"+fit_type+"/";}
                    if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(samplename) ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<samplename<<"' not found ! Skip...")<<endl; continue;}
                    // cout<<"dir_hist/samplename "<<dir_hist<<samplename<<endl;
                    h_tmp->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+samplename))->GetBinContent(bin_to_read));
                    file_input->Close();
                    //--------------------------------------------

                    //-- Read bin error
                    //--------------------------------------------
                    dir_hist = fit_type+"/"; //Reminder: read per-year histograms to get individual contributions from processes, *but* read total-sum histogram to get full postfit error
                    inputfilename = path_dir_shapes_binError + "shapes_"+fit_type+"_datacard_bin"+Convert_Number_To_TString(ibin)+"_"+total_var_list[ivar]+".root";
                    if(template_name.Contains("SM")) {inputfilename = path_dir_shapes_binError + "shapes_"+fit_type+"_COMBINED_Datacard_TemplateFit_"+total_var_list[ivar]+"_Run2.root";} //Special case: for prefit NN-SM plots, consider full templates directly
                    if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                    file_input = TFile::Open(inputfilename, "READ");
                    // cout<<"inputfilename "<<inputfilename<<endl;
                    if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains("TotalProcs") ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<"TotalProcs' not found ! Skip...")<<endl; continue;}
                    float binerror = ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinError(bin_to_read);
                    float bincontent = ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinContent(bin_to_read); //Now read total yield (for error vector)
                    file_input->Close();
                    //--------------------------------------------

                    //-- Errors
                    if(v_y[ibin-1] < 0) //Need to fill total error only once (not for each process)
                    {
                        float bin_width = (xmax_tmp-xmin_tmp) / nIndivBins;
                        v_x[ibin-1] = xmin_tmp + ((ibin-1)*bin_width) + (bin_width/2.);
                        v_y[ibin-1] = bincontent;
                        v_eyl[ibin-1] = binerror;
                        v_eyh[ibin-1] = binerror;
                        v_exl[ibin-1] = bin_width / 2; v_exh[ibin-1] = bin_width / 2;
                        // cout<<"bin "<<ibin<<" / "<<v_y[ibin-1]<<endl;
                    }
                } //nbins

				//-- Set histo style (use color vector filled in main) //NB: run for all sub-histos... for simplicity
                //---------------------------------------------------
				h_tmp->SetFillStyle(1001);
				if(samplename == "Fakes") {h_tmp->SetFillStyle(3005);}
		        else if(samplename == "QFlip" ) {h_tmp->SetFillStyle(3006);}

				h_tmp->SetFillColor(color_list[isample]);
				h_tmp->SetLineColor(kBlack);
                h_tmp->SetLineColor(color_list[isample]);

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
                //---------------------------------------------------

                //-- Fill vector of MC histograms
                if(v_vector_MC_histo[ivar].size() <=  index_MC_sample) {v_vector_MC_histo[ivar].push_back((TH1F*) h_tmp->Clone());}
                else if(!v_vector_MC_histo[ivar][index_MC_sample] && h_tmp) {v_vector_MC_histo[ivar][index_MC_sample] = (TH1F*) h_tmp->Clone();}
                else {v_vector_MC_histo[ivar][index_MC_sample]->Add((TH1F*) h_tmp->Clone());}
                if(v_vector_MC_histo[ivar][index_MC_sample]) {v_vector_MC_histo[ivar][index_MC_sample]->SetDirectory(0);} //Dis-associate histo from TFile //https://root.cern.ch/root/htmldoc/guides/users-guide/ObjectOwnership.html

				delete h_tmp; h_tmp = NULL; //No crash ? (else only delete if new)
			} //end sample loop


// #####    ##   #####   ##
// #    #  #  #    #    #  #
// #    # #    #   #   #    #
// #    # ######   #   ######
// #    # #    #   #   #    #
// #####  #    #   #   #    #

            TString dataname = "data_obs";
            hdata_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
            for(int ibin=1; ibin<nIndivBins+1; ibin++)
            {
                int bin_to_read = 1; //Default (for split-per-bin postfit plots, always need to read bin1); but for prefit NN_SM templates, need to read current ibin (because we are reading full templates directly, not 1 bin at a time)
                if(template_name.Contains("SM")) {bin_to_read=ibin;}

                //NB: sum per-year histograms, but could as well read directly the total summed histo !
                //-- Read bin content
                //--------------------------------------------
                inputfilename = path_dir_shapes_binContent + "shapes_"+fit_type+"_datacard_bin"+Convert_Number_To_TString(ibin)+"_"+total_var_list[ivar]+".root";
                if(template_name.Contains("SM")) {inputfilename = path_dir_shapes_binContent + "shapes_"+fit_type+"_COMBINED_Datacard_TemplateFit_"+total_var_list[ivar]+"_Run2.root";} //Special case: for prefit NN-SM plots, consider full templates directly
                if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                file_input = TFile::Open(inputfilename, "READ");
                // cout<<"inputfilename "<<inputfilename<<endl;
                TString dir_hist = "bin" + Convert_Number_To_TString(ibin) + "_" + total_var_list[ivar] + "_" + v_lumiYears[iyear] + "_"+fit_type+"/";
                if(template_name.Contains("SM")) {dir_hist = total_var_list[ivar] + "_" + v_lumiYears[iyear] + "_"+fit_type+"/";}
                if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(dataname) ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<dataname<<"' not found ! Skip...")<<endl; continue;}
                // cout<<"dir_hist/dataname "<<dir_hist<<dataname<<endl;
                float bincontent = ((TH1F*) file_input->Get(dir_hist+dataname))->GetBinContent(bin_to_read);
                float binerror = ((TH1F*) file_input->Get(dir_hist+dataname))->GetBinError(bin_to_read);
                file_input->Close();
                //--------------------------------------------

                if(isnan(bincontent)) {cout<<FRED("ERROR: NaN data (dir_hist="<<dir_hist<<") !")<<endl; continue;} //Can happen when input data is 0 (?)
                hdata_tmp->SetBinContent(ibin, bincontent); //Get content/error from individual bin
                hdata_tmp->SetBinError(ibin, binerror);

                // cout<<"dir_hist "<<dir_hist<<endl;
                // cout<<"hdata_tmp->GetBinContent(ibin) "<<hdata_tmp->GetBinContent(ibin)<<endl;
            } //nbins

            if(v_hdata[ivar] == NULL) {v_hdata[ivar] = (TH1F*) hdata_tmp->Clone();}
            else {v_hdata[ivar]->Add((TH1F*) hdata_tmp->Clone());}
            v_hdata[ivar]->SetMarkerStyle(20);

            v_hdata[ivar]->SetDirectory(0); //Dis-associate from TFile
            delete hdata_tmp; hdata_tmp = NULL; //No crash ? (else only delete if new)


// #####  #####  # #    #         ####    ##   #    # #####  #      ######
// #    # #    # # #    #        #       #  #  ##  ## #    # #      #
// #    # #    # # #    #         ####  #    # # ## # #    # #      #####
// #####  #####  # #    # ###         # ###### #    # #####  #      #
// #      #   #  #  #  #  ###    #    # #    # #    # #      #      #
// #      #    # #   ##   ###     ####  #    # #    # #      ###### ######

            if(!template_name.Contains("SM")) //If plotting SM+EFT/SM ratio (for postfit plots) --> Read TH1EFT objects in my own template files
            {
                //-- Protection: if want to plot private SMEFT samples, make sure they are included in the main sample list
                //-- NB: if using Combine file, treat private SMEFT samples like central samples (included in stack); can't rescale & superimpose on plot, because Combine does not store TH1EFT objects !
                bool PrivMC_sample_found = false;
                for(int isample=0; isample<sample_list.size(); isample++) {if(sample_list[isample].Contains("PrivMC")) {PrivMC_sample_found = true;}}

                inputfilename = "./outputs/Templates_"+template_name+"_EFT2_Run2.root";
                if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                file_input = TFile::Open(inputfilename, "READ");
                // cout<<"inputfilename "<<inputfilename<<endl;

                TH1EFT* th1eft_tmp = NULL;
                TH1F* th1f_new = NULL;

                TString samplename = "PrivMC_tZq";
                if(total_var_list[ivar].Contains("ttZ")) {samplename = "PrivMC_ttZ";}
                TString histo_name = total_var_list[ivar] + "_" + v_lumiYears[iyear] + "__" + samplename;
                // cout<<"histo_name "<<histo_name<<endl;

                if(!file_input->GetListOfKeys()->Contains(histo_name) ) {cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found ! Skip..."))<<endl; continue;}
                th1eft_tmp = (TH1EFT*) file_input->Get(histo_name);
                th1f_new = new TH1F("", "", th1eft_tmp->GetNbinsX(), th1eft_tmp->GetXaxis()->GetXmin(), th1eft_tmp->GetXaxis()->GetXmax());
                if(v_histo_ratio_sm[ivar] == NULL) {v_histo_ratio_sm[ivar] = (TH1F*) th1eft_tmp->Clone();}
                else {v_histo_ratio_sm[ivar]->Add((TH1F*) th1eft_tmp->Clone());}

                // cout<<"SM Integral = "<<th1eft_tmp->Integral()<<endl;

                for(int j=0; j<v3_EFTsig_binContents[ivar].size(); j++)
                {
                    if(v3_EFTsig_binContents[ivar][j].size()==0) {v3_EFTsig_binContents[ivar][j].resize(th1eft_tmp->GetNbinsX());}
                }
                for(int ibin=1; ibin<th1eft_tmp->GetNbinsX()+1; ibin++)
                {
                    v3_EFTsig_binContents[ivar][0][ibin-1]-= th1eft_tmp->GetBinContent(ibin);
                    v3_EFTsig_binContents[ivar][1][ibin-1]-= th1eft_tmp->GetBinContent(ibin);
                }

                EFTpoint_name1 = "", EFTpoint_name2 = "";
                if(total_var_list[ivar].Contains("ctz")) {EFTpoint_name1 = "rwgt_ctz_1.5";}
                else if(total_var_list[ivar].Contains("ctw")) {EFTpoint_name1 = "rwgt_ctw_1";}
                else if(total_var_list[ivar].Contains("cpq3"))
                {
                    EFTpoint_name1 = "rwgt_cpq3_1.5";
                }
                else if(total_var_list[ivar].Contains("5D"))
                {
                    EFTpoint_name1 = "rwgt_ctz_0.5_ctw_0.5_cpq3_1";
                    // if(total_var_list[ivar].Contains("SRttZ")) {EFTpoint_name1 = "rwgt_ctz_0.5_ctw_0.5";}
                    // if(total_var_list[ivar].Contains("SRttZ")) {EFTpoint_name1 = "rwgt_cpq3_2";} //Testing
                }
                else if(total_var_list[ivar].Contains("SM")) {EFTpoint_name1 = "rwgt_cpqm_5";}

                //-- Rescale TH1EFT accordingly to current reweight //Pay attention to operator exact names !
                WCPoint wcp = WCPoint(EFTpoint_name1, 1.);
                th1eft_tmp->Scale(wcp);
                // cout<<"EFTpoint_name1 "<<EFTpoint_name1<<" --> Integral = "<<th1eft_tmp->Integral()<<endl;
                for(int ibin=1; ibin<th1f_new->GetNbinsX()+1; ibin++)
                {
                    // cout<<"ibin "<<ibin<<" / content "<<th1eft_tmp->GetBinContent(ibin)<<" / error "<<th1eft_tmp->GetBinError(ibin)<<endl;
                    th1f_new->SetBinContent(ibin, th1eft_tmp->GetBinContent(ibin));
                    th1f_new->SetBinError(ibin, th1eft_tmp->GetBinError(ibin));
                }
                if(v2_histo_ratio_eft[ivar][0] == NULL) {v2_histo_ratio_eft[ivar][0] = (TH1F*) th1f_new->Clone();}
                else {v2_histo_ratio_eft[ivar][0]->Add((TH1F*) th1f_new->Clone());}

                //-- EFT scenario 1
                for(int ibin=1; ibin<th1eft_tmp->GetNbinsX()+1; ibin++)
                {
                    v3_EFTsig_binContents[ivar][0][ibin-1]+= th1eft_tmp->GetBinContent(ibin);
                }

                if(total_var_list[ivar].Contains("ctz")) {EFTpoint_name2 = "rwgt_ctz_3";}
                else if(total_var_list[ivar].Contains("ctw")) {EFTpoint_name2 = "rwgt_ctw_1.5";}
                // else if(total_var_list[ivar].Contains("ctw")) {EFTpoint_name2 = "rwgt_ctw_2";}
                else if(total_var_list[ivar].Contains("cpq3"))
                {
                    EFTpoint_name2 = "rwgt_cpq3_3";
                }
                else if(total_var_list[ivar].Contains("5D"))
                {
                    EFTpoint_name2 = "rwgt_ctz_1_ctw_1_cpq3_3";
                    // if(total_var_list[ivar].Contains("SRttZ")) {EFTpoint_name2 = "rwgt_ctz_1_ctw_1";}
                    // if(total_var_list[ivar].Contains("SRttZ")) {EFTpoint_name2 = "rwgt_cpq3_4";} //Testing
                }
                else if(total_var_list[ivar].Contains("SM")) {EFTpoint_name2 = "rwgt_cpt_5";}
                wcp = WCPoint(EFTpoint_name2, 1.);
                th1eft_tmp->Scale(wcp);
                // cout<<"EFTpoint_name2 "<<EFTpoint_name2<<" --> Integral = "<<th1eft_tmp->Integral()<<endl;
                for(int ibin=1; ibin<th1f_new->GetNbinsX()+1; ibin++)
                {
                    // cout<<"ibin "<<ibin<<" / content "<<th1eft_tmp->GetBinContent(ibin)<<" / error "<<th1eft_tmp->GetBinError(ibin)<<endl;
                    th1f_new->SetBinContent(ibin, th1eft_tmp->GetBinContent(ibin));
                    th1f_new->SetBinError(ibin, th1eft_tmp->GetBinError(ibin));
                }
                if(v2_histo_ratio_eft[ivar][1] == NULL) {v2_histo_ratio_eft[ivar][1] = (TH1F*) th1f_new->Clone();}
                else {v2_histo_ratio_eft[ivar][1]->Add((TH1F*) th1f_new->Clone());}

                //-- EFT scenario 2
                for(int ibin=1; ibin<th1eft_tmp->GetNbinsX()+1; ibin++)
                {
                    v3_EFTsig_binContents[ivar][1][ibin-1]+= th1eft_tmp->GetBinContent(ibin);
                }

                if(th1f_new) {delete th1f_new; th1f_new = NULL;}
            }

        } //years loop

        //-- Protection against nan/inf data points
        for(int ibin=1; ibin<v_hdata[ivar]->GetNbinsX()+1; ibin++)
        {
            // cout<<"histo_ratio_data["<<ibin<<"] = "<<histo_ratio_data->GetBinContent(ibin)<<endl;
            if(std::isnan(v_hdata[ivar]->GetBinContent(ibin)) || std::isinf(v_hdata[ivar]->GetBinContent(ibin))) {cout<<FRED("ERROR: v_hdata["<<ivar<<"]->GetBinContent("<<ibin<<")="<<v_hdata[ivar]->GetBinContent(ibin)<<" ! May cause plotting bugs !")<<endl; v_hdata[ivar]->SetBinContent(ibin, 0.);}
        }

        if(!template_name.Contains("SM"))
        {
            //-- Sum MC - (signal_SM) + (signal_SMEFT)
            v2_histo_ratio_eft_totalProcs_sm[ivar][0] = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
            v2_histo_ratio_eft_totalProcs_sm[ivar][1] = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
            v2_histo_ratio_eft_totalProcs_smeft[ivar][0] = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
            v2_histo_ratio_eft_totalProcs_smeft[ivar][1] = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
            for(int ibin=1; ibin<v2_histo_ratio_eft[ivar][0]->GetNbinsX()+1; ibin++)
            {
                inputfilename = path_dir_shapes_binContent + "shapes_"+fit_type+"_datacard_bin"+Convert_Number_To_TString(ibin)+"_"+total_var_list[ivar]+".root";
                if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                file_input = TFile::Open(inputfilename, "READ");
                TString dir_hist = fit_type+"/"; //Reminder: read per-year histograms to get individual contributions from processes, *but* read total-sum histogram to get full postfit error
                if(!file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains("TotalProcs") ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<"TotalProcs' not found ! Skip...")<<endl; continue;}
                v2_histo_ratio_eft_totalProcs_sm[ivar][0]->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinContent(1)); //Get content/error from individual bin
                v2_histo_ratio_eft_totalProcs_sm[ivar][1]->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinContent(1)); //Get content/error from individual bin
                v2_histo_ratio_eft_totalProcs_smeft[ivar][0]->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinContent(1)+v3_EFTsig_binContents[ivar][0][ibin-1]); //Get content/error from individual bin
                v2_histo_ratio_eft_totalProcs_smeft[ivar][1]->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinContent(1)+v3_EFTsig_binContents[ivar][1][ibin-1]); //Get content/error from individual bin
                file_input->Close();
            }

            // for(int ibin=1; ibin<v2_histo_ratio_eft[ivar][0]->GetNbinsX()+1; ibin++)
            // {
            //     cout<<"ibin "<<ibin<<" / v2_histo_ratio_eft_totalProcs_smeft[ivar][0]->GetBinContent(ibin) "<<v2_histo_ratio_eft_totalProcs_smeft[ivar][0]->GetBinContent(ibin)<<endl;
            //     cout<<"ibin "<<ibin<<" / v2_histo_ratio_eft_totalProcs_sm[ivar][0]->GetBinContent(ibin) "<<v2_histo_ratio_eft_totalProcs_sm[ivar][0]->GetBinContent(ibin)<<endl;
            // }
            v2_histo_ratio_eft_totalProcs_smeft[ivar][0]->Divide(v2_histo_ratio_eft_totalProcs_sm[ivar][0]);
            v2_histo_ratio_eft_totalProcs_smeft[ivar][1]->Divide(v2_histo_ratio_eft_totalProcs_sm[ivar][1]);
        }


// ##### #    #  ####  #####   ##    ####  #    #
//   #   #    # #        #    #  #  #    # #   #
//   #   ######  ####    #   #    # #      ####
//   #   #    #      #   #   ###### #      #  #
//   #   #    # #    #   #   #    # #    # #   #
//   #   #    #  ####    #   #    #  ####  #    #

    	//-- Add legend entries -- iterate backwards, so that last histo stacked is on top of legend
        v_stack[ivar] = new THStack;

		for(int i=v_vector_MC_histo[ivar].size()-1; i>=0; i--)
		{
			if(!v_vector_MC_histo[ivar][i]) {continue;} //Some templates may be null
			v_stack[ivar]->Add(v_vector_MC_histo[ivar][i]);

            if(v_histo_total_MC[ivar] == NULL) {v_histo_total_MC[ivar] = (TH1F*) v_vector_MC_histo[ivar][i]->Clone();}
            else {v_histo_total_MC[ivar]->Add(v_vector_MC_histo[ivar][i]);}

			// cout<<"Stacking sample "<<MC_samples_legend[i]<<" / integral "<<v_vector_MC_histo[ivar][i]->Integral()<<endl;
            // cout<<"stack bin 1 content = "<<((TH1*) v_stack[ivar]->GetStack()->Last())->GetBinContent(1)<<endl;
		}


// ###### #####  #####   ####  #####   ####      ####  #####   ##    ####  #    #
// #      #    # #    # #    # #    # #         #        #    #  #  #    # #   #
// #####  #    # #    # #    # #    #  ####      ####    #   #    # #      ####
// #      #####  #####  #    # #####       #         #   #   ###### #      #  #
// #      #   #  #   #  #    # #   #  #    #    #    #   #   #    # #    # #   #
// ###### #    # #    #  ####  #    #  ####      ####    #   #    #  ####  #    #

        //-- Use pointers to vectors : need to give the adress of first element (all other elements can then be accessed iteratively)
        double* eyl = &v_eyl[0];
        double* eyh = &v_eyh[0];
        double* exl = &v_exl[0];
        double* exh = &v_exh[0];
        double* xx = &v_x[0];
        double* yy = &v_y[0];

        v_gr_error[ivar] = new TGraphAsymmErrors(nIndivBins,xx,yy,exl,exh,eyl,eyh);
        v_gr_error[ivar]->SetFillStyle(3254); //3002 //3004
        v_gr_error[ivar]->SetFillColor(kBlack);

        //-- Debug printouts
        // for(int ibin=1; ibin<nIndivBins+1; ibin++)
        // {
        //     cout<<"-- ibin "<<ibin<<endl;
        //     cout<<"v_eyh[ibin] "<<v_eyh[ibin]<<" / v_exl[ibin] "<<v_exl[ibin]<<" / v_exh[ibin] "<<v_exh[ibin]<<" / v_x[ibin] "<<v_x[ibin]<<" / v_y[ibin] "<<v_y[ibin]<<endl;
        //     cout<<"v_hdata[ivar]->GetBinContent(ibin) "<<v_hdata[ivar]->GetBinContent(ibin)<<endl;
        //     cout<<"v_histo_total_MC[ivar]->GetBinContent(ibin) "<<v_histo_total_MC[ivar]->GetBinContent(ibin)<<endl;
        // }

        for(int ibin=1; ibin<nIndivBins+1; ibin++)
        {
            if(isnan(v_hdata[ivar]->GetBinContent(ibin))) {cout<<FRED("ERROR: v_hdata["<<ivar<<"]->GetBinContent("<<ibin<<") is nan ! May cause plotting bugs !")<<endl;}
        }


// #####  #####    ##   #    #
// #    # #    #  #  #  #    #
// #    # #    # #    # #    #
// #    # #####  ###### # ## #
// #    # #   #  #    # ##  ##
// #####  #    # #    # #    #

        c1->cd(ivar+1);

        //Draw stack
        v_stack[ivar]->Draw("hist");

        v_hdata[ivar]->Draw("e0p same");

        v_gr_error[ivar]->Draw("e2 same"); //Superimposes the uncertainties on stack


// #   # #    #   ##   #    #
//  # #  ##  ##  #  #   #  #
//   #   # ## # #    #   ##
//   #   #    # ######   ##
//   #   #    # #    #  #  #
//   #   #    # #    # #    #

        //-- Set minimum
        v_stack[ivar]->SetMinimum(1.5);
        // v_stack[ivar]->SetMinimum(2.5);

        if(total_var_list[ivar].Contains("ctw_SRtZq")) {v_stack[ivar]->SetMinimum(1.);} //Very low prediction in last bin

        //-- Set Yaxis maximum
        double ymax = 0;
        ymax = v_hdata[ivar]->GetMaximum(); //Data ymax
        if(ymax < v_stack[ivar]->GetMaximum()) {ymax = v_stack[ivar]->GetMaximum();} //MC ymax
        if(isLogY) {ymax*= 6.;}
        else {ymax*= 1.45;}
        v_stack[ivar]->SetMaximum(ymax);
        c1->Modified();


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
        v_tpad_ratio[ivar] = new TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 1.0);
        float middle_panel_size = 0.27; //NB: must adjust middle panel size (margins computed w.r.t. pad height)
        if(template_name.Contains("SM"))
        {
            v_tpad_ratio[ivar]->SetTopMargin(0.70);
            v_tpad_ratio[ivar]->SetBottomMargin(0.13);
        }
        else
        {
            v_tpad_ratio[ivar]->SetTopMargin(top_panel_size);
            v_tpad_ratio[ivar]->SetBottomMargin(middle_panel_size);
        }
        v_tpad_ratio[ivar]->SetFillColor(0);
    	v_tpad_ratio[ivar]->SetFillStyle(0);
    	v_tpad_ratio[ivar]->SetGridy(1);
    	v_tpad_ratio[ivar]->Draw();
    	v_tpad_ratio[ivar]->cd(0);

        if(rightmargin>0) {v_tpad_ratio[ivar]->SetRightMargin(rightmargin);}
        // if(template_name.Contains("SM")) {v_tpad_ratio[ivar]->SetBottomMargin(0.10);}

        v_histo_ratio_data[ivar] = (TH1F*) v_hdata[ivar]->Clone();

    	if(!show_pulls_ratio) //Compute ratios (with error bars)
    	{
    		//To get correct error bars in ratio plot, must only account for errors from data, not MC ! (MC error shown as separate band)
    		for(int ibin=1; ibin<v_histo_total_MC[ivar]->GetNbinsX()+1; ibin++)
    		{
    			v_histo_total_MC[ivar]->SetBinError(ibin, 0.);
    		}

    		v_histo_ratio_data[ivar]->Divide(v_histo_total_MC[ivar]);
    	}
     	else //-- Compute pulls (no error bars)
    	{
    		for(int ibin=1; ibin<v_histo_ratio_data[ivar]->GetNbinsX()+1; ibin++)
    		{
    			//Add error on signal strength (since we rescale signal manually)
    			// double bin_error_mu = v_vector_MC_histo[ivar].at(index_tZq_sample)->GetBinError(ibin) * sig_strength_err;
    			// cout<<"bin_error_mu = "<<bin_error_mu<<endl;

    			double bin_error_mu = 0; //No sig strength uncert. for prefit ! //-- postfit -> ?

    			//Quadratic sum of systs, stat error, and sig strength error
    			double bin_error = pow(pow(v_histo_total_MC[ivar]->GetBinError(ibin), 2) + pow(v_histo_ratio_data[ivar]->GetBinError(ibin), 2) + pow(bin_error_mu, 2), 0.5);

    			if(!v_histo_total_MC[ivar]->GetBinError(ibin)) {v_histo_ratio_data[ivar]->SetBinContent(ibin,-99);} //Don't draw null markers
    			else{v_histo_ratio_data[ivar]->SetBinContent(ibin, (v_histo_ratio_data[ivar]->GetBinContent(ibin) - v_histo_total_MC[ivar]->GetBinContent(ibin)) / bin_error );}
    		}

    		//-- Don't draw null data
    		for(int ibin=1; ibin<v_histo_ratio_data[ivar]->GetNbinsX()+1; ibin++)
    		{
                if(std::isnan(v_histo_ratio_data[ivar]->GetBinContent(ibin)) || std::isinf(v_histo_ratio_data[ivar]->GetBinContent(ibin)) || v_histo_ratio_data[ivar]->GetBinContent(ibin) == 0) {v_histo_ratio_data[ivar]->SetBinContent(ibin, -99);}
    		}
    	}

    	if(show_pulls_ratio) {v_histo_ratio_data[ivar]->GetYaxis()->SetTitle("Pulls");}
        else {v_histo_ratio_data[ivar]->GetYaxis()->SetTitle("#frac{Data}{Pred.}");}
        // else {v_histo_ratio_data[ivar]->GetYaxis()->SetTitle("Data/MC");}
        v_histo_ratio_data[ivar]->GetYaxis()->CenterTitle(true);
    	// v_histo_ratio_data[ivar]->GetYaxis()->SetTickLength(0.);
        v_histo_ratio_data[ivar]->GetYaxis()->SetLabelFont(42);
        v_histo_ratio_data[ivar]->GetYaxis()->SetLabelSize(0.04);
    	v_histo_ratio_data[ivar]->GetXaxis()->SetTitleFont(42);
    	v_histo_ratio_data[ivar]->GetYaxis()->SetTitleFont(42);
        v_histo_ratio_data[ivar]->GetYaxis()->SetNdivisions(503); //grid draw on primary tick marks only
        v_histo_ratio_data[ivar]->GetXaxis()->SetNdivisions(505); //'-' to force Ndivisions
        v_histo_ratio_data[ivar]->GetYaxis()->SetTitleSize(0.04);
        v_histo_ratio_data[ivar]->GetYaxis()->SetTitleOffset(1.9); //1.5
        v_histo_ratio_data[ivar]->GetXaxis()->SetTickLength(0.01); //0.4
    	v_histo_ratio_data[ivar]->SetMarkerStyle(20);
    	v_histo_ratio_data[ivar]->SetMarkerSize(1.2);
        v_histo_ratio_data[ivar]->GetYaxis()->SetTickLength(0.15);
        v_histo_ratio_data[ivar]->GetYaxis()->SetNdivisions(303);

        if(template_name.Contains("SM")) //SM vs SM --> Display title in this panel
        {
            //-- SET X_AXIS TITLES
            v_histo_ratio_data[ivar]->GetXaxis()->SetTitle(Get_Template_XaxisTitle(total_var_list[ivar], true, this->use_NN_cpq3_SRttZ));
            v_histo_ratio_data[ivar]->GetYaxis()->SetTitleSize(0.05);
            v_histo_ratio_data[ivar]->GetYaxis()->SetTitleOffset(1.30);
            v_histo_ratio_data[ivar]->GetYaxis()->SetLabelSize(0.05);

            //-- Arbitrary bin names
            v_histo_ratio_data[ivar]->GetXaxis()->SetLabelOffset(0.02); //Add some x-label offset
            for(int i=1;i<v_histo_ratio_data[ivar]->GetNbinsX()+1;i++)
            {
                TString label = Convert_Number_To_TString(i);
                v_histo_ratio_data[ivar]->GetXaxis()->SetBinLabel(i, label);
            }

            v_histo_ratio_data[ivar]->GetXaxis()->SetLabelSize(0.09);
            v_histo_ratio_data[ivar]->GetXaxis()->SetTickLength(0.02);
            v_histo_ratio_data[ivar]->GetXaxis()->SetNdivisions(-v_histo_ratio_data[ivar]->GetNbinsX()); //'-' to force Ndivisions
        }
        else //SM vs EFT --> Will display title via additional bottom panel
        {
            v_histo_ratio_data[ivar]->GetXaxis()->SetNdivisions(-v2_histo_ratio_eft[ivar][0]->GetNbinsX()); //'-' to force Ndivisions //NB: must use same pattern as bottom TPad (if present) !
            v_histo_ratio_data[ivar]->GetXaxis()->SetTitleSize(0.);
            v_histo_ratio_data[ivar]->GetXaxis()->SetLabelSize(0.);
        }

        //-- If a point is outside the y-range of the ratio pad defined by SetMaximum/SetMinimum(), it disappears with its error
        //-- Trick: fill 2 histos with points either above/below y-range, to plot some markers indicating missing points (cleaner)
        //NB: only for ratio plot, not pulls
        float ratiopadmin = 0.4, ratiopadmax = 1.6; //Define ymin/ymax for ratio plot
        TH1F* h_pointsAboveY = (TH1F*) v_histo_ratio_data[ivar]->Clone();
        h_pointsAboveY->SetMarkerStyle(26); //Open triangle pointing up
        h_pointsAboveY->SetMarkerSize(1.5);
        TH1F* h_pointsBelowY = (TH1F*) v_histo_ratio_data[ivar]->Clone();
        h_pointsBelowY->SetMarkerStyle(32); //Open triangle pointing down
        h_pointsBelowY->SetMarkerSize(1.5);
        if(show_pulls_ratio)
    	{
    		v_histo_ratio_data[ivar]->SetMinimum(-2.99);
    		v_histo_ratio_data[ivar]->SetMaximum(2.99);
    	}
    	else
    	{
            //-- Default
            v_histo_ratio_data[ivar]->SetMinimum(ratiopadmin); //NB: removes error bars if data point is below ymin...?
            v_histo_ratio_data[ivar]->SetMaximum(ratiopadmax);

            //-- Fill histos with points outside yrange
            for(int ibin=1; ibin<v_histo_ratio_data[ivar]->GetNbinsX()+1; ibin++)
            {
                //-- Default: make point invisible
                h_pointsAboveY->SetBinContent(ibin, -999);
                h_pointsBelowY->SetBinContent(ibin, -999);

                if(v_histo_ratio_data[ivar]->GetBinContent(ibin) > ratiopadmax && v_hdata[ivar]->GetBinContent(ibin) >= 1)
                {
                    //Adjust error
                    float initial_y = v_histo_ratio_data[ivar]->GetBinContent(ibin);
                    float initial_err = v_histo_ratio_data[ivar]->GetBinError(ibin);
                    float new_err = initial_err - (initial_y-ratiopadmax);
                    if(new_err<0) {new_err=0.;}

                    h_pointsAboveY->SetBinContent(ibin, ratiopadmax-0.05); //Add some padding
                    h_pointsAboveY->SetBinError(ibin, new_err);
                }
                else if(v_histo_ratio_data[ivar]->GetBinContent(ibin) < ratiopadmin && v_hdata[ivar]->GetBinContent(ibin) >= 1)
                {
                    //Adjust error
                    float initial_y = v_histo_ratio_data[ivar]->GetBinContent(ibin);
                    float initial_err = v_histo_ratio_data[ivar]->GetBinError(ibin);
                    float new_err = initial_err - (ratiopadmin-initial_y);
                    if(new_err<0) {new_err=0.;}

                    h_pointsBelowY->SetBinContent(ibin, ratiopadmin+(ratiopadmin/10.)); //Add some padding
                    h_pointsBelowY->SetBinError(ibin, new_err);
                }
            }
    	}

    	if(show_pulls_ratio) {v_histo_ratio_data[ivar]->Draw("HIST P");} //Draw ratio points
        else
        {
            v_histo_ratio_data[ivar]->Draw("E1 X0 P"); //Draw ratio points ; E1 : perpendicular lines at end ; X0 : suppress x errors

            h_pointsAboveY->Draw("E1 X0 P same");
            h_pointsBelowY->Draw("E1 X0 P same");
        }


// ###### #####  #####   ####  #####   ####     #####    ##   ##### #  ####
// #      #    # #    # #    # #    # #         #    #  #  #    #   # #    #
// #####  #    # #    # #    # #    #  ####     #    # #    #   #   # #    #
// #      #####  #####  #    # #####       #    #####  ######   #   # #    #
// #      #   #  #   #  #    # #   #  #    #    #   #  #    #   #   # #    #
// ###### #    # #    #  ####  #    #  ####     #    # #    #   #   #  ####


		//Copy previous TGraphAsymmErrors, then modify it -> error TGraph for ratio plot
		TGraphAsymmErrors *thegraph_tmp = NULL;
		double *theErrorX_h;
		double *theErrorY_h;
		double *theErrorX_l;
		double *theErrorY_l;
		double *theY;
		double *theX;

		thegraph_tmp = (TGraphAsymmErrors*) v_gr_error[ivar]->Clone();
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

		v_gr_ratio_error[ivar] = new TGraphAsymmErrors(thegraph_tmp->GetN(), theX , theY ,  theErrorX_l, theErrorX_h, theErrorY_l, theErrorY_h);
        v_gr_ratio_error[ivar]->SetFillStyle(3254); //3002 //3004
        v_gr_ratio_error[ivar]->SetFillColor(kBlack); //kBlue+2 //kCyan

		if(!show_pulls_ratio) {v_gr_ratio_error[ivar]->Draw("e2 same");} //Draw error bands in ratio plot


// ###### ###### #####    #####    ##   ##### #  ####
// #      #        #      #    #  #  #    #   # #    #
// #####  #####    #      #    # #    #   #   # #    #
// #      #        #      #####  ######   #   # #    #
// #      #        #      #   #  #    #   #   # #    #
// ###### #        #      #    # #    #   #   #  ####

        if(!template_name.Contains("SM")) //Add second bottom TPad for SMEFT/SM ratio (not for SM prefit plots)
        {
            //-- Create subpad to plot ratio
            v_tpad_ratio_eft[ivar] = new TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 1.0);
            v_tpad_ratio_eft[ivar]->SetTopMargin(1-middle_panel_size); //+0.02 to add space
            v_tpad_ratio_eft[ivar]->SetFillColor(0);
            v_tpad_ratio_eft[ivar]->SetFillStyle(0);
            // v_tpad_ratio_eft[ivar]->SetGridy(1);
            v_tpad_ratio_eft[ivar]->Draw();
            v_tpad_ratio_eft[ivar]->cd(0);
            // v_tpad_ratio_eft[ivar]->SetLogy();
            if(rightmargin>0) {v_tpad_ratio_eft[ivar]->SetRightMargin(rightmargin);}
            v_tpad_ratio_eft[ivar]->SetBottomMargin(0.10); //Decrease from default 0.13

            //-- Debug printouts
            // cout<<"v2_histo_ratio_eft[ivar][0]->GetBinContent(1) "<<v2_histo_ratio_eft[ivar][0]->GetBinContent(1)<<endl;
            // cout<<"v_histo_ratio_sm[ivar]->GetBinContent(1) "<<v_histo_ratio_sm[ivar]->GetBinContent(1)<<endl;

            v2_histo_ratio_eft[ivar][0]->Divide(v_histo_ratio_sm[ivar]);
            v2_histo_ratio_eft[ivar][1]->Divide(v_histo_ratio_sm[ivar]);

            //-- Debug printouts
            // for(int ibin=1; ibin<v2_histo_ratio_eft[ivar][0]->GetNbinsX()+1; ibin++)
            // {
            //     cout<<"bin "<<ibin<<" EFT1 --> "<<v2_histo_ratio_eft[ivar][0]->GetBinContent(ibin)<<endl;
            //     cout<<"bin "<<ibin<<" EFT2 --> "<<v2_histo_ratio_eft[ivar][1]->GetBinContent(ibin)<<endl;
            // }

            v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetTitle("#frac{SM+EFT}{SM}");

            v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetTickLength(0.);
            v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetTitleOffset(1.9); //1.8
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetTitleOffset(1.05);
            if(total_var_list[ivar].Contains("cpq3_SRtZq")) {v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetTitleOffset(1.04);}
            if(total_var_list[ivar].Contains("ttZ")) {v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetTitleOffset(0.95);} //Takes more space
            v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetLabelSize(0.04);
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetLabelSize(0.07);
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetLabelFont(42);
            v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetLabelFont(42);
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetTitleFont(42);
            v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetTitleFont(42);
            v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetNdivisions(505); //grid drawn on primary tick marks only
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetNdivisions(505);
            v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetTitleSize(0.04);
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetTickLength(0.01); //0.4
            v2_histo_ratio_eft[ivar][0]->SetMarkerStyle(20);
            v2_histo_ratio_eft[ivar][0]->SetMarkerSize(1.2); //changed from 1.4
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetTitleSize(0.05); //changed from 0.06
            v2_histo_ratio_eft[ivar][0]->GetYaxis()->CenterTitle(true);
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetTitleSize(0.045); //changed from 0.06 //NB: when using 0.05, got extra space before tbar in cpq3_SRttZ, for no reason...
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetNdivisions(-v2_histo_ratio_eft[ivar][0]->GetNbinsX()); //'-' to force Ndivisions

            v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetTickLength(0.1);
            // v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetNdivisions(503); //grid drawn on primary tick marks only

            //-- SET X_AXIS TITLES
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetTitle(Get_Template_XaxisTitle(total_var_list[ivar], true, this->use_NN_cpq3_SRttZ));

            //-- Arbitrary bin names
            // v2_histo_ratio_eft[ivar][0]->GetYaxis()->SetMoreLogLabels();
            v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetLabelOffset(0.02); //Add some x-label offset
            for(int i=1;i<v2_histo_ratio_eft[ivar][0]->GetNbinsX()+1;i++)
            {
                // TString label = "Bin " + Convert_Number_To_TString(i);
                TString label = Convert_Number_To_TString(i);
                v2_histo_ratio_eft[ivar][0]->GetXaxis()->SetBinLabel(i, label);
                // v2_histo_ratio_eft[ivar][0]->LabelsOption("v", "X"); //X labels vertical
            }

            v2_histo_ratio_eft[ivar][0]->SetMinimum(0.8); //Default
            // v2_histo_ratio_eft[ivar][0]->SetMinimum(1.);
            // v2_histo_ratio_eft[ivar][0]->SetMinimum(1.01);
            v2_histo_ratio_eft[ivar][0]->SetMaximum(9.99); //Default

            if(total_var_list[ivar].Contains("ctz_SRttZ"))
            {
                // v2_histo_ratio_eft[ivar][0]->SetMinimum(1.);
                v2_histo_ratio_eft[ivar][0]->SetMaximum(12.99);
            }
            else if(total_var_list[ivar].Contains("ctz_SRtZq"))
            {
                v2_histo_ratio_eft[ivar][0]->SetMaximum(4.99);
            }
            else if(total_var_list[ivar].Contains("ctw_SRtZq"))
            {
                v2_histo_ratio_eft[ivar][0]->SetMaximum(6.69);
                // v2_histo_ratio_eft[ivar][0]->SetMaximum(11.99);
            }
            else if(total_var_list[ivar].Contains("ctw_SRttZ"))
            {
                v2_histo_ratio_eft[ivar][0]->SetMaximum(2.49);
                // v2_histo_ratio_eft[ivar][0]->SetMaximum(3.99);
            }
            else if(total_var_list[ivar].Contains("cpq3_SRtZq"))
            {
                v2_histo_ratio_eft[ivar][0]->SetMaximum(5.99);
            }
            else if(total_var_list[ivar].Contains("cpq3_SRttZ"))
            {
                if(!this->use_NN_cpq3_SRttZ) {v2_histo_ratio_eft[ivar][0]->SetMaximum(1.069);}
                else {v2_histo_ratio_eft[ivar][0]->SetMaximum(1.129);}

                // v2_histo_ratio_eft[ivar][0]->SetMaximum(1.139);
                v2_histo_ratio_eft[ivar][0]->SetMinimum(1.);
            }
            else if(total_var_list[ivar].Contains("5D_SRttZ"))
            {
                // v2_histo_ratio_eft[ivar][0]->SetMaximum(1.99);
                if(!this->use_NN_cpq3_SRttZ) {v2_histo_ratio_eft[ivar][0]->SetMaximum(1.79);}
                else {v2_histo_ratio_eft[ivar][0]->SetMaximum(1.99);}
                v2_histo_ratio_eft[ivar][0]->SetMinimum(0.95);
            }
            else if(total_var_list[ivar].Contains("5D_SRtZq"))
            {
                v2_histo_ratio_eft[ivar][0]->SetMaximum(9.99);
            }
            else if(total_var_list[ivar].Contains("SM"))
            {
                if(total_var_list[ivar].Contains("SRtZq"))
                {
                    v2_histo_ratio_eft[ivar][0]->SetMaximum(1.19);
                }
                else
                {
                    v2_histo_ratio_eft[ivar][0]->SetMinimum(0.4);
                    v2_histo_ratio_eft[ivar][0]->SetMaximum(1.69);
                }
            }

            v2_histo_ratio_eft[ivar][0]->SetLineWidth(3.);
            v2_histo_ratio_eft[ivar][1]->SetLineWidth(3.);

            // v2_histo_ratio_eft[ivar][0]->SetLineColor(kRed);
            // v2_histo_ratio_eft[ivar][1]->SetLineColor(kBlue);
            // v2_histo_ratio_eft[ivar][0]->SetLineColor(v_custom_colors[10]->GetNumber());
            // v2_histo_ratio_eft[ivar][1]->SetLineColor(v_custom_colors[15]->GetNumber());
            // v2_histo_ratio_eft[ivar][0]->SetLineColor(v_custom_colors[12]->GetNumber());
            // v2_histo_ratio_eft[ivar][0]->SetLineColor(v_custom_colors[4]->GetNumber());
            // v2_histo_ratio_eft[ivar][0]->SetLineColor(v_custom_colors[8]->GetNumber());
            v2_histo_ratio_eft[ivar][1]->SetLineColor(v_custom_colors[11]->GetNumber());
            v2_histo_ratio_eft[ivar][0]->SetLineColor(kMagenta+2);
            v2_histo_ratio_eft[ivar][1]->SetLineColor(kPink+1);

            v2_histo_ratio_eft[ivar][0]->Draw("hist");
            v2_histo_ratio_eft[ivar][1]->Draw("hist same");

            TString EFTpointlabel1, EFTpointlabel2;
            vector<pair<TString,float>> v = Parse_EFTreweight_ID(EFTpoint_name1);
            for(int i=0; i<v.size(); i++)
            {
                if(v[i].second != 0)
                {
                    if(EFTpointlabel1 != "") {EFTpointlabel1+= ",";}
                    if(total_var_list[ivar].Contains("5D")) {EFTpointlabel1+= Convert_Number_To_TString(v[i].second);}
                    else {EFTpointlabel1+= Convert_Number_To_TString(v[i].second);}
                    // else {EFTpointlabel1+= Get_EFToperator_label((TString) v[i].first) + "=" + Convert_Number_To_TString(v[i].second);}
                }
            }
            v = Parse_EFTreweight_ID(EFTpoint_name2);
            for(int i=0; i<v.size(); i++)
            {
                if(v[i].second != 0)
                {
                    if(EFTpointlabel2 != "") {EFTpointlabel2+= ",";}
                    if(total_var_list[ivar].Contains("5D")) {EFTpointlabel2+= Convert_Number_To_TString(v[i].second);}
                    else {EFTpointlabel2+= Convert_Number_To_TString(v[i].second);}
                    // else {EFTpointlabel2+= Get_EFToperator_label((TString) v[i].first) + "=" + Convert_Number_To_TString(v[i].second);}
                }
            }
            if(total_var_list[ivar].Contains("5D"))
            {
                EFTpointlabel1 = "("+EFTpointlabel1+")";
                EFTpointlabel2 = "("+EFTpointlabel2+")";
            }

            TString proc_tmp = "";
            if(total_var_list[ivar].Contains("SRtZq")) {proc_tmp = "tZq";}
            if(total_var_list[ivar].Contains("SRttZ")) {proc_tmp = "t#bar{t}Z";}
            // if(total_var_list[ivar].Contains("5D")) {header+= "("+Get_EFToperator_label("ctz")+","+Get_EFToperator_label("ctw")+","+Get_EFToperator_label("cpq3")+")";}

            if(total_var_list[ivar].Contains("5D")) //5D: more text = more space
            {
                if(total_var_list[ivar].Contains("SRttZ") && !this->use_NN_cpq3_SRttZ)
                {
                    v_tlegend_ratio_eft[ivar] = new TLegend(0.18,0.16,0.80,0.27); //Smaller x
                }
                else
                {
                    v_tlegend_ratio_eft[ivar] = new TLegend(0.18,0.16,0.86,0.27);
                }
                // v_tlegend_ratio_eft[ivar]->SetTextSize(0.03);

                // TString header = "("+Get_EFToperator_label("ctz")+", "+Get_EFToperator_label("ctw")+", "+Get_EFToperator_label("cpq3")+")";
                // if(total_var_list[ivar].Contains("ttZ")) {header = "("+Get_EFToperator_label("ctz")+", "+Get_EFToperator_label("ctw")+")";}

                TString header = "("+Get_EFToperator_label("ctz")+" /#Lambda^{2}, "+Get_EFToperator_label("ctw")+" /#Lambda^{2}, "+Get_EFToperator_label("cpq3")+" /#Lambda^{2})";
                if(total_var_list[ivar].Contains("ttZ") && !this->use_NN_cpq3_SRttZ) {header = "("+Get_EFToperator_label("ctz")+" /#Lambda^{2}, "+Get_EFToperator_label("ctw")+" /#Lambda^{2})";}
                header+= " [TeV^{-2}]";

                v_tlegend_ratio_eft[ivar]->SetHeader(header);
                v_tlegend_ratio_eft[ivar]->SetNColumns(4);
                v_tlegend_ratio_eft[ivar]->SetTextSize(0.04);
                // v_tlegend_ratio_eft[ivar]->SetTextSize(0.035);
            }
            else //Default
            {
                v_tlegend_ratio_eft[ivar] = new TLegend(0.18,0.17,0.80,0.27);
                // v_tlegend_ratio_eft[ivar] = new TLegend(0.18,0.17,0.65,0.26);
                v_tlegend_ratio_eft[ivar]->SetTextSize(0.04);
                // v_tlegend_ratio_eft[ivar]->SetNColumns(2);
                v_tlegend_ratio_eft[ivar]->SetNColumns(4);

                TString header = Get_EFToperator_label("ctz")+" /#Lambda^{2}";
                if(total_var_list[ivar].Contains("ctw")) {header = Get_EFToperator_label("ctw")+" /#Lambda^{2}";}
                else if(total_var_list[ivar].Contains("cpq3")) {header = Get_EFToperator_label("cpq3")+" /#Lambda^{2}";}
                header+= " [TeV^{-2}]";
                v_tlegend_ratio_eft[ivar]->SetHeader(header);
            }
            // v_tlegend_ratio_eft[ivar]->SetTextSize(0.035);
            // v_tlegend_ratio_eft[ivar]->SetMargin(0.2); //x-axis fractional size of legend entry symbol //Default 0.25
            v_tlegend_ratio_eft[ivar]->SetBorderSize(0);
            v_tlegend_ratio_eft[ivar]->SetFillStyle(0); //transparent
            v_tlegend_ratio_eft[ivar]->SetTextAlign(12); //align = 10*HorizontalAlign + VerticalAlign //Horiz: 1=left adjusted, 2=centered, 3=right adjusted //Vert: 1=bottom adjusted, 2=centered, 3=top adjusted

            //-- Dummy histograms to add 'solid/dashed lines' in legend
            v_tlegend_dummy[ivar][0] = new TH1F("", "", 1, 0, 1);
            v_tlegend_dummy[ivar][1] = new TH1F("", "", 1, 0, 1);
            v_tlegend_dummy[ivar][0]->SetFillColor(kBlack);
            v_tlegend_dummy[ivar][1]->SetFillColor(kBlack);
            v_tlegend_dummy[ivar][1]->SetLineStyle(2.);
            v_tlegend_dummy[ivar][0]->SetLineWidth(2);
            v_tlegend_dummy[ivar][1]->SetLineWidth(2);

            v_tlegend_ratio_eft[ivar]->AddEntry(v2_histo_ratio_eft[0][0], EFTpointlabel1, "L");
            v_tlegend_ratio_eft[ivar]->AddEntry(v2_histo_ratio_eft[0][1], EFTpointlabel2, "L");
            // v_tlegend_ratio_eft[ivar]->AddEntry((TObject*)0, "", "");

            v_tlegend_ratio_eft[ivar]->AddEntry(v_tlegend_dummy[ivar][0], proc_tmp, "L");
            // v_tlegend_ratio_eft[ivar]->AddEntry(v_tlegend_dummy[ivar][0], proc_tmp+" only", "L");
            if(total_var_list[ivar].Contains("5D")) {v_tlegend_ratio_eft[ivar]->AddEntry(v_tlegend_dummy[ivar][1], "Total pred.", "L");}
            else {v_tlegend_ratio_eft[ivar]->AddEntry(v_tlegend_dummy[ivar][1], "Total prediction", "L");}

            v_tlegend_ratio_eft[ivar]->Draw("same");

            //-- Obsolete
            // if(total_var_list[ivar].Contains("SM") || total_var_list[ivar].Contains("cpq3_SRttZ"))
            // {
            //     v_tlatex_legend_eft[ivar] = new TLatex();
            //     v_tlatex_legend_eft[ivar]->SetNDC();
            //     v_tlatex_legend_eft[ivar]->SetTextFont(42);
            //     v_tlatex_legend_eft[ivar]->SetTextAlign(31);
            //     v_tlatex_legend_eft[ivar]->SetTextSize(0.04);
            //     v_tlatex_legend_eft[ivar]->DrawLatex(0.24, 0.16, header);
            // }

            v2_histo_ratio_eft_totalProcs_smeft[ivar][0]->SetLineColor(kMagenta+2);
            v2_histo_ratio_eft_totalProcs_smeft[ivar][1]->SetLineColor(kPink+1);
            // v2_histo_ratio_eft_totalProcs_smeft[ivar][1]->SetLineColor(v_custom_colors[11]->GetNumber());
            // v2_histo_ratio_eft_totalProcs_smeft[ivar][0]->SetLineColor(v_custom_colors[8]->GetNumber());
            // v2_histo_ratio_eft_totalProcs_smeft[ivar][1]->SetLineColor(v_custom_colors[11]->GetNumber());
            v2_histo_ratio_eft_totalProcs_smeft[ivar][0]->SetLineStyle(2);
            v2_histo_ratio_eft_totalProcs_smeft[ivar][1]->SetLineStyle(2);
            v2_histo_ratio_eft_totalProcs_smeft[ivar][0]->SetLineWidth(3.);
            v2_histo_ratio_eft_totalProcs_smeft[ivar][1]->SetLineWidth(3.);
            v2_histo_ratio_eft_totalProcs_smeft[ivar][0]->Draw("hist same");
            v2_histo_ratio_eft_totalProcs_smeft[ivar][1]->Draw("hist same");
        }


//  ####   ####   ####  #    # ###### ##### #  ####   ####
// #    # #    # #      ##  ## #        #   # #    # #
// #      #    #  ####  # ## # #####    #   # #       ####
// #      #    #      # #    # #        #   # #           #
// #    # #    # #    # #    # #        #   # #    # #    #
//  ####   ####   ####  #    # ######   #   #  ####   ####

    	//-- Draw ratio y-lines manually
        v_tpad_ratio[ivar]->cd();
    	v_hlines1[ivar] = new TH1F("","",this->nbins, xmin_tmp, xmax_tmp);
    	v_hlines2[ivar] = new TH1F("","",this->nbins, xmin_tmp, xmax_tmp);
    	for(int ibin=1; ibin<this->nbins +1; ibin++)
    	{
    		if(show_pulls_ratio)
    		{
    			v_hlines1[ivar]->SetBinContent(ibin, -1);
    			v_hlines2[ivar]->SetBinContent(ibin, 1);
    		}
    		else
    		{
                v_hlines1[ivar]->SetBinContent(ibin, 0.75);
                v_hlines2[ivar]->SetBinContent(ibin, 1.25);
    		}
    	}
    	v_hlines1[ivar]->SetLineStyle(6); v_hlines2[ivar]->SetLineStyle(6);
    	// v_hlines1[ivar]->Draw("hist same"); v_hlines2[ivar]->Draw("hist same"); //Removed extra lines

        TString Y_label = "Events / bin";
        if(v_stack[ivar]) //Must be drawn first
    	{
    		v_stack[ivar]->GetYaxis()->SetLabelFont(42);
    		v_stack[ivar]->GetYaxis()->SetTitleFont(42);
    		v_stack[ivar]->GetYaxis()->SetTitleSize(0.06);
            v_stack[ivar]->GetYaxis()->SetTickLength(0.04);
    		v_stack[ivar]->GetYaxis()->SetLabelSize(0.048);
    		v_stack[ivar]->GetYaxis()->SetNdivisions(506);
            v_stack[ivar]->GetYaxis()->SetTitleOffset(1.20); //1.15
    		v_stack[ivar]->GetYaxis()->SetTitle(Y_label);
            v_stack[ivar]->GetXaxis()->SetLabelSize(0.);
            v_stack[ivar]->GetXaxis()->SetTickLength(0.);
    	}

    	//----------------
    	// CAPTIONS //
    	//----------------
    	// -- using https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines
        // -- About fonts: https://root.cern.ch/doc/master/classTAttText.html#T5

    	float l = pad->GetLeftMargin();
    	float t = pad->GetTopMargin();

    	TString cmsText = "CMS";
    	TLatex latex;
    	latex.SetNDC();
    	latex.SetTextColor(kBlack);
        latex.SetTextFont(62); //Changed
    	latex.SetTextAlign(11);
    	latex.SetTextSize(0.06);
        // latex.DrawLatex(l + 0.04, 0.87, cmsText);
        if(use_paperStyle) {latex.DrawLatex(l + 0.04, 0.87, cmsText);} //CMS guideline: within frame
        else {latex.DrawLatex(l + 0.01, 0.94, cmsText);} //Default: outside frame

    	float lumi = lumiValue;
    	TString lumi_ts = Convert_Number_To_TString(lumi);
    	lumi_ts += " fb^{-1} (13 TeV)";
    	latex.SetTextFont(42);
    	latex.SetTextAlign(31);
    	latex.SetTextSize(0.04);
        latex.DrawLatex(0.96, 0.94,lumi_ts);

        TString extraText = "Preliminary";
        if(!use_paperStyle) //Default is without (for paper)
        {
            latex.SetTextFont(52);
            latex.SetTextSize(0.05);
            if(total_var_list[ivar] == "NN_SM_SRother") {latex.DrawLatex(l + 0.40, 0.94, extraText);} //Hardcoded: due to different pad size (?), need to add a bit of extra space for this one
            else {latex.DrawLatex(l + 0.38, 0.94, extraText);}
        }

        TString info_data = Get_Region_Label(region, total_var_list[ivar]);
        TLatex text2;
        text2.SetNDC();
        text2.SetTextAlign(13);
        text2.SetTextSize(0.04);
        // text2.SetTextSize(0.045);
        text2.SetTextFont(42);
        if(info_data != "")
        {
            float left = l+0.04;
            if(total_var_list[ivar].Contains("other")) {left = l+0.03;} //Need more space
            if(use_paperStyle) {text2.DrawLatex(left,0.85,info_data);} //0.83
            else {text2.DrawLatex(left,0.90,info_data);} //Default: outside frame //0.90
        }

        pad = NULL;
    } //Var loop
//--------------------------------------------


// ####### #
//    #    #       ######  ####  ###### #    # #####
//    #    #       #      #    # #      ##   # #    #
//    #    #       #####  #      #####  # #  # #    #
//    #    #       #      #  ### #      #  # # #    #
//    #    #       #      #    # #      #   ## #    #
//    #    ####### ######  ####  ###### #    # #####

    int ivar = 0; //Read first element by default

    int n_columns = ceil(nSampleGroups/2.) > 6 ? 6 : ceil(nSampleGroups/2.); //ceil = upper int
    float x_left = 0.94-n_columns*0.12; //Each column allocated same x-space //0.12 needed for most crowded plots
    if(x_left < 0.4) {x_left = 0.4;} //Leave some space for region label

    TLegend* qw = new TLegend(x_left-0.05,0.80,0.94,0.92); //Default /
    qw->SetTextSize(0.04);
    qw->SetNColumns(n_columns);
    qw->SetBorderSize(0);
    qw->SetFillStyle(0); //transparent
    qw->SetTextAlign(12); //align = 10*HorizontalAlign + VerticalAlign //Horiz: 1=left adjusted, 2=centered, 3=right adjusted //Vert: 1=bottom adjusted, 2=centered, 3=top adjusted
    // cout<<"x_left "<<x_left<<endl;
    // cout<<"ceil(nSampleGroups/2.) "<<ceil(nSampleGroups/2.)<<endl;

    //-- Dummy object, only used to display uncertainty band also in legend
    TH1F* h_uncert = new TH1F("h_uncert", "h_uncert", 1, 0, 1);
    h_uncert->SetFillStyle(3254); //3002 //3004
    h_uncert->SetFillColor(kBlack);
    h_uncert->SetLineWidth(0.);
    qw->AddEntry(h_uncert, "Unc.", "F");

	//-- Data on top of legend
    qw->AddEntry(v_hdata[0], "Data" , "ep");

	for(int i=0; i<v_vector_MC_histo[ivar].size(); i++)
	{
        if(MC_samples_legend[i].Contains("tZq")) {qw->AddEntry(v_vector_MC_histo[ivar][i], "tZq", "f");}
        else if(MC_samples_legend[i].EndsWith("ttZ") ) {qw->AddEntry(v_vector_MC_histo[ivar][i], "t#bar{t}Z", "f");}
        else if(MC_samples_legend[i].EndsWith("tWZ") ) {qw->AddEntry(v_vector_MC_histo[ivar][i], "tWZ", "f");}
        else if(MC_samples_legend[i] == "ttW" || MC_samples_legend[i] == "tX") {qw->AddEntry(v_vector_MC_histo[ivar][i], "t(#bar{t})X", "f");}
        else if(MC_samples_legend[i] == "WZ") {qw->AddEntry(v_vector_MC_histo[ivar][i], "WZ", "f");}
        else if(MC_samples_legend[i] == "WWZ" || MC_samples_legend[i] == "VVV") {qw->AddEntry(v_vector_MC_histo[ivar][i], "VV(V)", "f");}
        else if(MC_samples_legend[i] == "TTGamma_Dilep" || MC_samples_legend[i] == "XG") {qw->AddEntry(v_vector_MC_histo[ivar][i], "X#gamma", "f");}
        else if(MC_samples_legend[i] == "TTbar_DiLep" || MC_samples_legend[i] == "NPL" || MC_samples_legend[i] == "NPL_DATA") {qw->AddEntry(v_vector_MC_histo[ivar][i], "NPL", "f");}
	}

    for(int ivar=0; ivar<nvar; ivar++)
    {
        c1->cd(ivar+1);
        qw->Draw("same");
    }


// #    # #####  # ##### ######     ####  #    # ##### #####  #    # #####
// #    # #    # #   #   #         #    # #    #   #   #    # #    #   #
// #    # #    # #   #   #####     #    # #    #   #   #    # #    #   #
// # ## # #####  #   #   #         #    # #    #   #   #####  #    #   #
// ##  ## #   #  #   #   #         #    # #    #   #   #      #    #   #
// #    # #    # #   #   ######     ####   ####    #   #       ####    #

    TString outdir = "plots/paperPlots/";
    mkdir(outdir.Data(), 0777);
    TString output_plot_name = outdir + ""+fit_type+"Templates_signalRegions_";
    output_plot_name+= template_name;
    if(!use_paperStyle) {output_plot_name+= "_prelim";}
    c1->SaveAs(output_plot_name + ".png");
    c1->SaveAs(output_plot_name + ".eps");
    c1->SaveAs(output_plot_name + ".pdf");

    delete c1; c1 = NULL;
    delete qw; qw = NULL;
    if(h_uncert) {delete h_uncert; h_uncert = NULL;}

    for(int ivar=0; ivar<nvar; ivar++)
    {
        if(v_gr_error[ivar]) {delete v_gr_error[ivar]; v_gr_error[ivar] = NULL;}
        if(v_stack[ivar]) {delete v_stack[ivar]; v_stack[ivar] = NULL;}
        if(v_gr_ratio_error[ivar]) {delete v_gr_ratio_error[ivar]; v_gr_ratio_error[ivar] = NULL;}
        if(v_hlines1[ivar]) {delete v_hlines1[ivar]; v_hlines1[ivar] = NULL;}
        if(v_hlines2[ivar]) {delete v_hlines2[ivar]; v_hlines2[ivar] = NULL;}
        if(v_tlatex_legend_eft[ivar]) {delete v_tlatex_legend_eft[ivar]; v_tlatex_legend_eft[ivar] = NULL;}
        if(v_tlegend_ratio_eft[ivar]) {delete v_tlegend_ratio_eft[ivar]; v_tlegend_ratio_eft[ivar] = NULL;}

        if(v2_histo_ratio_eft_totalProcs_sm[ivar][0]) {delete v2_histo_ratio_eft_totalProcs_sm[ivar][0]; v2_histo_ratio_eft_totalProcs_sm[ivar][0] = NULL;}
        if(v2_histo_ratio_eft_totalProcs_sm[ivar][1]) {delete v2_histo_ratio_eft_totalProcs_sm[ivar][1]; v2_histo_ratio_eft_totalProcs_sm[ivar][1] = NULL;}
        if(v2_histo_ratio_eft[ivar][0]) {delete v2_histo_ratio_eft[ivar][0]; v2_histo_ratio_eft[ivar][0] = NULL;}
        if(v2_histo_ratio_eft[ivar][1]) {delete v2_histo_ratio_eft[ivar][1]; v2_histo_ratio_eft[ivar][1] = NULL;}
        if(v_tlegend_dummy[ivar][0]) {delete v_tlegend_dummy[ivar][0]; v_tlegend_dummy[ivar][0] = NULL;}
        if(v_tlegend_dummy[ivar][1]) {delete v_tlegend_dummy[ivar][1]; v_tlegend_dummy[ivar][1] = NULL;}
    }

    return;
}










/**
 * Hard-coded function to produce control plots for the paper.
 */
void TopEFT_analysis::Make_PaperPlot_ControlPlots()
{
//--------------------------------------------

    cout<<endl<<BYEL("                          ")<<endl<<endl;
	cout<<FYEL("--- Producing paper figure [CONTROL PLOTS] ---")<<endl;
    cout<<endl<<BYEL("                          ")<<endl<<endl;


//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

    //--------------------------
    // DEFINE VARIABLES TO PLOT
    //--------------------------

    //-- Hardcoded list of variables
    vector<TString> total_var_list;
    total_var_list.clear();
    total_var_list.push_back("recoZ_dPhill");
    total_var_list.push_back("maxDeepJet");
    total_var_list.push_back("jPrimeAbsEta");
    total_var_list.push_back("nbjets");
    total_var_list.push_back("lAsymmetry");
    total_var_list.push_back("metEt");

    // total_var_list.push_back("jprime_Pt");
    // total_var_list.push_back("mTW");
    // total_var_list.push_back("nbjets");
    // total_var_list.push_back("jPrimeAbsEta");
    // total_var_list.push_back("jprime_Pt");
    // total_var_list.push_back("metEt");
    // total_var_list.push_back("njets");
    // total_var_list.push_back("mTW");
    // total_var_list.push_back("lAsymmetry");
    // total_var_list.push_back("recoZ_Pt");
    // total_var_list.push_back("recoZ_dPhill");
    // total_var_list.push_back("dR_tZ");
    int nvar = total_var_list.size();

    //--------------------------
    // DIVIDE CANVAS IN 4 PADS
    //--------------------------

    //-- Canvas definition
    Load_Canvas_Style(); //Default top/bottom/left/right margins: 0.07/0.13/0.16/0.03
    TCanvas* c1 = new TCanvas("c1","c1", 1200, 800);
    // TCanvas* c1 = new TCanvas("c1","c1", 1400, 800);
    c1->Divide(3, 2, 1E-11, 1E-11); //(x,y)

    //--------------------------
    // CREATE VECTORS OF OBJECTS
    //--------------------------

    vector<TH1F*> v_hdata(nvar);
    vector<THStack*> v_stack(nvar);
    vector<TH1F*> v_histo_total_MC(nvar);
    vector<vector<TH1F*> > v_vector_MC_histo(nvar); //Store separately the histos for each MC sample --> stack them after loops
    vector<TH1F*> v_histo_ratio_data(nvar);
    vector<TH1F*> v_hlines1(nvar), v_hlines2(nvar); //Draw TLinesin ratio plots

    vector<TGraphAsymmErrors*> v_gr_error(nvar);
    vector<TGraphAsymmErrors*> v_gr_ratio_error(nvar);

    vector<TPad*> v_tpad_ratio(nvar);

    vector<TString> MC_samples_legend; //List the MC samples to mention in legend

    int nIndivBins; float xmin_tmp, xmax_tmp;


// #       ####   ####  #####   ####
// #      #    # #    # #    # #
// #      #    # #    # #    #  ####
// #      #    # #    # #####       #
// #      #    # #    # #      #    #
// ######  ####   ####  #       ####

//--------------------------------------------
	for(int ivar=0; ivar<nvar; ivar++)
	{
        cout<<endl<<FBLU("== VARIABLE: "<<total_var_list[ivar]<<"")<<endl;

    	TString primitive_name = "c1_" + Convert_Number_To_TString(ivar+1);
        TPad* pad = (TPad*) c1->GetPrimitive(primitive_name);
        TString inputfilename = "";
        float rightmargin = -1; //Default (use default right margin)
        if(ivar==2 || ivar==5) {rightmargin = 0.04;} //Add some margin for rightmost plots

        // if(total_var_list[ivar].Contains("SRtZq")) //Left
        if(ivar==0) //Left
    	{
    		pad->SetPad(0, 0.5, 0.33, 1);
    	}
    	else if(ivar==1) //Right
    	{
    		pad->SetPad(0.33, 0.5, 0.66, 1);
        }
        else if(ivar==2)
        {
            pad->SetPad(0.66, 0.5, 1., 1);
        }
        else if(ivar==3)
        {
            pad->SetPad(0., 0., 0.33, 0.5);
        }
        else if(ivar==4)
        {
            pad->SetPad(0.33, 0., 0.66, 0.5);
        }
        else if(ivar==5)
        {
            pad->SetPad(0.66, 0., 1., 0.5);
        }
        if(rightmargin>0) {pad->SetRightMargin(rightmargin);}
        pad->SetBottomMargin(0.30);

        TFile* file_input = NULL;
        TH1F* h_tmp = NULL; //Tmp storing histo
        TH1F* hdata_tmp = NULL; //Tmp storing data histo

		//-- Init error vectors
		double x, y, errory_low, errory_high;

		vector<double> v_eyl, v_eyh, v_exl, v_exh, v_x, v_y; //Contain the systematic errors (used to create the TGraphError)

        float bin_width = -1; //Get bin width of histograms for current variable

        string EFTpoint_name1 = "", EFTpoint_name2 = "";

        //-- All histos are for given lumiYears and sub-channels --> Need to sum them all for plots
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {


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

                // cout<<"sample_list[isample] "<<sample_list[isample]<<endl;
                // cout<<"index_MC_sample "<<index_MC_sample<<endl;

				//-- In Combine, some individual contributions are merged as "Rares"/"EWK", etc.
				//-- If using Combine file, change the names of the samples we look for, and look only once for histogram of each "group"
				TString samplename = sample_list[isample];
				if(isample > 0 && sample_groups[isample] == sample_groups[isample-1]) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //if same group as previous sample, skip it
                else if(make_SMvsEFT_templates_plots && (sample_groups[isample] == "tZq" || sample_groups[isample] == "ttZ" || sample_groups[isample] == "tWZ")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //SM vs EFT --> use private signal samples
                else {samplename = sample_groups[isample];}

				//-- Protections, special cases
				if(sample_list[isample] == "DATA") {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;}
                else if(!make_SMvsEFT_templates_plots && sample_list[isample].Contains("PrivMC")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //SM configuration --> only stack central samples (not private samples)
                else if(make_SMvsEFT_templates_plots && (sample_list[isample] == "tZq" || sample_list[isample] == "ttZ" || sample_list[isample] == "tWZ")) {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //EFT configuration --> only stack private samples (at SM point), not central samples
                else if(sample_list[isample] == "NPL_DATA")  {samplename = "NPL";} //Instead of 'NPL_DATA' and 'NPL_MC', we only want to read the merged histo 'NPL'
                else if(sample_list[isample] == "NPL_MC")  {v_isSkippedSample[isample] = true; nof_skipped_samples++; continue;} //NPL_MC gets substracted from NPL histograms and deleted --> Ignore this vector element //Remove ?

                //-- Add sample name to list (used for legend) //NB: add even if histo was not found and skipped, because expect that it will be found for some other year/channel/... But if not found at all, legend will be wrong
                if(iyear==0 && ivar==0 && samplename != "DATA")
                {
                    if(v_vector_MC_histo[ivar].size() <=  index_MC_sample) {MC_samples_legend.push_back(samplename);}
                }
                if(v_isSkippedSample[isample] == true) {continue;} //Skip this sample

				// cout<<endl<<UNDL(FBLU("-- Sample : "<<sample_list[isample]<<" : "))<<endl;

				h_tmp = NULL;
                Get_Variable_Range(total_var_list[ivar], nIndivBins, xmin_tmp, xmax_tmp);

                //HARDCODED -- same as in Produce_Templates()
                nIndivBins = 10;
                if(total_var_list[ivar] == "njets") {xmin_tmp = 2; xmax_tmp = 8; nIndivBins = 6;}
                else if(total_var_list[ivar] == "nbjets") {xmin_tmp = 1; xmax_tmp = 4; nIndivBins = 3;}
                else if(total_var_list[ivar] == "jPrimeAbsEta" || total_var_list[ivar] == "jprime_Pt" || total_var_list[ivar] == "dR_tZ" || total_var_list[ivar] == "maxDeepJet") {nIndivBins = 8;}
                // cout<<"total_var_list[ivar] "<<total_var_list[ivar]<<" / nIndivBins = "<<nIndivBins<<" / xmin_tmp = "<<xmin_tmp<<" / xmax_tmp = "<<xmax_tmp<<endl;

                h_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
                v_eyl.resize(nIndivBins); v_eyh.resize(nIndivBins); v_exl.resize(nIndivBins); v_exh.resize(nIndivBins); v_y.resize(nIndivBins); v_x.resize(nIndivBins);
                std::fill(v_y.begin(), v_y.end(), -1); //Init errors positions to <0 (invisible)

                inputfilename = "./outputs/dir_shapes_tmp/shapes_prefit_COMBINED_Datacard_TemplateFit_"+total_var_list[ivar]+"__Run2.root";
                if(!Check_File_Existence(inputfilename)) {cout<<FRED("File "<<inputfilename<<" not found !")<<endl; continue;}
                file_input = TFile::Open(inputfilename, "READ");
                // cout<<"inputfilename "<<inputfilename<<endl;

                for(int ibin=1; ibin<nIndivBins+1; ibin++)
                {
                    TString dir_hist = total_var_list[ivar] + "_" + v_lumiYears[iyear] + "_prefit/";
                    if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(samplename) ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<samplename<<"' not found ! Skip...")<<endl; continue;}
                    // cout<<"dir_hist/samplename "<<dir_hist<<samplename<<endl;

                    h_tmp->SetBinContent(ibin, ((TH1F*) file_input->Get(dir_hist+samplename))->GetBinContent(ibin)); //Get content/error from individual bin
                    h_tmp->SetBinError(ibin, ((TH1F*) file_input->Get(dir_hist+samplename))->GetBinError(ibin));

                    //-- Errors
                    if(v_y[ibin-1] < 0) //Need to fill total error only once (not for each process)
                    {
                        dir_hist = "prefit/"; //Reminder: read per-year histograms to get individual contributions from processes, *but* read total-sum histogram to get full postfit error
                        if(!file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains("TotalProcs") ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<"TotalProcs' not found ! Skip...")<<endl; continue;}
                        float bin_width = (xmax_tmp-xmin_tmp) / nIndivBins;
                        v_x[ibin-1] = xmin_tmp + ((ibin-1)*bin_width) + (bin_width/2.);
                        v_y[ibin-1] = ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinContent(ibin);
                        v_eyl[ibin-1] = ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinError(ibin);
                        v_eyh[ibin-1] = ((TH1F*) file_input->Get(dir_hist+"TotalProcs"))->GetBinError(ibin);
                        v_exl[ibin-1] = bin_width / 2; v_exh[ibin-1] = bin_width / 2;
                        // cout<<"bin "<<ibin<<" / error = "<<v_eyh[ibin-1]<<endl;
                    }
                } //nbins

				//-- Set histo style (use color vector filled in main) //NB: run for all sub-histos... for simplicity
                //---------------------------------------------------
				h_tmp->SetFillStyle(1001);
				if(samplename == "Fakes") {h_tmp->SetFillStyle(3005);}
		        else if(samplename == "QFlip" ) {h_tmp->SetFillStyle(3006);}

				h_tmp->SetFillColor(color_list[isample]);
				h_tmp->SetLineColor(kBlack);
                h_tmp->SetLineColor(color_list[isample]);

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
                //---------------------------------------------------

                //-- Fill vector of MC histograms
                if(v_vector_MC_histo[ivar].size() <=  index_MC_sample) {v_vector_MC_histo[ivar].push_back((TH1F*) h_tmp->Clone());}
                else if(!v_vector_MC_histo[ivar][index_MC_sample] && h_tmp) {v_vector_MC_histo[ivar][index_MC_sample] = (TH1F*) h_tmp->Clone();}
                else {v_vector_MC_histo[ivar][index_MC_sample]->Add((TH1F*) h_tmp->Clone());}
                if(v_vector_MC_histo[ivar][index_MC_sample]) {v_vector_MC_histo[ivar][index_MC_sample]->SetDirectory(0);} //Dis-associate histo from TFile //https://root.cern.ch/root/htmldoc/guides/users-guide/ObjectOwnership.html

				delete h_tmp; h_tmp = NULL; //No crash ? (else only delete if new)
			} //end sample loop


// #####    ##   #####   ##
// #    #  #  #    #    #  #
// #    # #    #   #   #    #
// #    # ######   #   ######
// #    # #    #   #   #    #
// #####  #    #   #   #    #

            TString dataname = "data_obs";
            hdata_tmp = new TH1F("", "", nIndivBins, xmin_tmp, xmax_tmp);
            hdata_tmp->SetDirectory(0); //Dis-associate from TFile

            for(int ibin=1; ibin<nIndivBins+1; ibin++)
            {
                TString dir_hist = total_var_list[ivar] + "_" + v_lumiYears[iyear] + "_prefit/";
                // cout<<"dir_hist/dataname "<<dir_hist<<dataname<<endl;
                if(!file_input->GetDirectory(dir_hist) || !file_input->GetDirectory(dir_hist)->GetListOfKeys()->Contains(dataname) ) {cout<<FRED("Directory '"<<dir_hist<<"' or histogram '"<<dir_hist<<dataname<<"' not found ! Skip...")<<endl; continue;}

                int bin_to_read = 1; //Default (for split-per-bin postfit plots, always need to read bin1); but for prefit NN_SM templates, need to read current ibin (because we are reading full templates directly, not 1 bin at a time)
                float bin_content = ((TH1F*) file_input->Get(dir_hist+dataname))->GetBinContent(ibin);
                float bin_error = ((TH1F*) file_input->Get(dir_hist+dataname))->GetBinError(ibin);

                if(isnan(bin_content)) {cout<<FRED("ERROR: NaN data (dir_hist="<<dir_hist<<") !")<<endl; continue;} //Can happen when input data is 0 (?)
                hdata_tmp->SetBinContent(ibin, bin_content); //Get content/error from individual bin
                hdata_tmp->SetBinError(ibin, bin_error);

                // cout<<"dir_hist "<<dir_hist<<endl;
                // cout<<"hdata_tmp->GetBinContent(ibin) "<<hdata_tmp->GetBinContent(ibin)<<endl;
            } //nbins

            file_input->Close();

            if(v_hdata[ivar] == NULL) {v_hdata[ivar] = (TH1F*) hdata_tmp->Clone();}
            else {v_hdata[ivar]->Add((TH1F*) hdata_tmp->Clone());}
            v_hdata[ivar]->SetMarkerStyle(20);

            v_hdata[ivar]->SetDirectory(0); //Dis-associate from TFile
            delete hdata_tmp; hdata_tmp = NULL; //No crash ? (else only delete if new)
        } //Years loop


// ##### #    #  ####  #####   ##    ####  #    #
//   #   #    # #        #    #  #  #    # #   #
//   #   ######  ####    #   #    # #      ####
//   #   #    #      #   #   ###### #      #  #
//   #   #    # #    #   #   #    # #    # #   #
//   #   #    #  ####    #   #    #  ####  #    #

    	//-- Add legend entries -- iterate backwards, so that last histo stacked is on top of legend
        v_stack[ivar] = new THStack;

		for(int i=v_vector_MC_histo[ivar].size()-1; i>=0; i--)
		{
			if(!v_vector_MC_histo[ivar][i]) {continue;} //Some templates may be null
			v_stack[ivar]->Add(v_vector_MC_histo[ivar][i]);

            if(v_histo_total_MC[ivar] == NULL) {v_histo_total_MC[ivar] = (TH1F*) v_vector_MC_histo[ivar][i]->Clone();}
            else {v_histo_total_MC[ivar]->Add(v_vector_MC_histo[ivar][i]);}

			// cout<<"Stacking sample "<<MC_samples_legend[i]<<" / integral "<<v_vector_MC_histo[ivar][i]->Integral()<<endl;
            // cout<<"stack bin 1 content = "<<((TH1*) v_stack[ivar]->GetStack()->Last())->GetBinContent(1)<<endl;
		}


// ###### #####  #####   ####  #####   ####      ####  #####   ##    ####  #    #
// #      #    # #    # #    # #    # #         #        #    #  #  #    # #   #
// #####  #    # #    # #    # #    #  ####      ####    #   #    # #      ####
// #      #####  #####  #    # #####       #         #   #   ###### #      #  #
// #      #   #  #   #  #    # #   #  #    #    #    #   #   #    # #    # #   #
// ###### #    # #    #  ####  #    #  ####      ####    #   #    #  ####  #    #

        //-- Use pointers to vectors : need to give the adress of first element (all other elements can then be accessed iteratively)
        double* eyl = &v_eyl[0];
        double* eyh = &v_eyh[0];
        double* exl = &v_exl[0];
        double* exh = &v_exh[0];
        double* xx = &v_x[0];
        double* yy = &v_y[0];

        v_gr_error[ivar] = new TGraphAsymmErrors(nIndivBins,xx,yy,exl,exh,eyl,eyh);
        v_gr_error[ivar]->SetFillStyle(3254); //3002 //3004
        v_gr_error[ivar]->SetFillColor(kBlack);

        //-- Debug printouts
        // for(int ibin=1; ibin<nIndivBins+1; ibin++)
        // {
        //     cout<<"-- ibin "<<ibin<<endl;
        //     cout<<"v_eyh[ibin] "<<v_eyh[ibin]<<" / v_exl[ibin] "<<v_exl[ibin]<<" / v_exh[ibin] "<<v_exh[ibin]<<" / v_x[ibin] "<<v_x[ibin]<<" / v_y[ibin] "<<v_y[ibin]<<endl;
        //     cout<<"v_hdata[ivar]->GetBinContent(ibin) "<<v_hdata[ivar]->GetBinContent(ibin)<<endl;
        //     cout<<"v_histo_total_MC[ivar]->GetBinContent(ibin) "<<v_histo_total_MC[ivar]->GetBinContent(ibin)<<endl;
        // }

        for(int ibin=1; ibin<nIndivBins+1; ibin++)
        {
            if(isnan(v_hdata[ivar]->GetBinContent(ibin))) {cout<<FRED("ERROR: v_hdata["<<ivar<<"]->GetBinContent("<<ibin<<") is nan ! May cause plotting bugs !")<<endl;}
        }


// #####  #####    ##   #    #
// #    # #    #  #  #  #    #
// #    # #    # #    # #    #
// #    # #####  ###### # ## #
// #    # #   #  #    # ##  ##
// #####  #    # #    # #    #

        c1->cd(ivar+1);

        //Draw stack
        v_stack[ivar]->Draw("hist");

        v_hdata[ivar]->Draw("e0p same");

        v_gr_error[ivar]->Draw("e2 same"); //Superimposes the uncertainties on stack


// #   # #    #   ##   #    #
//  # #  ##  ##  #  #   #  #
//   #   # ## # #    #   ##
//   #   #    # ######   ##
//   #   #    # #    #  #  #
//   #   #    # #    # #    #

        //-- Set minimum
        v_stack[ivar]->SetMinimum(1.5);
        // v_stack[ivar]->SetMinimum(2.5);

        if(total_var_list[ivar].Contains("ctw_SRtZq")) {v_stack[ivar]->SetMinimum(1.);} //Very low prediction in last bin

        //-- Set Yaxis maximum
        double ymax = 0;
        ymax = v_hdata[ivar]->GetMaximum(); //Data ymax
        if(ymax < v_stack[ivar]->GetMaximum()) {ymax = v_stack[ivar]->GetMaximum();} //MC ymax
        if(total_var_list[ivar] == "lAsymmetry") {ymax*= 1.1;}
        else {ymax*= 1.3;}
        // ymax*= 1.5; //Default
        v_stack[ivar]->SetMaximum(ymax);
        c1->Modified();


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
        v_tpad_ratio[ivar] = new TPad("pad_ratio", "pad_ratio", 0.0, 0.0, 1.0, 1.0);
        v_tpad_ratio[ivar]->SetTopMargin(0.70);
        // v_tpad_ratio[ivar]->SetBottomMargin(0.13);
        v_tpad_ratio[ivar]->SetFillColor(0);
    	v_tpad_ratio[ivar]->SetFillStyle(0);
    	v_tpad_ratio[ivar]->SetGridy(1);
    	v_tpad_ratio[ivar]->Draw();
    	v_tpad_ratio[ivar]->cd(0);

        if(rightmargin>0) {v_tpad_ratio[ivar]->SetRightMargin(rightmargin);}

        v_histo_ratio_data[ivar] = (TH1F*) v_hdata[ivar]->Clone();

    	if(!show_pulls_ratio) //Compute ratios (with error bars)
    	{
    		//To get correct error bars in ratio plot, must only account for errors from data, not MC ! (MC error shown as separate band)
    		for(int ibin=1; ibin<v_histo_total_MC[ivar]->GetNbinsX()+1; ibin++)
    		{
    			v_histo_total_MC[ivar]->SetBinError(ibin, 0.);
    		}

    		v_histo_ratio_data[ivar]->Divide(v_histo_total_MC[ivar]);
    	}
     	else //-- Compute pulls (no error bars)
    	{
    		for(int ibin=1; ibin<v_histo_ratio_data[ivar]->GetNbinsX()+1; ibin++)
    		{
    			//Add error on signal strength (since we rescale signal manually)
    			// double bin_error_mu = v_vector_MC_histo[ivar].at(index_tZq_sample)->GetBinError(ibin) * sig_strength_err;
    			// cout<<"bin_error_mu = "<<bin_error_mu<<endl;

    			double bin_error_mu = 0; //No sig strength uncert. for prefit ! //-- postfit -> ?

    			//Quadratic sum of systs, stat error, and sig strength error
    			double bin_error = pow(pow(v_histo_total_MC[ivar]->GetBinError(ibin), 2) + pow(v_histo_ratio_data[ivar]->GetBinError(ibin), 2) + pow(bin_error_mu, 2), 0.5);

    			if(!v_histo_total_MC[ivar]->GetBinError(ibin)) {v_histo_ratio_data[ivar]->SetBinContent(ibin,-99);} //Don't draw null markers
    			else{v_histo_ratio_data[ivar]->SetBinContent(ibin, (v_histo_ratio_data[ivar]->GetBinContent(ibin) - v_histo_total_MC[ivar]->GetBinContent(ibin)) / bin_error );}
    		}

    		//-- Don't draw null data
    		for(int ibin=1; ibin<v_histo_ratio_data[ivar]->GetNbinsX()+1; ibin++)
    		{
                if(std::isnan(v_histo_ratio_data[ivar]->GetBinContent(ibin)) || std::isinf(v_histo_ratio_data[ivar]->GetBinContent(ibin)) || v_histo_ratio_data[ivar]->GetBinContent(ibin) == 0) {v_histo_ratio_data[ivar]->SetBinContent(ibin, -99);}
    		}
    	}

        //-- SET X_AXIS TITLES
        //-- NB: function Get_Variable_Name() (uses 'l' instead of '\ell') --> Hardcode var names here if needed
        TString xtitle = Get_Variable_Name(total_var_list[ivar]);
        if(total_var_list[ivar] == "lAsymmetry") {xtitle = "Lepton asymmetry";}
        else if(total_var_list[ivar] == "recoZ_dPhill") {xtitle = "\\Delta\\phi(\\ell_{1}^{\\text{Z}},\\ell_{2}^{\\text{Z}})";}
        else if(total_var_list[ivar] == "maxDeepJet") {xtitle = "Maximum DeepJet discriminant";}
        // else if(total_var_list[ivar] == "maxDeepJet") {xtitle = "Max. DeepJet discriminant";}
        // if(total_var_list[ivar] == "lAsymmetry") {xtitle = "\\text{q}_{\\ell} \\centerdot \\left|\\eta(\\ell)\\right|";}
        // else if(total_var_list[ivar] == "recoZ_dPhill") {xtitle = "\\Delta\\varphi(\\ell_{1}^{\\text{Z}},\\ell_{2}^{\\text{Z}})";}

        v_histo_ratio_data[ivar]->GetXaxis()->SetTitle(xtitle);

        //-- AXES OPTIONS
    	if(show_pulls_ratio) {v_histo_ratio_data[ivar]->GetYaxis()->SetTitle("Pulls");}
        else {v_histo_ratio_data[ivar]->GetYaxis()->SetTitle("#frac{Data}{Pred.}");}
        // else {v_histo_ratio_data[ivar]->GetYaxis()->SetTitle("Data/MC");}
        v_histo_ratio_data[ivar]->GetYaxis()->CenterTitle(true);
        v_histo_ratio_data[ivar]->GetYaxis()->SetLabelFont(42);
        v_histo_ratio_data[ivar]->GetYaxis()->SetLabelSize(0.05);
    	v_histo_ratio_data[ivar]->GetXaxis()->SetTitleFont(42);
    	v_histo_ratio_data[ivar]->GetYaxis()->SetTitleFont(42);
        v_histo_ratio_data[ivar]->GetYaxis()->SetNdivisions(303); //grid draw on primary tick marks only
        v_histo_ratio_data[ivar]->GetXaxis()->SetNdivisions(505); //'-' to force Ndivisions //NB: based on previous CMS publications, it does not matter whether ticks are placed at bin boundaries
        v_histo_ratio_data[ivar]->GetYaxis()->SetTitleSize(0.05);
        v_histo_ratio_data[ivar]->GetYaxis()->SetTitleOffset(1.40); //1.35
        v_histo_ratio_data[ivar]->GetXaxis()->SetTickLength(0.04);
        v_histo_ratio_data[ivar]->GetYaxis()->SetTickLength(0.15);
    	v_histo_ratio_data[ivar]->SetMarkerStyle(20);
    	v_histo_ratio_data[ivar]->SetMarkerSize(1.2);

        if(total_var_list[ivar] == "nbjets" || total_var_list[ivar] == "njets")
        {
            v_histo_ratio_data[ivar]->GetXaxis()->CenterLabels(true); //Special case: want to have label (nof jets) at center of bin
            v_histo_ratio_data[ivar]->GetXaxis()->SetNdivisions(-v_histo_ratio_data[ivar]->GetNbinsX()); //'-' to force Ndivisions //NB: must use same pattern as bottom TPad (if present) !
        }

        //-- OPEN-HEADED TRIANGLES
        //-- If a point is outside the y-range of the ratio pad defined by SetMaximum/SetMinimum(), it disappears with its error
        //-- Trick: fill 2 histos with points either above/below y-range, to plot some markers indicating missing points (cleaner)
        //NB: only for ratio plot, not pulls
        float ratiopadmin = 0.4, ratiopadmax = 1.6; //Define ymin/ymax for ratio plot
        TH1F* h_pointsAboveY = (TH1F*) v_histo_ratio_data[ivar]->Clone();
        h_pointsAboveY->SetMarkerStyle(26); //Open triangle pointing up
        h_pointsAboveY->SetMarkerSize(1.5);
        TH1F* h_pointsBelowY = (TH1F*) v_histo_ratio_data[ivar]->Clone();
        h_pointsBelowY->SetMarkerStyle(32); //Open triangle pointing down
        h_pointsBelowY->SetMarkerSize(1.5);
        if(show_pulls_ratio)
    	{
    		v_histo_ratio_data[ivar]->SetMinimum(-2.99);
    		v_histo_ratio_data[ivar]->SetMaximum(2.99);
    	}
    	else
    	{
            //-- Default
            v_histo_ratio_data[ivar]->SetMinimum(ratiopadmin); //NB: removes error bars if data point is below ymin...?
            v_histo_ratio_data[ivar]->SetMaximum(ratiopadmax);

            //-- Fill histos with points outside yrange
            for(int ibin=1; ibin<v_histo_ratio_data[ivar]->GetNbinsX()+1; ibin++)
            {
                //-- Default: make point invisible
                h_pointsAboveY->SetBinContent(ibin, -999);
                h_pointsBelowY->SetBinContent(ibin, -999);

                if(v_histo_ratio_data[ivar]->GetBinContent(ibin) > ratiopadmax && v_hdata[ivar]->GetBinContent(ibin) >= 1)
                {
                    //Adjust error
                    float initial_y = v_histo_ratio_data[ivar]->GetBinContent(ibin);
                    float initial_err = v_histo_ratio_data[ivar]->GetBinError(ibin);
                    float new_err = initial_err - (initial_y-ratiopadmax);
                    if(new_err<0) {new_err=0.;}

                    h_pointsAboveY->SetBinContent(ibin, ratiopadmax-0.05); //Add some padding
                    h_pointsAboveY->SetBinError(ibin, new_err);
                }
                else if(v_histo_ratio_data[ivar]->GetBinContent(ibin) < ratiopadmin && v_hdata[ivar]->GetBinContent(ibin) >= 1)
                {
                    //Adjust error
                    float initial_y = v_histo_ratio_data[ivar]->GetBinContent(ibin);
                    float initial_err = v_histo_ratio_data[ivar]->GetBinError(ibin);
                    float new_err = initial_err - (ratiopadmin-initial_y);
                    if(new_err<0) {new_err=0.;}

                    h_pointsBelowY->SetBinContent(ibin, ratiopadmin+(ratiopadmin/10.)); //Add some padding
                    h_pointsBelowY->SetBinError(ibin, new_err);
                }
            }
    	}

    	if(show_pulls_ratio) {v_histo_ratio_data[ivar]->Draw("HIST P");} //Draw ratio points
        else
        {
            v_histo_ratio_data[ivar]->Draw("E1 X0 P"); //Draw ratio points ; E1 : perpendicular lines at end ; X0 : suppress x errors

            h_pointsAboveY->Draw("E1 X0 P same");
            h_pointsBelowY->Draw("E1 X0 P same");
        }


// ###### #####  #####   ####  #####   ####     #####    ##   ##### #  ####
// #      #    # #    # #    # #    # #         #    #  #  #    #   # #    #
// #####  #    # #    # #    # #    #  ####     #    # #    #   #   # #    #
// #      #####  #####  #    # #####       #    #####  ######   #   # #    #
// #      #   #  #   #  #    # #   #  #    #    #   #  #    #   #   # #    #
// ###### #    # #    #  ####  #    #  ####     #    # #    #   #   #  ####


		//Copy previous TGraphAsymmErrors, then modify it -> error TGraph for ratio plot
		TGraphAsymmErrors *thegraph_tmp = NULL;
		double *theErrorX_h;
		double *theErrorY_h;
		double *theErrorX_l;
		double *theErrorY_l;
		double *theY;
		double *theX;

		thegraph_tmp = (TGraphAsymmErrors*) v_gr_error[ivar]->Clone();
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

		v_gr_ratio_error[ivar] = new TGraphAsymmErrors(thegraph_tmp->GetN(), theX , theY ,  theErrorX_l, theErrorX_h, theErrorY_l, theErrorY_h);
        v_gr_ratio_error[ivar]->SetFillStyle(3254); //3002 //3004
        v_gr_ratio_error[ivar]->SetFillColor(kBlack); //kBlue+2 //kCyan
        // v_gr_ratio_error[ivar]->SetFillColorAlpha(kBlack, 0.1); //kBlue+2 //kCyan //Only works for pnf ?
        // gStyle->SetHatchesSpacing(0.1);

		if(!show_pulls_ratio) {v_gr_ratio_error[ivar]->Draw("e2 same");} //Draw error bands in ratio plot


//  ####   ####   ####  #    # ###### ##### #  ####   ####
// #    # #    # #      ##  ## #        #   # #    # #
// #      #    #  ####  # ## # #####    #   # #       ####
// #      #    #      # #    # #        #   # #           #
// #    # #    # #    # #    # #        #   # #    # #    #
//  ####   ####   ####  #    # ######   #   #  ####   ####

    	//-- Draw ratio y-lines manually
        v_tpad_ratio[ivar]->cd();
    	v_hlines1[ivar] = new TH1F("","",this->nbins, xmin_tmp, xmax_tmp);
    	v_hlines2[ivar] = new TH1F("","",this->nbins, xmin_tmp, xmax_tmp);
    	for(int ibin=1; ibin<this->nbins +1; ibin++)
    	{
    		if(show_pulls_ratio)
    		{
    			v_hlines1[ivar]->SetBinContent(ibin, -1);
    			v_hlines2[ivar]->SetBinContent(ibin, 1);
    		}
    		else
    		{
                v_hlines1[ivar]->SetBinContent(ibin, 0.75);
                v_hlines2[ivar]->SetBinContent(ibin, 1.25);
    		}
    	}
    	v_hlines1[ivar]->SetLineStyle(6); v_hlines2[ivar]->SetLineStyle(6);
    	// v_hlines1[ivar]->Draw("hist same"); v_hlines2[ivar]->Draw("hist same"); //Removed extra lines

        TString Y_label = "Events / bin";
        Y_label = "Events / " + Convert_Number_To_TString( (xmax_tmp - xmin_tmp) / nIndivBins, 3); //Automatically get the Y label depending on binning
        Y_label+= Get_Unit_Variable(total_var_list[ivar]);
        if(total_var_list[ivar] == "nbjets") {Y_label = "Events / 1 unit";}

        if(v_stack[ivar]) //Must be drawn first
    	{
    		v_stack[ivar]->GetYaxis()->SetLabelFont(42);
    		v_stack[ivar]->GetYaxis()->SetTitleFont(42);
    		v_stack[ivar]->GetYaxis()->SetTitleSize(0.06);
            v_stack[ivar]->GetYaxis()->SetTickLength(0.04);
    		v_stack[ivar]->GetYaxis()->SetLabelSize(0.048);
    		v_stack[ivar]->GetYaxis()->SetNdivisions(506);
            v_stack[ivar]->GetYaxis()->SetTitleOffset(1.40); //1.30
    		v_stack[ivar]->GetYaxis()->SetTitle(Y_label);
            v_stack[ivar]->GetXaxis()->SetLabelSize(0.);
            v_stack[ivar]->GetXaxis()->SetTickLength(0.);
    	}

    	//----------------
    	// CAPTIONS //
    	//----------------
    	// -- using https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines
        // -- About fonts: https://root.cern.ch/doc/master/classTAttText.html#T5

    	float l = pad->GetLeftMargin();
    	float t = pad->GetTopMargin();

    	TString cmsText = "CMS";
    	TLatex latex;
    	latex.SetNDC();
    	latex.SetTextColor(kBlack);
        latex.SetTextFont(62); //Changed
    	latex.SetTextAlign(11);
    	latex.SetTextSize(0.06);
        // latex.DrawLatex(l + 0.04, 0.87, cmsText);
        if(use_paperStyle) {latex.DrawLatex(l + 0.04, 0.87, cmsText);} //CMS guideline: within frame
        else {latex.DrawLatex(l + 0.01, 0.94, cmsText);} //Default: outside frame

    	float lumi = lumiValue;
    	TString lumi_ts = Convert_Number_To_TString(lumi);
    	lumi_ts += " fb^{-1} (13 TeV)";
    	latex.SetTextFont(42);
    	latex.SetTextAlign(31);
    	latex.SetTextSize(0.04);
        latex.DrawLatex(0.96, 0.94,lumi_ts);

        TString extraText = "Preliminary";
        if(!use_paperStyle) //Default is without (for paper)
        {
            latex.SetTextFont(52);
            latex.SetTextSize(0.05);
            if(total_var_list[ivar] == "NN_SM_SRother") {latex.DrawLatex(l + 0.40, 0.94, extraText);} //Hardcoded: due to different pad size (?), need to add a bit of extra space for this one
            else {latex.DrawLatex(l + 0.38, 0.94, extraText);}
        }

        TString info_data = Get_Region_Label(region, total_var_list[ivar]);
        TLatex text2;
        text2.SetNDC();
        text2.SetTextAlign(13);
        text2.SetTextSize(0.04);
        // text2.SetTextSize(0.045);
        text2.SetTextFont(42);
        if(info_data != "")
        {
            float left = l+0.04;
            if(total_var_list[ivar].Contains("other")) {left = l+0.03;} //Need more space
            if(use_paperStyle) {text2.DrawLatex(left,0.85,info_data);} //0.83
            else {text2.DrawLatex(left,0.90,info_data);} //Default: outside frame //0.90
        }

        pad = NULL;
    } //Var loop
//--------------------------------------------


// ####### #
//    #    #       ######  ####  ###### #    # #####
//    #    #       #      #    # #      ##   # #    #
//    #    #       #####  #      #####  # #  # #    #
//    #    #       #      #  ### #      #  # # #    #
//    #    #       #      #    # #      #   ## #    #
//    #    ####### ######  ####  ###### #    # #####

    int n_columns = ceil(nSampleGroups/2.) > 6 ? 6 : ceil(nSampleGroups/2.); //ceil = upper int
    float x_left = 0.94-n_columns*0.12; //Each column allocated same x-space //0.12 needed for most crowded plots
    if(x_left < 0.4) {x_left = 0.4;} //Leave some space for region label

    TLegend* qw = new TLegend(0.20,0.58,0.87,0.77);
    // TLegend* qw = new TLegend(x_left-0.05,0.80,0.94,0.92); //Default
    qw->SetTextSize(0.04);
    qw->SetNColumns(n_columns);
    qw->SetBorderSize(0);
    qw->SetFillStyle(0); //transparent
    qw->SetTextAlign(12); //align = 10*HorizontalAlign + VerticalAlign //Horiz: 1=left adjusted, 2=centered, 3=right adjusted //Vert: 1=bottom adjusted, 2=centered, 3=top adjusted
    // cout<<"x_left "<<x_left<<endl;
    // cout<<"ceil(nSampleGroups/2.) "<<ceil(nSampleGroups/2.)<<endl;

    //-- Dummy object, only used to display uncertainty band also in legend
    TH1F* h_uncert = new TH1F("h_uncert", "h_uncert", 1, 0, 1);
    h_uncert->SetFillStyle(3254); //3002 //3004
    h_uncert->SetFillColor(kBlack);
    h_uncert->SetLineWidth(0.);
    qw->AddEntry(h_uncert, "Unc.", "F");

	//-- Data on top of legend
    qw->AddEntry(v_hdata[0], "Data" , "ep");

	for(int i=0; i<v_vector_MC_histo[0].size(); i++)
	{
        if(MC_samples_legend[i].Contains("tZq")) {qw->AddEntry(v_vector_MC_histo[0][i], "tZq", "f");}
        else if(MC_samples_legend[i].EndsWith("ttZ") ) {qw->AddEntry(v_vector_MC_histo[0][i], "t#bar{t}Z", "f");}
        else if(MC_samples_legend[i].EndsWith("tWZ") ) {qw->AddEntry(v_vector_MC_histo[0][i], "tWZ", "f");}
        else if(MC_samples_legend[i] == "ttW" || MC_samples_legend[i] == "tX") {qw->AddEntry(v_vector_MC_histo[0][i], "t(#bar{t})X", "f");}
        else if(MC_samples_legend[i] == "WZ") {qw->AddEntry(v_vector_MC_histo[0][i], "WZ", "f");}
        else if(MC_samples_legend[i] == "WWZ" || MC_samples_legend[i] == "VVV") {qw->AddEntry(v_vector_MC_histo[0][i], "VV(V)", "f");}
        else if(MC_samples_legend[i] == "TTGamma_Dilep" || MC_samples_legend[i] == "XG") {qw->AddEntry(v_vector_MC_histo[0][i], "X#gamma", "f");}
        else if(MC_samples_legend[i] == "TTbar_DiLep" || MC_samples_legend[i] == "NPL" || MC_samples_legend[i] == "NPL_DATA") {qw->AddEntry(v_vector_MC_histo[0][i], "NPL", "f");}
	}

    for(int ivar=0; ivar<nvar; ivar++)
    {
        if(ivar!=1) {continue;} //Only display legend for maxDeepJet (upper middle)

        c1->cd(ivar+1);
        qw->Draw("same");
    }


// #    # #####  # ##### ######     ####  #    # ##### #####  #    # #####
// #    # #    # #   #   #         #    # #    #   #   #    # #    #   #
// #    # #    # #   #   #####     #    # #    #   #   #    # #    #   #
// # ## # #####  #   #   #         #    # #    #   #   #####  #    #   #
// ##  ## #   #  #   #   #         #    # #    #   #   #      #    #   #
// #    # #    # #   #   ######     ####   ####    #   #       ####    #

    TString outdir = "plots/paperPlots/";
    mkdir(outdir.Data(), 0777);
    TString output_plot_name = outdir + "controlPlots_prefit";
    if(!use_paperStyle) {output_plot_name+= "_prelim";}
    c1->SaveAs(output_plot_name + ".png");
    c1->SaveAs(output_plot_name + ".eps");
    c1->SaveAs(output_plot_name + ".pdf");

    delete c1; c1 = NULL;
    delete qw; qw = NULL;
    if(h_uncert) {delete h_uncert; h_uncert = NULL;}

    for(int ivar=0; ivar<nvar; ivar++)
    {
        if(v_gr_error[ivar]) {delete v_gr_error[ivar]; v_gr_error[ivar] = NULL;}
        if(v_stack[ivar]) {delete v_stack[ivar]; v_stack[ivar] = NULL;}
        if(v_gr_ratio_error[ivar]) {delete v_gr_ratio_error[ivar]; v_gr_ratio_error[ivar] = NULL;}
        if(v_hlines1[ivar]) {delete v_hlines1[ivar]; v_hlines1[ivar] = NULL;}
        if(v_hlines2[ivar]) {delete v_hlines2[ivar]; v_hlines2[ivar] = NULL;}
    }

    return;
}
















//--------------------------------------------
//  #######  ######## ##     ## ######## ########   ######
// ##     ##    ##    ##     ## ##       ##     ## ##    ##
// ##     ##    ##    ##     ## ##       ##     ## ##
// ##     ##    ##    ######### ######   ########   ######
// ##     ##    ##    ##     ## ##       ##   ##         ##
// ##     ##    ##    ##     ## ##       ##    ##  ##    ##
//  #######     ##    ##     ## ######## ##     ##  ######
//--------------------------------------------

/**
 * Create output file containing the MVA scores for different NNs for a given process (e.g. to compute correlations between different NN scores)
 */
void TopEFT_analysis::Dump_Scores_allNNs()
{
//--------------------------------------------
    TString process = "PrivMC_tZq";
    // TString process = "PrivMC_ttZ";
//--------------------------------------------

    cout<<endl<<BYEL("                          ")<<endl<<endl;
    cout<<FYEL("--- TESTING ---")<<endl;
    cout<<endl<<BYEL("                          ")<<endl<<endl;

    vector<TString> v_templates;
    v_templates.push_back("NN_SM");
    v_templates.push_back("NN_ctz");
    v_templates.push_back("NN_ctw");
    if(process == "PrivMC_tZq") {v_templates.push_back("NN_cpq3");}
    v_templates.push_back("NN_5D");

    TString outname = "./NN_scores_"+process+".root";
    TFile* f = new TFile(outname, "RECREATE");
    TTree* t = new TTree("NN_scores", "NN_scores");

    for(int itemplate=0; itemplate<v_templates.size(); itemplate++)
    {
        cout<<"-- v_templates[itemplate] "<<v_templates[itemplate]<<endl;

        t->ResetBranchAddresses();

        vector<TString> var_list_NN, v_NN_nodeLabels; int NN_iMaxNode = -1; TString NN_strategy = "", NN_inputLayerName = "", NN_outputLayerName = ""; int NN_nNodes = -1; std::vector<float> minmax_bounds;

        TString NNinfo_input_path = Get_MVAFile_InputPath(v_templates[itemplate], "tZq", "Run2", false, true, true, 2, false, this->use_NN_cpq3_SRttZ);
        TString MVA_input_path = Get_MVAFile_InputPath(v_templates[itemplate], "tZq", "Run2", false, true, false, 2, false, this->use_NN_cpq3_SRttZ);
        TFModel* clfy_tmp = NULL;

        //-- Load neural network model
        if(Extract_Values_From_NNInfoFile(NNinfo_input_path, var_list_NN, v_NN_nodeLabels, NN_inputLayerName, NN_outputLayerName, NN_iMaxNode, NN_nNodes, minmax_bounds))
        {
            // cout<<"clfy_tmp "<<clfy_tmp<<endl;
            // cout<<"MVA_input_path "<<MVA_input_path<<endl;
            // cout<<"var_list_NN.size() "<<var_list_NN.size()<<endl;
            // cout<<"NN_inputLayerName "<<NN_inputLayerName<<endl;
            // cout<<"NN_outputLayerName "<<NN_outputLayerName<<endl;
            // cout<<"NN_nNodes "<<NN_nNodes<<endl;
            clfy_tmp = new TFModel(MVA_input_path.Data(), var_list_NN.size(), NN_inputLayerName.Data(), NN_nNodes, NN_outputLayerName.Data());
        }
        else {cout<<"Missing NN information ! "<<endl; continue;} //Error: missing NN infos

        //-- Create an input tensor
        long long int n_inputs = var_list_NN.size() > 0? var_list_NN.size():1;
        tensorflow::Tensor input(tensorflow::DT_FLOAT, { 1, n_inputs }); // single batch of dimension 10
        std::vector<tensorflow::Tensor> outputs; //Store outputs
        std::vector<float> clfy_outputs(NN_nNodes);

        if(v_templates[itemplate] == "NN_SM")
        {
            t->Branch("NN_SM_tZq", &clfy_outputs[0], "NN_SM_tZq/F");
            t->Branch("NN_SM_ttZ", &clfy_outputs[1], "NN_SM_ttZ/F");
            t->Branch("NN_SM_other", &clfy_outputs[2], "NN_SM_other/F");

            // cout<<"&clfy_outputs[0] "<<&clfy_outputs[0]<<endl;
            // cout<<"&clfy_outputs[1] "<<&clfy_outputs[1]<<endl;
        }
        else
        {
            TString branch_name = v_templates[itemplate];
            // if(v_templates[itemplate] == "NN_ctz")

            t->Branch(branch_name, &clfy_outputs[0], (branch_name+"/F").Data());
        }

        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
            cout<<"-- v_lumiYears[iyear] "<<v_lumiYears[iyear]<<endl;

            TString inputfile_name = "./input_ntuples/"+v_lumiYears[iyear]+"/"+process+".root";
            if(!Check_File_Existence(inputfile_name)) {cout<<FRED("File "<<inputfile_name<<" not found !")<<endl; continue;}
            cout<<"-- Opening file : "<<inputfile_name<<endl;
            TFile* f_input = TFile::Open(inputfile_name);
            TTree* t_input = (TTree*) f_input->Get("result");
            vector<float> v_float_inputs(var_list_NN.size());

            t_input->SetBranchStatus("*", 0);
            for(int i=0; i<var_list_NN.size(); i++)
            {
                t_input->SetBranchStatus(var_list_NN[i], 1);
                t_input->SetBranchAddress(var_list_NN[i], &v_float_inputs[i]);
            }

            int nmax = t_input->GetEntries();
            nmax = 1000;
            for(int ientry=0; ientry<nmax; ientry++)
            {
                t_input->GetEntry(ientry);
                if(ientry%5000==0) {cout<<DIM("Entry "<<ientry<<"")<<endl;}

                // for(int i=0; i<var_list_NN.size(); i++)
                // {
                //     cout<<"var "<<i<<" --> "<<v_float_inputs[i]<<endl;
                // }

                // for(int ioutput=0; ioutput<clfy_outputs.size(); ioutput++)
                // {
                //     cout<<"output "<<ioutput<<" --> "<<clfy_tmp->evaluate(v_float_inputs)[ioutput]<<" / "<<&clfy_tmp->evaluate(v_float_inputs)[ioutput]<<endl;
                // }

                for(int inode=0; inode<clfy_outputs.size(); inode++)
                {
                    clfy_outputs[inode] = clfy_tmp->evaluate(v_float_inputs)[inode];
                }

                t->Fill();
            } //Tree entries

            f_input->Close();
        } //Years

        delete clfy_tmp; clfy_tmp = NULL;
    } //Templates

    f->cd();
    t->Write();
    f->Close();

    cout<<endl<<FYEL("==> Created root file: ")<<f->GetName()<<endl;

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
