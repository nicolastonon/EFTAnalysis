//by Nicolas Tonon (DESY)

//--- LIST OF FUNCTIONS (for quick search) :
//--------------------------------------------
// Train_BDT

// Produce_Templates

// Draw_Templates

// SetBranchAddress_SystVariationArray()
// Merge_Templates_ByProcess()

//--------------------------------------------

/** NOTES **

# Include private MC samples in code :
- only select central/private via samplename => update all name checks
- what if want to use both ? give priority in plots for central...?

*/

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
TopEFT_analysis::TopEFT_analysis(vector<TString> thesamplelist, vector<TString> thesamplegroups, vector<TString> thesystlist, vector<TString> thesystTreelist, vector<TString> thechannellist, vector<TString> thevarlist, vector<TString> set_v_cut_name, vector<TString> set_v_cut_def, vector<bool> set_v_cut_IsUsedForBDT, vector<TString> set_v_add_var_names, TString theplotextension, vector<TString> set_lumi_years, bool show_pulls, TString region, TString signal_process, TString classifier_name, bool use_custom_colorPalette)
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

    this->signal_process = signal_process;

	// this->is_blind = is_blind;

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


 // #      #    # #    # #
 // #      #    # ##  ## #
 // #      #    # # ## # #
 // #      #    # #    # #
 // #      #    # #    # #
 // ######  ####  #    # #

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
    // dir_ntuples = "./input_ntuples/"+lumiYear+"/";
	// cout<<"dir_ntuples : "<<dir_ntuples<<endl;

	//-- Get colors
    this->use_custom_colorPalette = use_custom_colorPalette;
	color_list.resize(sample_list.size());
	Get_Samples_Colors(color_list, v_custom_colors, sample_list, 0, use_custom_colorPalette); //Read hard-coded sample colors
    // if(use_custom_colorPalette) {Set_Custom_ColorPalette(v_custom_colors, color_list, sample_groups);} //Replace colors with custom color list

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
		v_cut_char.push_back(0);

		//NOTE : it is a problem if a variable is present in more than 1 list, because it will cause SetBranchAddress conflicts (only the last SetBranchAddress to a branch will work)
		//---> If a variable is present in 2 lists, erase it from other lists !
		for(int ivar=0; ivar<thevarlist.size(); ivar++)
		{
			if(thevarlist[ivar].BeginsWith("is_") || thevarlist[ivar].BeginsWith("passed") )
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
			if(set_v_add_var_names[ivar].BeginsWith("is_") || set_v_add_var_names[ivar].BeginsWith("passed") )
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
    array_LepEff_mu = NULL;
    array_LepEff_el = NULL;

	//Store the "cut name" that will be written as a suffix in the name of each output file
	this->filename_suffix = "";
	TString tmp = "";
	for(int ivar=0; ivar<v_cut_name.size(); ivar++)
	{
		if(v_cut_name[ivar].BeginsWith("is_") || v_cut_name[ivar].BeginsWith("passed") ) {continue;} //No need to appear in filename
		else if(v_cut_name[ivar] == "nLightJets_Fwd40") {this->filename_suffix+= "_fwdCut"; continue;} //No need to appear in filename

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

    //Read the list of input variables directly from .txt file generated at DNN training
    //Also read the mean/variance of each var, in order to rescale the input values correspondingly
    v_inputs_rescaling.resize(0); var_list_DNN.resize(0);
    if(classifier_name == "DNN")
    {
        TString file_var_path = "./weights/DNN/"+lumiName+"/DNN_infos.txt";
        if(!Check_File_Existence(file_var_path) ) {cout<<BOLD(FRED("Error ! DNN infos not found ("<<file_var_path<<")"))<<endl; return;}
        else{cout<<"Reading list of DNN input variables from : "<<file_var_path<<endl;}
        ifstream file_in(file_var_path);
        string line;
        while(!file_in.eof( ))
        {
            getline(file_in, line);
            // TString ts_line(line);
            stringstream ss(line);
            TString varname; float tmp1, tmp2; //Values tmp1 and tmp2 could be the mean and variance, or min and max, etc. depending on the rescaling
            ss >> varname >> tmp1 >> tmp2;
            if(varname != "") //Last line may be empty
            {
                if(tmp1 == -1 && tmp2 == -1) {DNN_inputLayerName = varname;} //Name of input layer
                else if(tmp1 == -2 && tmp2 == -2) {DNN_outputLayerName = varname;} //Name of output layer
                else if(tmp1 == -3 && tmp2 == -3) {nNodes = Convert_TString_To_Number(varname);} //Number of output nodes
                else
                {
                    var_list_DNN.push_back(varname);
                    std::pair <float,float> pair_tmp = std::make_pair(tmp1, tmp2);
                    v_inputs_rescaling.push_back(pair_tmp);
                }
                // cout<<"-->  "<<varname<<endl;
            }
        }
        // cout<<"-->  "<<DNN_inputLayerName<<endl;
        // cout<<"-->  "<<DNN_outputLayerName<<endl;
        // cout<<"-->  "<<nNodes<<endl;
        // cout<<"v_inputs_rescaling.first "<<v_inputs_rescaling[0].first<<" / v_inputs_rescaling.second "<<v_inputs_rescaling[0].second<<endl;
    }

    cout<<endl<<endl<<BLINK(BOLD(FBLU("[Region : "<<region<<"]")))<<endl;
    cout<<endl<<BLINK(BOLD(FBLU("[Luminosity : "<<lumiName<<"]")))<<endl<<endl<<endl;

    usleep(1000000); //Pause for 1s (in microsec)
}

TopEFT_analysis::~TopEFT_analysis()
{
    // cout<<"~TopEFT_analysis "<<endl;

    if(use_custom_colorPalette)
    {
        for(int icol=0; icol<v_custom_colors.size(); icol++)
        {
            if(v_custom_colors[icol] != 0) {delete v_custom_colors[icol];}
        }
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
void TopEFT_analysis::Train_BDT(TString channel, bool write_ranking_info)
{
//--- Options ---------------------------------

	bool use_relative_weights = false; //false <-> use abs(weight), much faster if there are many negative weights

//--------------------------------------------

    cout<<endl<<YELBKG("                          ")<<endl<<endl;
    cout<<FYEL("--- TRAINING ---")<<endl;
    cout<<endl<<YELBKG("                          ")<<endl<<endl;

	if(use_relative_weights) {cout<<"-- Using "<<BOLD(FGRN("*RELATIVE weights*"))<<" --"<<endl<<endl<<endl;}
	else {cout<<"-- Using "<<BOLD(FGRN("*ABSOLUTE weights*"))<<" --"<<endl<<endl<<endl;}

    if(classifier_name != "BDT") {cout<<BOLD(FRED("ERROR ! Can only train BDTs within TMVA for now... Abort ! (classifier_name = "<<classifier_name<<")"))<<endl; return;}

	mkdir("weights", 0777);
    mkdir("weights/BDT", 0777);
    mkdir(("weights/BDT/"+lumiName).Data(), 0777);

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
	// TString cat_tmp = "";
	// cat_tmp = Get_Category_Boolean_Name(nLep_cat, region, analysis_type, "", scheme);

	//Even if ask templates in the SR, need to use training (looser) category for training !
	// if(cat_tmp.Contains("_SR") )
	// {
	// 	int i = cat_tmp.Index("_SR"); //Find index of substring
	// 	cat_tmp.Remove(i); //Remove substring
	// }
    // tmp+= cat_tmp + "==1";

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

	TString weights_dir = "weights";

    //The TMVA::DataLoader object will hold the training and test data, and is later passed to the TMVA::Factory
	TMVA::DataLoader* myDataLoader = new TMVA::DataLoader(weights_dir); //If no TString given in arg, path of weightdir *in TTree* will be : default/weights/...

	//--- Could modify here the name of local dir. storing the BDT weights (default = "weights")
	//By setting it to "", weight files will be stored directly at the path given to myDataLoader
	//Complete path for weight files is : [path_given_toDataloader]/[fWeightFileDir]
	//Apparently, TMVAGui can't handle nested repos in path given to myDataLoader... so split path in 2 here
	TMVA::gConfig().GetIONames().fWeightFileDir = "BDT/"+lumiName;
	// if(classifier_name == "DNN" && DNN_type == "TMVA") {TMVA::gConfig().GetIONames().fWeightFileDir = "DNN/TMVA/"+lumiName;}

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
	for(int i=0; i<v_cut_name.size(); i++)
	{
		// cout<<"Is "<<v_cut_name[i]<<" used ? "<<(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("=="))<<endl;

		// if we ask "var == x", all the selected events will be equal to x, so can't use it as discriminant variable !
		if(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("==")) {myDataLoader->AddVariable(v_cut_name[i].Data(), 'F');}
		// else {myDataLoader->AddSpectator(v_cut_name[i].Data(), v_cut_name[i].Data(), 'F');}
	}
	// for(int i=0; i<v_add_var_names.size(); i++) //Don't add spectator variables anymore ; would allow to get their correlations with training variables, etc.
	// {
		// myDataLoader->AddSpectator(v_add_var_names[i].Data(), v_add_var_names[i].Data(), 'F');
	// }

	double nEvents_sig = 0;
	double nEvents_bkg = 0;


//--------------------------------------------
 //                          #
 //  ####  #  ####          #     #####  #    #  ####
 // #      # #    #        #      #    # #   #  #    #
 //  ####  # #            #       #####  ####   #
 //      # # #  ###      #        #    # #  #   #  ###
 // #    # # #    #     #         #    # #   #  #    #
 //  ####  #  ####     #          #####  #    #  ####
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
            if(signal_process == "tZq")
            {
                if(!samplename_tmp.Contains("tZq") && !samplename_tmp.EndsWith("ttZ") && !samplename_tmp.Contains("ttH") && samplename_tmp.Contains("ttW") && samplename_tmp.Contains("WZ") && samplename_tmp.Contains("ZZ4l") && samplename_tmp.Contains("TTbar_DiLep") ) {continue;}
            }

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

	if(mycuts != mycutb) {cout<<__LINE__<<FRED("PROBLEM : cuts are different for signal and background ! If this is normal, modify code -- Abort")<<endl; return;}

	// If nTraining_Events=nTesting_Events="0", half of the events in the tree are used for training, and the other half for testing
	//NB : when converting nEvents to TString, make sure to ask for sufficient precision !

	float trainingEv_proportion = 0.5;
	// float trainingEv_proportion = 0.7;

	//-- Choose dataset splitting
	TString nTraining_Events_sig = "", nTraining_Events_bkg = "", nTesting_Events_sig = "", nTesting_Events_bkg = "";

    int nmaxEv = 50000; //max nof events for train or test
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
    TString output_file_name = "outputs/" + classifier_name + "_" + signal_process;
	if(channel != "") {output_file_name+= "_" + channel;}
    output_file_name+= "_" + lumiName;
	output_file_name+= this->filename_suffix + ".root";

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

    mkdir("outputs/Rankings", 0777);
    mkdir(("outputs/Rankings/"+lumiName).Data(), 0777); //Dir. containing variable ranking infos

	TString ranking_file_path = "outputs/Rankings/"+lumiName+"/rank_"+classifier_name+"_"+signal_process+".txt";

	if(write_ranking_info) cout<<endl<<endl<<endl<<FBLU("NB : Temporarily redirecting standard output to file '"<<ranking_file_path<<"' in order to save Ranking Info !!")<<endl<<endl<<endl;

	std::ofstream out("ranking_info_tmp.txt"); //Temporary name
    out<<endl;
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	if(write_ranking_info) std::cout.rdbuf(out.rdbuf()); //redirect std::cout to text file --> Ranking info will be saved !

    // Train MVAs using the set of training events
    factory->TrainAllMethods();

	if(write_ranking_info) std::cout.rdbuf(coutbuf); //reset to standard output again

	//-- NB : Test & Evaluation recap in the output files
    factory->TestAllMethods(); // ---- Evaluate all MVAs using the set of test events
    factory->EvaluateAllMethods(); // ----- Evaluate and compare performance of all configured MVAs

	//Could retrieve ROC graph directly
	// TMultiGraph* rocgraph = f.GetROCCurveAsMultiGraph("<datasetname>");

    // --------------------------------------------------------------
    // Save the output
    output_file->Close();
    std::cout << "==> Wrote root file: " << output_file->GetName() << std::endl;
    std::cout << "==> TMVA is done!" << std::endl;

	if(write_ranking_info)
	{
		MoveFile("./ranking_info_tmp.txt", ranking_file_path);
		Extract_Ranking_Info(ranking_file_path, channel); //Extract only ranking info from TMVA output
	}
	else {int tmp = system("rm ./ranking_info_tmp.txt");} //Else remove the temporary ranking file

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


void TopEFT_analysis::Produce_Templates(TString template_name, bool makeHisto_inputVars)
{
//--------------------------------------------
    bool noSysts_inputVars = true; //true <-> don't compute syst weights for histos of input variables (not worth the CPU)

    bool use_specificMVA_eachYear = false;
//--------------------------------------------

    cout<<endl<<YELBKG("                          ")<<endl<<endl;
	if(makeHisto_inputVars) {cout<<FYEL("--- Producing Input variables histograms ---")<<endl;}
	else if(template_name == "") {cout<<FYEL("--- Producing "<<template_name<<" Templates ---")<<endl;}
	else {cout<<BOLD(FRED("--- ERROR : invalid arguments ! Exit !"))<<endl; cout<<"Valid template names are : ttbar / ttV / 2D / 2Dlin !"<<endl; return;}
    cout<<endl<<YELBKG("                          ")<<endl<<endl;

	if(classifier_name != "BDT" && classifier_name != "DNN") {cout<<BOLD(FRED("Error : classifier_name value not supported !"))<<endl; return;}

    TString restore_classifier_name = classifier_name;
	if(makeHisto_inputVars) {classifier_name = "";} //For naming conventions

    vector<TString> restore_syst_list = syst_list;
    vector<TString> restore_systTree_list = systTree_list;
    //Don't make systematics shifted histos for input vars (too long)
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
        array_Btag = new double[4];
        array_LepEff_mu = new double[4];
        array_LepEff_el = new double[4];
    }

//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

	TH1::SetDefaultSumw2();

	//Output file name
	//-- For BDT templates
	TString output_file_name = "outputs/Templates_" + classifier_name + template_name + "_" + region + "_" + lumiName + filename_suffix + ".root";

	//-- For input vars
	if(makeHisto_inputVars) {output_file_name = "outputs/ControlHistograms_" + region + "_" + lumiName + filename_suffix +".root";}

    //Create output file
	TFile* file_output = TFile::Open(output_file_name, "RECREATE");

    // vector<pair<float,float>> v_inputs_rescaling;
    // TString DNN_inputLayerName = ""; TString DNN_outputLayerName = ""; int nNodes = 1;
    if(!makeHisto_inputVars && classifier_name == "BDT")
    {
        //NB : TMVA requires floats, and nothing else, to ensure reproducibility of results (training done with floats) => Need to recast e.g. doubles as flots
        //See : https://sourceforge.net/p/tmva/mailman/message/836453/
        reader = new TMVA::Reader("!Color:!Silent");

        // Name & adress of local variables which carry the updated input values during the event loop
        // NB : the variable names MUST corresponds in name and type to those given in the weight file(s) used -- same order
        // NB : if booking 2 BDTs, must make sure that they use the same input variables... or else, find some way to make it work in the code)
        for(int i=0; i<var_list.size(); i++)
        {
            //cout<<"Added variable "<<var_list[i]<<endl;
            reader->AddVariable(var_list[i].Data(), &var_list_floats[i]);
        }

        for(int i=0; i<v_cut_name.size(); i++)
        {
            if(v_cut_IsUsedForBDT[i] && !v_cut_def[i].Contains("=="))
            {
                reader->AddVariable(v_cut_name[i].Data(), &v_cut_float[i]);
            }
        }
    }
    else if(!makeHisto_inputVars && classifier_name == "DNN") //DNN
    {
        var_list = var_list_DNN; //Use DNN input features
        var_list_floats.resize(var_list.size());

        //--- Load DNN model
        TString DNNmodel_path = "./weights/DNN/"+lumiName+"/model.pb";
        if(!Check_File_Existence(DNNmodel_path) ) {cout<<BOLD(FRED("Model "<<DNNmodel_path<<" not found ! Abort"))<<endl; return;}

        clfy1 = new TFModel(DNNmodel_path.Data(), var_list.size(), DNN_inputLayerName.Data(), nNodes, DNN_outputLayerName.Data()); //Specify names of I/O layers, and nof I/O nodes //These names can be read from the 'model.pbtxt' file produced at DNN training : look for name of first and last nodes *WHICH IS NOT A TRAINING NODE*
    }

	//Input TFile and TTree, called for each sample
	TFile* file_input;
	TTree* tree(0);

	//Template binning
	double xmin = -1, xmax = 1;
	nbins = 10;
    if(classifier_name == "DNN") {xmin = 0;}

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
        if(classifier_name == "BDT") {total_var_list.push_back(classifier_name);}
		else
        {
            for(int inode=0; inode<nNodes; inode++)
            {
                total_var_list.push_back(classifier_name + Convert_Number_To_TString(inode));
            }
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

	// cout<<endl<<ITAL("-- Ntuples directory : "<<dir_ntuples<<"")<<endl<<endl;

    // float tmp_compare = 0;
    float total_nentries_toProcess = Count_Total_Nof_Entries(dir_ntuples, t_name, sample_list, systTree_list, v_cut_name, v_cut_def, v_lumiYears, makeHisto_inputVars, noSysts_inputVars);

    cout<<endl<<FBLU(OVERLINE("                           "))<<endl;
    cout<<FBLU(BOLD("Will process "<<std::setprecision(12)<<total_nentries_toProcess<<" entries..."))<<endl;
    cout<<FBLU(UNDL("                           "))<<endl<<endl<<endl;

    //Draw progress bar
    bool draw_progress_bar = true;
    if(total_nentries_toProcess < 200000) {draw_progress_bar = false;}
    Int_t ibar = 0; //event counter
    TMVA::Timer timer(total_nentries_toProcess, "", true);
    TMVA::gConfig().SetDrawProgressBar(1);
    TMVA::gConfig().SetUseColor(1);

    bool MVA_already_booked = false; //If reading same MVA for multiple years, need to book it only once
    for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
    {
        cout<<endl<<UNDL(FMAG("=== YEAR : "<<v_lumiYears[iyear]<<""))<<endl<<endl;

    	// --- Book the MVA method(s)
    	TString dir = "weights/" + classifier_name + "/" + lumiName;
        if(use_specificMVA_eachYear) {dir = "weights/" + classifier_name + "/" + v_lumiYears[iyear];} //Use year-specific MVA

    	TString MVA_method_name = "";
    	TString weightfile = "";
    	TString template_name_MVA = "";
    	if(!makeHisto_inputVars && classifier_name == "BDT")
    	{
            //Method name does not matter, it is just a way to identify a given method from a given weight file, to be able to call it later
            if(use_specificMVA_eachYear) {template_name_MVA = "BDT_"+signal_process+"_"+v_lumiYears[iyear];} //1 method for each year
            else {template_name_MVA = "BDT_"+signal_process;} //single method
    		MVA_method_name = template_name_MVA + " method";

            weightfile = dir + "/" + classifier_name + "_" + signal_process + ".weights.xml";

            //If weightfile not found for a specific year, try to use instead the weightfile of a BDT training with full Run 2 samples
            if(!Check_File_Existence(weightfile) )
            {
                dir = "weights/" + classifier_name + "/Run2";
                weightfile = dir + "/" + classifier_name + "_" + signal_process + ".weights.xml";

                cout<<BOLD(FRED("Warning : Weight file not found for '"<<lumiName<<"' ! Using '"<<weightfile<<"' (trained on full Run 2) instead ! About to proceed..."))<<endl;
                if(use_specificMVA_eachYear) {usleep(5000000);} //Pause
                else {usleep(1000000);} //Pause
            }

    		if(!Check_File_Existence(weightfile) ) {cout<<BOLD(FRED("Weight file "<<weightfile<<" not found ! Abort"))<<endl; return;}

            // cout<<"MVA_method_name "<<MVA_method_name<<endl;
    		if(use_specificMVA_eachYear || !MVA_already_booked) {reader->BookMVA(MVA_method_name, weightfile); MVA_already_booked = true;}
    	}

    	//SAMPLE LOOP
    	for(int isample=0; isample<sample_list.size(); isample++)
    	{
    		// cout<<endl<<endl<<UNDL(FBLU("Sample : "<<sample_list[isample]<<""))<<endl;

    		//Open input TFile
    		TString inputfile = dir_ntuples + v_lumiYears[iyear] + "/" + sample_list[isample] + ".root";

    		// cout<<"inputfile "<<inputfile<<endl;
    		if(!Check_File_Existence(inputfile))
    		{
    			cout<<endl<<"File "<<inputfile<<FRED(" not found!")<<endl;
    			continue;
    		}

    		file_input = TFile::Open(inputfile, "READ");

            bool isPrivMC = false;
            if(sample_list[isample].Contains("PrivMC")) {isPrivMC = true;}

            vector<float> v_SWE; //Store Sums of Weights (SWE) for all reweight points -- for private MC samples only
            // if(isPrivMC) //Not available yet
            // {
            //     //Read and store sums of weights (SWE)
            //     TH1F* h_SWE = (TH1F*) f->Get("h_SWE");
            //     for(int ibin=0; ibin<h_SWE->GetNbinsX(); ibin++)
            //     {
            //         v_SWE.push_back(h_SWE->GetBinContent(ibin+1)); //1 SWE stored for each stored weight
            //         // cout<<"v_SWE[ibin] = "<<v_SWE[ibin]<<endl;
            //     }
            //     delete h_SWE;
            // }

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
    			tree->SetBranchStatus("*", 0); //disable all branches, speed up
    			// tree->SetBranchStatus("xxx", 1);

    			if(makeHisto_inputVars)
    			{
    				for(int i=0; i<total_var_list.size(); i++)
    				{
                        tree->SetBranchStatus(total_var_list[i], 1);
                        tree->SetBranchAddress(total_var_list[i], &total_var_floats[i]);
    				}
    			}
    			else //Book input variables in same order as trained BDT
    			{
                    for(int i=0; i<var_list.size(); i++)
                    {
                        tree->SetBranchStatus(var_list[i], 1);
                        tree->SetBranchAddress(var_list[i], &var_list_floats[i]);
                        // cout<<"Activate var '"<<var_list[i]<<"'"<<endl;
                    }
    			}

    			for(int i=0; i<v_cut_name.size(); i++)
    			{
                    // cout<<"v_cut_name[i] "<<v_cut_name[i]<<endl;

    				tree->SetBranchStatus(v_cut_name[i], 1);
                    if(v_cut_name[i].BeginsWith("is_") || v_cut_name[i].BeginsWith("passed") ) //Categories are encoded into Char_t, not float
    				{
    					tree->SetBranchAddress(v_cut_name[i], &v_cut_char[i]);
    				}
    				else //All others are floats
    				{
    					tree->SetBranchAddress(v_cut_name[i], &v_cut_float[i]);
    				}
    			}

    			//--- Cut on relevant event selection (e.g. 3l SR, ttZ CR, etc.) -- stored as Char_t
    			Char_t is_goodCategory; //Categ. of event
    			// TString cat_name = Get_Category_Boolean_Name(nLep_cat, region, analysis_type, sample_list[isample], scheme);
    			// tree->SetBranchStatus(cat_name, 1);
    			// tree->SetBranchAddress(cat_name, &is_goodCategory);
    			// cout<<"Categ <=> "<<cat_name<<endl;

    			//--- Cut on relevant categorization (lepton flavour, btagging, charge)
    			float channel;
    			tree->SetBranchStatus("channel", 1);
    			tree->SetBranchAddress("channel", &channel);

                //--- Event weights
                double eventWeight; //Corresponds to 'nominalMEWeight_' in TopAnalysis code
                float eventMCFactor; //Sample-specific SF (xsec*lumi/SWE)
                tree->SetBranchStatus("eventWeight", 1);
    			tree->SetBranchAddress("eventWeight", &eventWeight);
                tree->SetBranchStatus("eventMCFactor", 1);
    			tree->SetBranchAddress("eventMCFactor", &eventMCFactor);
                // float mc_weight_originalValue;
    			// tree->SetBranchStatus("mc_weight_originalValue", 1);
    			// tree->SetBranchAddress("mc_weight_originalValue", &mc_weight_originalValue);

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

                            v3_histo_chan_syst_var[ichan][isyst][ivar] = new TH1F("", "", nbins, xmin, xmax);

                            if(isPrivMC) {v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = new TH1EFT("", "", nbins, xmin, xmax);}
                            else {v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = NULL;}
    					}
    				} //syst
    			} //chan

// ###### #    # ###### #    # #####    #       ####   ####  #####
// #      #    # #      ##   #   #      #      #    # #    # #    #
// #####  #    # #####  # #  #   #      #      #    # #    # #    #
// #      #    # #      #  # #   #      #      #    # #    # #####
// #       #  #  #      #   ##   #      #      #    # #    # #
// ######   ##   ###### #    #   #      ######  ####   ####  #

    			// cout<<"* Tree '"<<systTree_list[itree]<<"' :"<<endl;

    			// int nentries = 10;
    			int nentries = tree->GetEntries();

                // if(!draw_progress_bar) {cout<<endl<< "--- "<<sample_list[isample]<<" : Processing: " << tree->GetEntries() << " events" << std::endl;}
                cout<<endl<< "--- "<<sample_list[isample]<<" : Processing " << nentries << " events" << std::endl;

    			for(int ientry=0; ientry<nentries; ientry++)
    			{
    				// cout<<FGRN("ientry "<<ientry<<"")<<endl;

                    if(isPrivMC && ientry%1000==0) {cout<<DIM("Entry "<<ientry<<"")<<endl;} //Very slow, print progress

    				std::fill(var_list_floats.begin(), var_list_floats.end(), 0); //Reset vectors reading inputs to 0

    				tree->GetEntry(ientry);

    				if(isnan(eventWeight*eventMCFactor) || isinf(eventWeight*eventMCFactor))
    				{
    					cout<<BOLD(FRED("* Found event with eventWeight*eventMCFactor = "<<eventWeight*eventMCFactor<<" ; remove it..."))<<endl; continue;
    				}

    				//--- Cut on category value
    				// if(!is_goodCategory) {continue;}

//---- APPLY CUTS HERE (defined in main)  ----
    				bool pass_all_cuts = true;
    				for(int icut=0; icut<v_cut_name.size(); icut++)
    				{
    					if(v_cut_def[icut] == "") {continue;}

    					//Categories are encoded into Char_t. Convert them to float for code automation
    					if(v_cut_name[icut].BeginsWith("is_") || v_cut_name[icut].BeginsWith("passed") ) {v_cut_float[icut] = (float) v_cut_char[icut];}
    					// cout<<"Cut : name="<<v_cut_name[icut]<<" / def="<<v_cut_def[icut]<<" / value="<<v_cut_float[icut]<<" / pass ? "<<Is_Event_Passing_Cut(v_cut_def[icut], v_cut_float[icut])<<endl;
    					if(!Is_Event_Passing_Cut(v_cut_def[icut], v_cut_float[icut])) {pass_all_cuts = false; break;}
    				}
    				if(!pass_all_cuts) {continue;}
//--------------------------------------------

                    ibar++;
                    if(draw_progress_bar && ibar%50000==0) {timer.DrawProgressBar(ibar, ""); cout<<ibar<<" / "<<total_nentries_toProcess<<endl; }

                    // cout<<"eventWeight "<<eventWeight<<endl;
                    // cout<<"eventMCFactor "<<eventMCFactor<<endl;

    				//Get MVA value to make template
                    if(!makeHisto_inputVars)
                    {
                        if(classifier_name == "BDT") {total_var_floats[0] = reader->EvaluateMVA(MVA_method_name);}

                        //NB -- slow evaluation ! ==> Don't rescale inputs, add lambda layer in model to rescale inputs !
                        else //DNN
                        {
                            //Convert my vector storing input values into array
                            // float clfy1_inputs[var_list_floats.size()];
                            // std::copy(var_list_floats.begin(), var_list_floats.end(), clfy1_inputs);
                            // cout<<"//-------------------------------------------- "<<endl;
                            // for(int i=0; i<var_list_floats.size(); i++) {cout<<"clfy1_inputs["<<i<<"] "<<clfy1_inputs[i]<<std::endl;}
                            // for(int i=0; i<var_list.size(); i++)
                            // {
                            //     var_list_floats[i] = Rescale_Input_Variable(var_list_floats[i], v_inputs_rescaling[i].first, v_inputs_rescaling[i].second);
                            // }
                            // for(int i=0; i<var_list_floats.size(); i++) {cout<<"clfy1_inputs["<<i<<"] "<<clfy1_inputs[i]<<std::endl;}
                            // std::vector<float> clfy1_outputs = clfy1->evaluate(clfy1_inputs);
                            // for(unsigned i=0; i<clfy1_outputs.size(); i++) {cout<<"clfy1_outputs["<<i<<"] "<<clfy1_outputs[i]<<endl;}

                            //Evaluate output nodes values
                            std::vector<float> clfy1_outputs = clfy1->evaluate(var_list_floats);

                            for(int ivar=0; ivar<total_var_list.size(); ivar++) {total_var_floats[ivar] = clfy1_outputs[ivar];}
                            // cout<<"ientry "<<ientry<<" ==> "<<clfy1_outputs[0]<<endl;
                        }
                    }
                    // else
                    // {
                    //     for(int ivar=0; ivar<total_var_list.size(); ivar++)
                    //     {
                    //         //-- Rescale input variables here, to plot them and make sure they are properly rescaled
                    //         // cout<<"//--------------------------------------------"<<endl;
                    //         // cout<<"sample_list["<<isample<<"] "<<sample_list[isample]<<endl;
                    //         // cout<<"total_var_list["<<ivar<<"] "<<total_var_list[ivar]<<endl;
                    //         // cout<<"total_var_floats["<<ivar<<"] "<<total_var_floats[ivar]<<endl;
                    //         // cout<<"v_inputs_rescaling["<<ivar<<"].first "<<v_inputs_rescaling[ivar].first<<endl;
                    //         // cout<<"v_inputs_rescaling["<<ivar<<"].second "<<v_inputs_rescaling[ivar].second<<endl;
                    //         // total_var_floats[ivar] = Rescale_Input_Variable(total_var_floats[ivar], v_inputs_rescaling[ivar].first, v_inputs_rescaling[ivar].second);
                    //         // cout<<"===> total_var_floats["<<ivar<<"] "<<total_var_floats[ivar]<<endl;
                    //
                    //         //Some special variables are already read, must get their proper values
                    //         // if(total_var_list[ivar] == "channel") {total_var_floats[ivar] = channel;}
                    //     }
                    // }

    				//-- Fill histos for all subcategories
    				for(int ichan=0; ichan<channel_list.size(); ichan++)
    				{
                        //Fill histos for sub-channels with corresponding events
                        if(channel_list[ichan] == "uuu" && channel != 0) {continue;}
                        if(channel_list[ichan] == "uue" && channel != 1) {continue;}
                        if(channel_list[ichan] == "eeu" && channel != 2) {continue;}
                        if(channel_list[ichan] == "eee" && channel != 3) {continue;}

    					for(int isyst=0; isyst<syst_list.size(); isyst++)
    					{
                            //-- Protections : not all syst weights apply to all samples, etc.
                            if(sample_list[isample] == "DATA" && syst_list[isyst] != "") {break;}
                            else if(systTree_list[itree] != "") {break;} //Syst event weights only stored in nominal TTree
                            else if(v_lumiYears[iyear] == "2018" && syst_list[isyst].BeginsWith("prefiring") ) {continue;} //no prefire in 2018

    						double weight_tmp = 0; //Fill histo with this weight ; manipulate differently depending on syst

    						// cout<<"-- sample "<<sample_list[isample]<<" / channel "<<channel_list[ichan]<<" / syst "<<syst_list[isyst]<<endl;

    						weight_tmp = eventWeight*eventMCFactor; //Nominal (no syst)

                            //--- add protection for prefire/2018 (set to 1) => need to run separately on each year !
                            // if(syst_list[isyst] != "" && isnan(*(v_double_systWeights[isyst])) ) {*(v_double_systWeights[isyst]) = 1;}

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

                            for(int ivar=0; ivar<total_var_list.size(); ivar++)
                            {
                                Fill_TH1F_UnderOverflow(v3_histo_chan_syst_var[ichan][isyst][ivar], total_var_floats[ivar], weight_tmp);

                                if(isPrivMC)
                                {
                                    Fill_TH1EFT(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar], total_var_floats[ivar], v_ids, v_wgts, v_SWE, weight_tmp);
                                }
                            }
    					} //syst loop

    					if(channel_list[ichan] != "") {break;} //subcategories are orthogonal ; if already found, can break subcat. loop
    				} //subcat/chan loop

    			} //end TTree entries loop
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
    						// if(makeHisto_inputVars)
    						{
    							output_histo_name = total_var_list[ivar] + "_" + region;
    							if(channel_list[ichan] != "") {output_histo_name+= "_" + channel_list[ichan];}
                                output_histo_name+= "_" + v_lumiYears[iyear];
    							output_histo_name+= "__" + samplename;
    							if(syst_list[isyst] != "") {output_histo_name+= "__" + Get_Modified_SystName(syst_list[isyst], v_lumiYears[iyear]);}
    							else if(systTree_list[itree] != "") {output_histo_name+= "__" + systTree_list[itree];}
    						}
    						// else
    						// {
    						// 	output_histo_name = classifier_name + template_name + "_" + region;
    						// 	if(channel_list[ichan] != "") {output_histo_name+= "_" + channel_list[ichan];}
                            //     output_histo_name+= "_" + v_lumiYears[iyear];
    						// 	output_histo_name+= "__" + samplename;
    						// 	if(syst_list[isyst] != "") {output_histo_name+= "__" + Get_Modified_SystName(syst_list[isyst], v_lumiYears[iyear]);}
    						// 	else if(systTree_list[itree] != "") {output_histo_name+= "__" + systTree_list[itree];}
    						// }

    						file_output->cd();

    						v3_histo_chan_syst_var[ichan][isyst][ivar]->Write(output_histo_name);
                            // cout<<"Wrote histo : "<<output_histo_name<<endl;

                            if(isPrivMC)
                            {
                                output_histo_name = "TH1EFT_"+ output_histo_name;
                                v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]->Write(output_histo_name);
                            }

    						delete v3_histo_chan_syst_var[ichan][isyst][ivar]; v3_histo_chan_syst_var[ichan][isyst][ivar] = NULL;
                            if(v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]) {delete v3_TH1EFT_chan_syst_var[ichan][isyst][ivar]; v3_TH1EFT_chan_syst_var[ichan][isyst][ivar] = NULL;}
    					} //var loop
    				} //syst loop
    			} //chan loop

                if(isPrivMC) {delete v_wgts; delete v_ids;}

    			// cout<<"Done with "<<sample_list[isample]<<" sample"<<endl;

    			tree->ResetBranchAddresses(); //Detach tree from local variables (safe)
    			delete tree; tree = NULL;
    		} //end tree loop

    		file_input->Close(); file_input = NULL;
    	} //end sample loop

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

    if(!makeHisto_inputVars)
    {
        if(classifier_name == "BDT")  {delete reader; reader = NULL;}
        else {delete clfy1; clfy1 = NULL;}
    }

    //-- Can verify that the total nof processed entries computed from Count_Total_Nof_Entries() was effectively the nof processed entries
    // cout<<"total_nentries_toProcess --> "<<total_nentries_toProcess<<endl;
    // cout<<"tmp_compare --> "<<tmp_compare<<endl;

    //Restore potentially modified variables
    classifier_name = restore_classifier_name;
    syst_list = restore_syst_list;
    systTree_list = restore_systTree_list;

    if(!makeHisto_inputVars || !noSysts_inputVars) //free memory
    {
        delete array_PU; array_PU = NULL;
        delete array_prefiringWeight; array_prefiringWeight = NULL;
        delete array_Btag; array_Btag = NULL;
        delete array_LepEff_mu; array_LepEff_mu = NULL;
        delete array_LepEff_el; array_LepEff_el = NULL;
    }

// #    # ###### #####   ####  ######
// ##  ## #      #    # #    # #
// # ## # #####  #    # #      #####
// #    # #      #####  #  ### #
// #    # #      #   #  #    # #
// #    # ###### #    #  ####  ######

    if(!makeHisto_inputVars) //For COMBINE fit, want to directly merge contributions from different processes into single histograms
    {
        Merge_Templates_ByProcess(output_file_name, template_name, total_var_list);
    }

    // Test_TH1EFT();

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

void TopEFT_analysis::Draw_Templates(bool drawInputVars, TString channel, TString template_name, bool prefit, bool use_combine_file)
{
//--------------------------------------------
	// bool doNot_stack_signal = false;

	bool draw_errors = true; //true <-> superimpose error bands on plot/ratio plot

	bool draw_logarithm = false;

    bool superimpose_GENhisto = false; //Superimpose corresponding GEN-level EFT histogram, for shape comparison...?

//--------------------------------------------

    cout<<endl<<YELBKG("                          ")<<endl<<endl;
	if(drawInputVars) {cout<<FYEL("--- Producing Input Variables Plots / channel : "<<channel<<" ---")<<endl;}
	else if(template_name == "") {cout<<FYEL("--- Producing "<<classifier_name<<" Template Plots ---")<<endl;}
	else {cout<<FRED("--- ERROR : invalid args !")<<endl;}
    cout<<endl<<YELBKG("                          ")<<endl<<endl;

	if(drawInputVars)
	{
		classifier_name = ""; //For naming conventions
		use_combine_file = false;
		if(drawInputVars && !prefit) {cout<<"Error ! Can not draw postfit input vars yet !"<<endl; return;}
		if(template_name == "categ" && !prefit) {cout<<"Can not plot yields per subcategory using the Combine output file ! Will plot [PREFIT] instead of [POSTFIT] !"<<endl; prefit = true;}
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

	//Get input TFile
	if(!prefit)
	{
		use_combine_file = true;
	}
	if(drawInputVars && use_combine_file)
	{
		cout<<"-- Setting 'use_combine_file = false' !"<<endl;
		use_combine_file = false;
	}

	TString input_name = "";

	if(use_combine_file)
	{
        //TRY 1 : look for Combine file
		input_name = "./outputs/fitDiagnostics_";
		input_name+= classifier_name + template_name + "_" + region + filename_suffix + ".root";
		if(!Check_File_Existence(input_name)) {input_name = "./outputs/fitDiagnostics.root";} //Try another name

        //Combine file not found !
		if(!Check_File_Existence(input_name))
		{
			cout<<FBLU("-- NB : File ")<<input_name<<FBLU(" (produced by COMBINE) not found !")<<endl;
			if(!prefit) {cout<<FRED("===> Can not produce postfit plots ! Abort !")<<endl; return;}
			use_combine_file = false;

            //TRY 2 : look for my own file containing prefit templates/control histos
			if(drawInputVars) {input_name = "outputs/ControlHistograms_" + region + "_" + lumiName + filename_suffix + ".root";}
			else //Templates
			{
				input_name = "outputs/Templates_" + classifier_name + template_name + "_" + region + "_" + lumiName + filename_suffix + ".root";
			}
            if(!Check_File_Existence(input_name) && lumiName != "Run2") //If did not find year-specific file, also try to look for full Run 2 file (contains contributions from each year)
            {
                if(drawInputVars) {input_name = "outputs/ControlHistograms_" + region + "_Run2" + filename_suffix + ".root";}
    			else //Templates
    			{
    				input_name = "outputs/Templates_" + classifier_name + template_name + "_" + region + "_Run2" + filename_suffix + ".root";
    			}
            }

            //Did not even find my own file --> Can not plot anything !
			if(!Check_File_Existence(input_name))
			{
				cout<<FRED("Did not find any file containing histos to plot ! Either the files do not exist, or naming conventions are not respected... (check the code) ! Abort !")<<endl;
				return;
			}
			else {cout<<FBLU("--> Using file ")<<input_name<<FBLU(" instead !")<<endl; usleep(3000000);}
		}
		else {cout<<FBLU("--> Using Combine output file : ")<<input_name<<FBLU(" (NB : total error included)")<<endl; use_combine_file = true;}
	}

	else //Using my own template file
	{
		if(drawInputVars) {input_name = "outputs/ControlHistograms_" + region + "_" + lumiName + filename_suffix + ".root";}
		else //Templates
		{
			input_name = "outputs/Templates_" + classifier_name + template_name + "_" + region + "_" + lumiName + filename_suffix + ".root";
		}

        if(!Check_File_Existence(input_name) && lumiName != "Run2") //If did not find year-specific file, also try to look for full Run 2 file (contains contributions from each year)
        {
            if(drawInputVars) {input_name = "outputs/ControlHistograms_" + region + "_Run2" + filename_suffix + ".root";}
            else //Templates
            {
                input_name = "outputs/Templates_" + classifier_name + template_name + "_" + region + "_Run2" + filename_suffix + ".root";
            }
        }

		if(!Check_File_Existence(input_name))
		{
            cout<<FRED("Did not find any file containing histos to plot ! Either the files do not exist, or their naming convention is wrong (check code) ! Abort !")<<endl;
			return;
		}
	}
	cout<<endl<<endl<<endl;
	usleep(1000000); //Pause for 1s (in microsec)

	//Input file containing histos
	TFile* file_input = TFile::Open(input_name);

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

	//Want to plot ALL selected variables
	vector<TString> total_var_list;
	if(drawInputVars)
	{
		// for(int i=0; i<v_cut_name.size(); i++)
		// {
		// 	if(v_cut_name[icut].BeginsWith("is_") || v_cut_name[icut].BeginsWith("passed")) {continue;} //Don't care about plotting the categories
		//
		// 	total_var_list.push_back(v_cut_name[i]);
		// }
		for(int i=0; i<var_list.size(); i++)
		{
			// if(!v_var_tmp->at(i).Contains("mem") || var_list[i].Contains("e-")) {continue;}

			total_var_list.push_back(var_list[i]);
		}
		for(int i=0; i<v_add_var_names.size(); i++)
		{
			total_var_list.push_back(v_add_var_names[i]);
		}
	}
	else
	{
        if(classifier_name == "DNN")
        {
            for(int inode=0; inode<nNodes; inode++) {total_var_list.push_back(classifier_name + Convert_Number_To_TString(inode));} //1 template per output node
        }
        else {total_var_list.push_back(classifier_name);} //Single BDT template
	}


// #       ####   ####  #####   ####
// #      #    # #    # #    # #
// #      #    # #    # #    #  ####
// #      #    # #    # #####       #
// #      #    # #    # #      #    #
// ######  ####   ####  #       ####

	for(int ivar=0; ivar<total_var_list.size(); ivar++)
	{
		if(drawInputVars) {cout<<endl<<FBLU("* Variable : "<<total_var_list[ivar]<<" ")<<endl<<endl;}

		//TH1F* to retrieve distributions
		TH1F* h_tmp = 0; //Tmp storing histo

		TH1F* h_tzq = 0; //Store tZq shape
		TH1F* h_ttz = 0; //Store ttZ shape
		TH1F* h_sum_data = 0; //Will store data histogram
		vector<TH1F*> v_MC_histo; //Will store all MC histograms (1 TH1F* per MC sample)

		TGraphAsymmErrors* g_data = 0; //If using Combine file, data are stored in TGAE
		TGraphAsymmErrors* g_tmp = 0; //Tmp storing graph

		vector<TString> MC_samples_legend; //List the MC samples which are actually used (to get correct legend)

		//-- Init error vectors
		double x, y, errory_low, errory_high;

		vector<double> v_eyl, v_eyh, v_exl, v_exh, v_x, v_y; //Contain the systematic errors (used to create the TGraphError)
		int nofbins=-1;

		vector<float> v_yield_sig, v_yield_bkg;

//-- All histos are for given lumiYears and sub-channels --> Need to sum them all for plots
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
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
    			dir_hist+= classifier_name + template_name + "_" + region;
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

    				// cout<<endl<<UNDL(FBLU("-- Sample : "<<sample_list[isample]<<" : "))<<endl;

    				h_tmp = 0;

    				TString histo_name = "";
    				if(use_combine_file) {histo_name = dir_hist + samplename;}
    				else
    				{
    					histo_name = total_var_list[ivar] + "_" + region;
    					if(channel != "") {histo_name+= "_" + channel;}
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

    					//-- NEW, draw all errors
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

    				if(!samplename.Contains("DATA") )
    				{
    					if(v_MC_histo.size() <=  index_MC_sample) {MC_samples_legend.push_back(samplename);}
    					// MC_samples_legend.push_back(samplename); //Fill vector containing existing MC samples names
    					// cout<<"ADD samplename "<<samplename<<endl;
    					// cout<<"MC_samples_legend.size() "<<MC_samples_legend.size()<<endl;
    				}

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
    						if(v_isSkippedSample[isample-k]) {continue;}
    						else if(color_list[isample] == color_list[isample-k]) {h_tmp->SetLineColor(color_list[isample]); break;}
    					}
    					else {break;}
    				}

    				// v_MC_histo.push_back((TH1F*) h_tmp->Clone());
    				if(v_MC_histo.size() <=  index_MC_sample) {v_MC_histo.push_back((TH1F*) h_tmp->Clone());}
    				else if(!v_MC_histo[index_MC_sample]) {v_MC_histo[index_MC_sample] = (TH1F*) h_tmp->Clone();} //For FakeEle and FakeMu
    				else {v_MC_histo[index_MC_sample]->Add((TH1F*) h_tmp->Clone());}

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

    				delete h_tmp; h_tmp = 0;
    			} //end sample loop
    		} //channel loop
        } //years loop

		//Printout S/B for each bin
		// for(int ibin=0; ibin<v_yield_sig.size(); ibin++)
		// {
		// 	cout<<endl<<"ibin : S/B = "<<v_yield_sig[ibin]/v_yield_bkg[ibin]<<" // bkg = "<<v_yield_bkg[ibin]<<" // sig = "<<v_yield_sig[ibin]<<endl;
		// }

// #####    ##   #####   ##
// #    #  #  #    #    #  #
// #    # #    #   #   #    #
// #    # ######   #   ######
// #    # #    #   #   #    #
// #####  #    #   #   #    #


        //--- Retrieve DATA histo(s)
        TString data_histo_name = "";
        for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
        {
    		h_tmp = 0;
    		if(use_combine_file) {data_histo_name = "data";}
    		else
    		{
                // data_histo_name = classifier_name + total_var_list[ivar] + "_" + region;
                data_histo_name = total_var_list[ivar] + "_" + region;
    			if(channel != "") {data_histo_name+= "_" + channel;}
    			data_histo_name+= "_" + v_lumiYears[iyear];
                data_histo_name+= "__data_obs";
    		}

    		if(use_combine_file)
    		{
    			for(int ichan=0; ichan<channel_list.size(); ichan++)
    			{
    				if(channel != "" && channel_list[ichan] != channel) {continue;}

    				//Combine file : histos stored in subdirs -- define dir name
    				TString dir_hist = "";
    				if(prefit) {dir_hist = "shapes_prefit/";}
    				else {dir_hist = "shapes_fit_s/";}
    				dir_hist+= classifier_name + template_name + "_" + region;
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
    				if(h_sum_data == 0) {h_sum_data = (TH1F*) h_tmp->Clone();}
    				else {h_sum_data->Add((TH1F*) h_tmp->Clone());} //not needed anymore (1 channel only)

    				delete h_tmp; h_tmp = NULL;
    			}
    		}
        } //years loop

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
				}
			}
			else
			{
				for(int ibin = 1; ibin<h_sum_data->GetNbinsX()+1; ibin++)
				{
					if(h_sum_data->GetBinContent(ibin) < 0) {h_sum_data->SetBinContent(ibin, 0);}
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
		int index_tZq_sample = -99;
		int index_ttZ_sample = -99;
		int index_NPL_sample = -99;

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


// ####### #
//    #    #       ######  ####  ###### #    # #####
//    #    #       #      #    # #      ##   # #    #
//    #    #       #####  #      #####  # #  # #    #
//    #    #       #      #  ### #      #  # # #    #
//    #    #       #      #    # #      #   ## #    #
//    #    ####### ######  ####  ###### #    # #####

		TLegend* qw = 0;
        float x_left = 0.94-ceil(nSampleGroups/2.)*0.10; //each column allocated same x-space
        // qw = new TLegend(x_left,0.88,0.99,0.99);
        // qw = new TLegend(x_left,0.75,0.95,0.89); //y-space for 2 rows
        qw = new TLegend(x_left,0.78,0.94,0.87); //y-space for 2 rows
        qw->SetNColumns(ceil(nSampleGroups/2.));
        // qw->SetLineColor(1);
        qw->SetBorderSize(0);
		qw->SetTextSize(0.04);

		//--Data on top of legend
        if(use_combine_file && g_data != 0) {qw->AddEntry(g_data, "Data" , "ep");}
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
            else if(MC_samples_legend[i] == "ttW") {qw->AddEntry(v_MC_histo[i], "t#bar{t}X", "f");}
            else if(MC_samples_legend[i] == "tHq") {qw->AddEntry(v_MC_histo[i], "tX", "f");}
			else if(MC_samples_legend[i] == "WZ") {qw->AddEntry(v_MC_histo[i], "VV(V)", "f");}
			else if(MC_samples_legend[i] == "Fakes") {qw->AddEntry(v_MC_histo[i], "Non-prompt", "f");}
			else if(MC_samples_legend[i] == "QFlip") {qw->AddEntry(v_MC_histo[i], "Flip", "f");}
            else if(MC_samples_legend[i] == "GammaConv") {qw->AddEntry(v_MC_histo[i], "#gamma-conv.", "f");}
            else if(MC_samples_legend[i] == "DY" ) {qw->AddEntry(v_MC_histo[i], "V+jets", "f");}
            else if(MC_samples_legend[i] == "TTbar_DiLep" || MC_samples_legend[i] == "TTbar") {qw->AddEntry(v_MC_histo[i], "t#bar{t}", "f");}

            //group names
            else if(MC_samples_legend[i] == "ttX" ) {qw->AddEntry(v_MC_histo[i], "t#bar{t}X", "f");}
            else if(MC_samples_legend[i] == "tX" ) {qw->AddEntry(v_MC_histo[i], "tX", "f");}
            else if(MC_samples_legend[i] == "VV" ) {qw->AddEntry(v_MC_histo[i], "VV", "f");}
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

			stack_MC->Add(v_MC_histo[i]);
			// cout<<"Stacking sample "<<MC_samples_legend[i]<<" / integral "<<v_MC_histo[i]->Integral()<<endl;
		}

		//Set Yaxis maximum & minimum
		if(use_combine_file && data_notEmpty)
		{
			Long64_t locmax = TMath::LocMax(g_data->GetN(), g_data->GetY()); //the corresponding x value can be obtained with double xmax = gr->GetX()[locmax];
			double ymax = g_data->GetY()[locmax];

			if(ymax > stack_MC->GetMaximum() ) {stack_MC->SetMaximum(ymax*1.3);}
			else stack_MC->SetMaximum(stack_MC->GetMaximum()*1.5);
		}
		else if(!use_combine_file && data_notEmpty)
		{
			if(h_sum_data->GetMaximum() > stack_MC->GetMaximum() ) {stack_MC->SetMaximum(h_sum_data->GetMaximum()+0.3*h_sum_data->GetMaximum());}
			else {stack_MC->SetMaximum(stack_MC->GetMaximum()*1.5);}
		}
        else if(!data_notEmpty) {stack_MC->SetMaximum(stack_MC->GetMaximum()*1.5);}

		stack_MC->SetMinimum(0.0001); //Remove '0' label

		if(draw_logarithm)
		{
			stack_MC->SetMinimum(0.5);
			stack_MC->SetMaximum(stack_MC->GetMaximum()*6);
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

		//Superimpose shape of signal
		// if(doNot_stack_signal)
		// {
		// 	if(h_tzq != 0) {h_tzq->Draw("same hist");}
		// 	if(h_ttz != 0) {h_ttz->Draw("same hist");}
		// }

		qw->Draw("same"); //Draw legend


// ###### ###### #####     ####  ###### #    #
// #      #        #      #    # #      ##   #
// #####  #####    #      #      #####  # #  #
// #      #        #      #  ### #      #  # #
// #      #        #      #    # #      #   ##
// ###### #        #       ####  ###### #    #

        TH1F* h_GEN = 0;
        if(superimpose_GENhisto && drawInputVars)
        {
            Get_Pointer_GENHisto(h_GEN, var_list[ivar]);
        }

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
		// histo_ratio_data->GetXaxis()->SetTitleOffset(1); //FIXME
		histo_ratio_data->GetYaxis()->SetTitleOffset(1.2);
		histo_ratio_data->GetYaxis()->SetLabelSize(0.048);
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
			histo_ratio_data->SetMinimum(-0.2);
			histo_ratio_data->SetMaximum(2.2);
		}

		if(drawInputVars) {histo_ratio_data->GetXaxis()->SetTitle(Get_Variable_Name(total_var_list[ivar]));}
		else
		{
            // histo_ratio_data->GetXaxis()->SetTitle(classifier_name+" (vs "+template_name + ")");
            histo_ratio_data->GetXaxis()->SetTitle(total_var_list[ivar]);

            //Hardcode DNN output nodes names...?
            if(total_var_list[ivar] == "DNN0") {histo_ratio_data->GetXaxis()->SetTitle("DNN (tZq node)");}
            else if(total_var_list[ivar] == "DNN1" && nNodes == 3) {histo_ratio_data->GetXaxis()->SetTitle("DNN (ttZ node)");}
            else if(total_var_list[ivar] == "DNN2" && nNodes == 3) {histo_ratio_data->GetXaxis()->SetTitle("DNN (Bkgs node)");}

			if(template_name == "categ") //Vertical text X labels (categories names)
			{
				histo_ratio_data->GetXaxis()->SetTitle("Categ.");
				histo_ratio_data->GetXaxis()->SetLabelSize(0.08);
				histo_ratio_data->GetXaxis()->SetLabelOffset(0.02);

				// int nx = 2;
				// if(nLep_cat == "2l")
				// {
				// 	const char *labels[5]  = {"ee", "e#mu bl", "e#mu bt", "#mu#mu bl", "#mu#mu bt"};
				// 	for(int i=1;i<=5;i++) {histo_ratio_data->GetXaxis()->SetBinLabel(i,labels[i-1]);}
				// }
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
					h_line1->SetBinContent(ibin, 0.5);
					h_line2->SetBinContent(ibin, 1.5);
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
			stack_MC->GetYaxis()->SetTitleOffset(1.2);
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
		text2.DrawLatex(0.23,0.86,info_data);


// #    # #####  # ##### ######     ####  #    # ##### #####  #    # #####
// #    # #    # #   #   #         #    # #    #   #   #    # #    #   #
// #    # #    # #   #   #####     #    # #    #   #   #    # #    #   #
// # ## # #####  #   #   #         #    # #    #   #   #####  #    #   #
// ##  ## #   #  #   #   #         #    # #    #   #   #      #    #   #
// #    # #    # #   #   ######     ####   ####    #   #       ####    #

		if(drawInputVars)
		{
			mkdir("plots/input_vars", 0777);
            mkdir(("plots/input_vars/"+region).Data(), 0777);
            mkdir(("plots/input_vars/"+region+"/"+lumiName).Data(), 0777);
		}
		else
		{
			mkdir("plots/templates", 0777);
            mkdir( ("plots/templates/"+region).Data(), 0777);
            mkdir( ("plots/templates/"+region+"/"+lumiName).Data(), 0777);
			if(prefit) {mkdir( ("plots/templates/"+region+"/"+lumiName+"/prefit").Data(), 0777);}
			else {mkdir( ("plots/templates/"+region+"/"+lumiName+"/postfit").Data(), 0777);}
		}

		//Output
		TString output_plot_name;

		if(drawInputVars)
		{
			output_plot_name = "plots/input_vars/" + region + "/" + lumiName + "/" + total_var_list[ivar];
		}
		else
		{
			output_plot_name = "plots/templates/" + region + "/" + lumiName;
			if(prefit) {output_plot_name+= "/prefit/";}
			else {output_plot_name+= "/postfit/";}
            // output_plot_name+= classifier_name + template_name + "_template_" + signal_process;
            output_plot_name+= total_var_list[ivar] + "_template_" + signal_process;
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

void TopEFT_analysis::Compare_TemplateShapes_Processes(TString template_name, TString channel)
{
	bool drawInputVars = true;

	bool normalize = true;

	if(drawInputVars)
	{
		template_name = "metEt"; //Can hardcode here another variable
		classifier_name = ""; //For naming convention
	}

//--------------------------------------------
    //-- Hardcode samples here... or could filter the main sample list
	vector<TString> v_samples; vector<TString> v_groups; vector<int> v_colors;
	v_samples.push_back("tZq"); v_groups.push_back("tZq"); v_colors.push_back(kBlack);

    vector<TString> v_syst;
    v_syst.push_back("");
    // v_syst.push_back("JESUp");

//--------------------------------------------

    cout<<endl<<YELBKG("                          ")<<endl<<endl;
	if(drawInputVars) {cout<<FYEL("--- Producing Input Vars Plots / channel : "<<channel<<" ---")<<endl;}
    cout<<endl<<YELBKG("                          ")<<endl<<endl;

//  ####  ###### ##### #    # #####
// #      #        #   #    # #    #
//  ####  #####    #   #    # #    #
//      # #        #   #    # #####
// #    # #        #   #    # #
//  ####  ######   #    ####  #

	//Get input TFile
	TString input_name;
	if(drawInputVars)
	{
		input_name = "outputs/ControlHistograms_" + region + filename_suffix + ".root";
	}
	else
	{
		input_name = "outputs/Templates_" + classifier_name + template_name + "_" + region + filename_suffix + ".root";
	}

	if(!Check_File_Existence(input_name))
	{
		cout<<FRED("File "<<input_name<<" (prefit templates) not found ! Did you specify the region/background ? Abort")<<endl;
		return;
	}
	else {cout<<FBLU("--> Using file ")<<input_name<<FBLU(" instead (NB : only stat. error will be included)")<<endl;}
	cout<<endl<<endl<<endl;

	//Input file containing histos
	TFile* file_input = 0;
	file_input = TFile::Open(input_name);

	//TH1F* to retrieve distributions
	TH1F* h_tmp = 0; //Tmp storing histo

	vector<vector<TH1F*>> v2_histos(v_samples.size()); //store histos, for each sample/syst


// #       ####   ####  #####   ####
// #      #    # #    # #    # #
// #      #    # #    # #    #  ####
// #      #    # #    # #####       #
// #      #    # #    # #      #    #
// ######  ####   ####  #       ####

	//Combine output : all histos are given for subcategories --> Need to sum them all
	for(int ichan=0; ichan<channel_list.size(); ichan++)
	{
		if(channel_list[ichan] != channel) {continue;}

		for(int isample = 0; isample<v_samples.size(); isample++)
		{
			cout<<endl<<UNDL(FBLU("-- Sample : "<<v_samples[isample]<<" : "))<<endl;

			v2_histos[isample].resize(v_syst.size());

            TFile* f;
            f = file_input;

			TString samplename = v_samples[isample];

            for(int isyst=0; isyst<v_syst.size(); isyst++)
            {
				if(v_samples[isample].Contains("Fake") && !v_syst[isyst].Contains("Clos") && !v_syst[isyst].Contains("FR") && v_syst[isyst] != "") {continue;}

				// cout<<"syst "<<v_syst[isyst]<<endl;

                h_tmp = 0;

    			TString histo_name = "";
    			histo_name = classifier_name + template_name + "_" + region;
    			if(channel != "") {histo_name+= "_" + channel;}
    			histo_name+= + "__" + samplename;
                if(v_syst[isyst] != "") {histo_name+= "__" + v_syst[isyst];}

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

                v2_histos[isample][isyst] = (TH1F*) h_tmp->Clone();

				// cout<<"v2_histos["<<isample<<"]["<<isyst<<"]->Integral() "<<v2_histos[isample][isyst]->Integral()<<endl;

    			delete h_tmp; h_tmp = 0;
            } //end syst loop

			// f->Close();
		} //end sample loop
	} //subcat loop


// #####  #####    ##   #    #
// #    # #    #  #  #  #    #
// #    # #    # #    # #    #
// #    # #####  ###### # ## #
// #    # #   #  #    # ##  ##
// #####  #    # #    # #    #

	//Canvas definition
	Load_Canvas_Style();
	gStyle->SetOptTitle(1);
	TCanvas* c = new TCanvas("", "", 1000, 800);
	c->SetTopMargin(0.1);

	// c->SetLogy();

	TLegend* qw;
	qw = new TLegend(0.75,.70,1.,1.);

	c->cd();

	for(int isample=0; isample<v2_histos.size(); isample++)
	{
		TString systlist = "";

		for(int isyst=0; isyst<v_syst.size(); isyst++)
		{
			if(v_samples[isample].Contains("Fake") && !v_syst[isyst].Contains("Clos") && !v_syst[isyst].Contains("FR") && v_syst[isyst] != "") {continue;}

			if(v2_histos[isample][isyst] == 0) {cout<<"Null histo ! Skip"<<endl; continue;}

			if(isample == 0)
			{
				if(v_syst[isyst] != "")
				{
					systlist+= " / " + v_syst[isyst];
				}

				v2_histos[isample][isyst]->SetTitle(systlist);
				// if(normalize) {v2_histos[isample][isyst]->SetMaximum(0.45);}
				// else {v2_histos[isample][isyst]->SetMaximum(6.);}
				// else {v2_histos[isample][isyst]->SetMaximum(11.);}

				if(drawInputVars) {v2_histos[isample][isyst]->GetXaxis()->SetTitle(template_name);}
                // else {v2_histos[isample][isyst]->GetXaxis()->SetTitle(classifier_name+" (vs "+template_name + ")");}
                else {v2_histos[isample][isyst]->GetXaxis()->SetTitle(classifier_name);}
			}

			if(isample < v2_histos.size() - 1)
			{
				if(v_groups[isample] == v_groups[isample+1])
				{
					v2_histos[isample+1][isyst]->Add(v2_histos[isample][isyst]); //Merge with next sample
					continue; //Will draw merged histo, not this single one
				}
			}

			if(normalize) {v2_histos[isample][isyst]->SetMaximum(0.5);}
			v2_histos[isample][isyst]->Draw("hist same");

			if(v_syst[isyst] == "")
			{
                if(v_groups[isample].Contains("tZq")) {qw->AddEntry(v2_histos[isample][isyst], "tZq", "f");}
                else if(v_groups[isample].EndsWith("ttZ") ) {qw->AddEntry(v2_histos[isample][isyst], "t#bar{t}Z", "f");}
                else if(v_groups[isample] == "ttW") {qw->AddEntry(v2_histos[isample][isyst], "t#bar{t}X", "f");}
                else if(v_groups[isample] == "tHq") {qw->AddEntry(v2_histos[isample][isyst], "tX", "f");}
    			else if(v_groups[isample] == "WZ") {qw->AddEntry(v2_histos[isample][isyst], "VV(V)", "f");}
    			else if(v_groups[isample] == "Fakes") {qw->AddEntry(v2_histos[isample][isyst], "Non-prompt", "f");}
    			else if(v_groups[isample] == "QFlip") {qw->AddEntry(v2_histos[isample][isyst], "Flip", "f");}
                else if(v_groups[isample] == "GammaConv") {qw->AddEntry(v2_histos[isample][isyst], "#gamma-conv.", "f");}
                else if(v_groups[isample] == "DY" ) {qw->AddEntry(v2_histos[isample][isyst], "V+jets", "f");}
                else if(v_groups[isample] == "TTbar_DiLep" ) {qw->AddEntry(v2_histos[isample][isyst], "t#bar{t}", "f");}
			}

			//HARDCODED
			if(v_syst[isyst] == "JESUp") {qw->AddEntry(v2_histos[isample][isyst], "JES Up", "L");}
			if(v_syst[isyst] == "JESDown") {qw->AddEntry(v2_histos[isample][isyst], "JES Down", "L");}
		}
	} //sample loop

	qw->Draw("same");


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

    mkdir("plots/templates_shapes", 0777);
    mkdir(("plots/templates_shapes/"+lumiName).Data(), 0777);

	//Output
	TString output_plot_name = "plots/templates_shapes/";
	output_plot_name+= classifier_name + template_name + "_" + region +"_templatesShapes";
	if(channel != "") {output_plot_name+= "_" + channel;}
	output_plot_name+= this->filename_suffix + this->plot_extension;

	c->SaveAs(output_plot_name);

	for(int isample=0; isample<v2_histos.size(); isample++)
	{
		for(int isyst=0; isyst<v_syst.size(); isyst++)
		{
			delete v2_histos[isample][isyst];
		}
	}

	delete c; c = NULL;
	delete qw; qw = NULL;

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
//NB : some inconsistent indices : down variation may be element 0 or 1... hence, hard-coded !
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
        if(systname.EndsWith("HUp")) {index = 0;}
        else if(systname.EndsWith("HDown")) {index = 1;}
        else if(systname.EndsWith("LUp")) {index = 2;}
        else if(systname.EndsWith("LDown")) {index = 3;}
    }
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
    else{cout<<FRED("ERROR ! Systematic '"<<systname<<"' not included in function SetBranchAddress_SystVariation() from Helper.cxx ! Can *not* compute it !")<<endl; return;}

    t->SetBranchStatus(array_name, 1);
    t->SetBranchAddress(array_name, address_memberArray);
    v_doubles[isyst] = &address_memberArray[index];

    return;
}









//--------------------------------------------
 // #    # ###### #####   ####  ######    ##### ###### #    # #####  #        ##   ##### ######  ####
 // ##  ## #      #    # #    # #           #   #      ##  ## #    # #       #  #    #   #      #
 // # ## # #####  #    # #      #####       #   #####  # ## # #    # #      #    #   #   #####   ####
 // #    # #      #####  #  ### #           #   #      #    # #####  #      ######   #   #           #
 // #    # #      #   #  #    # #           #   #      #    # #      #      #    #   #   #      #    #
 // #    # ###### #    #  ####  ######      #   ###### #    # #      ###### #    #   #   ######  ####
//--------------------------------------------

/**
 * The main code is producing/plotting template histograms for each subprocess independently
 * However, in COMBINE, may want to group some processes (e.g. "Rares", "EWK", ...)
 * ===> In addition to individual histos, also merge the relevant subprocesses together and store the merged histos
 * NB : here the order of loops is important because we sum histograms recursively, and the 'sample_list' loop must be the most nested one !
 */
void TopEFT_analysis::Merge_Templates_ByProcess(TString filename, TString template_name, vector<TString> total_var_list, bool force_normTemplate_positive/*=false*/)
{
	cout<<FYEL("==> Merging some templates in file : ")<<filename<<endl;

	if(!Check_File_Existence(filename) ) {cout<<endl<<FRED("File "<<filename<<" not found! Abort template merging !")<<endl; return;}
	TFile* f = TFile::Open(filename, "UPDATE");

	//NB :here the order of loops is important because we sum histograms recursively ! The 'sample_list' loop *must be the most nested one* !
    for(int ivar=0; ivar<total_var_list.size(); ivar++) //There may be more than 1 template, e.g. several DNN output nodes
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

        				TH1F* h_merging = 0; //Merged histogram

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

        					TString histoname = total_var_list[ivar] + "_" + region;
        					if(channel_list[ichan] != "") {histoname+= "_" + channel_list[ichan];}
                            histoname+= "_" + v_lumiYears[iyear];
                            histoname+= "__" + samplename;
        					if(syst_list[isyst] != "") {histoname+= "__" + Get_Modified_SystName(syst_list[isyst], v_lumiYears[iyear]);}
        					else if(systTree_list[itree] != "") {histoname+= "__" + systTree_list[itree];}
        					// cout<<"histoname "<<histoname<<endl;

        					if(!f->GetListOfKeys()->Contains(histoname) && systTree_list[itree] == "" && syst_list[isyst] == "")
        					{
        						cout<<FRED("Histo "<<histoname<<" not found in file "<<filename<<" !")<<endl;
        					 	continue;
        					}

        					TH1F* h_tmp = (TH1F*) f->Get(histoname); //Get individual histograms
        					// cout<<"h_tmp->Integral() = "<<h_tmp->Integral()<<endl;

        					int factor = +1; //Addition
        					// if(sample_list[isample] == "Fakes_MC") {factor = -1;} //Substraction of 'MC Fakes' (prompt contribution to fakes)

        					if(h_tmp != 0)
        					{
        						if(!h_merging) {h_merging = (TH1F*) h_tmp->Clone();}
        						else {h_merging->Add(h_tmp, factor);}
        					}
        					else {cout<<"h_tmp null !"<<endl;}

        					// cout<<"h_merging->Integral() = "<<h_merging->Integral()<<endl;

        					delete h_tmp; h_tmp = 0;
        					if(!h_merging) {cout<<"Syst "<<syst_list[isyst]<<systTree_list[itree]<<" / chan "<<channel_list[ichan]<<" / sample "<<sample_list[isample]<<endl; cout<<"h_merging is null ! Fix this first"<<endl; return;}

        					//Check if next sample will be merged with this one, or else if must write the histogram
        					if(isample < sample_list.size()-1 && sample_groups[isample+1] == sample_groups[isample]) {continue;}
        					else
        					{
        						// if(force_normTemplate_positive)
        						// {
        						// 	//If integral of histo is negative, set to 0 (else COMBINE crashes) -- must mean that norm is close to 0 anyway
        						// 	if(h_merging->Integral() <= 0)
        						// 	{
        						// 		// cout<<endl<<"While merging processes by groups ('Rares'/...) :"<<endl<<FRED(" h_merging->Integral() = "<<h_merging->Integral()<<" (<= 0) ! Distribution set to ~>0 (flat), to avoid crashes in COMBINE !")<<endl;
        						// 		Set_Histogram_FlatZero(h_merging, true, "h_merging");
        						// 		cout<<"(Syst "<<syst_list[isyst]<<systTree_list[itree]<<" / chan "<<channel_list[ichan]<<" / sample "<<sample_list[isample]<<")"<<endl;
        						// 	}
        						// }
        						// cout<<"h_merging->Integral() = "<<h_merging->Integral()<<endl;

        						TString histoname_new = total_var_list[ivar] + "_" + region;
        						if(channel_list[ichan] != "") {histoname_new+="_"  + channel_list[ichan];}
                                histoname_new+= "_" + v_lumiYears[iyear];
        						histoname_new+= "__" + sample_groups[isample];
        						if(syst_list[isyst] != "") {histoname_new+= "__" + Get_Modified_SystName(syst_list[isyst], v_lumiYears[iyear]);}
        						else if(systTree_list[itree] != "") {histoname_new+= "__" + systTree_list[itree];}

        						f->cd();
        						h_merging->Write(histoname_new, TObject::kOverwrite);

        						// cout<<"- Writing merged histo "<<histoname_new<<" with integral "<<h_merging->Integral()<<endl;

        						delete h_merging; h_merging = 0;
        					} //write histo
        				} //sample loop
        			} //syst loop
        		} //tree loop
        	} //channel loop
        } //years loop
    } //vars loop

	f->Close();

	cout<<endl<<FYEL("... Done")<<endl<<endl<<endl;

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

//FIXME -- check weight, ...
void TopEFT_analysis::Fill_TH1EFT(TH1EFT*& h, float x, vector<string>* v_ids, vector<float>* v_wgts, vector<float> v_SWE, float wgt_nominal)
{
    bool debug = false;

    WCFit eft_fit("myfit");

    //May only loop on minimal required number of points for WCFit (depends on n.of WCs) -- will get fit warning otherwise
    //=== WARNING : ASSUMES THAT FIRST WEIGHT is 'rwgt_1', and second is 'rwgt_sm' ===
    for(int iwgt=2; iwgt<25; iwgt++) //just the necessary nof weights to overconstrain fit with 5 WCs
    // for(int iwgt=2; iwgt<v_ids->size(); iwgt++)
    {
        WCPoint wc_pt(v_ids->at(iwgt), wgt_nominal*v_wgts->at(iwgt)/v_wgts->at(1));
        eft_fit.points.push_back(wc_pt);
    }

    //-- Include 'manually' the SM point as first element (not included automatically because named 'SM' and not via its operator values)
    eft_fit.points.insert(eft_fit.points.begin(), eft_fit.points[0]); //Duplicate the first element
    eft_fit.points[0].setSMPoint(); //Set (new) first element to SM coeffs (all null)
    eft_fit.points[0].wgt = v_wgts->at(1); //Set (new) first element to SM weight

    eft_fit.fitPoints();

    if(debug) //Printout WC values, compare true weight to corresponding fit result
    {
        cout<<"//-------------------------------------------- "<<endl;
        eft_fit.dump(); //Printout all names and coefficients of WC pairs
        for (uint i=0; i < eft_fit.points.size(); i++)
        {
            WCPoint wc_pt = eft_fit.points.at(i);
            double fit_val = eft_fit.evalPoint(&wc_pt);
            wc_pt.dump(); //Printout names and values of all WCs for this point
            std::cout << "===> " << std::setw(3) << i << ": " << std::setw(12) << wc_pt.wgt << " | " << std::setw(12) << fit_val << " | " << std::setw(12) << (wc_pt.wgt-fit_val) << std::endl; //Printout : i / true weight / evaluated weight / diff
        }
        cout<<endl<<endl<<endl;
    }

    //Fill TH1EFT with SM weight by default
    h->Fill(x, wgt_nominal , eft_fit);

    return;
}

void TopEFT_analysis::Test_TH1EFT()
{
    cout<<FMAG("Enter Test_TH1EFT() ")<<endl;

    TString filename = "/home/ntonon/Postdoc/EFTSimulations/myAnalysis/outputs/Templates_DNN_SR_2017.root";
    if(!Check_File_Existence(filename) ) {cout<<BOLD(FRED("File not found ("<<filename<<")"))<<endl; return;}
    TFile* f = TFile::Open(filename, "READ");

    TString th1eft_name = "TH1EFT_DNN0_SR_2017__PrivMC_tZq";
    if(!f->GetListOfKeys()->Contains(th1eft_name)) {cout<<BOLD(FRED("Histo not found ("<<th1eft_name<<")"))<<endl; return;}
    TH1EFT* h = (TH1EFT*) f->Get(th1eft_name);

    // cout<<"h->IsA()->InheritsFrom(TH1EFT::Class()) "<<h->IsA()->InheritsFrom(TH1EFT::Class())<<endl;
    // cout<<"h->GetEntries() "<<h->GetEntries()<<endl;
    // cout<<"h->GetNbinsX() "<<h->GetNbinsX()<<endl;

    h->DumpFits();

    WCFit fit = h->GetSumFit();
    vector<string> vname = fit.getNames();
    vector<std::pair<int,int>> vpair = fit.getPairs();
    vector<double> vcoeff = fit.getCoefficients();
    for(int i=0; i<vname.size(); i++)
    {
        cout<<"vname[i] "<<vname[i]<<endl;
    }
    for(int i=0; i<vpair.size(); i++)
    {
        cout<<"vpair[i] first / second "<<vpair[i].first<<" / "<<vpair[i].second<<endl;
    }
    for(int i=0; i<vcoeff.size(); i++)
    {
        cout<<"vcoeff[i] "<<vcoeff[i]<<endl<<endl;
    }

    // WCPoint wcp("smpt", 1.);
    // cout<<"fit.evalPoint(SM) "<<fit.evalPoint(&wcp)<<endl;
    // wcp = WCPoint("rwgt_ctZ_2.0_ctW_2.0_cpQM_2.0_cpQ3_2.0_cpt_2.0", 1.);
    // cout<<"fit.evalPoint(xx) "<<fit.evalPoint(&wcp)<<endl;
    // wcp = WCPoint("rwgt_ctZ_3.0_ctW_3.0_cpQM_3.0_cpQ3_3.0_cpt_3.0", 1.);
    // cout<<"fit.evalPoint(xx) "<<fit.evalPoint(&wcp)<<endl;

    return;
}