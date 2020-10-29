#include "Helper.h"

using namespace std;

//-- HARDCODE THE BASE DIR. CONTAINING PRIVATE NTUPLES (used by all codes)
TString NTUPLEDIR = "./input_ntuples/"; //LOCAL
// TString NTUPLEDIR = "/nfs/dust/cms/user/ntonon/CMSSW_10_2_20/src/potato_nicolas/potato-nicolas/nicolas/output/Analyzer3l-V10-AllSamples-d20201023-t151834/merged_ntuples/"; //CMSSW


//--------------------------------------------
// ##        #######  ##      ##         ##       ######## ##     ## ######## ##
// ##       ##     ## ##  ##  ##         ##       ##       ##     ## ##       ##
// ##       ##     ## ##  ##  ##         ##       ##       ##     ## ##       ##
// ##       ##     ## ##  ##  ## ####### ##       ######   ##     ## ######   ##
// ##       ##     ## ##  ##  ##         ##       ##        ##   ##  ##       ##
// ##       ##     ## ##  ##  ##         ##       ##         ## ##   ##       ##
// ########  #######   ###  ###          ######## ########    ###    ######## ########

// ##     ## ######## ##       ########  ######## ########
// ##     ## ##       ##       ##     ## ##       ##     ##
// ##     ## ##       ##       ##     ## ##       ##     ##
// ######### ######   ##       ########  ######   ########
// ##     ## ##       ##       ##        ##       ##   ##
// ##     ## ##       ##       ##        ##       ##    ##
// ##     ## ######## ######## ##        ######## ##     ##

// ######## ##     ## ##    ##  ######   ######
// ##       ##     ## ###   ## ##    ## ##    ##
// ##       ##     ## ####  ## ##       ##
// ######   ##     ## ## ## ## ##        ######
// ##       ##     ## ##  #### ##             ##
// ##       ##     ## ##   ### ##    ## ##    ##
// ##        #######  ##    ##  ######   ######
//--------------------------------------------

//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Existence(const TString& name, bool printError_missingFile)
{
    struct stat buffer;
    bool found = (stat (name.Data(), &buffer) == 0); //true if file exists
    if(!found && printError_missingFile) {cout<<endl<<BRED(BOLD(FWHT("ERROR : file "<<name<<" not found !")))<<endl<<endl;}
    return found;
}

//Move file with bash command 'mv'
int MoveFile(TString origin_path, TString dest_path)
{
	TString command = "mv "+origin_path + " " + dest_path;

	return system(command.Data() );
}

//Copy file with bash command 'cp'
int CopyFile(TString origin_path, TString dest_path)
{
	TString command = "cp "+origin_path + " " + dest_path;

	return system(command.Data() );;
}

//Convert a double into a TString
//@ precision --> can choose if TString how many digits the TString should display
//NB: alternatively, can simply use TString::Form(format,var). Ex: Form("%d", (int) a), Form("%f", (float) f), ...
TString Convert_Number_To_TString(double number, int precision/*=3*/)
{
	stringstream ss;
	ss << std::setprecision(precision) << number;
	TString ts = ss.str();
	return ts;
}

//Convert a TString into a float
float Convert_TString_To_Number(TString ts)
{
	float number = 0;
	string s = ts.Data();
	stringstream ss(s);

    //Bad characters may produce errors. Read 'out' these chars
    while(ss >> number || !ss.eof())
    {
        if(ss.fail()) //Found bad char. Remove it, read the remaining string
        {
            ss.clear();
            char dummy;
            ss >> dummy;
            continue;
        }
    }

	return number;
}

//Find a number into a TString, and returns it as a float
float Find_Number_In_TString(TString ts)
{
	TString tmp = ""; int number = 0;
	string s = ts.Data(); int ts_size = s.size(); //Only to get TString size

	for (int i=0; i < ts_size; i++)
	{
		if( isdigit(ts[i]) )
		{
			do
			{
				tmp += ts[i];
				i++;
			} while(isdigit(ts[i]) || ts[i] == '.'); //NB : Pay attention to quotes : "a" creates a 2-char array (letter+terminator) -> string. 'a' identifies a single character !

		  break;
		}
	}

	return Convert_TString_To_Number(tmp);
}

//Translates signs into words
TString Convert_Sign_To_Word(TString ts_in)
{
	TString sign = "";

	if(ts_in.Contains("<=")) {sign = "MaxEq";}
	else if(ts_in.Contains(">=")) {sign = "MinEq";}
	else if(ts_in.Contains(">")) {sign = "Min";}
	else if(ts_in.Contains("<")) {sign = "Max";}
	else if(ts_in.Contains("==")) {sign = "Eq";}
	else if(ts_in.Contains("!=")) {sign = "Diff";}

	return sign;
}

//Used when a "cut" is composed of 2 conditions -> breaks it into a pair of TStrings <cut1,cut2>
pair<TString,TString> Break_Cuts_In_Two(TString multiple_cut)
{
	TString cut1 = "", cut2 = "";

	string s = multiple_cut.Data(); int size = s.size(); //To get TString size

	for (int i=0; i < size; i++) //Extract condition 1
	{
		if(multiple_cut[i] == '&' || multiple_cut[i] == '|') {break;} //stop when we find '&&' or '||'
		if(multiple_cut[i] != ' ') {cut1 += multiple_cut[i];} //else, store characters in TString
	}

	TString tmp = "";
	for (int i=0; i < size; i++) //Extract condition 2
	{
		tmp+= multiple_cut[i];
		if( (!tmp.Contains("&&") && !tmp.Contains("||")) || multiple_cut[i] == '&' || multiple_cut[i] == '|') {continue;}
		if(multiple_cut[i] != ' ') {cut2 += multiple_cut[i];}
	}

	pair<TString,TString> the_cuts; the_cuts.first = cut1; the_cuts.second = cut2;

	return the_cuts;
}

//Check if a given value passes a given cut definition
bool Is_Event_Passing_Cut(TString cut_def, double value)
{
	double cut_tmp;

	if(!cut_def.Contains("&&") && !cut_def.Contains("||")) //If cut contains only 1 condition
	{
		cut_tmp = Find_Number_In_TString(cut_def);
		if(cut_def.Contains(">=") && value < cut_tmp)		 {return false;}
		else if(cut_def.Contains("<=") && value > cut_tmp)	 {return false;}
		else if(cut_def.Contains(">") && !cut_def.Contains(">=") && value <= cut_tmp)	 {return false;}
		else if(cut_def.Contains("<") && !cut_def.Contains("<=") && value >= cut_tmp) 	 {return false;}
		else if(cut_def.Contains("==") && value != cut_tmp)  {return false;}
		else if(cut_def.Contains("!=") && value == cut_tmp)  {return false;}
	}
	else if(cut_def.Contains("&&") && cut_def.Contains("||"))
	{
		cout<<BOLD(FRED("ERROR ! Cut definition too complex ; break it !"))<<endl;
	}
	else if(cut_def.Contains("&&") )//If '&&' in the cut def, break it in 2
	{
		TString cut1 = Break_Cuts_In_Two(cut_def).first;
		TString cut2 = Break_Cuts_In_Two(cut_def).second;
		//CUT 1
		cut_tmp = Find_Number_In_TString(cut1);
		if(cut1.Contains(">=") && value < cut_tmp)			 {return false;}
		else if(cut1.Contains("<=") && value > cut_tmp)		 {return false;}
		else if(cut1.Contains(">") && value <= cut_tmp) 	{return false;}
		else if(cut1.Contains("<") && value >= cut_tmp) 	{return false;}
		else if(cut1.Contains("==") && value != cut_tmp) 	 {return false;}
		else if(cut1.Contains("!=") && value == cut_tmp) 	 {return false;}
		//CUT 2
		cut_tmp = Find_Number_In_TString(cut2);
		if(cut2.Contains(">=") && value < cut_tmp)			 {return false;}
		else if(cut2.Contains("<=") && value > cut_tmp)		 {return false;}
		else if(cut2.Contains(">") && value <= cut_tmp)	{return false;}
		else if(cut2.Contains("<") && value >= cut_tmp)  {return false;}
		else if(cut2.Contains("==") && value != cut_tmp) 	 {return false;}
		else if(cut2.Contains("!=") && value == cut_tmp) 	 {return false;}
	}
	else if(cut_def.Contains("||") )//If '||' in the cut def, break it in 2
	{
		TString cut1 = Break_Cuts_In_Two(cut_def).first;
		TString cut2 = Break_Cuts_In_Two(cut_def).second;

		bool pass_cut1 = true; bool pass_cut2 = true; //Need to pass at least 1 cut

		//CUT 1
		cut_tmp = Find_Number_In_TString(cut1);
		if(cut1.Contains(">=") && value < cut_tmp)			 {pass_cut1 = false;}
		else if(cut1.Contains("<=") && value > cut_tmp)		 {pass_cut1 = false;}
		else if(cut1.Contains(">") && value <= cut_tmp)	{pass_cut1 = false;}
		else if(cut1.Contains("<") && value >= cut_tmp) 	{pass_cut1 = false;}
		else if(cut1.Contains("==") && value != cut_tmp) 	 {pass_cut1 = false;}
		else if(cut1.Contains("!=") && value == cut_tmp) 	 {pass_cut1 = false;}
		//CUT 2
		cut_tmp = Find_Number_In_TString(cut2);
		if(cut2.Contains(">=") && value < cut_tmp)			 {pass_cut2 = false;}
		else if(cut2.Contains("<=") && value > cut_tmp)		 {pass_cut2 = false;}
		else if(cut2.Contains(">") && value <= cut_tmp)		 	{pass_cut2 = false;}
		else if(cut2.Contains("<") && value >= cut_tmp) 	{pass_cut2 = false;}
		else if(cut2.Contains("==") && value != cut_tmp) 	 {pass_cut2 = false;}
		else if(cut2.Contains("!=") && value == cut_tmp) 	 {pass_cut2 = false;}

		if(!pass_cut1 && !pass_cut2) {return false;}
	}

	//If no condition was matched, it means the values passes the cut
	return true;
}

//Rescales the input value so that it retains its relative position in the new scale
double Compute_RangeScaled_Value(double value, double min_old, double max_old, double min_new, double max_new)
{
	return ((value - min_old) / (max_old - min_old)) * (max_new - min_new) + min_new;
}

//Returns the nof std deviations w.r.t. central value of the input value
double Compute_StdDevScaled_Value(double value, double mean, double dev)
{
	return (value - mean) / dev;
}

//Pass the names of the dir contents as references
bool Get_Dir_Content(string dir, vector<TString> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL)
    {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return false;
    }

    while ((dirp = readdir(dp)) != NULL)
    {
        files.push_back(string(dirp->d_name));
    }

    closedir(dp);

    return true;
}

//Split a long TString in smaller TStrings which are all separated by a common delimiter
//Only to remember the logic, no useful return for now (could return the array, but then loose control of the dynamic pointer...)
TString Split_TString_Into_Keys(TString ts, TString delim)
{
    TObjArray *tx = ts.Tokenize("_");
    // tx->Print();

    for(Int_t i = 0; i < tx->GetEntries(); i++)
    {
        TString tmp = ((TObjString *)(tx->At(i)))->String();
        std::cout << tmp << std::endl;
    }

    delete tx;

    return ts;
}

/**
 * Count the number of histograms found in a TFile
 * NB: very slow for huge files !
 */
int Count_nofHistos_inTFile(TFile* f)
{
    TKey* key = NULL;
    Int_t total = 0;

    TIter next((TList*) f->GetListOfKeys());
    while(key = (TKey*)next())
    {
        TClass* cl = gROOT->GetClass(key->GetClassName());
        if(cl->InheritsFrom("TH1"))
        {
            TH1* h = (TH1*) key->ReadObj();
            total++;
            // cout << "Histo found: " << h->GetName() << " - " << h->GetTitle() << endl;
        }
    }

    // cout << "Found " << total << " histograms" << endl;
    return total;
}


/**
 * Overloaded function
 */
int Count_nofHistos_inTFile(TString fname)
{
    TFile *f = TFile::Open(fname, "READ"); //Open
    if(!f || f->IsZombie()) {cout << "Unable to open " << fname << " for reading..." <<endl; return -1;}

    int total = Count_nofHistos_inTFile(f); //Count

    f->Close(); //Close

    return total; //Return count
}

/**
 * Round float to bottom/upper value, with precision determined by 'roundLevel' (e.g. 0.1)
 */
float Round_Float_toDecimal(float x, float roundLevel, bool roundDown)
{
    if(roundDown) {return int(x/roundLevel)*roundLevel;}
    return int(x/roundLevel + 1)*roundLevel;; //Round up
}








//--------------------------------------------
// ########     ###     ######  ####  ######
// ##     ##   ## ##   ##    ##  ##  ##    ##
// ##     ##  ##   ##  ##        ##  ##
// ########  ##     ##  ######   ##  ##
// ##     ## #########       ##  ##  ##
// ##     ## ##     ## ##    ##  ##  ##    ##
// ########  ##     ##  ######  ####  ######

// ##     ## ######## ##       ########  ######## ########
// ##     ## ##       ##       ##     ## ##       ##     ##
// ##     ## ##       ##       ##     ## ##       ##     ##
// ######### ######   ##       ########  ######   ########
// ##     ## ##       ##       ##        ##       ##   ##
// ##     ## ##       ##       ##        ##       ##    ##
// ##     ## ######## ######## ##        ######## ##     ##

// ######## ##     ## ##    ##  ######
// ##       ##     ## ###   ## ##    ##
// ##       ##     ## ####  ## ##
// ######   ##     ## ## ## ## ##
// ##       ##     ## ##  #### ##
// ##       ##     ## ##   ### ##    ##
// ##        #######  ##    ##  ######
//--------------------------------------------

//Increment weight of last bin by 'weight'
void Fill_Last_Bin_TH1F(TH1F* h, double weight)
{
	h->SetBinContent(h->GetNbinsX(), h->GetBinContent(h->GetNbinsX()) + weight);

	return;
}

//Increment weight of first bin by 'weight'
void Fill_First_Bin_TH1F(TH1F* h, double weight)
{
	h->SetBinContent(1, h->GetBinContent(1) + weight);

	return;
}

//Enforce desired canvas style
void Load_Canvas_Style()
{
	// For the canvas:
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(0); // must be kWhite but I dunno how to do that in PyROOT
	gStyle->SetCanvasDefH(600); //Height of canvas
	gStyle->SetCanvasDefW(600); //Width of canvas
	gStyle->SetCanvasDefX(0);   //POsition on screen
	gStyle->SetCanvasDefY(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0); // kWhite
	gStyle->SetPadGridX(0); //false
	gStyle->SetPadGridY(0); //false
	gStyle->SetGridColor(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(1);
	gStyle->SetFrameFillColor(0);
	gStyle->SetFrameFillStyle(0);
	gStyle->SetFrameLineColor(1);
	gStyle->SetFrameLineStyle(1);
	gStyle->SetFrameLineWidth(1);
	gStyle->SetHistLineColor(1);
	gStyle->SetHistLineStyle(0);
	gStyle->SetHistLineWidth(1);
	gStyle->SetEndErrorSize(2);
	gStyle->SetOptFit(1011);
	gStyle->SetFitFormat("5.4g");
	gStyle->SetFuncColor(2);
	gStyle->SetFuncStyle(1);
	gStyle->SetFuncWidth(1);
	gStyle->SetOptDate(0);
	gStyle->SetOptFile(0);
	gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
	gStyle->SetStatColor(0); // kWhite
	gStyle->SetStatFont(42);
	gStyle->SetStatFontSize(0.04);
	gStyle->SetStatTextColor(1);
	gStyle->SetStatFormat("6.4g");
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.15);
	gStyle->SetPadTopMargin(0.07);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadRightMargin(0.03);
	gStyle->SetOptTitle(0);
	gStyle->SetTitleFont(42);
	gStyle->SetTitleColor(1);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetTitleColor(1, "XYZ");
	gStyle->SetTitleFont(42, "XYZ");
	gStyle->SetTitleSize(0.06, "XYZ");
	gStyle->SetTitleXOffset(0.9);
	gStyle->SetTitleYOffset(1.25);
	gStyle->SetLabelColor(1, "XYZ");
	gStyle->SetLabelFont(42, "XYZ");
	gStyle->SetLabelOffset(0.007, "XYZ");
	gStyle->SetLabelSize(0.05, "XYZ");
	gStyle->SetAxisColor(1, "XYZ");
	gStyle->SetStripDecimals(1); // kTRUE
	gStyle->SetTickLength(0.03, "XYZ");
	gStyle->SetNdivisions(510, "XYZ");
	gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	gStyle->SetPadTickY(1);
	gStyle->SetOptLogx(0);
	gStyle->SetOptLogy(0);
	gStyle->SetOptLogz(0);
	gStyle->SetPaperSize(20.,20.);
}

//Enforce the Combine name conventions
TString Combine_Naming_Convention(TString name)
{
	if(name.Contains("Up") || name.Contains("Down") ) {return name;} //If systematics is already following Combine convention

	else if(name.Contains("__plus"))
	{
		int i = name.Index("__plus"); //Find index of substring
		name.Remove(i); //Remove substring
		name+= "Up"; //Add Combine syst. suffix
	}

	else if(name.Contains("__minus"))
	{
		int i = name.Index("__minus"); //Find index of substring
		name.Remove(i); //Remove substring
		name+= "Down"; //Add Combine syst. suffix
	}

	return name;
}

//Modifies the TMVA standard output to extract only the useful ranking information (hard-coded way)
void Extract_Ranking_Info(TString TMVA_output_file, TString channel)
{
	ifstream file_in(TMVA_output_file.Data() );

	ofstream file_out("ranking_tmp.txt");

	file_out<<"--- "<<channel<<" : Variable Ranking ---"<<endl<<endl;

	string string_start = "Ranking input variables (method specific)...";
	string string_stop = "----------------------------------------------";

	string line;

	while(!file_in.eof( ) )
	{
		getline(file_in, line);

		if(line.find(string_start) != std::string::npos) {break;} //When 'string_start' is found, start copying lines into output file
	}

	int k=0;
	while(!file_in.eof( ) )
	{
		getline(file_in, line);

		if(k>3) //We know that the first 3 lines are useless
		{
			if(line.find(string_stop) != std::string::npos) {break;} //Don't stop until all variables are read

			TString string_to_write = line;
			int index = string_to_write.First(":"); //Find first occurence of ':' in the lines, and remove the text before it (useless)
		    string_to_write.Remove(0,index);

		    string s = string_to_write.Data(); int size = s.size(); //Compute size of the TString this way (ugly)
		    for(int i=0; i<size; i++)
		    {
		        if(isdigit(string_to_write[i]) ) //Only copy what is after the ranking index, remove the rest
		        {
		            string_to_write.Remove(0,i);
		            break;
		        }
		    }

			file_out<<string_to_write.Data()<<endl;
		}

		k++;
	}

	TString mv_command = "mv ranking_tmp.txt " + TMVA_output_file;
	int tmp = system(mv_command.Data() ); //Replace TMVA output textfile by modified textfile !

	return;
}

//Read modified TMVA ranking file in order to fill 2 vectors (passed by reference) ordered by decreasing rank (var names & importances)
//NB : Hard-coded !
void Get_Ranking_Vectors(TString channel, vector<TString> &v_BDTvar_name, vector<double> &v_BDTvar_importance)
{
	TString file_rank_path = "./outputs/Rankings/RANKING_"+channel+".txt";
	if(!Check_File_Existence(file_rank_path.Data()) ) {std::cout<<file_rank_path<<" not found -- Abort"<<endl; return;}
	ifstream file_in(file_rank_path);

	string line;
	getline(file_in, line); getline(file_in, line);

	while(!file_in.eof())
	{
		getline(file_in, line);

		TString ts = line; //Convert to TString for easy modification

		int index = ts.First(":"); //Find first occurence of ':' in the lines, and remove the text before it (useless)
		ts.Remove(0,index+1); //Remove first useless characters

		TString var_name = "";
		string ts_tmp = ts.Data(); int size = ts_tmp.size();
		for(int ichar=0; ichar<size; ichar++)
		{
			if(var_name != "" && ts[ichar] == ' ') {break;} //end of var name
			else if(ts[ichar] == ' ') {continue;} //Begin at var name

			var_name+=ts[ichar];
		}

		v_BDTvar_name.push_back(var_name);


		index = ts.First(":"); //Find other occurence of ':' (placed before variable importance number)
		ts.Remove(0,index+1); //After that, only the importance number remains in the TString

 		TString var_importance_string = "";
 		ts_tmp = ts.Data(); size = ts_tmp.size();
		for(int ichar=0; ichar<size; ichar++)
		{
			if(ts[ichar] == ' ') {continue;}
			var_importance_string+=ts[ichar];
		}

		//Convert tstring to number (importance of the var)
		double var_importance = 0;
		if(var_importance_string.Contains("e-01") )
		{
			index = var_importance_string.Index('e');
			var_importance_string.Remove(index);
			var_importance = Convert_TString_To_Number(var_importance_string);
			var_importance*= 0.1;
		}
		else if(var_importance_string.Contains("e-02") )
		{
			index = var_importance_string.Index('e');
			var_importance_string.Remove(index);
			var_importance = Convert_TString_To_Number(var_importance_string);
			var_importance*= 0.01;
		}

		v_BDTvar_importance.push_back(var_importance);
	}


	//Because while loop on 'eof' --> one line too many at the end --> erase last entry
	v_BDTvar_name.erase(v_BDTvar_name.begin() + v_BDTvar_name.size()-1);
	v_BDTvar_importance.erase(v_BDTvar_importance.begin() + v_BDTvar_importance.size()-1);

	return;
}

//Loof for differences between 2 different histos stored in 2 different TFiles
void Compare_Histograms(TString filepath1, TString filepath2, TString histname1, TString histname2)
{
    TFile* f1 = TFile::Open(filepath1, "READ");
    TFile* f2 = TFile::Open(filepath2, "READ");

    if(!f1->GetListOfKeys()->Contains(histname1) ) {cout<<FRED("Error : "<<histname1<<" not found in file "<<filepath1<<"")<<endl; return;}
    if(!f2->GetListOfKeys()->Contains(histname2) ) {cout<<FRED("Error : "<<histname2<<" not found in file "<<filepath2<<"")<<endl; return;}

	TH1F* h1 = (TH1F*) f1->Get(histname1);
	TH1F* h2 = (TH1F*) f2->Get(histname2);

	for(int ibin=0; ibin<h1->GetNbinsX(); ibin++)
	{
		double diff = fabs((h1->GetBinContent(ibin+1) - h2->GetBinContent(ibin+1))) / h1->GetBinContent(ibin+1);

		if(diff > 1./100)
		{
			cout<<"WARNING : bin "<<ibin+1<<" different ! ("<<h1->GetBinContent(ibin+1)<<" / "<<h2->GetBinContent(ibin+1)<<")"<<endl;
		}
	}

	delete h1;
	delete h2;

	f1->Close();
	f2->Close();

	return;
}

//Rescale input value, depending on the mean and variance of this variable
float Rescale_Input_Variable(float value, float val1, float val2)
{
    return (value * val2) + val1; //sklearn.preprocessing.MinMaxScaler (-1;+1)

    // return (value - val1) / sqrt(val2); //sklearn.preprocessing.StandardScaler
}

//Create WCFit object by fitting a sufficient number of WCPoint objects, corresponding to EFT points.
void Get_WCFit(WCFit*& eft_fit, vector<string>* v_ids, vector<float>* v_wgts, const vector<float>& v_SWE, float wgt_baseline, float weightMENominal, float w_SMpoint, int idx_sm, int nweights)
{
    bool debug = false;


    vector<WCPoint> v_wcpts; //Vector storing the WCPoint objects, to be passed to the WCFit object
    // eft_fit->clear(); //Re-initialize object

    for(int iwgt=1; iwgt<nweights; iwgt++) //Only loop on the minimal nof weights required in my case to overconstrain fit //Skip first point, expected to correspond to SM in my private samples (treated separately)
    {
        if(ToLower(v_ids->at(iwgt)) == "rwgt_sm") {continue;} //SM point added manually

        //-- Formula for a given point P: w = (baseline_event_weight / w_MEnominal) * reweight_P / SWE_baseline
        //-- NB: account for 'baseline_event_weight / w_MEnominal' to account for scale factors (included in the default event weight 'baseline_event_weight')
        float w = wgt_baseline * v_wgts->at(iwgt) / (weightMENominal * v_SWE[idx_sm]);
        // float w = wgt_nominal * v_wgts->at(iwgt) / v_SWE[iwgt];

        //-- Printouts
        // cout<<"wgt_baseline "<<wgt_baseline<<endl;
        // cout<<"v_wgts->at(iwgt) "<<v_wgts->at(iwgt)<<endl;
        // cout<<"weightMENominal "<<weightMENominal<<endl;
        // cout<<"v_SWE[idx_sm] "<<v_SWE[idx_sm]<<endl;
        // cout<<"w "<<w<<endl;

        //-- Create WCPoint and copy it into WCFit object
        // WCPoint wc_pt(v_ids->at(iwgt), w); //Create WCPoint object
        // eft_fit->points.push_back(wc_pt); //Add WCPoint to WCFit object directly

        v_wcpts.push_back(WCPoint(v_ids->at(iwgt), w));
    }

    //-- Include 'manually' the SM point as first element (not included automatically because named 'SM' and not via its operator values)
    if(idx_sm >= 0)
    {
        // eft_fit->points.insert(eft_fit->points.begin(), eft_fit->points[0]); //Duplicate the first element
        // eft_fit->points[0].setSMPoint(); //Set (new) first element to SM coeffs (all null)
        // eft_fit->points[0].wgt = w_SMpoint; //Set (new) first element to SM weight

        v_wcpts.insert(v_wcpts.begin(), v_wcpts[0]); //Duplicate the first element
        v_wcpts[0].setSMPoint(); //Set (new) first element to SM coeffs (all null)
        v_wcpts[0].wgt = w_SMpoint; //Set (new) first element to SM weight
    }

    // eft_fit->fitPoints(); //Fit the WCPoint objects
    eft_fit->fitPoints(v_wcpts); //Fit the WCPoint objects

    if(debug) //Printout WC values, compare true weight to corresponding fit result
    {
        cout<<"//-------------------------------------------- "<<endl;
        eft_fit->dump(); //Printout all names and coefficients of WC pairs
        for (uint i=0; i < eft_fit->points.size(); i++)
        {
            WCPoint wc_pt = eft_fit->points.at(i);
            double fit_val = eft_fit->evalPoint(&wc_pt);
            wc_pt.dump(); //Printout names and values of all WCs for this point
            std::cout << "===> " << std::setw(3) << i << ": " << std::setw(12) << wc_pt.wgt << " | " << std::setw(12) << fit_val << " | " << std::setw(12) << (wc_pt.wgt-fit_val) << std::endl; //Printout : i / true weight / evaluated weight / diff
        }
        cout<<endl<<endl<<endl;
    }

    return;
}

//Get mirror histo (rescaled to nominal if necessary)
//As done by Benjamin/ttH : https://github.com/stiegerb/cmgtools-lite/blob/80X_M17_tHqJan30_bbcombination/TTHAnalysis/python/plotter/makeShapeCardsTHQ.py#L609-L633
//Directly modify the histos passed in args
void Get_Mirror_Histogram(TH1F*& hcentral, TH1F*& hvariation, TH1F*& hmirror, bool is_ShapeOnly)
{
	hmirror = 0; //empty at beginning
	hmirror = (TH1F*) hcentral->Clone(); //Get binning, etc.

	//Rescale variation to nominal (shape-only)
    if(is_ShapeOnly) {hvariation->Scale(hcentral->Integral() / hvariation->Integral());}

	for(int ibin=1; ibin<=hcentral->GetNbinsX(); ibin++)
	{
		double yA = hvariation->GetBinContent(ibin); //'up = alternate'
		double y0 = hcentral->GetBinContent(ibin); //'central'
		double yM = y0; //'down = mirror'

		if(y0 > 0 && yA > 0) {yM = y0 * y0/yA;}
		else if(yA == 0) {yM = 2*y0;}

		hmirror->SetBinContent(ibin, yM);
	}

    if(is_ShapeOnly) {hmirror->Scale(hcentral->Integral() / hmirror->Integral());} //Rescale mirror to nominal (shape-only)
	else //Mirror normalization, cf github
	{
    	double mirror_norm = pow(hcentral->Integral(), 2) / hvariation->Integral();
		hmirror->Scale(mirror_norm / hvariation->Integral());
	}

	return;
}

//Get mirror histo, as done by ttH2017 :  https://github.com/CERN-PH-CMG/cmgtools-lite/blob/759bdc4213c50db48cb695ae498f7a97794a1410/TTHAnalysis/python/plotter/uncertaintyFile.py#L106
//Directly modify the hdown histo
void Get_TemplateSymm_Histogram(TH1F*& hcentral, TH1F*& hvariation, TH1F*& hmirror, bool is_ShapeOnly)
{
	hmirror = 0; //empty at beginning
	hmirror = (TH1F*) hcentral->Clone();

    if(hcentral->Integral() <= 0 || hvariation->Integral() <= 0) {return;}

	//Protection -- check if I artificially set the bin to a float ~null distribution
	//If that is the case, just do mirror = central
	if(hcentral->Integral() == hcentral->GetNbinsX() * hcentral->GetBinContent(1) || hvariation->Integral() == hvariation->GetNbinsX() * hvariation->GetBinContent(1))
	{
		return;
	}

	hmirror->Multiply(hmirror);
	hmirror->Divide(hvariation);

	// cout<<"hmirror->Integral() "<<hmirror->Integral()<<endl;

	return;
}

/**
 * In FCNC analysis, we want to reuse the closure shape systematics from ttH, even though they were derived in different regions
 * Thus we want to inflate these systematics by 20%. Can be done in datacard parser for norm syst. For shape syst, check the different (syst - nominal) in each bin, and inflate by 20%
 */
void Inflate_Syst_inShapeTemplate(TH1F*& h_toInflate, TH1F* h_nominal, float inflation_factor)
{
	if(!h_toInflate || !h_nominal || !h_toInflate->GetEntries() || !h_nominal->GetEntries()) {cout<<"Error in Inflate_Syst_inShapeTemplate() !"<<endl; return;}

	for(int ibin=1; ibin<h_toInflate->GetNbinsX()+1; ibin++)
	{
		//syst - nominal
		float diff_bin = h_toInflate->GetBinContent(ibin) - h_nominal->GetBinContent(ibin);

		h_toInflate->SetBinContent(ibin, h_toInflate->GetBinContent(ibin) + inflation_factor * diff_bin); //+- X%, depending on direction of variation
	}

	return;
}


//--------------------------------------------
//    ###    ##    ##    ###    ##       ##    ##  ######  ####  ######
//   ## ##   ###   ##   ## ##   ##        ##  ##  ##    ##  ##  ##    ##
//  ##   ##  ####  ##  ##   ##  ##         ####   ##        ##  ##
// ##     ## ## ## ## ##     ## ##          ##     ######   ##   ######
// ######### ##  #### ######### ##          ##          ##  ##        ##
// ##     ## ##   ### ##     ## ##          ##    ##    ##  ##  ##    ##
// ##     ## ##    ## ##     ## ########    ##     ######  ####  ######

//  ######  ########  ########  ######  #### ######## ####  ######
// ##    ## ##     ## ##       ##    ##  ##  ##        ##  ##    ##
// ##       ##     ## ##       ##        ##  ##        ##  ##
//  ######  ########  ######   ##        ##  ######    ##  ##
//       ## ##        ##       ##        ##  ##        ##  ##
// ##    ## ##        ##       ##    ##  ##  ##        ##  ##    ##
//  ######  ##        ########  ######  #### ##       ####  ######

// ######## ##     ## ##    ##  ######   ######
// ##       ##     ## ###   ## ##    ## ##    ##
// ##       ##     ## ####  ## ##       ##
// ######   ##     ## ## ## ## ##        ######
// ##       ##     ## ##  #### ##             ##
// ##       ##     ## ##   ### ##    ## ##    ##
// ##        #######  ##    ##  ######   ######
//--------------------------------------------

//Modifies arguments passed by reference according to command args
bool Apply_CommandArgs_Choices(int argc, char **argv, vector<TString>& v_lumiYear, TString& region_choice)
{
	if(argc >= 2)
	{
        TString arg1 = (TString) argv[1]; arg1.ToLower();

        if(arg1 == "2016") {v_lumiYear.resize(0); v_lumiYear.push_back("2016");}
        else if(arg1 == "2017") {v_lumiYear.resize(0); v_lumiYear.push_back("2017");}
        else if(arg1 == "2018") {v_lumiYear.resize(0); v_lumiYear.push_back("2018");}
        else if(arg1 == "201617") {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2017");}
        else if(arg1 == "201618") {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2018");}
        else if(arg1 == "201718") {v_lumiYear.resize(0); v_lumiYear.push_back("2017"); v_lumiYear.push_back("2018");}
        else if(arg1 == "run2") {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2017"); v_lumiYear.push_back("2018");}

        else if(arg1 == "tzq") {region_choice = "tZq";}
        else if(arg1 == "ttz") {region_choice = "ttZ";}
        else if(arg1 == "ttz4l") {region_choice = "ttZ4l";}
        else if(arg1 == "twz") {region_choice = "tWZ";}
        else if(arg1 == "signal") {region_choice = "signal";}
        else if(arg1 == "vg") {region_choice = "Vg";}
        else if(arg1 == "zz") {region_choice = "zz";}
        else if(arg1 == "tt") {region_choice = "tt";}
        else if(arg1 == "wz") {region_choice = "wz";}
        else if(arg1 == "dy") {region_choice = "dy";}
        else if(arg1 == "tx") {region_choice = "tX";} //Obsolete

        else
        {
            cout<<BOLD(FRED("Error : wrong arguments at execution !"))<<endl;
            cout<<"argc = "<<argc<<endl;
            cout<<"argv[1] = '"<<argv[1]<<"'"<<endl;
            if(argc >= 3) {cout<<"argv[2] = '"<<argv[2]<<"'"<<endl;}
            cout<<UNDL("--> Syntax : ./analysis_main.exe [lumi] [region]")<<endl;
            cout<<UNDL("(Please check conventions for 'lumi' arg. in code)")<<endl;

            return 0;
        }

		if(argc >= 3)
		{
            TString arg2 = (TString) argv[2]; arg2.ToLower();

            if(arg2 == "2016") {v_lumiYear.resize(0); v_lumiYear.push_back("2016");}
            else if(arg2 == "2017") {v_lumiYear.resize(0); v_lumiYear.push_back("2017");}
            else if(arg2 == "2018") {v_lumiYear.resize(0); v_lumiYear.push_back("2018");}
            else if(arg2 == "201617") {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2017");}
            else if(arg2 == "201618") {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2018");}
            else if(arg2 == "201718") {v_lumiYear.resize(0); v_lumiYear.push_back("2017"); v_lumiYear.push_back("2018");}
            else if(arg2 == "run2") {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2017"); v_lumiYear.push_back("2018");}

            else if(arg2 == "tzq") {region_choice = "tZq";}
            else if(arg2 == "ttz") {region_choice = "ttZ";}
            else if(arg2 == "ttz4l") {region_choice = "ttZ4l";}
            else if(arg2 == "twz") {region_choice = "tWZ";}
            else if(arg2 == "signal") {region_choice = "signal";}
            else if(arg2 == "vg") {region_choice = "Vg";}
            else if(arg2 == "zz") {region_choice = "zz";}
            else if(arg2 == "tt") {region_choice = "tt";}
            else if(arg2 == "wz") {region_choice = "wz";}
            else if(arg2 == "dy") {region_choice = "dy";}
            else if(arg2 == "tx") {region_choice = "tX";} //Obsolete

			else
			{
				cout<<BOLD(FRED("Error : wrong arguments at execution !"))<<endl;
				cout<<"argc = "<<argc<<endl;
				cout<<"argv[1] = '"<<argv[1]<<"'"<<endl;
				cout<<"argv[2] = '"<<argv[2]<<"'"<<endl;
                cout<<UNDL("--> Syntax : ./analysis_main.exe [lumi] [region]")<<endl;
                cout<<UNDL("(Please check conventions for 'lumi' arg. in code)")<<endl;

				return 0;
			}
		}
	}

	return 1;
}

//Select a color for each activated sample
//NB : also account for the DATA sample in color vector, but keep its element empty !
void Get_Samples_Colors(vector<int>& v_colors, std::vector<TColor*>& v_custom_colors, vector<TString> v_samples, vector<TString> v_groups, int color_scheme)
{
    //Define some TColors, which we can then use to define the default color vector. NB : must only delete them after all plots have been created.
    v_custom_colors.resize(10, 0);

    //-- Colors definition : Number + RGB value, from 0 to 1
    v_custom_colors[0] = new TColor(9000, 31./255., 120./255., 180./255.); //Dark blue
    v_custom_colors[1] = new TColor(9001, 166./255., 206./255., 227./255.); //Light blue
    v_custom_colors[2] = new TColor(9002, 51./255., 160./255., 44./255.); //Dark green 2
    v_custom_colors[3] = new TColor(9003, 178./255., 223./255., 138./255.); //Light green
    v_custom_colors[4] = new TColor(9004, 255./255., 129./255., 123./255.); //Salmon
    v_custom_colors[5] = new TColor(9005, 53./255., 184./255., 99./255.); //Mild green
    v_custom_colors[6] = new TColor(9006, 227./255., 65./255., 61./255.); //Dark red
    v_custom_colors[7] = new TColor(9007, 255./255., 210./255., 47./255.); //Mild yellow
    v_custom_colors[8] = new TColor(9008, 152./255., 78./255., 163./255.); //Dark purple

    //Also nice
    //kOrange-2,

//--------------------------------------------
    if(color_scheme == 0) //Custom color scheme
    {
        for(int isample=0; isample<v_groups.size(); isample++)
		{
            //Signals
            if(v_groups[isample] == "tZq" || v_groups[isample] == "PrivMC_tZq") {v_colors[isample] = v_custom_colors[6]->GetNumber();}
            else if(v_groups[isample] == "ttZ" || v_groups[isample] == "PrivMC_ttZ") {v_colors[isample] = v_custom_colors[0]->GetNumber();}
            else if(v_groups[isample] == "tWZ" || v_groups[isample] == "PrivMC_tWZ") {v_colors[isample] = v_custom_colors[7]->GetNumber();} //FIXME

            //t(t)X
            else if(v_groups[isample] == "tX") {v_colors[isample] = v_custom_colors[1]->GetNumber();}

            //WZ
            else if(v_groups[isample] == "WZ") {v_colors[isample] = v_custom_colors[3]->GetNumber();}

            //ZZ

            //VVV
            else if(v_groups[isample] == "VVV") {v_colors[isample] = v_custom_colors[5]->GetNumber();}
            // else if(v_groups[isample] == "ZZ") {v_colors[isample] = v_custom_colors[5]->GetNumber();}
            // else if(v_groups[isample] == "VVV") {v_colors[isample] = kViolet+2;}

            //X+g
            else if(v_groups[isample] == "Vg" || v_groups[isample] == "Xg" || v_groups[isample] == "VG" || v_groups[isample] == "XG") {v_colors[isample] = kPink+2;}

            //NPL
            else if(v_groups[isample] == "NPL") {v_colors[isample] = kGray;}
		}
    }

    else if(color_scheme == 1) //'Default' color scheme
	{
		for(int isample=0; isample<v_samples.size(); isample++)
		{
            //Signals
            if(v_samples[isample] == "tZq" || v_samples[isample] == "PrivMC_tZq") {v_colors[isample] = kOrange+10;}
            // if(v_samples[isample] == "tZq") {v_colors[isample] = test.GetNumber();}
            else if(v_samples[isample] == "ttZ" || v_samples[isample] == "PrivMC_ttZ") {v_colors[isample] = kOrange+6;}
            else if(v_samples[isample] == "tWZ" || v_samples[isample] == "PrivMC_tWZ") {v_colors[isample] = kOrange+2;}

            //ttX
            else if(v_samples[isample] == "ttH") {v_colors[isample] = kGreen-5;}
            else if(v_samples[isample] == "ttW") {v_colors[isample] = kGreen-5;}
            else if(v_samples[isample] == "tttt") {v_colors[isample] = kGreen-5;}
            else if(v_samples[isample] == "ttZZ") {v_colors[isample] = kGreen-5;}
            else if(v_samples[isample] == "ttWW") {v_colors[isample] = kGreen-5;}
            else if(v_samples[isample] == "ttWZ") {v_colors[isample] = kGreen-5;}
            else if(v_samples[isample] == "ttZH") {v_colors[isample] = kGreen-5;}
            else if(v_samples[isample] == "ttWH") {v_colors[isample] = kGreen-5;}

            //tX
            else if(v_samples[isample] == "tHq") {v_colors[isample] = kSpring+2;}
            else if(v_samples[isample] == "tHW") {v_colors[isample] = kSpring+2;}
            else if(v_samples[isample] == "ST") {v_colors[isample] = kSpring+2;}

            //VV(V)
            else if(v_samples[isample] == "WZ") {v_colors[isample] = kViolet-6;}
            else if(v_samples[isample] == "ZZ4l") {v_colors[isample] = kViolet-6;}
            else if(v_samples[isample] == "ZZZ") {v_colors[isample] = kViolet-6;}
            else if(v_samples[isample] == "WZZ") {v_colors[isample] = kViolet-6;}
            else if(v_samples[isample] == "WWW") {v_colors[isample] = kViolet-6;}
            else if(v_samples[isample] == "WWZ") {v_colors[isample] = kViolet-6;}
            else if(v_samples[isample] == "WZ2l2q") {v_colors[isample] = kViolet-6;}
            else if(v_samples[isample] == "ZZ2l2q") {v_colors[isample] = kViolet-6;}
            else if(v_samples[isample] == "ZG2l2g") {v_colors[isample] = kViolet-6;}

            //Fakes
            else if(v_samples[isample] == "DY") {v_colors[isample] = kGray;} //kAzure-7

            else if(v_samples[isample] == "TTbar_DiLep") {v_colors[isample] = kGray+2;} //kPink-4, kCyan-6;
            else if(v_samples[isample] == "TTbar_SemiLep") {v_colors[isample] = kGray+2;}
            else if(v_samples[isample].Contains("NPL")) {v_colors[isample] = kGray+2;}
		}
	}

    else if(color_scheme == 2) //David
	{
		for(int isample=0; isample<v_samples.size(); isample++)
		{
            //Signals
            if(v_samples[isample] == "tZq" || v_samples[isample] == "PrivMC_tZq") {v_colors[isample] = v_custom_colors[0]->GetNumber();}
            else if(v_samples[isample] == "ttZ" || v_samples[isample] == "PrivMC_ttZ") {v_colors[isample] = v_custom_colors[1]->GetNumber();}
            else if(v_samples[isample] == "tWZ") {v_colors[isample] = kOrange+2;}

            //ttX
            else if(v_samples[isample] == "ttH") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttW") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "tttt") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttZZ") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttWW") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttWZ") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttZH") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttWH") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttHH") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "tHq") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "tHW") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            // else if(v_samples[isample] == "tWZ") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ST") {v_colors[isample] = v_custom_colors[2]->GetNumber();}

            //WZ
            else if(v_samples[isample] == "WZ") {v_colors[isample] = kViolet-6;}

            //ZZ
            else if(v_samples[isample] == "ZZ4l") {v_colors[isample] = kTeal+2;}
            else if(v_samples[isample] == "ggToZZTo4l") {v_colors[isample] = kTeal+2;}

            //VVV
            else if(v_samples[isample] == "ZZZ") {v_colors[isample] = kGreen-3;}
            else if(v_samples[isample] == "WZZ") {v_colors[isample] = kGreen-3;}
            else if(v_samples[isample] == "WWW") {v_colors[isample] = kGreen-3;}
            else if(v_samples[isample] == "WWZ") {v_colors[isample] = kGreen-3;}
            else if(v_samples[isample] == "WZ2l2q") {v_colors[isample] = kGreen-3;}
            else if(v_samples[isample] == "ZZ2l2q") {v_colors[isample] = kGreen-3;}
            else if(v_samples[isample] == "ZG2l2g") {v_colors[isample] = kGreen-3;}

            //X+g
            else if(v_samples[isample] == "TTGamma_Dilep") {v_colors[isample] = kPink+2;}
            else if(v_samples[isample] == "tGJets") {v_colors[isample] = kPink+2;}
            else if(v_samples[isample] == "WGToLNuG") {v_colors[isample] = kPink+2;}
            else if(v_samples[isample] == "ZGToLLG_01J") {v_colors[isample] = kPink+2;}

            //Fakes
            else if(v_samples[isample] == "DY") {v_colors[isample] = kGray;}
            else if(v_samples[isample] == "TTbar_DiLep") {v_colors[isample] = kGray;}
            else if(v_samples[isample] == "TTbar_SemiLep") {v_colors[isample] = kGray;}
            else if(v_samples[isample].Contains("NPL")) {v_colors[isample] = kGray;}
		}
	}

	return;
}


/*
//Use custom color palette
//-- Idea : take good-looking/efficient color palettes from web and apply it manually
void Set_Custom_ColorPalette(vector<TColor*> &v_custom_colors, vector<int> &v, vector<TString> v_sampleGroups)
{ TColor* col = new TColor(1700, 141./255., 211./255., 199./255.);
    // col.SetRGB(141./255., 211./255., 199./255.);
    // cout<<col->GetNumber()<<endl;

    // Int_t dark   = TColor::GetColorDark(color_index);
    // Int_t bright = TColor::GetColorBright(color_index);

    // Int_t ci = TColor::GetFreeColorIndex(); //returns a free color index which can be used to define a user custom color.
    // TColor *color = new TColor(ci, 0.1, 0.2, 0.3);

	v.resize(10);
	v_custom_colors.resize(10, 0);

	//-- Colors definition : Number + RGB value, from 0 to 1
    v_custom_colors[0] = 0;
    v_custom_colors[1] = new TColor(9002, 31./255., 120./255., 180./255.);
    v_custom_colors[2] = new TColor(9001, 166./255., 206./255., 227./255.);
    v_custom_colors[3] = new TColor(9004, 51./255., 160./255., 44./255.);
    v_custom_colors[4] = new TColor(9003, 178./255., 223./255., 138./255.);
    // v_custom_colors[5] = new TColor(9005, 251./255., 251./255., 153./255.);
	// v_custom_colors[6] = new TColor(9006, 227./255., 26./255., 28./255.);

    // int index_customColor = 0;
    // for(int isample=0; isample<v_sampleGroups.size(); isample++)
    // {
    //     if(isample >=1 && v_sampleGroups[isample] != v_sampleGroups[isample-1]) {index_customColor++;}
    //     if(v_custom_colors[index_customColor] != 0) {v[isample] = v_custom_colors[index_customColor]->GetNumber();}
    // }

    return;
}
*/

//Store here the binnings and ranges of all the variables to be plotted via ControlHistograms files
bool Get_Variable_Range(TString var, int& nbins, double& xmin, double& xmax)
{
    // {nbins = 20; xmin = -1; xmax = 1; return true;}

    //Categories are either 0 or 1 (NB : encoded in Char_t!)
    if(var.BeginsWith("is_") || var.BeginsWith("passed") ) {nbins = 2; xmin = 0; xmax = 2;}

    else if(var == "metEt") {nbins = 20; xmin = 0; xmax = 250;}
	else if(var == "dEtaFwdJetBJet") {nbins = 20; xmin = 0; xmax = 4;}
	else if(var == "dEtaFwdJetClosestLep") {nbins = 20; xmin = 0; xmax = 4;}
	else if(var == "minDRll") {nbins = 15; xmin = 0; xmax = 3.5;}
    else if(var == "maxEtaJet") {nbins = 10; xmin = 0; xmax = 5;}
    else if(var == "forwardJetAbsEta") {nbins = 10; xmin = 0; xmax = 5;}
    else if(var == "jPrimeAbsEta") {nbins = 10; xmin = 0; xmax = 5;}
	else if(var == "dPhiHighestPtSSPair") {nbins = 20; xmin = 0; xmax = 3.2;}
	else if(var == "nJetEta1") {nbins = 7; xmin = 0.5; xmax = 7.5;}
	else if(var == "lepCharge") {nbins = 3; xmin = -1.5; xmax = 1.5;}
	else if(var == "hardestBjetPt") {nbins = 20; xmin = 20; xmax = 300;}
    else if(var == "hardestBjetEta") {nbins = 40; xmin = 0; xmax = 2.5;}
	else if(var == "FwdJetPt") {nbins = 20; xmin = 20; xmax = 200;}
	else if(var == "inv_mll") {nbins = 20; xmin = 50; xmax = 130;}
    else if(var == "mHT") {nbins = 10; xmin = 200; xmax = 1200;}
    else if(var == "top_mass") {nbins = 15; xmin = 100; xmax = 300;}
    else if(var == "dPhijj_max") {nbins = 10; xmin = 0; xmax = 3;}
    else if(var == "maxDijet_M") {nbins = 10; xmin = 100; xmax = 1000;}
    else if(var == "maxDiJet_Pt") {nbins = 20; xmin = 100; xmax = 400;}
    else if(var == "maxDijet_dPhi") {nbins = 20; xmin = 0.; xmax = 3.5;}
    else if(var == "maxDelPhiLL") {nbins = 10; xmin = 0; xmax = 3.5;}
    else if(var == "m3l" || var == "Mass_3l") {nbins = 20; xmin = 50; xmax = 500;}
    else if(var == "leptonCharge") {nbins = 3; xmin = -1.5; xmax = 1.5;}
    else if(var == "mTW") {nbins = 20; xmin = 0.; xmax = 200;}
    else if(var == "recoZ_Mass") {nbins = 20; xmin = 75.; xmax = 110;}
    else if(var == "lAsymmetry") {nbins = 20; xmin = -3.; xmax = 3.;}
    else if(var == "Mass_tZ") {nbins = 10; xmin = 200; xmax = 1000;}
    else if(var == "maxDeepCSV" || var == "maxDeepJet" || var == "deepCSV_2nd" || var == "deepJet_2nd") {nbins = 20; xmin = 0.3; xmax = 1.1;}
    else if(var == "maxDelRbL") {nbins = 20; xmin = 1.; xmax = 5;}
    else if(var == "Top_delRbl" || var == "Top_delRbW") {nbins = 20; xmin = 0.; xmax = 4.5;}
    else if(var == "channel") {nbins = 4; xmin = 0.; xmax = 4;}
    else if(var == "mbjMax") {nbins = 20; xmin = 0.; xmax = 8.;}
    else if(var == "recoZ_Eta") {nbins = 20; xmin = 0.; xmax = 1.;}
    else if(var == "recoLepTop_Eta") {nbins = 20; xmin = -3.; xmax = 3.;}
    else if(var == "recoLepTop_Pt") {nbins = 20; xmin = 0.; xmax = 300.;}
    else if(var == "dR_blW" || var == "dR_bW") {nbins = 20; xmin = 0.; xmax = 4.;}

    else if(var == "TopZsystem_M") {nbins = 20; xmin = 150.; xmax = 1000.;}
    else if(var == "jet1_pt") {nbins = 20; xmin = 50.; xmax = 600.;}
    else if(var == "lep1_pt") {nbins = 20; xmin = 50.; xmax = 300.;}
    else if(var == "recoZ_Pt") {nbins = 20; xmin = 0.; xmax = 250.;}
    else if(var == "recoZ_Eta") {nbins = 20; xmin = -3.; xmax = 3.;}

    else if(var.Contains("CSV", TString::kIgnoreCase)) {nbins = 20; xmin = 0.; xmax = 1.1;}
    else if(var.Contains("dR") || var.Contains("DelR", TString::kIgnoreCase) ) {nbins = 20; xmin = 0; xmax = 7.;}
    else if(var.Contains("deta", TString::kIgnoreCase) ) {nbins = 20; xmin = 0; xmax = 3.;}
    else if(var.Contains("eta", TString::kIgnoreCase) ) {nbins = 20; xmin = -3.; xmax = 3.;}
    else if(var.Contains("phi", TString::kIgnoreCase) ) {nbins = 20; xmin = 0; xmax = 3.;}
    else if(var.BeginsWith("cos", TString::kIgnoreCase)) {nbins = 20; xmin = -1.; xmax = 1.;}

    else if(var == "njets" || var == "nbjets") {nbins = 6; xmin = 0.; xmax = 6;}

    // else if(var == "Lep3Pt") {nbins = 20; xmin = 0; xmax = 80;}

    else {nbins = 20; xmin = -5.; xmax = 150.;}
    // else {return false;}

	return true;
}

//Return the binning of the variable passed in arg
//NB: use Contains() function because I am passing full variable names as arguments (e.g. 'NN_SRtZq')
void Get_Template_Range(int& nbins, float& xmin, float& xmax, TString template_name, bool use_SManalysis_strategy, bool make_SMvsEFT_templates_plots, int categorization_strategy, bool plot_onlyMaxNodeEvents, int& nbjets_min, int& nbjets_max, int& njets_min, int& njets_max)
{
    nbins = 15; //Default

    xmin = -1; xmax = 1; //BDT: [-1,1]

    if(template_name.Contains("NN")) //NN: [0,1]
    {
        xmin = 0;
        if(!make_SMvsEFT_templates_plots && categorization_strategy==2 && plot_onlyMaxNodeEvents) {xmin = 0.3; nbins = 14;} //Special case: if we plot SM vs SM multiclass NN and only plot events in their max. node, then by construction there can be no events with x<1/3 --> Adapt axis
    }
    else if(template_name.Contains("categ")) {nbins = (nbjets_max-nbjets_min+1)*(njets_max-njets_min+1); xmin = 0; xmax = nbins;} //1 bin per sub-category
    else if(template_name.Contains("ZptCos")) {nbins = 12; xmin = 0; xmax = 12;} //2D Zpt-cosThetaStarPolZ (as in TOP-18-009) //4bins in Zpt, 3 in cosTheta
    else if(template_name.Contains("Zpt")) {nbins = 5; xmin = 0; xmax = 400;} //1D Zpt //1 bin per sub-category

    if(use_SManalysis_strategy) {xmin = 0;}

    if(template_name.Contains("NN") || template_name.Contains("BDT")) {nbins = 10;}
    if(template_name.Contains("mTW"))
    {
        nbins=15; xmin=0; xmax=150;
        if(use_SManalysis_strategy) {nbins=10; xmax=150;} //Keep my binning for now
    }
    if(template_name.Contains("countExp")) {nbins=1; xmax=1;}
    if(template_name.Contains("channel")) {nbins=2; xmax=2;}

    return;
}

//Return beautified variable name
TString Get_Variable_Name(TString var)
{
    // l -> italic l, different font
    // add [GeV] units for some vars ?
    // 'l' for LateX curvy l ; but belongs to TMathText, and not compatible with TLatex '#' chars?
    //'#leftrightarrow' : <->
    //Big delimiters : #leftX, #rightX, with X=|,(,{,[ (and opposite for right)

    if(var == "maxDijetDelR") {return "max. #DeltaR(j,j)";}
    if(var == "dEtaFwdJetBJet") {return "#Delta#eta#left(j^{fwd},b#right)";}
    if(var == "dEtaFwdJetClosestLep") {return "#Delta#eta#left(j^{fwd},l^{closest}#right)";}
    if(var == "mHT") {return "m_{HT}";}
    if(var == "mTW") {return "m_{T}^{W}";}
    if(var == "Mass_3l") {return "m_{3l}";}
    if(var == "forwardJetAbsEta") {return "#left|#eta#left(j^{fwd}#right)#right|";}
    if(var == "jPrimeAbsEta") {return "#left|#eta(j')#right|";}
    if(var == "maxDeepCSV") {return "max. DeepCSV";}
    if(var == "deepCSV_2nd") {return "2nd max. DeepCSV";}
    if(var == "maxDeepJet") {return "max. DeepJet";}
    if(var == "deepJet_2nd") {return "2nd max. DeepJet";}
    if(var == "delRljPrime") {return "#DeltaR(j',l)";}
    if(var == "lAsymmetry") {return "q_{l} #upoint #left|#eta(l)#right|";}
    if(var == "maxDijetMass") {return "max. m_{j,j}";}
    if(var == "maxDelPhiLL") {return "max. #Delta#phi(l,l)";}
    if(var == "metEt") {return "E_{T}^{miss}";}
    if(var == "recoZ_Mass") {return "m_{Z}";}
    if(var == "Mass_tZ") {return "m_{tZ}";}
    if(var == "maxDijetPt") {return "max. p_{T}(j,j)";}
    if(var == "maxDijetDelPhi") {return "max #Delta#phi(j,j)";}
    if(var == "minDelRbL") {return "min. #DeltaR(b,l)";}
    if(var == "Top_delRbl") {return "#DeltaR#left(l^{t},b^{t}#right)";}
    if(var == "Top_delRbW") {return "#DeltaR#left(W^{t},b^{t}#right)";}
    if(var == "cosThetaStarPolTop") {return "cos#left(#theta^{* t}_{pol}#right)";}
    if(var == "cosThetaStarPolZ") {return "cos#left(#theta^{* Z}_{pol}#right)";}

    //GenPlotter variables
    if(var == "Z_pt") {return "p_{T}(Z)";}
    if(var == "Z_eta") {return "#eta(Z)";}
    if(var == "Z_m") {return "m(Z)";}
    if(var == "Zreco_m") {return "m(Z_{reco})";}
    if(var == "Top_pt") {return "p_{T}(t)";}
    if(var == "Top_eta") {return "#eta(t)";}
    if(var == "Top_m") {return "m(t)";}
    if(var == "TopZsystem_m") {return "m(tZ)";}
    if(var == "LeadingTop_pt") {return "p_{T}(t^{lead})";}
    if(var == "LeadingTop_eta") {return "#eta(t^{lead})";}
    if(var == "Zreco_dPhill") {return "#Delta#Phi_{ll}";}
    if(var == "cosThetaStarPol_Z") {return "cos(#theta^{*}_{Z})";}
    if(var == "cosThetaStarPol_Top") {return "cos(#theta^{*}_{t})";}

    return var;
}

//Get var name corresponding to category
//isFake <-> use different flags for prompt/fake events
//NB: always call this function to automate cuts on flags, since they differ for NPL sample, etc.
TString Get_Category_Boolean_Name(TString region, bool isFake)
{
    TString name = "";

    // if(region=="tZq") {return "is_tZq_3l_SR";}
    // else if(region=="ttZ") {return "is_ttZ_3l_SR";}
    // else if(region=="tWZ") {return "is_tWZ_3l_SR";}

    if(region=="signal") {name = "is_signal_SR";}
    else if(region=="tZq") {name = "is_tzq_SR";}
    else if(region=="ttZ") {name = "is_ttz_SR";}
    else if(region=="ttZ4l") {name = "is_ttz4l_SR";}

    else if(region=="Vg") {name = "is_Vg_CR";}
    else if(region=="zz") {name = "is_zz_CR";}
    else if(region=="tX") {name = "is_tX_CR";}
    else if(region=="tt") {name = "is_tt_CR";}
    else if(region=="wz") {name = "is_wz_CR";}
    else if(region=="dy") {name = "is_dy_CR";}

    else {cout<<"Warning: empty category name..."<<endl;}

    if(isFake && name != "") {name+= "Fake";}

    return name;
}

//Computes total nof entries which will be processed by the Produce_Templates() function, so that the Timebar is correct
//NB1 : don't account for nof syst weights, since they are computed all at once when reading an event
//Obsolete -- NB2 : this returns the nof entries which will get read, not processed (else, should take cuts into account, etc.)
float Count_Total_Nof_Entries(TString dir_ntuples, TString t_name_nominal, vector<TString> v_samples, vector<TString> v_systTrees, vector<TString> v_cut_name, vector<TString> v_cut_def, vector<TString> v_lumiYears, bool makeHisto_inputVars, bool noSysts_inputVars)
{
    float total_nentries = 0;

    if(makeHisto_inputVars && noSysts_inputVars) //syst weights not computed
    {
        v_systTrees.resize(1);
    }

    //--- Define the cut chain to apply to events
    TString cuts_chain = "";
    /*
    for(int ivar=0; ivar<v_cut_name.size(); ivar++)
    {
        if(v_cut_def[ivar] != "")
        {
            if(cuts_chain != "") {cuts_chain+= " && ";}

            if(!v_cut_def[ivar].Contains("&&") && !v_cut_def[ivar].Contains("||")) {cuts_chain+= v_cut_name[ivar] + v_cut_def[ivar];} //If cut contains only 1 condition
            else if(v_cut_def[ivar].Contains("&&") && v_cut_def[ivar].Contains("||")) {cout<<BOLD(FRED("ERROR ! Wrong cut definition !"))<<endl;}
            else if(v_cut_def[ivar].Contains("&&") )//If '&&' in the cut, break it in 2
            {
                cuts_chain+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).first;
                cuts_chain+= " && ";
                cuts_chain+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).second;
            }
            else if(v_cut_def[ivar].Contains("||") )//If '||' in the cut, break it in 2
            {
                cuts_chain+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).first;
                cuts_chain+= " || ";
                cuts_chain+= v_cut_name[ivar] + Break_Cuts_In_Two(v_cut_def[ivar]).second;
            }
        }
    }*/

    for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
    {
        TString dir_tmp = dir_ntuples + v_lumiYears[iyear] + "/";

        for(int isample=0; isample<v_samples.size(); isample++)
        {
            // cout<<"v_samples[isample] "<<v_samples[isample]<<endl;

            //Open input TFile
            TString inputfile = dir_tmp + v_samples[isample] + ".root";
            // cout<<"inputfile "<<inputfile<<endl;
            if(!Check_File_Existence(inputfile))
            {
                // cout<<endl<<"File "<<inputfile<<FRED(" not found!")<<endl;
                continue;
            }
            TFile* file_input = TFile::Open(inputfile, "READ");

            for(int itree=0; itree<v_systTrees.size(); itree++)
    		{
                TTree* tree = 0;
                TString tmp = v_systTrees[itree];
    			if(tmp == "") {tmp = t_name_nominal;}
                tree = (TTree*) file_input->Get(tmp);
    			if(!tree)
    			{
    				// cout<<BOLD(FRED("ERROR : tree '"<<tmp<<"' not found in file : "<<inputfile<<" ! Skip !"))<<endl;
    				continue; //Skip sample
    			}

                // float nentries_tmp = tree->GetEntries("passedBJets");
                float nentries_tmp = tree->GetEntries(cuts_chain);

                total_nentries+= nentries_tmp; //Multiply nof entries by nof loops

                // cout<<"total_nentries = "<<total_nentries<<endl;
            } //trees
        } //samples
    } //years

    // cout<<"total_nentries "<<total_nentries<<endl;

    return total_nentries;
}

//-- For systs uncorrelated between years, need to add the year in the systName (to make it unique)
//-- NB: for first arg. can pass the concatenation of 'syst_list[isyst]+systTree_list[itree]' because only 1 can be non-empty at a time !
//-- NB: could also use the 'nuisance edit rename process channel oldname newname' directive directly in Combine datacards
TString Get_Modified_SystName(TString systname, TString lumiYear, TString samplename)
{
    //-- Only list the systematics which are *not* to be correlated in-between years
    //Need to add year as suffix (<-> unique)
    if(systname.BeginsWith("BtagHFstats")
    || systname.BeginsWith("BtagLFstats")
    || systname.BeginsWith("BtagCF")
    || systname.BeginsWith("JER")
    || systname.BeginsWith("UnclEn")
    || systname.BeginsWith("MET") //New name
    )
    {
        if(systname.EndsWith("Up"))
        {
            int i = systname.Index("Up"); //Find position of suffix
            systname.Remove(i); //Remove suffix
            systname+= lumiYear+"Up"; //Add year + suffix
        }
        else if(systname.EndsWith("Down"))
        {
            int i = systname.Index("Down"); //Find position of suffix
            systname.Remove(i); //Remove suffix
            systname+= lumiYear+"Down"; //Add year + suffix
        }
    }
    //-- Only list the systematics which are *not* to be correlated in-between samples
    //Add sample name as suffix <-> unique
    else if(systname.BeginsWith("ISR")
    || systname.BeginsWith("ME"))
    {
        if(systname.EndsWith("Up"))
        {
            int i = systname.Index("Up"); //Find position of suffix
            systname.Remove(i); //Remove suffix
            systname+= samplename+"Up"; //Add year + suffix
        }
        else if(systname.EndsWith("Down"))
        {
            int i = systname.Index("Down"); //Find position of suffix
            systname.Remove(i); //Remove suffix
            systname+= samplename+"Down"; //Add year + suffix
        }
    }

    return systname;
}

//Get the address of the relevant GEN-level histogram, so that it can be superimposed on data/MC comparison plots
void Get_Pointer_GENHisto(TH1F*& h, TString variable)
{
    TString input_file_name = "../Pheno/Plotter/GENhisto.root";
    if(!Check_File_Existence(input_file_name)) {cout<<BOLD(FRED("Input histo file "<<input_file_name<<" not found ! Can't superimpose GEN-level EFT histo..."))<<endl; return;}

    TFile* f = TFile::Open(input_file_name);

    TString hname = "Z_pt_tzq_sm"; //-- hardcoded for testing
    if(!f->GetListOfKeys()->Contains(hname) ) {cout<<ITAL("Histogram '"<<hname<<"' : not found...")<<endl; return;}
    h = (TH1F*) f->Get(hname);

    return;
}

//Parse the name of an EFT reweight point to extract the names and strengths of each operator
vector<pair<TString,float>> Parse_EFTreweight_ID(TString id)
{
    vector<pair<TString,float>> v;

    std::vector<std::string> words;
    split_string((string) id, words, "_");
    // if(words.at(0).compare("rwgt") != 0) {cout<<"Error: EFT ID should start with rwgt_ prefix ! "<<endl; return v;}
    if(words.at(0).compare("rwgt") != 0 && words.at(0).find("EFTrwgt") == std::string::npos) {cout<<"Error: EFT ID should start with rwgt_ prefix ! "<<endl; return v;} //Add TOP19001 convention

    for (uint i=1; i < words.size(); i+= 2)
    {
        if(i+1 >= words.size()) {std::cout<<"Error : i > words.size() ! Full string : "<<id<<std::endl; return v;}

        pair<TString, float> pair_tmp;
        pair_tmp.first = words.at(i);
        pair_tmp.second = stod(words.at(i+1));

        v.push_back(pair_tmp);
    }

    return v;
}

//Get the arbitrary x-axis value corresponding to the jet and bjet multiplicities of the event
//Logic: check if njets<=min, >=max, or ==x for any other multiplicity in-between to determine event bin
float Get_x_jetCategory(float njets, float nbjets, int nbjets_min, int nbjets_max, int njets_min, int njets_max)
{
    // int nbjets_min = 1, nbjets_max = 2;
    // int njets_min = 2, njets_max = 6;

    int bin = 0; //will return bin x-axis value, assuming fix width of 1

    bool stop_loop = false;

    {
        for(int ibjet=nbjets_min; ibjet<=nbjets_max; ibjet++) //Loop on possible nbjets values
        {
            if(stop_loop) {break;} //Bin found

            for(int ijet=njets_min; ijet<=njets_max; ijet++) //Loop on possible njets values
            {
                if(nbjets==ibjet || (nbjets<=nbjets_min && ibjet==nbjets_min) || (nbjets>=nbjets_max && ibjet==nbjets_max)) //nbjets categ found
                {
                    if(njets==ijet || (njets<=njets_min && ijet==njets_min) || (njets>=njets_max && ijet==njets_max)) //njets categ found
                    {
                        stop_loop = true; break; //Bin found
                    }
                }

                bin++; //Bin not found, increment
            }
        }
    }
    if(!stop_loop) {cout<<BOLD(FRED("Warning: jet category not found ! Return -1"))<<endl; return -1.;} //Category not found

    // cout<<"njets "<<njets<<" / nbjets "<<nbjets<<" --> categ "<<bin + 0.5<<endl;

    return bin + 0.5; //Assume bin width of 1. --> add event in center of bin
}

//Get x center position of bin corresponding to the Zpt and cosThetaStarPolZ values given in arg (following TOP-18-009 conventions)
float Get_x_ZptCosCategory(float Zpt, float cosThetaStarPolZ)
{
    if(Zpt < 100)
    {
        if(cosThetaStarPolZ < -0.6) {return 1-0.5;} //Corresponds to SR1 --> Fill at 0.5 (first bin)
        else if(cosThetaStarPolZ < 0.6) {return 2-0.5;}
        else {return 3-0.5;}
    }
    else if(Zpt < 200)
    {
        if(cosThetaStarPolZ < -0.6) {return 4-0.5;} //Corresponds to SR1 --> Fill at 0.5 (first bin)
        else if(cosThetaStarPolZ < 0.6) {return 5-0.5;}
        else {return 6-0.5;}
    }
    else if(Zpt < 400)
    {
        if(cosThetaStarPolZ < -0.6) {return 7-0.5;} //Corresponds to SR1 --> Fill at 0.5 (first bin)
        else if(cosThetaStarPolZ < 0.6) {return 8-0.5;}
        else {return 9-0.5;}
    }
    else
    {
        if(cosThetaStarPolZ < -0.6) {return 10-0.5;} //Corresponds to SR1 --> Fill at 0.5 (first bin)
        else if(cosThetaStarPolZ < 0.6) {return 11-0.5;}
        else {return 12-0.5;}
    }

    cout<<FRED("ERROR: relevant bin not found in Get_x_ZptCosCategory() ! Zpt = "<<Zpt<<" / cosThetaStarPolZ = "<<cosThetaStarPolZ<<"")<<endl;

    return -1;
}

//Get the path of the relevant MVA input file (.xml for BDT, .pb for NN), depending on specific analysis options. Intended for use in Produce_Templates() function
TString Get_MVAFile_InputPath(TString MVA_type, TString signal_process, TString year, bool use_specificMVA_eachYear, bool MVA_EFT, bool load_NN_info, int categorization_strategy)
{
    if(MVA_type != "BDT" && MVA_type != "NN") {cout<<BOLD(FRED("ERROR: wrong [MVA_type] argument !"))<<endl; return "";}
    // if(signal_process != "tZq" && signal_process != "ttZ" && signal_process != "Multiclass") {cout<<BOLD(FRED("ERROR: wrong [signal_process] argument !"))<<endl; return "";}
    if(signal_process != "tZq" && signal_process != "ttZ" && signal_process != "Multiclass") {cout<<"Argument [signal_process] not recognized... Setting it to 'Other'"<<endl; signal_process = "Other";}

    TString fullpath = ""; //Full input path
    TString basepath = "./weightsMVA/" + MVA_type + "/"; //Base path
    TString path_suffix = ""; //Suffix to append to base path; in-between the 2 we also need the year information, which we may want to modify if the nominal file is not found (see below)

    if(MVA_type == "BDT") {path_suffix = "BDT_" + signal_process + ".weights.xml";}
    else
    {
        if(categorization_strategy==0) {path_suffix = "tmp/";} //TEST dir.
        else //SM or EFT
        {
            if(!MVA_EFT) //SM vs SM
            {
                path_suffix = "SM/";
                if(categorization_strategy==1) {path_suffix+= signal_process + "/";} //Read MVA-tZq or MVA-ttZ
                else {path_suffix+= "Multiclass/";} //Read MVA_multiclass
            }
            else {path_suffix = "EFT/" + signal_process + "/";} //SM vs EFT
        }
        if(load_NN_info) {path_suffix+= "NN_info.txt";} //Read NN info file (input features, etc.)
        else {path_suffix+= "model.pb";} //Load NN model
    }

    if(use_specificMVA_eachYear)
    {
        fullpath = basepath + year + "/" + path_suffix;
        cout<<DIM("Trying to open MVA file "<<fullpath<<" ... ");

        if(!Check_File_Existence(fullpath)) //If did not find year-specific file, also try to look for full Run 2 file
        {
            fullpath = basepath + "Run2/" + path_suffix; //Try to read Run2 file (i.e. not year-specific) file instead
            if(year == "Run2" || !Check_File_Existence(fullpath)) {cout<<BOLD(FRED("ERROR: no MVA input file found !"))<<endl; return "";}
            else {cout<<DIM("NOT FOUND !")<<endl; cout<<FBLU("--> Using file ")<<fullpath<<FBLU(" instead !")<<endl; usleep(3000000);}
        }
        else {cout<<DIM("FOUND !")<<endl;}
    }
    else
    {
        fullpath = basepath + "Run2/" + path_suffix;
        cout<<DIM("Trying to open MVA file "<<fullpath<<" ... ");
        if(Check_File_Existence(fullpath)) {cout<<DIM("FOUND !")<<endl;}
        else {cout<<BOLD(FRED("ERROR: no MVA input file found !"))<<endl; return "";}
    }

    // cout<<"fullpath "<<fullpath<<endl;
    return fullpath;
}

//Get the path of the file containing the relevant histograms (either templates or input variables), depending on specific analysis options. Intended for use in Draw_Templates() function
TString Get_HistoFile_InputPath(bool is_templateFile, TString template_type, TString region, TString year, bool use_CombineFile, TString filename_suffix, bool MVA_EFT, int categorization_strategy, bool make_fixedRegions_templates, bool parametrized)
{
    TString fullpath = ""; //Full input path
    if(region != "") {region = "_" + region;}

	if(use_CombineFile) //Reading a template file produced by Combine (specific naming conventions)
	{
        //TRY 1 : look for Combine file
		fullpath = "./outputs/fitDiagnostics_";
		fullpath+= template_type + region + filename_suffix + ".root";

        cout<<DIM("Trying to open input file "<<fullpath<<" ... ");
        if(Check_File_Existence(fullpath)) {cout<<DIM("FOUND !")<<endl; return fullpath;}
		else
        {
            fullpath = "./outputs/fitDiagnostics.root"; //Try a generic name

            cout<<DIM("Trying to open input file "<<fullpath<<" ... ");
            if(Check_File_Existence(fullpath)) {cout<<DIM("FOUND !")<<endl; return fullpath;}
            else {cout<<BOLD(FRED("ERROR: fitDiagnostics file from Combine not found !"))<<endl; return "";}
        }
    }
	else //Reading my own file
	{
        TString MVA_type = "";
        TString region_tmp = region;
        if(is_templateFile && categorization_strategy>0 && !make_fixedRegions_templates)
        {
            region_tmp = "";
            if(MVA_EFT)
            {
                MVA_type+= "EFT" + Convert_Number_To_TString(categorization_strategy);
                if(parametrized) {MVA_type+= "param";}
            }
            else {MVA_type+= "SM";}
        }
        if(MVA_type != "") {MVA_type = "_" + MVA_type;}

		if(!is_templateFile) {fullpath = "outputs/ControlHistograms" + region + "_" + year + filename_suffix + ".root";} //Input variables
		else {fullpath = "outputs/Templates_" + template_type + MVA_type + region_tmp + "_" + year + filename_suffix + ".root";} //Templates

        cout<<DIM("Trying to open input file "<<fullpath<<" ... ");
        if(Check_File_Existence(fullpath)) {cout<<DIM("FOUND !")<<endl; return fullpath;}
        else
        {
            cout<<BOLD(FRED("ERROR: file "<<fullpath<<" not found !"))<<endl;
            if(year != "Run2") //Retry with Run2 file
            {
                if(!is_templateFile) {fullpath = "outputs/ControlHistograms" + region + "_Run2" + filename_suffix + ".root";} //Input variables
                else {fullpath = "outputs/Templates_" + template_type + MVA_type + region_tmp + "_Run2" + filename_suffix + ".root";} //Templates
                cout<<DIM("Trying to open input file "<<fullpath<<" ... ");
                if(Check_File_Existence(fullpath)) {cout<<DIM("FOUND !")<<endl; return fullpath;}
                else
                {
                    cout<<BOLD(FRED("ERROR: file "<<fullpath<<" not found ! Can not plot any histogram !"))<<endl;
                    return "";
                }
            }
            else {return "";}
        }
    }

    // cout<<"fullpath "<<fullpath<<endl;
    return fullpath;
}

//Read a NN info file and fill values via argument
bool Extract_Values_From_NNInfoFile(TString NNinfo_input_path, vector<TString>& var_list_NN, vector<TString>& v_NN_nodeLabels, TString& NN_inputLayerName, TString& NN_outputLayerName, int& NN_iMaxNode, int& NN_nNodes, vector<float>& minmax_bounds, TString* NN_strategy)
{
    cout<<DIM("Retrieving MVA information from : "<<NNinfo_input_path<<"")<<endl;
    ifstream file_in(NNinfo_input_path);
    string line;
    while(!file_in.eof())
    {
        getline(file_in, line);
        // TString ts_line(line);
        stringstream ss(line);
        TString varname; float tmp1, tmp2; //Values tmp1 and tmp2 could be the mean and variance, or min and max, etc. depending on the rescaling
        ss >> varname >> tmp1 >> tmp2;
        // cout<<" varname "<<varname<<" tmp1 "<<tmp1<<" tmp2 "<<tmp2<<endl;
        if(varname != "") //Last line may be empty
        {
            if(tmp1 == -1 && tmp2 == -1) {NN_inputLayerName = varname;} //Name of input layer
            else if(tmp1 == -2 && tmp2 == -2) {NN_outputLayerName = varname;} //Name of output layer
            else if(tmp1 == -3 && tmp2 == -3) {NN_nNodes = Convert_TString_To_Number(varname);} //Number of output nodes
            else if(tmp1 == -4 && tmp2 == -4) {v_NN_nodeLabels.push_back(varname);} //Label(s) of the node(s), e.g. 'tZq'/'ttZ'/'Backgrounds'
            else if(NN_strategy && tmp1 == -5 && tmp2 == -5) {*NN_strategy = varname;} //NB: pass this argument as a pointer instead of reference so it can have a default value --> Must deference first
            else if(varname=="bounds")
            {
                minmax_bounds.push_back(tmp1); minmax_bounds.push_back(tmp2);
            }
            // else if(varname=="bounds") {minmax_bounds = std::make_pair(tmp1, tmp2);}
            else if(tmp1 != -5 && tmp2 != -5) //Input variables
            {
                var_list_NN.push_back(varname);
                // std::pair <float,float> pair_tmp = std::make_pair(tmp1, tmp2);
            }
        }
    }
    // cout<<DIM("-->  "<<NN_inputLayerName<<"")<<endl;
    // cout<<DIM("-->  "<<NN_outputLayerName<<"")<<endl;
    // cout<<DIM("-->  "<<NN_nNodes<<" nodes")<<endl;

    if(NN_inputLayerName == "" || NN_outputLayerName == "" || NN_nNodes == -1 || !minmax_bounds.size()) {cout<<endl<<FRED("Warning : NN input/output info not found !")<<endl; return false;} //Need this info for NN

    return true;
}

TString Get_Region_Label(TString region, TString variable)
{
    TString label = "";

    if(region=="signal") {label = "SR";}
    else if(region=="tZq") {label = "tZq SR";}
    else if(region=="ttZ") {label = "ttZ SR";}
    else if(region=="ttZ4l") {label = "ttZ 4l SR";}

    else if(region=="Vg") {label = "V#gamma CR";}
    else if(region=="zz") {label = "ZZ CR";}
    else if(region=="tX") {label = "tX CR";}
    else if(region=="tt") {label = "t#bar{t} CR";}
    else if(region=="wz") {label = "WZ CR";}
    else if(region=="dy") {label = "DY CR";}

    if(variable.Contains("SRtZq")) {label = "SR-tZq";}
    else if(variable.Contains("SRttZ")) {label = "SR-ttZ";}
    else if(variable.Contains("_SRother")) {label = "SR-Other";}

    return label;
}

//Set the list of variables based on arguments
void Fill_Variables_List(vector<TString>& variable_list, bool use_predefined_EFT_strategy, TString template_name, TString region, bool scanOperators_paramNN, int NN_nNodes, bool make_SMvsEFT_templates_plots, TString operator_scan1, TString operator_scan2, vector<float> v_WCs_operator_scan1, vector<float> v_WCs_operator_scan2, bool use_SManalysis_strategy, bool make_fixedRegions_templates)
{
    variable_list.push_back(template_name);

    //-- If using predefined SM vs SM categ. strategy, consider the variable in 3 hard-coded regions
    if(use_predefined_EFT_strategy)
    {
        vector<TString> var_list_tmp(variable_list); //Tmp copy of variable list
        variable_list.clear();
        for(int ivar=0; ivar<var_list_tmp.size(); ivar++)
        {
            if(region == "signal" || region == "tZq") variable_list.push_back(var_list_tmp[ivar] + "_SRtZq"); //Replaced: var_list_tmp[ivar] + "_" + MVA_type + "_SRtZq"
            if(region == "signal" || region == "ttZ") variable_list.push_back(var_list_tmp[ivar] + "_SRttZ");
            // if(region == "signal" && !scanOperators_paramNN) //Only include CR for non-parametrized MVA (simpler)
            if(region == "signal")
            {
                if(make_SMvsEFT_templates_plots) {variable_list.push_back("mTW_SRother");} //For SM vs EFT in CR, use mTW distribution for now
                else {variable_list.push_back(var_list_tmp[ivar] + "_SRother");}
            }
        }

        //-- EFT scan: for parameterized NN templates, need to consider 1 set of variable *per EFT point*
        if(scanOperators_paramNN)
        {
            vector<TString> var_list_tmp(variable_list); //Tmp copy of variable list
            variable_list.clear(); //Reset the actual list
            for(int ivar=0; ivar<var_list_tmp.size(); ivar++) //For each variable initially found in the list, duplicate it once for each considered EFT point
            {
                for(int iop1=0; iop1<v_WCs_operator_scan1.size(); iop1++)
                {
                    TString opname1 = operator_scan1 + "_" + Convert_Number_To_TString(v_WCs_operator_scan1[iop1]); //Default: scan 1 operator
                    TString opname2 = ""; //Optional: scan 2nd operator
                    for(int iop2=0; iop2<v_WCs_operator_scan2.size(); iop2++)
                    {
                        if(operator_scan2=="" && iop2>0) {break;}
                        else if(operator_scan2!="") {opname2 = "_" + operator_scan2 + "_" + Convert_Number_To_TString(v_WCs_operator_scan2[iop2]);}
                        variable_list.push_back(var_list_tmp[ivar] + "_" + opname1 + opname2);
                        // cout<<"Added variable: "<<variable_list[variable_list.size()-1]<<endl;
                    }
                }
            }
        } //MVA EFT scan
    } //Predefined strategy
    else if(!use_predefined_EFT_strategy && template_name == "NN")
    {
        variable_list.clear();
        for(int inode=0; inode<NN_nNodes; inode++) {variable_list.push_back("NN" + (NN_nNodes == 1? "" : Convert_Number_To_TString(inode)));} //1 template per output node
    }

    if(use_SManalysis_strategy) //Hard-coded, following SM analysis strategy
    {
        variable_list.clear();
        if(region == "tZq" || region == "ttZ") {variable_list.push_back("NN");}
        else if(region == "wz") {variable_list.push_back("mTW");}
        else if(region == "zz") {variable_list.push_back("countExp");}
        else if(region == "Vg") {variable_list.push_back("channel");}
        template_name = variable_list[0];
    }
    else if(make_fixedRegions_templates) //NB: order is important !
    {
        variable_list.clear();
        variable_list.push_back("countExp_SRttZ4l");
        variable_list.push_back("mTW_CRWZ");
        variable_list.push_back("countExp_CRZZ");
        variable_list.push_back("countExp_CRDY");
        template_name = "";
    }

    return;
}
