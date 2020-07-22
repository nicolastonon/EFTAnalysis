#include "Helper.h"

using namespace std;

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
bool Check_File_Existence(const TString& name)
{
  struct stat buffer;
  return (stat (name.Data(), &buffer) == 0); //true if file exists
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
void Get_WCFit(WCFit*& eft_fit, vector<string>* v_ids, vector<float>* v_wgts, const vector<float>& v_SWE, float wgt_baseline, float weightMENominal, float w_SMpoint, int idx_sm)
{
    bool debug = false;

    for(int iwgt=1; iwgt<25; iwgt++) //Only loop on the minimal nof weights required in my case to overconstrain fit //Skip first point, expected to correspond to SM (treated separately)
    // for(int iwgt=1; iwgt<v_ids->size(); iwgt++)
    {
        //Formula for a given point P: w = (baseline_event_weight / w_MEnominal) * reweight_P / SWE_baseline
        //NB: account for 'baseline_event_weight / w_MEnominal' to account for scale factors (included in the default event weight 'baseline_event_weight')
        float w = wgt_baseline * v_wgts->at(iwgt) / (weightMENominal * v_SWE[idx_sm]);

        // cout<<"wgt_baseline "<<wgt_baseline<<endl;
        // cout<<"v_wgts->at(iwgt) "<<v_wgts->at(iwgt)<<endl;
        // cout<<"weightMENominal "<<weightMENominal<<endl;
        // cout<<"v_SWE[idx_sm] "<<v_SWE[idx_sm]<<endl;
        // cout<<"w "<<w<<endl;

        if(ToLower(v_ids->at(iwgt)) == "rwgt_sm") {continue;} //SM point added manually

        // double w = wgt_nominal * v_wgts->at(iwgt) / v_SWE[iwgt];

        WCPoint wc_pt(v_ids->at(iwgt), w); //Create WCPoint object

        eft_fit->points.push_back(wc_pt); //Add WCPoint to WCFit object
    }

    //-- Include 'manually' the SM point as first element (not included automatically because named 'SM' and not via its operator values)
    eft_fit->points.insert(eft_fit->points.begin(), eft_fit->points[0]); //Duplicate the first element
    eft_fit->points[0].setSMPoint(); //Set (new) first element to SM coeffs (all null)
    // eft_fit->points[0].wgt = v_wgts->at(idx_sm); //Set (new) first element to SM weight
    eft_fit->points[0].wgt = w_SMpoint; //Set (new) first element to SM weight

    eft_fit->fitPoints(); //Fit the WCPoint objects

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
        if(!strcmp(argv[1],"2016") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2016");}
        else if(!strcmp(argv[1],"2017") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2017");}
        else if(!strcmp(argv[1],"2018") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2018");}
        else if(!strcmp(argv[1],"201617") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2017");}
        else if(!strcmp(argv[1],"201618") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2018");}
        else if(!strcmp(argv[1],"201718") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2017"); v_lumiYear.push_back("2018");}
        else if(!strcmp(argv[1],"Run2") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2017"); v_lumiYear.push_back("2018");}

        else if(!strcmp(argv[1],"tZq")) {region_choice = "tZq";}
        else if(!strcmp(argv[1],"ttZ")) {region_choice = "ttZ";}
        else if(!strcmp(argv[1],"tWZ")) {region_choice = "tWZ";}

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
            if(!strcmp(argv[2],"tZq")) {region_choice = "tZq";}
            else if(!strcmp(argv[2],"ttZ")) {region_choice = "ttZ";}
            else if(!strcmp(argv[2],"tWZ")) {region_choice = "tWZ";}

            else if(!strcmp(argv[2],"2016") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2016");}
            else if(!strcmp(argv[2],"2017") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2017");}
            else if(!strcmp(argv[2],"2018") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2018");}
            else if(!strcmp(argv[2],"201617") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2017");}
            else if(!strcmp(argv[2],"201618") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2018");}
            else if(!strcmp(argv[2],"201718") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2017"); v_lumiYear.push_back("2018");}
            else if(!strcmp(argv[2],"Run2") ) {v_lumiYear.resize(0); v_lumiYear.push_back("2016"); v_lumiYear.push_back("2017"); v_lumiYear.push_back("2018");}

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
void Get_Samples_Colors(vector<int>& v_colors, std::vector<TColor*>& v_custom_colors, vector<TString> v_samples, int color_scheme, bool use_custom_colorPalette)
{
	// int i_skipData = 0; //index //include data in color vector for now (empty element)

    //Define some TColors, which we can then use to define the default color vector. NB : must only delete them after all plots have been created.
    if(use_custom_colorPalette)
    {
        v_custom_colors.resize(10, 0);

        //-- Colors definition : Number + RGB value, from 0 to 1
        v_custom_colors[0] = new TColor(9000, 31./255., 120./255., 180./255.);
        v_custom_colors[1] = new TColor(9001, 166./255., 206./255., 227./255.);
        v_custom_colors[2] = new TColor(9003, 51./255., 160./255., 44./255.);
        v_custom_colors[3] = new TColor(9004, 178./255., 223./255., 138./255.);

    }

//--------------------------------------------
    if(use_custom_colorPalette) //Custom color scheme
    {
        for(int isample=0; isample<v_samples.size(); isample++)
		{
            //Signals
            if(v_samples[isample] == "tZq" || v_samples[isample] == "PrivMC_tZq") {v_colors[isample] = v_custom_colors[0]->GetNumber();}
            else if(v_samples[isample] == "ttZ" || v_samples[isample] == "PrivMC_ttZ") {v_colors[isample] = v_custom_colors[1]->GetNumber();}

            //ttX
            else if(v_samples[isample] == "ttH") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttW") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "tttt") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttZZ") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttWW") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttWZ") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttZH") {v_colors[isample] = v_custom_colors[2]->GetNumber();}
            else if(v_samples[isample] == "ttWH") {v_colors[isample] = v_custom_colors[2]->GetNumber();}

            //tX
            else if(v_samples[isample] == "tHq") {v_colors[isample] = v_custom_colors[3]->GetNumber();}
            else if(v_samples[isample] == "tHW") {v_colors[isample] = v_custom_colors[3]->GetNumber();}
            else if(v_samples[isample] == "tWZ") {v_colors[isample] = v_custom_colors[3]->GetNumber();}
            else if(v_samples[isample] == "tGJets") {v_colors[isample] = v_custom_colors[3]->GetNumber();}
            else if(v_samples[isample] == "ST") {v_colors[isample] = v_custom_colors[3]->GetNumber();}

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
		}
    }

	else if(color_scheme == 0) //'Default' color scheme
	{
		for(int isample=0; isample<v_samples.size(); isample++)
		{
            //Signals
            if(v_samples[isample] == "tZq" || v_samples[isample] == "PrivMC_tZq") {v_colors[isample] = kOrange+10;}
            // if(v_samples[isample] == "tZq") {v_colors[isample] = test.GetNumber();}
            else if(v_samples[isample] == "ttZ" || v_samples[isample] == "PrivMC_ttZ") {v_colors[isample] = kOrange+6;}

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
            else if(v_samples[isample] == "tWZ") {v_colors[isample] = kSpring+2;}
            else if(v_samples[isample] == "tGJets") {v_colors[isample] = kSpring+2;}
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
    else if(var == "maxDijetMass") {nbins = 10; xmin = 0; xmax = 1000;}
    else if(var == "maxDijetPt") {nbins = 20; xmin = 0; xmax = 400;}
    else if(var == "maxDelPhiLL") {nbins = 10; xmin = 0; xmax = 3.5;}
    else if(var == "m3l" || var == "Mass_3l") {nbins = 20; xmin = 100; xmax = 500;}
    else if(var == "leptonCharge") {nbins = 3; xmin = -1.5; xmax = 1.5;}
    else if(var == "mTW") {nbins = 20; xmin = 0.; xmax = 200;}
    else if(var == "recoZ_Mass") {nbins = 20; xmin = 75.; xmax = 110;}
    else if(var == "lAsymmetry") {nbins = 20; xmin = -3.; xmax = 3.;}
    else if(var == "Mass_tZ") {nbins = 10; xmin = 200; xmax = 1000;}
    else if(var == "maxDeepCSV") {nbins = 20; xmin = 0.3; xmax = 1.1;}
    else if(var == "maxDijetDelPhi") {nbins = 20; xmin = 0.; xmax = 3.5;}
    else if(var == "maxDelRbL") {nbins = 20; xmin = 1.; xmax = 5;}
    else if(var == "Top_delRbl" || var == "Top_delRbW") {nbins = 20; xmin = 0.; xmax = 4.5;}

    else if(var.Contains("CSV", TString::kIgnoreCase)) {nbins = 20; xmin = 0.; xmax = 1.1;}
    else if(var.Contains("dR") || var.Contains("DelR", TString::kIgnoreCase) ) {nbins = 20; xmin = 0; xmax = 7.;}
    else if(var.BeginsWith("cos", TString::kIgnoreCase)) {nbins = 20; xmin = -1.; xmax = 1.;}

    else if(var == "njets" || var == "nbjets") {nbins = 6; xmin = 0.; xmax = 6;}

    // else if(var == "Lep3Pt") {nbins = 20; xmin = 0; xmax = 80;}

    else {nbins = 20; xmin = -5.; xmax = 150.;}
    // else {return false;}

	return true;
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
    if(var == "cosThetaStarPol") {return "cos#left(#theta^{*}_{pol}#right)";}
    if(var == "cosThetaStar") {return "cos#left(#theta^{*}#right)";}

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
TString Get_Category_Boolean_Name(TString region)
{
    if(region=="tZq") {return "is_tZq_3l_SR";}
    else if(region=="ttZ") {return "is_ttZ_3l_SR";}
    else if(region=="tWZ") {return "is_tWZ_3l_SR";}

    return "";
}

//Computes total nof entries which will be processed by the Produce_Templates() function, so that the Timebar is correct
//NB1 : don't account for nof syst weights, since they are computed all at once when reading an event
//NB2 : this returns the nof entries which will get read, not processed (else, should take cuts into account, etc.)
float Count_Total_Nof_Entries(TString dir_ntuples, TString t_name, vector<TString> v_samples, vector<TString> v_systTrees, vector<TString> v_cut_name, vector<TString> v_cut_def, vector<TString> v_lumiYears, bool makeHisto_inputVars, bool noSysts_inputVars)
{
    float total_nentries = 0;

    if(makeHisto_inputVars && noSysts_inputVars) //syst weights not computed
    {
        v_systTrees.resize(1);
    }

    //--- Define the cut chain to apply to events
    TString cuts_chain = "";
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
    }

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
    			if(tmp == "") {tmp = t_name;}
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

//For systs uncorrelated between years, need to add the year in the systName (to make it unique)
TString Get_Modified_SystName(TString systname, TString lumiYear)
{
    //Only list the systematics which are *not* to be correlated in-between years (i.e., need to give them a unique name per year)
    if(systname.BeginsWith("prefiring")
    || systname.BeginsWith("Btag")
    || systname.BeginsWith("LepEff")
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
    if(words.at(0).compare("rwgt") != 0) {cout<<"Error: EFT ID should start with rwgt_ prefix ! "<<endl; return v;}

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
