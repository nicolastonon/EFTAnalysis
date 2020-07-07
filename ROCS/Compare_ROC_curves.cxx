//by Nicolas Tonon (DESY)

//This code produces ROC curves and superimpose them on a plot
//Mostly automated, use the main() to define all options and define the list of input files
//The histograms are either produced :
// 1) from a TMVA output TTree [BDT]
// 2) from a root file created when training the neural network, using the train/test datasets [DNN]
// 3) from my own root file containing my BDT or DNN templates, for validation [BDT/DNN] //Warning : this uses the full sample and not only the test set !

//--------------------------------------------
// == REMARKS ==
// # To get ROC from user-file containing histograms :
// -- File should be produced with sufficient amount of bins for precision (make it consistent in main)
// -- Set correctly the name of object=histogram to retrieve

// -- 'Custom' histogram mode does not work yet, because should list all bkg processes in order to get consistent bkg definition (or store a merged 'bkg' histo when creating templates)
//--------------------------------------------

#include "../Utils/Helper.h"

#define DEBUG false

using namespace std;

//--------------------------------------------
// ##     ## ######## ##       ########  ######## ########     ######## ##     ## ##    ##  ######   ######
// ##     ## ##       ##       ##     ## ##       ##     ##    ##       ##     ## ###   ## ##    ## ##    ##
// ##     ## ##       ##       ##     ## ##       ##     ##    ##       ##     ## ####  ## ##       ##
// ######### ######   ##       ########  ######   ########     ######   ##     ## ## ## ## ##        ######
// ##     ## ##       ##       ##        ##       ##   ##      ##       ##     ## ##  #### ##             ##
// ##     ## ##       ##       ##        ##       ##    ##     ##       ##     ## ##   ### ##    ## ##    ##
// ##     ## ######## ######## ##        ######## ##     ##    ##        #######  ##    ##  ######   ######
//--------------------------------------------

/**
 * Sum all histograms from a vector into a single histogram, passed via pointer reference in argument
 * Allows to sum several histos into a single 'signal' or 'bkg' histo
 */
void Sum_Histograms(TH1F* &h, vector<TH1F*> v_histos)
{
	if(DEBUG) {cout<<endl<<"--- Sum_Histograms()"<<endl<<endl;}

	if(v_histos.size() < 2) {cout<<FRED("No histogram to sum !")<<endl; return;}

	h = 0; //histogram to be filled

	for(int ihisto=0; ihisto<v_histos.size(); ihisto++)
	{
        if(!v_histos[ihisto]) {continue;}

		if(!h) {h = (TH1F*) v_histos[ihisto]->Clone();}
		else {h->Add(v_histos[ihisto]);}
	}

    if(DEBUG) {cout<<endl<<"--- Exit Sum_Histograms()"<<endl<<endl;}

	return;
}

/**
 * Returns a color depending on index given in arg
 */
int Get_Color(int index)
{
    int palette=1;

	vector<int> v_colors;

	//-- Can uncomment these 4 lines to get 2*2 same colors (4 ROCs total)
	// v_colors.push_back(kRed);
	// v_colors.push_back(kBlue);
	// v_colors.push_back(kRed);
	// v_colors.push_back(kBlue);

	//-- Can use this palette to get 4*2 same colors (8+ ROCs)
    if(palette==0)
    {
        v_colors.push_back(kRed);
        v_colors.push_back(kBlue);
        v_colors.push_back(kGreen+1);
        v_colors.push_back(kBlack);
        v_colors.push_back(kRed);
        v_colors.push_back(kBlue);
        v_colors.push_back(kGreen+1);
        v_colors.push_back(kBlack);
        v_colors.push_back(kViolet-1);
        v_colors.push_back(kOrange-1);
    }
    else
    {
        v_colors.push_back(kBlack);
        v_colors.push_back(kRed);
        v_colors.push_back(kRed-2);
        v_colors.push_back(kPink);
        v_colors.push_back(kBlue+4);
        v_colors.push_back(kBlue-4);
        v_colors.push_back(kBlue-8);
        v_colors.push_back(kViolet);
        v_colors.push_back(kYellow+2);
        v_colors.push_back(kGreen+2);
        v_colors.push_back(kGreen-3);
        v_colors.push_back(kCyan);
    }

	if(index < v_colors.size() ) {return v_colors[index];}
	else {return index;}
}

// -- using https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines
void Apply_Cosmetics(TCanvas* &c1)
{
	float l = c1->GetLeftMargin();
	float t = c1->GetTopMargin();

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

	c1->SetFrameLineWidth(2.);

    return;
}

/**
 * Change the numerical precision displayed
 */
double GetFloatPrecision(double value, double precision)
{
    return (floor((value * pow(10, precision) + 0.5)) / pow(10, precision));
}



//--------------------------------------------
//  ######   ######## ########       ###    ##     ##  ######
// ##    ##  ##          ##         ## ##   ##     ## ##    ##
// ##        ##          ##        ##   ##  ##     ## ##
// ##   #### ######      ##       ##     ## ##     ## ##
// ##    ##  ##          ##       ######### ##     ## ##
// ##    ##  ##          ##       ##     ## ##     ## ##    ##
//  ######   ########    ##       ##     ##  #######   ######
//--------------------------------------------



/**
 * Compute AUC (ROC integral) from histogram produced by TMVA
 */
double Get_AUC_From_TMVAfile(TString filename, TString variable, bool use_TrainSample)
{
	if(DEBUG) {cout<<FYEL("-- Get_AUC_From_TMVAfile() --")<<endl;}

	TH1F* h = 0; //Pointer to histogram passed as argument, filled inside function

	if(!Check_File_Existence(filename) ) {cout<<FRED("File "<<filename<<" not found ! Abort")<<endl; return 0;}
	TFile* f = TFile::Open(filename);

	TString dir_name = "weights/Method_BDT/" + variable; //Path of subdir, hard-coded
	TString h_name = "MVA_" + variable + "_"; //Name of ROC histo, hard-coded

	if(use_TrainSample) {h_name+= "trainingRejBvsS";} //Training ROC
	else {h_name+= "rejBvsS";} //Testing ROC

	h = (TH1F*) f->Get(dir_name + "/" + h_name);
	if(!h) {cout<<FRED("Histo "<<dir_name<<"/"<<h_name<<" not found ! Abort")<<endl; return 0;}

	// cout<<"h->Integral() = "<<h->Integral()<<endl;
	double AUC = roundf(h->Integral() ); //Change precision
	AUC/= 100.; //Normalize

	if(DEBUG) {cout<<"AUC = "<<AUC<<endl;}

	return AUC;
}



//--------------------------------------------
//  ######   ######## ########    ##     ## ####  ######  ########  #######
// ##    ##  ##          ##       ##     ##  ##  ##    ##    ##    ##     ##
// ##        ##          ##       ##     ##  ##  ##          ##    ##     ##
// ##   #### ######      ##       #########  ##   ######     ##    ##     ##
// ##    ##  ##          ##       ##     ##  ##        ##    ##    ##     ##
// ##    ##  ##          ##       ##     ##  ##  ##    ##    ##    ##     ##
//  ######   ########    ##       ##     ## ####  ######     ##     #######
//--------------------------------------------


// ##### #    # #    #   ##
//   #   ##  ## #    #  #  #
//   #   # ## # #    # #    #
//   #   #    # #    # ######
//   #   #    #  #  #  #    #
//   #   #    #   ##   #    #


/**
 * Get histogram from a TMVA output TTree, for given cuts, etc.
 * @param h        histogram to be filled
 * @param t        TTree to read
 * @param variable variable to plot
 * @param cuts     cuts to apply -- "1" <-> no cut
 */
bool Create_Histogram_From_TMVA_Tree(bool is_bkg, TH1F*& h, TString filename, TString variable, bool use_TrainSample, TString cuts="1")
{
	if(DEBUG) {cout<<FYEL("-- Create_Histogram_From_TMVA_Tree() --")<<endl;}

    double xmin = -1; double xmax = +1; int nbins = 100; //Parameters of histogram ; BDT output usually between -1;+1

	h = 0; //Pointer to histogram passed as argument, filled inside function

	if(!Check_File_Existence(filename) ) {cout<<FRED("File "<<filename<<" not found ! Abort")<<endl; return 0;}
	TFile* f = TFile::Open(filename);

	//Hard-coded tree name
	TString treename = "weights/TestTree";
	if(use_TrainSample) {treename = "weights/TrainTree";}

	TTree* t = (TTree*) f->Get(treename);
	if(!t) {cout<<FRED("Null TTree '"<<treename<<"' ! Are you sure you are indeed reading a TMVA-produced file (else change option) ? Abort")<<endl; return 0;}

    //Define root interactive command to draw desired histogram
	variable = variable + " >> h(" + Convert_Number_To_TString(nbins) + "," + Convert_Number_To_TString(xmin) + "," + Convert_Number_To_TString(xmax) + ")";
    if(DEBUG) {cout<<"variable = "<<variable<<endl;}

	//Option to interactively apply weights to histogram, and apply set of cuts
	//NB : in TMVA, sig <-> "classID==0" & bkg <-> "classID==1"
	// NB : in TrainTree, weight variable is wrong (don't' know what it correspond to...!) //But weight ok for testTree
	cuts = "weight * (" + cuts + " && classID==" + Convert_Number_To_TString(is_bkg) + ")";
    if(DEBUG) {cout<<"cuts = "<<cuts<<endl;}

	//Produce histogram (interactively)
	t->Draw(variable, cuts);

	//Retrieve generated histogram
	h = (TH1F*) gDirectory->Get("h")->Clone();
	h->SetDirectory(0); //NECESSARY so that histo is not associated with TFile, and doesn't get deleted when file closed !
	if(!h || h->GetEntries() == 0) {cout<<BOLD(FRED("Null or void histogram (made from TMVA TTree) ! Abort !"))<<endl; return 0;}

	//Retrieve underflow content
	double underflow = h->GetBinContent(0);
	if(underflow>0)
	{
		h->AddBinContent(1, underflow); //Add underflow
		h->SetBinContent(0, 0); //Remove underflow
	}

	//Retrieve overflow content
	double overflow = h->GetBinContent(nbins+1);
	if(overflow>0)
	{
		h->AddBinContent(nbins, overflow); //Add overflow
		h->SetBinContent(nbins+1, 0); //Remove overflow
	}

	f->Close();

	return true;
}

/**
 * Can directly retrieve ROC histogram from TMVA files
 * NB : not used for now, since implemented use of TGraph... keep it as before (Else should slightly modify main func, to plot 1 single hist directly instead of computing eff. from TGraphs...)
 */
bool Retrieve_Histogram_From_TMVA_File(TH1F*& h, TString filename, TString variable, bool is_bkg, bool use_TrainSample)
{
	if(DEBUG) {cout<<FYEL("-- Retrieve_Histogram_From_TMVA_File() --")<<endl;}

	h = 0; //Pointer to histogram passed as argument, filled inside function

	if(!Check_File_Existence(filename) ) {cout<<FRED("File "<<filename<<" not found ! Abort")<<endl; return 0;}
	TFile* f = TFile::Open(filename);

    //Hard-coded histogram name (naming conventions following those chosen in my TMVA training)
    TString dir_name = "weights/Method_BDT/" + variable; //Path of subdir, hard-coded
	TString h_name = "MVA_" + variable + "_"; //Name of ROC histo, hard-coded

	if(use_TrainSample) {h_name+= "trainingRejBvsS";} //Training dataset ROC curve
	else {h_name+= "rejBvsS";} //Testing dataset ROC curve

	h = (TH1F*) f->Get(dir_name + "/" + h_name);
	if(!h) {cout<<FRED("Histo "<<dir_name<<"/"<<h_name<<" not found ! Abort")<<endl; return 0;}

	if(DEBUG) {cout<<"Integral = "<<h->Integral()<<endl;}

	return true;
}


// #    # ###### #####    ##    ####
// #   #  #      #    #  #  #  #
// ####   #####  #    # #    #  ####
// #  #   #      #####  ######      #
// #   #  #      #   #  #    # #    #
// #    # ###### #    # #    #  ####

/**
 * Get histograms produced during DNN training, containing ROC curves for train/test datasets
 */
bool Get_Histogram_From_KerasFile(TH1F*& h, TString filename, TString histo_name)
{
	if(DEBUG) {cout<<FYEL("-- Get_Histogram_From_KerasFile() --")<<endl;}

	h = 0; //Pointer to histogram passed as argument, filled inside function

	if(!Check_File_Existence(filename) ) {cout<<FRED("File "<<filename<<" not found ! Abort")<<endl; return 0;}
	TFile* f = TFile::Open(filename);

    if(!f->GetListOfKeys()->Contains(histo_name) ) {cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found !"))<<endl; return false;}

	if(DEBUG) {cout<<"- Opening histo : "<<histo_name<<endl;}

	h = (TH1F*) f->Get(histo_name);
	h->SetDirectory(0); //NECESSARY so that histo is not associated with TFile, and doesn't get deleted when file closed !
	if(!h || h->GetEntries() == 0) {cout<<BOLD(FRED("Null or void histogram '"<<histo_name<<"' ! Abort !"))<<endl; return 0;}

	//Retrieve underflow content
	double underflow = h->GetBinContent(0);
	if(underflow>0)
	{
		h->AddBinContent(1, underflow); //Add underflow
		h->SetBinContent(0, 0); //Remove underflow
	}

	//Retrieve overflow content
	int nbins = h->GetNbinsX();
	double overflow = h->GetBinContent(nbins+1);
	if(overflow>0)
	{
		h->AddBinContent(nbins, overflow); //Add overflow
		h->SetBinContent(nbins+1, 0); //Remove overflow
	}

	f->Close();

    if(!h || h->GetEntries() == 0) {cout<<BOLD(FRED("Null or void histogram '"<<histo_name<<"' ! Abort !"))<<endl; return 0;}

    if(DEBUG) {cout<<FYEL("-- Exit Get_Histogram_From_KerasFile() --")<<endl;}

	return 1;
}


//  ####  #    #  ####  #####  ####  #    #
// #    # #    # #        #   #    # ##  ##
// #      #    #  ####    #   #    # # ## #
// #      #    #      #   #   #    # #    #
// #    # #    # #    #   #   #    # #    #
//  ####   ####   ####    #    ####  #    #


/**
 * Get histogram directly from any root file -- if not following my naming conventions, must hard-code histo names !
 * @param h          histogram to be filled
 * @param filename   path of TFile containing histogram
 * @param histo_name name of histogram to retrieve in file
 * NB : problem : THIS USES THE FULL SAMPLE, NOT JUST TEST SET ! Use other function to get test sample from Keras
 */
bool Get_Histogram_From_CustomFile(TH1F*& h, TString filename, TString histo_name)
{
	if(DEBUG) {cout<<FYEL("-- Get_Histogram_From_CustomFile() --")<<endl;}
    cout<<ITAL("Warning : function 'Get_Histogram_From_CustomFile()' used full sample, not only train/test dataset... !")<<endl;

	h = 0; //Init pointer to histogram passed as argument
	if(!Check_File_Existence(filename) ) {cout<<FRED("File "<<filename<<" not found ! Abort")<<endl; return 0;}
	TFile* f = TFile::Open(filename);

	if(DEBUG) {cout<<"-- Opening histo : "<<histo_name<<endl;}

    if(!f->GetListOfKeys()->Contains(histo_name) ) {cout<<ITAL(DIM("Histogram '"<<histo_name<<"' : not found !"))<<endl; return false;}

	h = (TH1F*) f->Get(histo_name);
	h->SetDirectory(0); //NECESSARY so that histo is not associated with TFile, and doesn't get deleted when file closed !
	if(!h || h->GetEntries() == 0) {cout<<BOLD(FRED("Null or void histogram '"<<histo_name<<"' (from user-file) ! Abort !"))<<endl; return 0;}

	//Retrieve underflow content
	double underflow = h->GetBinContent(0);
	if(underflow>0)
	{
		h->AddBinContent(1, underflow); //Add underflow
		h->SetBinContent(0, 0); //Remove underflow
	}

	//Retrieve overflow content
	int nbins = h->GetNbinsX();
	double overflow = h->GetBinContent(nbins+1);
	if(overflow>0)
	{
		h->AddBinContent(nbins, overflow); //Add overflow
		h->SetBinContent(nbins+1, 0); //Remove overflow
	}

	f->Close();

	return true;
}


//--------------------------------------------
//  ######   ######## ########    ########  ######   ########     ###    ########  ##     ##
// ##    ##  ##          ##          ##    ##    ##  ##     ##   ## ##   ##     ## ##     ##
// ##        ##          ##          ##    ##        ##     ##  ##   ##  ##     ## ##     ##
// ##   #### ######      ##          ##    ##   #### ########  ##     ## ########  #########
// ##    ##  ##          ##          ##    ##    ##  ##   ##   ######### ##        ##     ##
// ##    ##  ##          ##          ##    ##    ##  ##    ##  ##     ## ##        ##     ##
//  ######   ########    ##          ##     ######   ##     ## ##     ## ##        ##     ##


// ######## ######## ######## ####  ######  #### ######## ##    ##  ######  ##    ##
// ##       ##       ##        ##  ##    ##  ##  ##       ###   ## ##    ##  ##  ##
// ##       ##       ##        ##  ##        ##  ##       ####  ## ##         ####
// ######   ######   ######    ##  ##        ##  ######   ## ## ## ##          ##
// ##       ##       ##        ##  ##        ##  ##       ##  #### ##          ##
// ##       ##       ##        ##  ##    ##  ##  ##       ##   ### ##    ##    ##
// ######## ##       ##       ####  ######  #### ######## ##    ##  ######     ##
//

//  ######  ####  ######            ##    ########  ##    ##  ######
// ##    ##  ##  ##    ##          ##     ##     ## ##   ##  ##    ##
// ##        ##  ##               ##      ##     ## ##  ##   ##
//  ######   ##  ##   ####       ##       ########  #####    ##   ####
//       ##  ##  ##    ##       ##        ##     ## ##  ##   ##    ##
// ##    ##  ##  ##    ##      ##         ##     ## ##   ##  ##    ##
//  ######  ####  ######      ##          ########  ##    ##  ######
//--------------------------------------------



/**
 * For each bin of input signal/bkg histograms, will compute corresponding efficiency and store it in a TGraph (passed as arg)
 */
bool Produce_Efficiency_TGraph(TGraph* &g, double& AUC, TH1F* h_sig, TH1F* h_bkg)
{
	if(DEBUG) {cout<<FYEL("-- Produce_Efficiency_TGraph()")<<endl;}

	if(!h_sig || h_sig->GetEntries() == 0) {cout<<BOLD(FRED("Null or void signal TGraph ! Abort !"))<<endl; return 0;}
	if(!h_bkg || h_bkg->GetEntries() == 0) {cout<<BOLD(FRED("Null or void bkg TGraph ! Abort !"))<<endl; return 0;}

	int nbins = h_sig->GetNbinsX();

	g = new TGraph(nbins); //TGraph to be filled

	if(h_sig->GetNbinsX() != h_bkg->GetNbinsX()) {cout<<BOLD(FRED("Different nbins for sig & bkg histograms ! Abort"))<<endl;}
	if(nbins < 50) {cout<<"Warning : only "<<nbins<<" bins in histograms (Low precision) !"<<endl;}

    bool invert_ROC = true; //The convention used on the plot is such that we need the signal distributions to be 'more to the right' w.r.t. the bkg distribution. Otherwise, we need to take the symmetric curve w.r.t. the diagonal
    if(h_sig->GetMean() < h_bkg->GetMean()) {invert_ROC = false;}

    AUC = 0; //Also compute the area under curve
	double integral_sig = h_sig->Integral();
	double integral_bkg = h_bkg->Integral();
	// cout<<"Integral sig = "<<integral_sig<<endl;
	// cout<<"Integral bkg = "<<integral_bkg<<endl;

    //For each bin, compute efficiency and store in TGraph
	for (int ibin=1; ibin<nbins+1; ibin++)
	{
		//Compute efficiency for sig and bkg, from current bin to xmax
		double integral_sig_tmp = h_sig->Integral(ibin, nbins);
		double integral_bkg_tmp = h_bkg->Integral(ibin, nbins);

		double eff_sig_tmp = integral_sig_tmp / integral_sig;
		double eff_bkg_tmp = integral_bkg_tmp / integral_bkg;

		// eff_sig_tmp = eff_bkg_tmp; //set for xcheck, so that auc = 0.5
        // cout<<"ibin "<<ibin<<endl; cout<<"eff_sig_tmp "<<eff_sig_tmp<<endl; cout<<"eff_bkg_tmp "<<eff_bkg_tmp<<endl;

		//Fill efficiency graph
        if(invert_ROC) {g->SetPoint(ibin-1, eff_sig_tmp, 1-eff_bkg_tmp);}
        else {g->SetPoint(ibin-1, 1-eff_sig_tmp, eff_bkg_tmp);}
	}

    //ROOT's TGraph->Integral() computes the integral of the closed polygon defined by first/lasrt points
    //Since for ROCs, the coordinates of these points are always (0,1) and (1,0) respectively, it means it omits all the area under the "random guess" threshold => Include this area manually
    AUC = g->Integral() + 0.5;
    // cout<<"AUC "<<AUC<<endl;

	return true;
}



//--------------------------------------------
//  ######   ######## ########    ########   #######   ######   ######
// ##    ##  ##          ##       ##     ## ##     ## ##    ## ##    ##
// ##        ##          ##       ##     ## ##     ## ##       ##
// ##   #### ######      ##       ########  ##     ## ##        ######
// ##    ##  ##          ##       ##   ##   ##     ## ##             ##
// ##    ##  ##          ##       ##    ##  ##     ## ##    ## ##    ##
//  ######   ########    ##       ##     ##  #######   ######   ######
//--------------------------------------------



/**
 * Get all ROC curves histograms
 */
void Get_ROC_Curves(vector<TGraph*>& v_graph, vector<double>& v_AUC, vector<TString>& v_label, vector<TString> v_filepath, vector<TString> v_isTMVA_file, vector<TString> v_Filelabel, vector<bool> v_isTrainSample, TString region, vector<TString> v_processes, bool superimpose_allNodes_DNN, TString lumiYear, TString cuts="1")
{
	if(!v_filepath.size() || !v_isTMVA_file.size() ) {cout<<FRED("Passed void vector as argument ! Abort !")<<endl; return;}

    int index_tgraph = 0; //Keep track of TGraph indices

	//For each file, will get/produce sig & bkg histograms, compute efficiency, and plot it
	for(int ifile=0; ifile<v_filepath.size(); ifile++)
	{
		cout<<endl<<endl<<endl<<UNDL("* Open file "<<v_filepath[ifile]<<" :")<<endl;
        if(v_isTMVA_file[ifile] != "TMVA" && v_isTMVA_file[ifile] != "Keras" && v_isTMVA_file[ifile] != "Custom") {cout<<FRED("ERROR ! Wrong file type...")<<endl; return;}

        //-- For each signal (DNN : 1 different signal for each node), get signal histo and bkg histos (i.e. all processes which are not signal for this node)
        for(int isig=0; isig<v_processes.size(); isig++)
		{
            if(v_isTMVA_file[ifile] == "TMVA" && isig > 0) {break;} //TMVA file : don't support multiclass yet ; only consider first process (can define signal/bkg from there)
            if(!superimpose_allNodes_DNN && isig > 0) {break;} //Plot only first node

            cout<<"- Signal n°"<<isig+1<<" : "<<v_processes[isig]<<endl;

            //For each ROC curve, get 1 'signal' histo ; keep vector, only to allow easier automation in the future
            //For TMVA, consider single background source (i.e. support binary class only) ; for Keras DNN, consider all processes which are not signal as backgrounds
            vector<TH1F*> v_h_sig, v_h_bkg;

            v_h_sig.push_back(0); v_AUC.push_back(0);  //Consider 1 'signal' proc
            v_label.push_back(v_Filelabel[ifile] + " " + v_processes[isig]);
            if(v_isTMVA_file[ifile] == "TMVA") //Produce histo from TMVA TTree
            {
                v_h_bkg.resize(1); //Will consider 1 'background' proc
                // v_AUC[index_tgraph] = Get_AUC_From_TMVAfile(v_filepath[ifile], v_processes[isig], v_isTrainSample[ifile]);
                if(!Create_Histogram_From_TMVA_Tree(false, v_h_sig[0], v_filepath[ifile], v_processes[isig], v_isTrainSample[ifile], cuts) ) {return;}
            }
            else if(v_isTMVA_file[ifile] == "Keras") //Retrieve directly histogram from file produced by DNN training python script
            {
                v_h_bkg.resize(v_processes.size()-1); //Will consider as 'backgrounds' all the processes not considered as 'signal'
                TString hname_tmp = "hist_test_NODE_" + v_processes[isig] + "_CLASS_" + v_processes[isig]; //Performance for class A in node A
                if(!Get_Histogram_From_KerasFile(v_h_sig[0], v_filepath[ifile], hname_tmp)) {return;} //Check naming convention used in pure-Keras NN
            }
            else if(v_isTMVA_file[ifile] == "Custom") //Retrieve directly histogram from custom user file
            {
                v_h_bkg.resize(v_processes.size()-1); //Will consider as 'backgrounds' all the processes not considered as 'signal'
                TString hname_tmp = "DNN0_" + region + "_" + lumiYear + "__" + v_processes[isig]; //hardcode as needed
                if(!Get_Histogram_From_CustomFile(v_h_sig[0], v_filepath[ifile], hname_tmp)) {return;}
            }

            int idx_bkg = 0; //Even if process is second in process list, if it is the first bkg we want its index in vector to be 0, etc.

            //Get bkg histos (all processes which are not signal)
            for(int ibkg=0; ibkg<v_processes.size(); ibkg++)
    		{
                if(v_isTMVA_file[ifile] == "TMVA" && ibkg > 0) {continue;} //TMVA file : don't support multiclass yet ; only consider first process (can define signal/bkg from there)
                else if(v_isTMVA_file[ifile] == "Keras" && isig == ibkg) {continue;} //Keras file : 'signal' process can not be also considered as 'background'

                if(v_isTMVA_file[ifile] == "TMVA") //Produce histo from TMVA TTree
    			{
                    if(!Create_Histogram_From_TMVA_Tree(true, v_h_bkg[idx_bkg], v_filepath[ifile], v_processes[isig], v_isTrainSample[ifile], cuts) ) {return;}
    			}
    			else if(v_isTMVA_file[ifile] == "Keras") //Retrieve directly histogram from file produced by DNN training python script
    			{
                    TString hname_tmp = "hist_test_NODE_" + v_processes[isig] + "_CLASS_" + v_processes[ibkg]; //Performance for class A in node B
    				if(!Get_Histogram_From_KerasFile(v_h_bkg[idx_bkg], v_filepath[ifile], hname_tmp)) {return;} //Check naming convention used in pure-Keras NN
                }
                else if(v_isTMVA_file[ifile] == "Custom") //Retrieve directly histogram from custom user file
                {
                	TString hname_tmp = "DNN" + Convert_Number_To_TString(isig) + "_" + region + "_" + lumiYear + "__" + v_processes[ibkg]; //hardcode as needed
                	if(!Get_Histogram_From_CustomFile(v_h_bkg[idx_bkg], v_filepath[ifile], hname_tmp)) {return;}
                }

                idx_bkg++; //Increment bkg vector index
            }

            //Sum histograms (many processes --> signal/bkg)
            TH1F* h_sig = 0;
            if(v_h_sig.size() > 1) {Sum_Histograms(h_sig, v_h_sig);}
            else if(v_h_sig.size() == 1) {h_sig = (TH1F*) v_h_sig[0]->Clone();}
            else {cout<<BOLD(FRED("Wrong vector size !"))<<endl; return;}
            if(!h_sig) {cout<<"h_sig is null !"<<endl; return;}

            TH1F* h_bkg = 0;
            if(v_h_bkg.size() > 1) {Sum_Histograms(h_bkg, v_h_bkg);}
            else if(v_h_bkg.size() == 1) {h_bkg = (TH1F*) v_h_bkg[0]->Clone();}
            else {cout<<BOLD(FRED("Wrong vector size !"))<<endl; return;}

            //Produce efficiency TGraph from sig & bkg histograms
            v_graph.push_back(0);
            if(!Produce_Efficiency_TGraph(v_graph[index_tgraph], v_AUC[index_tgraph], h_sig, h_bkg) ) {return;}
            index_tgraph++;

            for(int isig=0; isig<v_h_sig.size(); isig++) {delete v_h_sig[isig]; v_h_sig[isig] = 0;}
            for(int ibkg=0; ibkg<v_h_bkg.size(); ibkg++) {delete v_h_bkg[ibkg]; v_h_bkg[ibkg] = 0;}
            delete h_sig; delete h_bkg;

        } //end signal loop

	} //end file loop

    return;
}



//--------------------------------------------
//  ######  ##     ## ########  ######## ########  #### ##     ## ########   #######   ######  ########    ########   #######   ######
// ##    ## ##     ## ##     ## ##       ##     ##  ##  ###   ### ##     ## ##     ## ##    ## ##          ##     ## ##     ## ##    ##
// ##       ##     ## ##     ## ##       ##     ##  ##  #### #### ##     ## ##     ## ##       ##          ##     ## ##     ## ##
//  ######  ##     ## ########  ######   ########   ##  ## ### ## ########  ##     ##  ######  ######      ########  ##     ## ##
//       ## ##     ## ##        ##       ##   ##    ##  ##     ## ##        ##     ##       ## ##          ##   ##   ##     ## ##
// ##    ## ##     ## ##        ##       ##    ##   ##  ##     ## ##        ##     ## ##    ## ##          ##    ##  ##     ## ##    ##
//  ######   #######  ##        ######## ##     ## #### ##     ## ##         #######   ######  ########    ##     ##  #######   ######
//--------------------------------------------



/**
 * Superimpose all ROC curves passed in arg on same plot
 */
void Superimpose_ROC_Curves(vector<TGraph*> v_graph, vector<TString> v_label, vector<double> v_AUC)
{
    TCanvas* c = new TCanvas("", "", 1000, 800);
    c->SetTopMargin(0.1);
    // c->SetGrid();
	c->cd();

	//Draw custom background (coordinates hard-coded...)
	TPad *p = new TPad("p","p",0.16,0.13,0.97,0.9);
	p->SetFillColorAlpha(kGray, 0.15);

	TLine* l_randGuess = new TLine(0, 1, 1, 0); //Display "random guess" ROC
	l_randGuess->SetLineStyle(2);

	vector<TLine*> v_gridlines_Y(9);
	vector<TLine*> v_gridlines_X(9);
	double ticklength = 0.03;

	for(int i=0; i<9; i++)
	{
		v_gridlines_Y[i] = new TLine((i+1)*0.1, 0+ticklength, (i+1)*0.1, 1-ticklength);

		v_gridlines_X[i] = new TLine(0+ticklength, (i+1)*0.1, 1-ticklength, (i+1)*0.1);

		v_gridlines_Y[i]->SetLineColor(0);
		v_gridlines_Y[i]->SetLineWidth(4);
		// v_gridlines_Y[i]->SetLineStyle(2);

		v_gridlines_X[i]->SetLineColor(0);
		v_gridlines_X[i]->SetLineWidth(4);
		// v_gridlines_X[i]->SetLineStyle(2);
	}

	//Need this to set the axis on plot
	TH1F* h_axis = new TH1F("", "", 10, 0, 1);
	h_axis->SetMinimum(0.0001); //Remove 0 from axis
	h_axis->SetMaximum(1);

	h_axis->GetXaxis()->SetTitle("Signal efficiency");
	h_axis->GetXaxis()->SetTitleOffset(1);
	h_axis->GetXaxis()->SetLabelFont(42);
	h_axis->GetXaxis()->SetTitleFont(42);
	h_axis->GetXaxis()->SetTickLength(0.04);

	h_axis->GetYaxis()->SetTitle("Background rejection");
	h_axis->GetYaxis()->SetTitleOffset(1);
	h_axis->GetYaxis()->SetLabelFont(42);
	h_axis->GetYaxis()->SetTitleFont(42);
	h_axis->GetYaxis()->SetTickLength(0.04);

    float y_min_leg = 0.4 - v_graph.size() * 0.05;
    if(y_min_leg < 0.01) {y_min_leg=0.01;}
	TLegend* legend = new TLegend(0.19, y_min_leg, 0.60, 0.45);
    legend->SetTextSize(0.03);
    if(v_graph.size()>8) {legend->SetTextSize(0.02);}
    legend->SetLineColor(kGray);

	h_axis->Draw(); //Draw axis first

	//Draw custom grid
	p->Draw("same");
	for(int i=0; i<v_gridlines_Y.size(); i++)
	{
		v_gridlines_Y[i]->Draw("same");
		v_gridlines_X[i]->Draw("same");
	}

	l_randGuess->Draw("same"); //random guess

	for(int igraph=0; igraph<v_graph.size(); igraph++)
	{
        // v_graph[igraph]->SetLineColor(Get_Color(igraph));
        v_graph[igraph]->SetLineColorAlpha(Get_Color(igraph), 0.75);
		v_graph[igraph]->SetLineWidth(4);

		//-- dashed lines for some ROCs ?
        // if(igraph > 1)
		if(igraph >= v_graph.size()/2)
		{
			// v_graph[igraph]->SetLineStyle(8);
            v_graph[igraph]->SetLineWidth(4);
        }

		// cout<<"v_AUC[igraph] "<<v_AUC[igraph]<<endl;
        // cout<<"v_label[igraph] "<<v_label[igraph]<<endl;

        if(v_AUC[igraph] != 0) {legend->AddEntry(v_graph[igraph], v_label[igraph] + " (AUC="+Convert_Number_To_TString(v_AUC[igraph]) + ")", "L");}
		else {legend->AddEntry(v_graph[igraph], v_label[igraph], "L");}

		v_graph[igraph]->Draw("same C"); //'C' <-> smooth curve
	}

	legend->Draw("same");
	Apply_Cosmetics(c);

	c->SaveAs("ROC_comparison.png");

	delete h_axis; h_axis = 0;
	delete p; p = NULL;
	delete legend; legend = 0;
	delete c; c = 0;

	for(int i=0; i<v_gridlines_X.size(); i++)
	{
		delete v_gridlines_X[i]; v_gridlines_X[i] = NULL;
		delete v_gridlines_Y[i]; v_gridlines_Y[i] = NULL;
	}
	delete l_randGuess; l_randGuess = NULL;

	return;
}



//--------------------------------------------
// ##     ##    ###    ##    ## ########    ########  ##        #######  ########
// ###   ###   ## ##   ##   ##  ##          ##     ## ##       ##     ##    ##
// #### ####  ##   ##  ##  ##   ##          ##     ## ##       ##     ##    ##
// ## ### ## ##     ## #####    ######      ########  ##       ##     ##    ##
// ##     ## ######### ##  ##   ##          ##        ##       ##     ##    ##
// ##     ## ##     ## ##   ##  ##          ##        ##       ##     ##    ##
// ##     ## ##     ## ##    ## ########    ##        ########  #######     ##
//--------------------------------------------



/**
 * Main function producing all plots required in main()
 */
void Make_Plot(vector<TString> v_filepath, vector<TString> v_Filelabel, vector<TString> v_isTMVA_file, vector<bool> v_isTrainSample, TString region, vector<TString> v_processes, bool superimpose_allNodes_DNN, TString lumiYear, TString cuts="1")
{
    cout<<endl<<YELBKG("                          ")<<endl<<endl;
	cout<<FYEL("--- Will superimpose all ROC curves on plot ---")<<endl;
	cout<<endl<<YELBKG("                          ")<<endl<<endl;

    vector<TGraph*> v_graph;

    //1 value of AUC (area under curve) per ROC -- NB : not available yet for PureKeras files (only TMVA)
	vector<double> v_AUC;
    vector<TString> v_label;

    Get_ROC_Curves(v_graph, v_AUC, v_label, v_filepath, v_isTMVA_file, v_Filelabel, v_isTrainSample, region, v_processes, superimpose_allNodes_DNN, lumiYear, cuts);

    Superimpose_ROC_Curves(v_graph, v_label, v_AUC);

    for(int igraph=0; igraph<v_graph.size(); igraph++) {delete v_graph[igraph]; v_graph[igraph] = 0;}

    return;
}



//--------------------------------------------
//--------------------------------------------
//--------------------------------------------



//--------------------------------------------
// ##     ##    ###    #### ##    ##
// ###   ###   ## ##    ##  ###   ##
// #### ####  ##   ##   ##  ####  ##
// ## ### ## ##     ##  ##  ## ## ##
// ##     ## #########  ##  ##  ####
// ##     ## ##     ##  ##  ##   ###
// ##     ## ##     ## #### ##    ##
//--------------------------------------------

int main(int argc, char **argv)
{
	Load_Canvas_Style(); //Styling

//-- Select options
//--------------------------------------------
	vector<TString> v_processes; //For DNNs, must declare all process classes to define 'background' properly //Not used for BDT (signal vs bkg)
    // v_processes.push_back("tZq");
    // v_processes.push_back("ttZ");
    // v_processes.push_back("Backgrounds");
    v_processes.push_back("SM");
    v_processes.push_back("EFT");

    TString lumiYear = "2017"; //'2016,'2017','2018','Run2'

    bool superimpose_allNodes_DNN = false; //true <-> will plot 1 ROC per process class (for DNN only)

    // TString channel = "";
    TString region = "SR"; //not used for now
	TString cuts = "1"; //1 <-> No cut
//--------------------------------------------

//-- Fill vectors
//--------------------------------------------
//Push_back 1 element to each vector, for each ROC curve
//NB : within tHq framework, objNames are expected to be the following :
//-- for regular BDT : "weights/TestTree"

	vector<TString> v_filepath; //Path of TFile containing TMVA TTree or histograms
	vector<TString> v_Filelabel; //File label to be displayed on plot
    vector<TString> v_isTMVA_file; //'TMVA' <-> looking for TMVA TTree ; 'Keras' <-> looking for histograms in file producing during Keras training ; 'Custom' <-> looking for my template histograms (or any hardcoded histo name...)
    vector<bool> v_isTrainSample; //True <-> looking for ROC from train sample ; else test sample
//--------------------------------------------

//--------------------------------------------
    // v_filepath.push_back("../outputs/BDT_"+v_processes[0]+"_"+lumiYear+".root");
    // v_Filelabel.push_back("BDT "+lumiYear);
    // v_isTMVA_file.push_back("TMVA"); v_isTrainSample.push_back(false);

    // v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+".root");
    // v_Filelabel.push_back("NN "+lumiYear);
    // v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);

    // v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+"_feature.root");
    // v_Filelabel.push_back("Z_pt");
    // v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);

    v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+"_ctw1.root");
    v_Filelabel.push_back("NN 5");
    v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);
    // v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+"_ctw-5.root");
    // v_Filelabel.push_back("NN -5");
    // v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);
    v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+"_recoZ_Pt_ctw1.root");
    v_Filelabel.push_back("Zpt 5");
    v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);

    // for(int i=-5; i<=5; i++)
    // {
    //     v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+"_ctw"+std::to_string(i)+".root");
    //     v_Filelabel.push_back("NN ctz_"+std::to_string(i));
    //     v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);
    // }

//--------------------------------------------

//== Run main function ==
//--------------------------------------------

	Make_Plot(v_filepath, v_Filelabel, v_isTMVA_file, v_isTrainSample, region, v_processes, superimpose_allNodes_DNN, lumiYear, cuts);

//--------------------------------------------

	return 0;
}
