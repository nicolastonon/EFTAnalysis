#include "Helper.h"

using namespace std;

//Hardcode nicer latex-compatible category names //Can use '\\mathbf{}' for bold, but does not work with greek letters ?
//NB: makes use of local definitions from AN/paper
TString Get_Category_LatexName(TString cat)
{
    TString name = cat;

    if(cat.Contains("signal")) {name = "\\sr";}
    else if(cat.Contains("ttZ4l")) {name = "\\srttzfour";}
    if(cat.Contains("wz")) {name = "\\WZ CR";}
    if(cat.Contains("zz")) {name = "\\ZZ CR";}
    if(cat.Contains("xg")) {name = "$X\\gamma$ CR";}
    if(cat.Contains("dy")) {name = "$\\dy$ CR";}
    if(cat.Contains("ttbar")) {name = "$t\\bar{t}}$ CR";}

    return name;
}


//--------------------------------------------
// ##    ## #### ######## ##       ########     ########    ###    ########  ##       ########
//  ##  ##   ##  ##       ##       ##     ##       ##      ## ##   ##     ## ##       ##
//   ####    ##  ##       ##       ##     ##       ##     ##   ##  ##     ## ##       ##
//    ##     ##  ######   ##       ##     ##       ##    ##     ## ########  ##       ######
//    ##     ##  ##       ##       ##     ##       ##    ######### ##     ## ##       ##
//    ##     ##  ##       ##       ##     ##       ##    ##     ## ##     ## ##       ##
//    ##    #### ######## ######## ########        ##    ##     ## ########  ######## ########
//--------------------------------------------

void Compute_Write_Yields(vector<TString> v_samples, vector<TString> v_label, TString category, TString signal, TString lumi, bool group_samples_together, bool remove_totalSF, TString channel, bool use_privSamples_asSignal)
{
    bool create_latex_table = true; //true <-> also output latex-format yield tables
        bool blind = false; //true <-> don't include DATA in latex tables (but still include in printouts)
        int precision = 1; //Nof decimals displayed after floating point

//--------------------------------------------

    cout<<endl<<BYEL("                          ")<<endl<<endl;
	cout<<FYEL("--- Will count the yields for each sample ---")<<endl;
	cout<<"(category : "<<category<<" / lumi : "<<lumi<<" / channel : "<<channel<<")"<<endl;
    cout<<endl<<BYEL("                          ")<<endl<<endl;

    mkdir("./outputs/", 0777);
    mkdir("./outputs/yields", 0777);
    TString outname = "./outputs/yields/Yields_"+category+"_"+lumi;
    if(channel != "") {outname+= "_" + channel;}
    outname+= ".txt";

    TString outname_latex = "./outputs/yields/latex/Yields_"+category+"_"+lumi;
    if(create_latex_table)
    {
        mkdir("./outputs/yields/latex", 0777);
        if(channel != "") {outname_latex+= "_" + channel;}
        outname_latex+= ".txt";
    }
    ofstream file_latex; //Declaration

    vector<TString> v_years; //'Run2' -> Sum all 3 years
    if(lumi == "Run2") {v_years.push_back("2016"); v_years.push_back("2017"); v_years.push_back("2018");}
    else {v_years.push_back(lumi);}

    if(group_samples_together == false)
    {
        for(int isample=0; isample<v_label.size(); isample++)
        {
            // cout<<"v_label[isample] "<<v_label[isample];
            v_label[isample] = v_samples[isample];
            // cout<<" ==> "<<v_label[isample]<<endl;
        }
    }

//--------------------------------------------

    if(create_latex_table)
    {
        file_latex.open(outname_latex.Data());

    	//NB : '\' character must be escaped as '\\' in C++!
        // file_latex<<"\\begin{table}[]"<<endl; //Horizontal
        file_latex<<"\\begin{sidewaystable}[]"<<endl; //Vertical //Requires package rotating
        file_latex<<"\\centering"<<endl;
        file_latex<<"\\small"<<endl; //normal, small, tiny
        file_latex<<"\\begin{tabular}{|c|"; //Init with 1 column (leftmost column, e.g. listing different event categories)
    	for(int isample=0; isample<v_label.size(); isample++)
    	{
            // cout<<"v_label[isample] "<<v_label[isample]<<endl; //Debug

            if(v_label[isample].Contains("TTbar") || v_label[isample].Contains("DY")) {continue;} //Consider DD NPL, not MC
            else if(v_samples[isample] == "DATA") {continue;} //Data appended manually
            else if(!use_privSamples_asSignal && v_samples[isample].Contains("PrivMC")) {continue;} //Yield tables: consider central samples
            else if(use_privSamples_asSignal && (v_samples[isample]=="tZq" || v_samples[isample]=="ttZ" || v_samples[isample]=="tWZ")) {continue;} //Yield tables: consider private samples

            // cout<<"Pass "<<(isample == v_label.size()-1 || v_label[isample] != v_label[isample+1])<<endl; //Debug

    		if(isample == v_label.size()-1 || v_label[isample] != v_label[isample+1])
    		{
    			file_latex<<"c|"; //1 column per process (-1 already printed above)
    		}
    	}
        file_latex<<"|c|"; //Total SM
        if(!blind) {file_latex<<"|c|";} //Data
    	file_latex<<"}"<<endl;
    	file_latex<<"\\hline"<<endl;
    	file_latex<<" & "; //Leave leftmost case empty
    	for(int isample=0; isample<v_label.size(); isample++) //Declare all process names
    	{
            if(v_label[isample].Contains("TTbar") || v_label[isample].Contains("DY")) {continue;}
            else if(v_samples[isample] == "DATA") {continue;}
            else if(!use_privSamples_asSignal && v_samples[isample].Contains("PrivMC")) {continue;} //Yield tables: consider central samples
            else if(use_privSamples_asSignal && (v_samples[isample]=="tZq" || v_samples[isample]=="ttZ" || v_samples[isample]=="tWZ")) {continue;} //Yield tables: consider private samples

            if(isample == v_label.size()-1 || v_label[isample] != v_label[isample+1])
    		{
                // if(v_label[isample] == "ttZ") {file_latex<<"$t\\bar{t}Z$ & ";}
                if(v_label[isample].Contains("tZq")) {file_latex<<"$\\tZq$ & ";}
                else if(v_label[isample].Contains("ttZ")) {file_latex<<"$\\ttZ$ & ";}
                else if(v_label[isample].Contains("tWZ")) {file_latex<<"$\\tWZ$ & ";}
                else if(v_label[isample] == "ttW") {file_latex<<"$t\\ttW$ & ";}
                else if(v_label[isample] == "ttH") {file_latex<<"$t\\ttH$ & ";}
                else if(v_label[isample] == "tX") {file_latex<<"$t(\\tX$ & ";}
                else if(v_label[isample] == "WZ") {file_latex<<"$\\WZ$ & ";}
                else if(v_label[isample] == "VVV") {file_latex<<"$\\VVV$ & ";}
                else if(v_label[isample] == "XG") {file_latex<<"$\\Xg$ & ";}
    			else {file_latex<<""<<v_label[isample]<<" & ";}
    		}
    	}
        file_latex<<"Total SM";
        if(!blind) {file_latex<<" & Data";}
    	file_latex<<"\\\\ \\hline"<<endl;

        // file_latex<<" & "; //Leftmost column case empty
        TString cat_latex_name = Get_Category_LatexName(category); //Leftmost column case = category name
        file_latex<<cat_latex_name + " & ";
    } //Latex

//--------------------------------------------

//--------------------------------------------
	ofstream file_out(outname.Data());
	file_out<<"## Yields  in "<<category<<" category, "<<channel<<" channel ("<<lumi<<") ##"<<endl;
	file_out<<"____________________________________________"<<endl;
	file_out<<"____________________________________________"<<endl;

	//NB : don't declare inside sample loop, cause might want to sum samples
    vector<double> v_yields_proc_allYears(v_samples.size(), 0), v_statErr_proc_allYears(v_samples.size(), 0); //Sum yields per process for all years
    double yield_tmp = 0;
	double yield_signals = 0;
	double yield_bkg = 0;
	double yield_DATA = 0;

	double statErr_tmp = 0;
	double statErr_signals = 0;

 // #   # ######   ##   #####     #       ####   ####  #####
 //  # #  #       #  #  #    #    #      #    # #    # #    #
 //   #   #####  #    # #    #    #      #    # #    # #    #
 //   #   #      ###### #####     #      #    # #    # #####
 //   #   #      #    # #   #     #      #    # #    # #
 //   #   ###### #    # #    #    ######  ####   ####  #

    for(int iyear=0; iyear<v_years.size(); iyear++)
    {
        TString dir_ntuples = NTUPLEDIR + v_years[iyear] + "/"; //NTUPLEDIR is defined in Utils/Helper.h
        // TString dir_ntuples = NTUPLEDIR + v_years[iyear] + "/";

    	//FIRST LOOP ON SAMPLES : check here if files are missing ; else, may interfer with summing of several processes (TTZ, Rares, ...)
    	for(int isample=0; isample<v_samples.size(); isample++)
    	{
    		TString filepath = dir_ntuples + v_samples[isample]+".root";
    		// cout<<"-- File "<<filepath<<endl;
    		if(!Check_File_Existence(filepath) )
    		{
    			//ERASE MISSING SAMPLES FROM VECTORS
    			v_samples.erase(v_samples.begin() + isample);
    			v_label.erase(v_label.begin() + isample);

    			cout<<FRED("File "<<filepath<<" not found ! Erased index '"<<isample<<"' from vectors")<<endl;
    		}
    	}


    //  ####    ##   #    # #####  #      ######    #       ####   ####  #####
    // #       #  #  ##  ## #    # #      #         #      #    # #    # #    #
    //  ####  #    # # ## # #    # #      #####     #      #    # #    # #    #
    //      # ###### #    # #####  #      #         #      #    # #    # #####
    // #    # #    # #    # #      #      #         #      #    # #    # #
    //  ####  #    # #    # #      ###### ######    ######  ####   ####  #

    	for(int isample=0; isample<v_samples.size(); isample++)
    	{
    		TString filepath = dir_ntuples + v_samples[isample]+".root";
    		cout<<"-- File "<<filepath<<endl;

    		if(!Check_File_Existence(filepath) )
    		{
    			cout<<FRED("File "<<filepath<<" not found !")<<endl;
    			continue;
    		}

    		// cout<<FBLU("Sample : "<<v_samples[isample]<<"")<<endl;

            TString treename = "result";
    		TFile* f = new TFile(filepath);
    		TTree* t = NULL;
            t = (TTree*) f->Get(treename);
            if(!t) {cout<<FRED("Tree '"<<treename<<"' not found ! Skip !")<<endl; continue;}
            t->SetBranchStatus("*", 0); //disable all branches, speed up reading

            TH1F* h_SWE = NULL;
            vector<float> v_SWE;
            vector<float>* v_reweights_floats = NULL;
            vector<string>* v_reweights_ids = NULL;
            int idx_sm = -1;
            //For private MC samples, read and store sums of weights (SWE), and read rese
            if(v_samples[isample].Contains("Priv"))
            {
                //Read branches
                v_reweights_floats = new vector<float>();
                v_reweights_ids = new vector<string>();
                t->SetBranchStatus("mc_EFTweightIDs", 1);
                t->SetBranchAddress("mc_EFTweightIDs", &v_reweights_ids);
                t->SetBranchStatus("mc_EFTweights", 1);
                t->SetBranchAddress("mc_EFTweights", &v_reweights_floats);

                //Find SM index in vectors
                t->GetEntry(0); //Read 1 entry
                for(int iwgt=0; iwgt<v_reweights_ids->size(); iwgt++)
                {
                    // cout<<"iwgt "<<iwgt<<" / "<<v_reweights_ids->at(iwgt)<<endl;

                    TString ts = v_reweights_ids->at(iwgt);
                    if(ts.Contains("_sm", TString::kIgnoreCase) ) {idx_sm = iwgt; t->SetBranchStatus("mc_EFTweightIDs", 0); break;} //Found relevant index, no need to read this branch anymore
                    else if(ts.Contains("EFTrwgt183_", TString::kIgnoreCase) ) {idx_sm = iwgt; t->SetBranchStatus("mc_EFTweightIDs", 0); break;} //TOP19001 convention
                }
                if(v_samples[isample].Contains("_c")) {idx_sm = -1;} //Pure-EFT sample
                if(idx_sm == -1) {cout<<FRED("SM point not found !")<<endl;}
                // cout<<"idx_sm "<<idx_sm<<endl;

                h_SWE = (TH1F*) f->Get("EFT_SumWeights");
                if(!h_SWE) {cout<<FRED("EFT_SumWeights not found ! ")<<endl;}
                for(int ibin=1; ibin<=h_SWE->GetNbinsX(); ibin++)
                {
                    v_SWE.push_back(h_SWE->GetBinContent(ibin)); //1 SWE stored for each stored weight
                    // v_SWE.push_back(h_SWE->GetBinContent(ibin+1));
                    // cout<<"v_SWE["<<ibin<<"] = "<<v_SWE[ibin]<<endl;
                }
            }

    		Double_t weight = 1., weight_avg = 0.; //Event weight (gen-level weight, smeared by systematics)
            Float_t eventMCFactor, weightMENominal; //Sample-dependent factor computed at Potato-level (lumi*xsec/SWE)
            t->SetBranchStatus("eventWeight", 1);
    		t->SetBranchAddress("eventWeight", &weight);
            t->SetBranchStatus("eventMCFactor", 1);
    		t->SetBranchAddress("eventMCFactor", &eventMCFactor);
            t->SetBranchStatus("weightMENominal", 1);
    		t->SetBranchAddress("weightMENominal", &weightMENominal);

            Float_t mTW; //Debugging
            t->SetBranchStatus("mTW", 1);
    		t->SetBranchAddress("mTW", &mTW);

            //--- Cut on relevant event selection (e.g. 3l SR, ttZ CR, etc.) -- stored as Char_t
            Char_t is_goodCategory; //Categ. of event
            TString category_tmp = category; //NB: could use function 'Get_Category_Boolean_Name' as well... but expects different argument (not the flag name itself)
            if(v_samples[isample].Contains("NPL", TString::kIgnoreCase) || v_samples[isample].Contains("DY", TString::kIgnoreCase) || v_samples[isample].Contains("ttbar", TString::kIgnoreCase)) {category_tmp+= "Fake";} //Different flags for fakes

            if(category_tmp != "")
            {
                t->SetBranchStatus(category_tmp, 1);
                t->SetBranchAddress(category_tmp, &is_goodCategory);
            }

            // UInt_t njets, nbjets;
            // t->SetBranchStatus("nJets", 1);
    		// t->SetBranchAddress("nJets", &njets);
            // t->SetBranchStatus("nBJets", 1);
    		// t->SetBranchAddress("nBJets", &nbjets);

            //Study impact of SFs
            Double_t weightPU, weightPrefire, weightMuonLoose, weightMuonTight, weightElectronLoose, weightElectronTight, btagEventWeight[5];
            if(remove_totalSF && v_samples[isample] != "DATA")
            {
                t->SetBranchStatus("weightPU", 1); t->SetBranchAddress("weightPU", &weightPU);
                t->SetBranchStatus("weightPrefire", 1); t->SetBranchAddress("weightPrefire", &weightPrefire);
                t->SetBranchStatus("weightMuonLoose", 1); t->SetBranchAddress("weightMuonLoose", &weightMuonLoose);
                t->SetBranchStatus("weightMuonTight", 1); t->SetBranchAddress("weightMuonTight", &weightMuonTight);
                t->SetBranchStatus("weightElectronLoose", 1); t->SetBranchAddress("weightElectronLoose", &weightElectronLoose);
                t->SetBranchStatus("weightElectronTight", 1); t->SetBranchAddress("weightElectronTight", &weightElectronTight);
                t->SetBranchStatus("btagEventWeight", 1); t->SetBranchAddress("btagEventWeight", &btagEventWeight);
            }

            float chan;
            t->SetBranchStatus("channel", 1);
            t->SetBranchAddress("channel", &chan);


     // ###### #    # ###### #    # #####    #       ####   ####  #####
     // #      #    # #      ##   #   #      #      #    # #    # #    #
     // #####  #    # #####  # #  #   #      #      #    # #    # #    #
     // #      #    # #      #  # #   #      #      #    # #    # #####
     // #       #  #  #      #   ##   #      #      #    # #    # #
     // ######   ##   ###### #    #   #      ######  ####   ####  #

    		int nentries = t->GetEntries();
    		for(int ientry=0; ientry<nentries; ientry++)
    		{
    			weight=1.; eventMCFactor = 1.;

    			t->GetEntry(ientry);

                if(channel == "uuu" && chan != 0) {continue;}
                if(channel == "uue" && chan != 1) {continue;}
                if(channel == "eeu" && chan != 2) {continue;}
                if(channel == "eee" && chan != 3) {continue;}

                //--- Cut on category value
                if(category != "" && !is_goodCategory) {continue;}

                // if(v_samples[isample] == "PrivMC_tWZ" && weight>10) {continue;}

                //Divide event weight by SF (<-> remove SFs), equivalent to using the MENominal weight instead of eventWeight
                if(v_years[iyear] == "2018") {weightPrefire=1;} //no prefire in 2018
                if(!weightPU) {weightPU=1;}
                if(remove_totalSF && v_samples[isample] != "DATA")
                {
                    weight = weightMENominal;

                    // float total_SF = weightPU*weightPrefire*weightMuonLoose*weightMuonTight*weightElectronLoose*weightElectronTight*btagEventWeight[0];
                    // weight/= total_SF;
                    // if(isnan(total_SF) || isinf (total_SF) || total_SF == 0)
                    // {
                    //     cout<<"weightPU "<<weightPU<<endl;
                    //     cout<<"weightPrefire "<<weightPrefire<<endl;
                    //     cout<<"weightMuonLoose "<<weightMuonLoose<<endl;
                    //     cout<<"weightMuonTight "<<weightMuonTight<<endl;
                    //     cout<<"weightElectronLoose "<<weightElectronLoose<<endl;
                    //     cout<<"weightElectronTight "<<weightElectronTight<<endl;
                    //     cout<<"btagEventWeight "<<btagEventWeight[0]<<endl;
                    // }
                }

                //Private MC sample : need to do some rescaling
                if(v_samples[isample].Contains("Priv") && idx_sm != -1)
                {
                    if(v_reweights_floats->size()==0) {cout<<"v_reweights_floats->size()==0 "<<endl; continue;} //Protection (should never happen)
                    // cout<<"v_reweights_floats.size() "<<v_reweights_floats.size()<<endl;

                    //--- SM reweight
                    //Factor (weight / weightMENominal) should account for the central systematics (applied to 'eventWeight')
                    weight*= v_reweights_floats->at(idx_sm) / (weightMENominal * v_SWE[idx_sm]); //with SFs

                    if(remove_totalSF) {weight = v_reweights_floats->at(idx_sm) / v_SWE[idx_sm];}  //no SF, basic formula

                    //-- Tmp fixes to xsec
                    // if(v_samples[isample] == "PrivMC_ttZ_TOP19001") {weight*= 2.482*20;} //wrong eventMCFactor and wrong SWEs
                    if(v_samples[isample] == "PrivMC_tZq_TOP19001") {weight*= 3.087*20;}
                }

                if(isnan(weight*eventMCFactor) || isinf(weight*eventMCFactor))
                {
                    cout<<BOLD(FRED("* Found event with weight*eventMCFactor = "<<weight<<"*"<<eventMCFactor<<" ; remove it..."))<<endl; break; //continue;
                }
                // else if(!weight*eventMCFactor) {cout<<"weight*eventMCFactor = "<<weight<<"*"<<eventMCFactor<<" ! Is it expected ?"<<endl;}

                //After sanity checks, can compute final event weight
                weight*= eventMCFactor;

                v_yields_proc_allYears[isample]+= weight;
                v_statErr_proc_allYears[isample]+= weight*weight;

                if(v_samples[isample] == "DATA") //DATA
    			{
    				yield_DATA+= weight;
    				// statErr_DATA+= weight;
    			}
                else //MC
                {
                    yield_tmp+= weight; statErr_tmp+= weight*weight;

                    if(v_label[isample] == signal || v_label[isample] == "signal" || (signal=="signal" && ( (use_privSamples_asSignal && v_label[isample].Contains("PrivMC")) || (!use_privSamples_asSignal && (v_label[isample]=="tZq" || v_label[isample]=="ttZ" || v_label[isample]=="tWZ")) ) ) ) //Signals, group together
                    {
                        yield_signals+= weight;
                        statErr_signals+= weight*weight;
                        // cout<<"yield_signals "<<yield_signals<<endl;
                    }
                    else if(v_samples[isample] != "DATA" && !v_samples[isample].Contains("tZq") && !v_samples[isample].Contains("ttZ") && !v_samples[isample].Contains("tWZ") && !v_samples[isample].Contains("TTbar") && !v_samples[isample].Contains("DY")) //Backgrounds //Don't consider: signals / private samples / MC fakes / ...
                    {
                        yield_bkg+= weight;
                    }
                }
            } //event loop

    		// cout<<"yield_bkg = "<<yield_bkg<<endl;
    		// cout<<"yield_tmp "<<yield_tmp<<endl;

    		//Check if restart counter, or merge processes
            if(lumi != "Run2")
            {
                if(v_samples[isample] != "DATA" && (isample == v_label.size()-1 || v_label[isample] != v_label[isample+1]))
                {
                    file_out<<"--------------------------------------------"<<endl;
                    file_out<<left<<setw(25)<<v_label[isample]<<setprecision(4)<<yield_tmp;
                    // file_out<<v_label[isample]<<"\\t"<<yield_tmp;

                    file_out<<" (+/- "<<sqrt(statErr_tmp)<<" stat.)"<<endl;
                    // cout<<left<<setw(25)<<v_label[isample]<<yield_tmp<<endl;

                    if(create_latex_table && !v_label[isample].Contains("TTbar") && !v_label[isample].Contains("DY") && !v_label[isample].Contains("DATA") && v_label[isample]!="tZq" && v_label[isample]!="ttZ" && v_label[isample]!="tWZ")
                    // if(create_latex_table && !v_label[isample].Contains("TTbar") && !v_label[isample].Contains("DY") && !v_label[isample].Contains("DATA")  && !v_label[isample].Contains("PrivMC"))
                    {
                        file_latex<<fixed<<setprecision(precision)<<abs(yield_tmp)<<" ($\\pm$"<<fixed<<setprecision(precision)<<sqrt(statErr_tmp)<<") & "; //Single process
                    }

                    yield_tmp = 0; //Reset after writing to file
                    statErr_tmp = 0;
                } //write result
            }

            if(v_samples[isample].Contains("Priv"))
            {
                delete h_SWE;
                delete v_reweights_ids; delete v_reweights_floats;
            }

        } //sample loop

    } //year loop

    if(lumi == "Run2") //Run 2: need to sum quadratically per-year errors
    {
        float yield_currentGroup = 0, statErr_currentGroup = 0;

        for(int isample=0; isample<v_samples.size(); isample++)
        {
            if(v_samples[isample] == "DATA") {continue;} //printed separately

            //Combine errors from all years quadratically
            yield_currentGroup+= v_yields_proc_allYears[isample], statErr_currentGroup+= v_statErr_proc_allYears[isample]; //NB: no need to square uncert. here, it is already ! (cf. above)

            if(group_samples_together && v_samples[isample] != "DATA" && isample < v_label.size()-1 && v_label[isample] == v_label[isample+1]) {continue;} //Sum processes from same group

            //Printout
            file_out<<"--------------------------------------------"<<endl;
            file_out<<left<<setw(25)<<v_label[isample]<<setprecision(4)<<yield_currentGroup;
            // file_out<<v_label[isample]<<"\\t"<<yield_tmp;

            file_out<<" (+/- "<<sqrt(statErr_currentGroup)<<" stat.)"<<endl;
            // cout<<left<<setw(25)<<v_label[isample]<<yield_tmp<<endl;

            if(create_latex_table && !v_label[isample].Contains("TTbar") && !v_label[isample].Contains("DY") && !v_label[isample].Contains("DATA") && v_label[isample]!="tZq" && v_label[isample]!="ttZ" && v_label[isample]!="tWZ")
            // if(create_latex_table && !v_label[isample].Contains("TTbar") && !v_label[isample].Contains("DY") && !v_label[isample].Contains("DATA")  && !v_label[isample].Contains("PrivMC"))
            {
                file_latex<<fixed<<setprecision(precision)<<abs(yield_currentGroup)<<" ($\\pm$"<<fixed<<setprecision(precision)<<sqrt(statErr_currentGroup)<<") & "; //Single process
            }

            yield_currentGroup = 0; statErr_currentGroup = 0; //Reset
        }
    }

	file_out<<endl<<"____________________________________________"<<endl;
	file_out<<"____________________________________________"<<endl;
	file_out<<left<<setw(25)<<"Signal"<<setprecision(5)<<yield_signals;
	file_out<<endl;

    file_out<<endl<<"____________________________________________"<<endl;
	file_out<<"____________________________________________"<<endl;
	file_out<<left<<setw(25)<<"Total background"<<setprecision(5)<<yield_bkg;
	file_out<<endl;

    file_out<<endl<<"____________________________________________"<<endl;
	file_out<<"____________________________________________"<<endl;
	file_out<<left<<setw(25)<<"Total MC"<<setprecision(5)<<yield_signals+yield_bkg;
	file_out<<endl;

	file_out<<endl<<"____________________________________________"<<endl;
	file_out<<"____________________________________________"<<endl;
	file_out<<left<<setw(25)<<"DATA"<<setprecision(5)<<yield_DATA<<endl;
	file_out<<"____________________________________________"<<endl;
	file_out<<"____________________________________________"<<endl;

    if(create_latex_table)
    {
        file_latex<<fixed<<setprecision(precision)<<yield_signals+yield_bkg; //Total SM
        if(!blind) {file_latex<<fixed<<setprecision(precision)<<" & "<<setprecision(precision)<<yield_DATA;} //Data
        file_latex<<" \\\\ \\hline"<<endl;
        file_latex<<"\\end{tabular}"<<endl;
        file_latex<<"\\caption{Event yields for the "<<lumi<<" data-taking period.}"<<endl;
        file_latex<<"\\label{tab:yields}"<<endl;
        // file_latex<<"\\end{table}"<<endl;
        file_latex<<"\\end{sidewaystable}"<<endl;
    }

	cout<<endl<<FYEL("-- Wrote file : "<<outname<<"")<<endl;
	if(create_latex_table) {cout<<FYEL("-- Wrote file : "<<outname_latex<<"")<<endl;}

	return;
}








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
//--------------------------------------------
//--------------------------------------------

int main(int argc, char **argv)
{
    cout<<FYEL("USAGE : ./Yield_Table.exe [region] [2016,2017,2018,Run2]")<<endl<<endl;

//== OPTIONS ==
//--------------------------------------------

    //-- Default args (can be over-riden via command line args)
    TString signal = "signal"; //signal=tZq/ttZ/tWZ

    //-- Category: '' <-> all events ; 'xxx' <-> only include events satisfying condition xxx //E.g.: 'is_signal_SR'
    TString category = "is_signal_SR";
    // TString category = "is_ttz4l_SR";

    TString lumi = "all"; //'2016','2017','2018','Run2,'all''
    TString channel = ""; //'',uuu,uue,eeu,eee
    bool group_samples_together = true; //true <-> group similar samples together
    bool remove_totalSF = false; //SFs are applied to default weights ; can divide weight by total SF again to get nominal weight
    bool use_privSamples_asSignal = true; //true <-> signals yields in the latex table correspond to private samples (else central)
    bool process_samples_byGroup = true; //true <-> read grouped samples (if already hadded together), else read individual samples and combine them when creating histograms if needed (default)

//--------------------------------------------

    TString region = ""; vector<TString> v_lumis(1);
    Apply_CommandArgs_Choices(argc, argv, v_lumis, region); //Get lumi/region via command line
    if(region != "") {category = Get_Category_Boolean_Name(region);}
    if(v_lumis.size() == 3) {lumi = "Run2";}
    else if(v_lumis[0] != "") {lumi = v_lumis[0];}

//--------------------------------------------

	//Sample names and labels //NB : labels must be latex-compatible
	vector<TString> v_samples; vector<TString> v_label;

    //-- Read ntuples merged by sample groups
    if(process_samples_byGroup)
    {
        v_samples.push_back("DATA"); v_label.push_back("DATA");
        v_samples.push_back("PrivMC_tZq"); v_label.push_back("PrivMC_tZq");
        v_samples.push_back("PrivMC_ttZ"); v_label.push_back("PrivMC_ttZ");
        v_samples.push_back("PrivMC_tWZ"); v_label.push_back("PrivMC_tWZ");
        v_samples.push_back("tZq"); v_label.push_back("tZq");
        v_samples.push_back("ttZ"); v_label.push_back("ttZ");
        v_samples.push_back("tWZ"); v_label.push_back("tWZ");
        v_samples.push_back("WZ"); v_label.push_back("WZ");
        v_samples.push_back("tX"); v_label.push_back("tX");
        v_samples.push_back("VVV"); v_label.push_back("VVV");
        v_samples.push_back("XG"); v_label.push_back("XG");
        v_samples.push_back("NPL"); v_label.push_back("NPL");
    }
    //-- Read individual ntuples
    else
    {
        v_samples.push_back("DATA"); v_label.push_back("DATA");

        v_samples.push_back("tZq"); v_label.push_back("tZq");
        v_samples.push_back("ttZ"); v_label.push_back("ttZ");
        v_samples.push_back("tWZ"); v_label.push_back("tWZ");

        v_samples.push_back("PrivMC_tZq"); v_label.push_back("PrivMC_tZq");
        v_samples.push_back("PrivMC_ttZ"); v_label.push_back("PrivMC_ttZ");
        v_samples.push_back("PrivMC_tWZ"); v_label.push_back("PrivMC_tWZ");

        v_samples.push_back("ttZ_M1to10"); v_label.push_back("tX");
        v_samples.push_back("tHq"); v_label.push_back("tX");
        v_samples.push_back("tHW"); v_label.push_back("tX");
        v_samples.push_back("ttH"); v_label.push_back("tX");
        v_samples.push_back("ttW"); v_label.push_back("tX");
        v_samples.push_back("ttZZ"); v_label.push_back("tX");
        v_samples.push_back("ttHH"); v_label.push_back("tX");
        v_samples.push_back("ttWW"); v_label.push_back("tX");
        v_samples.push_back("ttWZ"); v_label.push_back("tX");
        v_samples.push_back("ttZH"); v_label.push_back("tX");
        v_samples.push_back("ttWH"); v_label.push_back("tX");
        v_samples.push_back("tttt"); v_label.push_back("tX");

        v_samples.push_back("WZ"); v_label.push_back("WZ");

        v_samples.push_back("ZZ4l"); v_label.push_back("VVV");
        v_samples.push_back("ZZZ"); v_label.push_back("VVV");
        v_samples.push_back("WZZ"); v_label.push_back("VVV");
        v_samples.push_back("WWW"); v_label.push_back("VVV");
        v_samples.push_back("WWZ"); v_label.push_back("VVV");

        v_samples.push_back("TTGamma_Dilep"); v_label.push_back("XG");
        v_samples.push_back("tGJets"); v_label.push_back("XG");
        v_samples.push_back("WGToLNuG"); v_label.push_back("XG");
        v_samples.push_back("ZGToLLG_01J"); v_label.push_back("XG");

        //MC nonprompt fakes
        v_samples.push_back("TTbar_DiLep"); v_label.push_back("TTbar_DiLep");
        v_samples.push_back("TTbar_SemiLep"); v_label.push_back("TTbar_SemiLep");
        v_samples.push_back("DY"); v_label.push_back("DY");

        //DD NPL (substract MC prompt contribution)
        v_samples.push_back("NPL_DATA"); v_label.push_back("NPL");
        v_samples.push_back("NPL_MC"); v_label.push_back("NPL");
    }

    //TMP
    // v_samples.push_back("PrivMC_tZq_TOP19001"); v_label.push_back("PrivMC_tZq_TOP19001");
    // v_samples.push_back("PrivMC_ttZ_TOP19001"); v_label.push_back("PrivMC_ttZ_TOP19001");
    // v_samples.push_back("PrivMC_tWZ_PSweights"); v_label.push_back("PrivMC_tWZ_PSweights");
    // v_samples.push_back("PrivMC_ttZ_PSweights"); v_label.push_back("PrivMC_ttZ_PSweights");

//--------------------------------------------

    if(lumi == "all")
    {
        Compute_Write_Yields(v_samples, v_label, category, signal, "2016", group_samples_together, remove_totalSF, channel, use_privSamples_asSignal);
        Compute_Write_Yields(v_samples, v_label, category, signal, "2017", group_samples_together, remove_totalSF, channel, use_privSamples_asSignal);
        Compute_Write_Yields(v_samples, v_label, category, signal, "2018", group_samples_together, remove_totalSF, channel, use_privSamples_asSignal);
        Compute_Write_Yields(v_samples, v_label, category, signal, "Run2", group_samples_together, remove_totalSF, channel, use_privSamples_asSignal); //should sum all years
    }
    else {Compute_Write_Yields(v_samples, v_label, category, signal, lumi, group_samples_together, remove_totalSF, channel, use_privSamples_asSignal);}

	return 0;
}
