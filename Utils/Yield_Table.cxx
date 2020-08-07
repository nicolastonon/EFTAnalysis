#include "Helper.h"

using namespace std;

//--------------------------------------------
// ##    ## #### ######## ##       ########     ########    ###    ########  ##       ########
//  ##  ##   ##  ##       ##       ##     ##       ##      ## ##   ##     ## ##       ##
//   ####    ##  ##       ##       ##     ##       ##     ##   ##  ##     ## ##       ##
//    ##     ##  ######   ##       ##     ##       ##    ##     ## ########  ##       ######
//    ##     ##  ##       ##       ##     ##       ##    ######### ##     ## ##       ##
//    ##     ##  ##       ##       ##     ##       ##    ##     ## ##     ## ##       ##
//    ##    #### ######## ######## ########        ##    ##     ## ########  ######## ########
//--------------------------------------------

void Compute_Write_Yields(vector<TString> v_samples, vector<TString> v_label, TString category, TString signal, TString lumi, bool group_samples_together, bool remove_totalSF, TString channel)
{
    cout<<endl<<BYEL("                          ")<<endl<<endl;
	cout<<FYEL("--- Will count the yields for each sample ---")<<endl;
	cout<<"(category : "<<category<<" / lumi : "<<lumi<<" / channel : "<<channel<<")"<<endl;
    cout<<endl<<BYEL("                          ")<<endl<<endl;

    mkdir("./outputs/", 0777);
    mkdir("./outputs/yields", 0777);
    // mkdir("./outputs/yields/latex", 0777);

    TString outname = "./outputs/yields/Yields_"+category+"_"+lumi;
    if(channel != "") {outname+= "_" + channel;}
    outname+= ".txt";
    // TString outname_latex = "./outputs/yields/latex/Yields_"+category+"_"+lumi;
    // if(channel != "") {outname_latex+= "_" + channel;}
    // outname_latex+= ".txt";

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
/*
    ofstream file_latex(outname_latex.Data()); //NEW : directly write yields into latex table

	//NB : '\' character must be escaped as '\\' in C++!
    file_latex<<"\\begin{table}[]"<<endl;
    file_latex<<"\\centering"<<endl;
    file_latex<<"\\begin{tabular}{|c|";
	for(int isample=0; isample<v_label.size(); isample++)
	{
		// if(v_samples[isample] == "QFlip" && (nLep != "2l" || subcat == "mm" || subcat == "uu")) {continue;}
		// if(v_samples[isample].Contains("GammaConv") && (subcat == "mm" || subcat == "uu")) {continue;}

		if(v_samples[isample] != "DATA" && !v_label[isample].Contains("SM") && isample < v_label.size()-1 && v_label[isample] != v_label[isample+1])
		{
			file_latex<<"c|"; //1 column per process
		}
	}
	file_latex<<"c|"; //also add total bkg
	file_latex<<"}"<<endl;
	file_latex<<"\\hline"<<endl;
	file_latex<<" & ";
	for(int isample=0; isample<v_label.size(); isample++)
	{
		// if(v_samples[isample] == "QFlip" && (nLep != "2l" || subcat == "mm" || subcat == "uu")) {continue;}
		// if(v_samples[isample].Contains("GammaConv") && (subcat == "mm" || subcat == "uu")) {continue;}

		if(v_samples[isample] != "DATA" && !v_label[isample].Contains("SM") && isample < v_label.size()-1 && v_label[isample] != v_label[isample+1])
		{
			if(v_label[isample] == "TTZ") {file_latex<<"$\\mathbf{t\\bar{t}Z}$ & ";}
			else if(v_label[isample] == "ttW") {file_latex<<"$\\mathbf{t\\bar{t}W}$ & ";}
			else {file_latex<<"\\textbf{"<<v_label[isample]<<"} & ";}
		}
	}
	file_latex<<"Total SM"; //1 column per process
	file_latex<<"\\\\ \\hline"<<endl;
	if(subcat != "")
    {
        file_latex<<"\\textbf{";
        if(subcat == "uu") {file_latex<<"\\mumu";}
        else if(subcat == "ue") {file_latex<<"\\emu";}
        else if(subcat == "ee") {file_latex<<"ee";}
        else {file_latex<<nLep;}
        file_latex<<"} & ";
    }
	else {file_latex<<"\\textbf{"<<nLep<<"} & ";}
*/
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
        TString dir_ntuples = "./input_ntuples/" + v_years[iyear] + "/";
        // TString dir_ntuples = "./input_ntuples/top19001_3lSR/";

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
    		TTree* t = 0;
            t = (TTree*) f->Get(treename);
            if(!t) {cout<<FRED("Tree '"<<treename<<"' not found ! Skip !")<<endl; continue;}
            t->SetBranchStatus("*", 0); //disable all branches, speed up reading

            TH1F* h_SWE = 0;
            vector<float> v_SWE;
            vector<float>* v_reweights_floats = 0;
            vector<string>* v_reweights_ids = 0;
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
            // if(category != "")
            // {
            //     TString cat_name = Get_Category_Boolean_Name(category);
            //     t->SetBranchStatus(cat_name, 1);
            //     t->SetBranchAddress(cat_name, &is_goodCategory);
            // }

            if(category != "")
            {
                t->SetBranchStatus(category, 1);
                t->SetBranchAddress(category, &is_goodCategory);
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

                //-- test cuts on njets
                // if(njets != 2 || nbjets != 2) {continue;}

                //Private MC sample : need to do some rescaling
                if(v_samples[isample].Contains("Priv") && idx_sm != -1)
                {
                    //--- SM reweight
                    //Factor (weight / weightMENominal) should account for the central systematics (applied to 'eventWeight')
                    weight*= v_reweights_floats->at(idx_sm) / (weightMENominal * v_SWE[idx_sm]); //with SFs

                    if(remove_totalSF) {weight = v_reweights_floats->at(idx_sm) / v_SWE[idx_sm];}  //no SF, basic formula
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

                    if(v_label[isample] == signal || v_label[isample] == "signal") //Signals, group together
                    {
                        yield_signals+= weight;
                        statErr_signals+= weight*weight;
                        // cout<<"yield_signals "<<yield_signals<<endl;
                    }
                    else if(v_samples[isample] != "DATA" && !v_samples[isample].Contains("Priv")) //Backgrounds
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

                    file_out<<" (+/- "<<statErr_tmp<<" stat.)"<<endl;
                    // cout<<left<<setw(25)<<v_label[isample]<<yield_tmp<<endl;

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

    if(lumi == "Run2")
    {
        float yield_currentGroup = 0, statErr_currentGroup = 0;

        for(int isample=0; isample<v_samples.size(); isample++)
        {
            if(v_samples[isample] == "DATA") {continue;} //printed separately

            yield_currentGroup+= v_yields_proc_allYears[isample], statErr_currentGroup+= v_statErr_proc_allYears[isample]*v_statErr_proc_allYears[isample];

            if(group_samples_together && v_samples[isample] != "DATA" && isample < v_label.size()-1 && v_label[isample] == v_label[isample+1]) //Sum processes from same group
            {
                continue;
            }

            statErr_currentGroup = sqrt(statErr_currentGroup); //compute error

            //Printout
            file_out<<"--------------------------------------------"<<endl;
            file_out<<left<<setw(25)<<v_label[isample]<<setprecision(4)<<yield_currentGroup;
            // file_out<<v_label[isample]<<"\\t"<<yield_tmp;

            file_out<<" (+/- "<<statErr_currentGroup<<" stat.)"<<endl;
            // cout<<left<<setw(25)<<v_label[isample]<<yield_tmp<<endl;

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

	// file_latex<<setprecision(3)<<yield_bkg<<""; //Total bkg
	// file_latex<<" \\\\ \\hline"<<endl;
	// file_latex<<"\\end{tabular}"<<endl;
	// file_latex<<"\\caption{xxx}"<<endl;
	// file_latex<<"\\label{tab:my-table}"<<endl;
	// file_latex<<"\\end{table}"<<endl;

	cout<<endl<<FYEL("-- Wrote file : "<<outname<<"")<<endl;
	// cout<<FYEL("-- Wrote file : "<<outname_latex<<"")<<endl;

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
    cout<<FYEL("USAGE : ./Make_Yield_Table.exe [SR] [2016,2017,2018,all,Run2] [uuu,eeu,uue,eee]")<<endl<<endl;

//== OPTIONS ==
//--------------------------------------------

    //-- Default args (can be over-riden via command line args)
    TString signal = "tZq";
    // TString category = "is_tzq_SR"; //'' <-> all events ; 'xxx' <-> only include events satisfying condition xxx
    TString category = "is_signal_SR"; //'' <-> all events ; 'xxx' <-> only include events satisfying condition xxx
    // TString category = "is_tZq_3l_SR"; //'' <-> all events ; 'xxx' <-> only include events satisfying condition xxx
    // TString category = "is_ttz_SR"; //'' <-> all events ; 'xxx' <-> only include events satisfying condition xxx
    TString lumi = "all"; //'2016','2017','2018','Run2,'all''
    TString channel = ""; //'',uuu,uue,eeu,eee
    bool group_samples_together = true; //true <-> group similar samples together
    bool remove_totalSF = false; //SFs are applied to default weights ; can divide weight by total SF again to get nominal weight

//--------------------------------------------

    if(argc > 1)
	{
        if(!strcmp(argv[1],"2016") || !strcmp(argv[1],"2017") || !strcmp(argv[1],"2018") || !strcmp(argv[1],"Run2")) {lumi = argv[1];}
        // else if(!strcmp(argv[1],"tZq") || !strcmp(argv[1],"ttZ") || !strcmp(argv[1],"tWZ") ) {category = argv[1];}
        else if(!strcmp(argv[1],"uuu") || !strcmp(argv[1],"uue") || !strcmp(argv[1],"eeu") || !strcmp(argv[1],"eee")) {channel = argv[1];}
		else {cout<<"Wrong first arg !"<<endl; return 0;}

        if(argc > 2)
    	{
            if(!strcmp(argv[2],"2016") || !strcmp(argv[2],"2017") || !strcmp(argv[2],"2018") || !strcmp(argv[2],"Run2")) {lumi = argv[2];}
            // else if(!strcmp(argv[2],"tZq") || !strcmp(argv[2],"ttZ") || !strcmp(argv[2],"tWZ") ) {category = argv[2];}
            else if(!strcmp(argv[2],"") || !strcmp(argv[2],"uuu") || !strcmp(argv[2],"uue") || !strcmp(argv[2],"eeu") || !strcmp(argv[2],"eee")) {channel = argv[2];}
    		else {cout<<"Wrong second arg !"<<endl; return 0;}

            if(argc > 3)
        	{
                if(!strcmp(argv[3],"2016") || !strcmp(argv[3],"2017") || !strcmp(argv[3],"2018") || !strcmp(argv[3],"Run2")) {lumi = argv[3];}
                // else if(!strcmp(argv[3],"tZq") || !strcmp(argv[3],"ttZ") || !strcmp(argv[3],"tWZ") ) {category = argv[3];}
                else if(!strcmp(argv[3],"") || !strcmp(argv[3],"uuu") || !strcmp(argv[3],"uue") || !strcmp(argv[3],"eeu") || !strcmp(argv[3],"eee")) {channel = argv[3];}
        		else {cout<<"Wrong second arg !"<<endl; return 0;}
            }
    	}
	}

//--------------------------------------------

	//Sample names and labels //NB : labels must be latex-compatible
	vector<TString> v_samples; vector<TString> v_label;

    v_samples.push_back("DATA"); v_label.push_back("DATA");

    v_samples.push_back("tZq"); v_label.push_back("tZq");
    v_samples.push_back("ttZ"); v_label.push_back("ttZ");

    v_samples.push_back("PrivMC_tZq"); v_label.push_back("PrivMC_tZq");
    v_samples.push_back("PrivMC_ttZ"); v_label.push_back("PrivMC_ttZ");

    v_samples.push_back("tWZ"); v_label.push_back("tWZ");

    v_samples.push_back("tHq"); v_label.push_back("tX");
    v_samples.push_back("tHW"); v_label.push_back("tX");
    v_samples.push_back("tGJets"); v_label.push_back("tX");
    // v_samples.push_back("ST"); v_label.push_back("tX");

    v_samples.push_back("ttH"); v_label.push_back("ttH");
    v_samples.push_back("ttW"); v_label.push_back("ttW");
    v_samples.push_back("ttZZ"); v_label.push_back("ttX");
    v_samples.push_back("ttHH"); v_label.push_back("ttX");
    v_samples.push_back("ttWW"); v_label.push_back("ttX");
    v_samples.push_back("ttWZ"); v_label.push_back("ttX");
    v_samples.push_back("ttZH"); v_label.push_back("ttX");
    v_samples.push_back("ttWH"); v_label.push_back("ttX");
    v_samples.push_back("tttt"); v_label.push_back("ttX");

    v_samples.push_back("WZ"); v_label.push_back("WZ");

    v_samples.push_back("ZZ4l"); v_label.push_back("ZZ");
    v_samples.push_back("ZZZ"); v_label.push_back("VV");
    v_samples.push_back("WZZ"); v_label.push_back("VV");
    v_samples.push_back("WWW"); v_label.push_back("VV");
    v_samples.push_back("WWZ"); v_label.push_back("VV");

    v_samples.push_back("TTGamma_Dilep"); v_label.push_back("Xg");
    v_samples.push_back("tGJets"); v_label.push_back("Xg");
    v_samples.push_back("WGToLNuG"); v_label.push_back("Xg");
    v_samples.push_back("ZGToLLG_01J"); v_label.push_back("Xg");
    v_samples.push_back("ggToZZTo4l"); v_label.push_back("Xg");

    v_samples.push_back("DY"); v_label.push_back("DY");

    v_samples.push_back("TTbar_DiLep"); v_label.push_back("TTbar_DiLep");
    v_samples.push_back("TTbar_SemiLep"); v_label.push_back("TTbar_SemiLep");

//--------------------------------------------

    if(lumi == "all")
    {
        Compute_Write_Yields(v_samples, v_label, category, signal, "2016", group_samples_together, remove_totalSF, channel);
        Compute_Write_Yields(v_samples, v_label, category, signal, "2017", group_samples_together, remove_totalSF, channel);
        Compute_Write_Yields(v_samples, v_label, category, signal, "2018", group_samples_together, remove_totalSF, channel);
        Compute_Write_Yields(v_samples, v_label, category, signal, "Run2", group_samples_together, remove_totalSF, channel); //should sum all years
    }
    else {Compute_Write_Yields(v_samples, v_label, category, signal, lumi, group_samples_together, remove_totalSF, channel);}

	return 0;
}
