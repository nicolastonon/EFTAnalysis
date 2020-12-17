#include "ROC_Plotter.h"

int main(int argc, char **argv)
{
	Load_Canvas_Style(); //Styling

//-- Select options
//--------------------------------------------
	vector<TString> v_processes; //For DNNs, must declare all process classes to define 'background' properly //Not used for BDT (signal vs bkg)
    // v_processes.push_back("tZq");
    // v_processes.push_back("ttZ");
    // v_processes.push_back("Backgrounds");

    v_processes.push_back("EFT");
    v_processes.push_back("SM");
    // v_processes.push_back("refPoint");

    TString lumiYear = "2017"; //'2016,'2017','2018','Run2'
    // TString lumiYear = "Run2"; //'2016,'2017','2018','Run2'

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

    v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+".root");
    v_Filelabel.push_back("NN "+lumiYear);
    v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);

    v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+"_recoZ_Pt.root");
    v_Filelabel.push_back("Z_pt");
    v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);

    // v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+"_recoZ_Eta.root");
    // v_Filelabel.push_back("Z_eta");
    // v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);

    // v_filepath.push_back("../outputs/NN_"+v_processes[0]+"_"+lumiYear+"_dEta_tjprime.root");
    // v_Filelabel.push_back("dEta_tjprime");
    // v_isTMVA_file.push_back("Keras"); v_isTrainSample.push_back(false);

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
