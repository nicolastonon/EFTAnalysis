/**
 * CODE USED TO AUTOMATE LIMIT COMPUTATION, READ RESULTS AND COMPARE/SORT THEM FOR THQ ANALYSIS OPTIMIZATION
 * ADAPT TO YOUR NEEDS !
 */

/* BASH COLORS */
#define RST   "[0m"
#define KRED  "[31m"
#define KGRN  "[32m"
#define KYEL  "[33m"
#define KBLU  "[34m"
#define KMAG  "[35m"
#define KCYN  "[36m"
#define KWHT  "[37m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST
#define BOLD(x) "[1m" x RST
#define UNDL(x) "[4m" x RST

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TTree.h"
#include "TString.h"
#include "TColor.h"
#include "TCut.h"
#include "TLegend.h"
#include "TLine.h"
#include "THStack.h"
#include "TString.h"
#include "TObject.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>

//#include "TMVA/Tools.h"
//#include "TMVA/Factory.h"
//#include "TMVA/Reader.h"
//#include "TMVA/MethodCuts.h"
//#include "TMVA/Config.h"

#include <cassert>     //Can be used to terminate program if argument is not true.
//Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> // to be able to use mkdir

using namespace std;

//Use stat function (from library sys/stat) to check if a file exists
bool Check_File_Existence(const TString& name)
{
  struct stat buffer;
  return (stat (name.Data(), &buffer) == 0); //true if file exists
}

//Move file with bash command 'mv'
void MoveFile(TString origin_path, TString dest_path)
{
	TString command = "mv "+origin_path + " " + dest_path;
	system(command.Data() );

	return;
}

//Copy file with bash command 'cp'
void CopyFile(TString origin_path, TString dest_path)
{
	TString command = "cp "+origin_path + " " + dest_path;
	system(command.Data() );

	return;
}

//Convert a double into a TString
// precision --> can choose if TString how many digits the TString should display
TString Convert_Number_To_TString(double number, int precision/*=3*/)
{
	stringstream ss;
	ss << std::setprecision(precision) << number;
	TString ts = ss.str();
	return ts;
}

TString Compute_FitResult_From_TemplateFile(TString type, TString path_templatefile, bool expected, bool use_syst)
{
	CopyFile(path_templatefile, "../templates/Templates_BDT2Dlin_3l_SR.root"); //Copy file to templates dir.

	TString datacard_path = "../datacards/COMBINED_Datacard_TemplateFit_3l";

	if(!use_syst) datacard_path+= "_noShape";

	datacard_path+=".txt";

	cout<<FBLU("USING DATACARD : "<<datacard_path<<" ")<<endl;

	TString f_tmp_name = "limit_expected_info_tmp.txt";

	// if(expected) {system( ("combine -M ProfileLikelihood --significance -t -1 --expectSignal=1 "+datacard_path+" | tee "+f_tmp_name).Data() );} //Expected a-priori
	// else {system( ("combine -M ProfileLikelihood --signif --cminDefaultMinimizerType=Minuit2 "+datacard_path+" | tee "+f_tmp_name).Data() );} //Observed

	TString cmd = "";
	if(type == "limit")
	{
		if(expected) {cmd = "combine -M AsymptoticLimits --noFitAsimov --run blind " + datacard_path+" | tee "+f_tmp_name;}
		else {cmd = "combine -M AsymptoticLimits " + datacard_path+" | tee "+f_tmp_name;}
	}
	else if(type == "significance")
	{
		if(expected) {cmd = "combine -M Significance -t -1 --expectSignal 1 " + datacard_path+" | tee "+f_tmp_name;}
		else {cmd = "combine -M Significance " + datacard_path+" | tee "+f_tmp_name;}
	}
	else {cout<<FRED("ERROR : wrong 'type' arg : "<<type<<"")<<endl; return "XXX";}

	system(cmd.Data());

	ifstream file_in(f_tmp_name.Data()); //Store output from combine in text file, read it & extract value
	string line;
	TString ts;
	while(!file_in.eof())
	{
		getline(file_in, line);
		ts = line;

		if(type == "limit")
		{
			if(!ts.Contains("Expected 50.0")) {continue;}
		}
		else
		{
			if(!ts.Contains("Significance:")) {continue;}
		}

		int index = ts.First(": ");
		ts.Remove(0, index);
		break;
	}

	system( ("rm "+f_tmp_name).Data() ); //Remove tmp file

	return ts;
}

void Compare_AllConfigs(TString type, int nof_binConfigs, TString dir_path, bool expected, bool use_syst)
{
	ofstream file_out("Comparison_allLimits.txt");
	file_out<<"-- COMPARE LIMITS FROM COMBINE --"<<endl;
	file_out<<"//--------------------------------------------//"<<endl<<endl;

	// cout<<"nof_binConfigs "<<nof_binConfigs<<endl;
	for(int i=0; i<nof_binConfigs; i++)
	{
		cout<<FYEL("** CONFIG "<<i<<" : ")<<endl;

		TString template_path = dir_path + "Templates_BDT2Dlin_3l_SR_bin"+Convert_Number_To_TString(i, 5)+".root";
		if(!Check_File_Existence(template_path)) {cout<<FRED("Error : file "<<template_path<<" not found !")<<endl; continue;}
		cout<<"-- Will read template file : "<<template_path<<" ......"<<endl<<endl<<endl;

		file_out<<"Bin config "<<i<<" ---> "<<Compute_FitResult_From_TemplateFile(type, template_path, expected, use_syst)<<endl<<endl;
	}

	return;
}

int main()
{
	TString dir_path = "./test/";

	// int nof_configs = 16;
	int nof_configs = 150;
	// int nof_configs = 441;

    TString type = "significance"; //'significance' or 'limit'
	bool expected = true;
	bool use_syst = true;

	Compare_AllConfigs(type, nof_configs, dir_path, expected, use_syst);
}
