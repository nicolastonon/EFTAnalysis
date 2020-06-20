/* BASH CUSTOM */
#define RST   "\e[0m"
#define KRED  "\e[31m"
#define KGRN  "\e[32m"
#define KYEL  "\e[33m"
#define KBLU  "\e[34m"
#define KMAG  "\e[35m"
#define KCYN  "\e[36m"
#define KWHT  "\e[37m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST
#define BOLD(x) "\e[1m" x RST
#define ITAL(x) "\e[3m" x RST
#define UNDL(x) "\e[4m" x RST
#define STRIKE(x) "\e[9m" x RST
#define DIM(x) "\e[2m" x RST
#define DOUBLEUNDERLINE(x) "\e[21m" x RST
#define CURLYUNDERLINE(x) "\e[4:3m" x RST
#define BLINK(x) "\e[5m" x RST
#define REVERSE(x) "\e[7m" x RST
#define INVISIBLE(x) "\e[8m" x RST
#define OVERLINE(x) "\e[53m" x RST
#define TURQUOISE  "\e[38;5;42m"
#define SALMON  "\e[38;2;240;143;104m"
#define FTURQUOISE(x) TURQUOISE x RST
#define FSALMON(x) SALMON x RST
#define YELBKG(x) "\e[43m" x RST
#define CYANBKG(x) "\e[46m" x RST

#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"

using namespace std;



//--------------------------------------------
//  ######   ######  ########  #### ########  ########    ########  #######  ########
// ##    ## ##    ## ##     ##  ##  ##     ##    ##       ##       ##     ## ##     ##
// ##       ##       ##     ##  ##  ##     ##    ##       ##       ##     ## ##     ##
//  ######  ##       ########   ##  ########     ##       ######   ##     ## ########
//       ## ##       ##   ##    ##  ##           ##       ##       ##     ## ##   ##
// ##    ## ##    ## ##    ##   ##  ##           ##       ##       ##     ## ##    ##
//  ######   ######  ##     ## #### ##           ##       ##        #######  ##     ##

// ######## ######## ##     ## ########  ##          ###    ######## ########    ######## #### ########
//    ##    ##       ###   ### ##     ## ##         ## ##      ##    ##          ##        ##     ##
//    ##    ##       #### #### ##     ## ##        ##   ##     ##    ##          ##        ##     ##
//    ##    ######   ## ### ## ########  ##       ##     ##    ##    ######      ######    ##     ##
//    ##    ##       ##     ## ##        ##       #########    ##    ##          ##        ##     ##
//    ##    ##       ##     ## ##        ##       ##     ##    ##    ##          ##        ##     ##
//    ##    ######## ##     ## ##        ######## ##     ##    ##    ########    ##       ####    ##
//--------------------------------------------

/**
 * Produce script containing the commands to produce the datacards (single and combination) automatically
 */
void Script_Datacards_TemplateFit(char include_systematics, char include_statistical, TString template_name, TString region, vector<TString> v_templates, vector<TString> v_channel, vector<TString> v_regions, TString lumiName)
{
    //Check if use shape syst or not
	TString systChoice;
	if(include_systematics == 'y') {systChoice = "withShape";}
	else if(include_systematics == 'n') {systChoice = "noShape";}
    else {cout<<"Wrong arguments ! Abort !"<<endl; return;}

    //Check if use stat. uncert. or not
	TString statChoice;
	if(include_statistical == 'y') {statChoice = "withStat";}
	else if(include_statistical == 'n') {statChoice = "noStat";}
    else {cout<<"Wrong arguments ! Abort !"<<endl; return;}

    //If specific arguments were chosen at command line, modify the vectors defined in the main
    // if(template_name != "0") {v_templates.resize(0); v_templates.push_back(template_name);}
	// if(region != "0") //else if region = 0, use vector as it is defined in main
    // {
    //     if(region == "1")
    //     {
    //         v_regions.resize(0);
    //         v_regions.push_back("xxx");
    //     }
    //     else {v_regions.resize(0); v_regions.push_back(region);}
    // }

    vector<TString> v_lumiYears;
    if(lumiName.Contains("2016") ) {v_lumiYears.push_back("2016");}
    else if(lumiName.Contains("2017") ) {v_lumiYears.push_back("2017");}
    else if(lumiName.Contains("2018") ) {v_lumiYears.push_back("2018");}
    else if(lumiName.Contains("201617") ) {v_lumiYears.push_back("2016"); v_lumiYears.push_back("2017");}
    else if(lumiName.Contains("201618") ) {v_lumiYears.push_back("2016"); v_lumiYears.push_back("2018");}
    else if(lumiName.Contains("201718") ) {v_lumiYears.push_back("2017"); v_lumiYears.push_back("2018");}
    else if(lumiName.Contains("Run2") ) {v_lumiYears.push_back("2016"); v_lumiYears.push_back("2017"); v_lumiYears.push_back("2018");}

	ofstream file_out("makeDatacardsForTemplateFit.sh"); //output script

	TString dir = "./datacards_TemplateFit/";
	// file_out <<"mkdir "<<dir<<endl<<endl;

//--------------------------------------------
//--- First loop over years/regions/templates
//===> Commands to produce single datacards
    for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
    {
        for(int iregion=0; iregion<v_regions.size(); iregion++) //Loop over regions
        {
            // cout<<"v_regions[iregion] "<<v_regions[iregion]<<endl;

    		//Need to add "CR" prefix for CRs
    		TString region_tmp = v_regions[iregion];
    		if(region_tmp == "ttW" || region_tmp == "ttZ" || region_tmp == "WZ" || region_tmp == "ZZ") {region_tmp = "CR_" + region_tmp;}

            for(int itemplate=0; itemplate<v_templates.size(); itemplate++) //Loop over templates
            {
                TString var = v_templates[itemplate] + "_" + region_tmp;

    			//Make subdir for single datacards
    			file_out <<"mkdir "<<dir<<endl<<endl;
                // file_out <<"mkdir "<<dir+"/"+region_tmp<<endl<<endl;
                file_out <<"mkdir "<<dir+v_lumiYears[iyear]<<endl<<endl;

    		    // TString file_histos = "../templates/Combine_Input.root";
                // TString file_histos = "../../../templates/Templates_"+v_templates[itemplate]+"_"+region_tmp+"_"+lumiName+".root";
                TString tmp = (v_templates[itemplate].Contains("NN") ? "NN" : v_templates[itemplate]);
                TString file_histos = "../../../templates/Templates_"+tmp+"_"+region_tmp+"_"+lumiName+".root";

    			cout<<"---> Will use filepath : "<<file_histos<<endl<<endl;

    			for(int ilepchan=0; ilepchan<v_channel.size(); ilepchan++)
    			{
    				file_out<<"python Parser_Datacard_Template.py "
                    + var + " "
                    + v_channel[ilepchan] + " "
                    + v_lumiYears[iyear ] + " "
    				+ file_histos + " "
    				+ systChoice + " "
    				+ statChoice + " "
                    + dir + v_lumiYears[iyear] + " "
                    // + dir + region_tmp + " "
                    // + v_regions[iregion] + " "
    				<<endl;
    			}

    			file_out<<endl<<endl;
        	} //template loop
        } //region loop
    }

	//Command to execute the script which combines the datacards
	file_out<<"python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py ";

//--------------------------------------------
//--- Second loop over years/regions/templates
//===> Give all the single datacards as arguments to the script, for combination
    for(int iyear=0; iyear<v_lumiYears.size(); iyear++)
    {
        for(int iregion=0; iregion<v_regions.size(); iregion++) //Loop over regions
        {
    		//Need to add "CR" prefix for CRs
    		TString region_tmp = v_regions[iregion];
    		if(region_tmp == "ttW" || region_tmp == "ttZ" || region_tmp == "WZ" || region_tmp == "ZZ") {region_tmp = "CR_" + region_tmp;}

    	    for(int itemplate=0; itemplate<v_templates.size(); itemplate++) //Loop over templates
    	    {
                TString var = v_templates[itemplate] + "_" + region_tmp;

    			for(int ilepchan=0; ilepchan<v_channel.size(); ilepchan++) //Loop over channels
    			{
    				file_out<<var;
    				if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
                    file_out<<"_"+v_lumiYears[iyear];
                    file_out<<"=" + dir + v_lumiYears[iyear] + "/"
                    // file_out<<"=" + dir + region_tmp + "/"
    				+ "datacard_"+var;
    				if(v_channel[ilepchan] != "all") {file_out<<"_" + v_channel[ilepchan];}
    				file_out<<".txt ";
    			}
    		} //template loop
    	} //region loop
    } //years loop

	TString output_name = "COMBINED_Datacard_TemplateFit";
    if(systChoice == "noShape") output_name+= "_noShape";
    if(statChoice == "noStat") output_name+= "_noStat";
	output_name+= ".txt";

	file_out<<"> "<<output_name<<endl<<endl;

//datacard for single channels
/*
for(int ichan=0; ichan<v_channel.size(); ichan++)
{
	file_out<<"python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/scripts/combineCards.py ";

	for(int ivar=0; ivar<var_list.size(); ivar++)
	{
		file_out<<var_list[ivar] + "_" + v_channel[ichan] + "=datacard_"+v_channel[ichan]+"_"+var_list[ivar]+".txt ";
	}

	output_name = "COMBINED_Datacard_TemplateFit_";
	if(systChoice == "noShape") output_name+= "noShape_";
	output_name+= v_channel[ichan]+".txt";

	file_out<<"> "<<output_name<<endl<<endl;

	file_out<<"mv "<<output_name<<" datacards_TemplateFit/"<<endl<<endl;
}
*/

	// file_out<<"mv datacard_*.txt datacards_TemplateFit/"<<endl;

	file_out.close();

	system("chmod 755 makeDatacardsForTemplateFit.sh");

	cout<<FGRN("... Created script ./makeDatacardsForTemplateFit.sh !")<<endl;

	cout<<endl<<endl<<FGRN("-- Executing script ./makeDatacardsForTemplateFit.sh !")<<endl;

	system("./makeDatacardsForTemplateFit.sh");

	return;
}







//--------------------------------------------
//  ######   #######  ##     ## ##     ##    ###    ##    ## ########
// ##    ## ##     ## ###   ### ###   ###   ## ##   ###   ## ##     ##
// ##       ##     ## #### #### #### ####  ##   ##  ####  ## ##     ##
// ##       ##     ## ## ### ## ## ### ## ##     ## ## ## ## ##     ##
// ##       ##     ## ##     ## ##     ## ######### ##  #### ##     ##
// ##    ## ##     ## ##     ## ##     ## ##     ## ##   ### ##     ##
//  ######   #######  ##     ## ##     ## ##     ## ##    ## ########

// ##       #### ##    ## ########       ###    ########   ######    ######
// ##        ##  ###   ## ##            ## ##   ##     ## ##    ##  ##    ##
// ##        ##  ####  ## ##           ##   ##  ##     ## ##        ##
// ##        ##  ## ## ## ######      ##     ## ########  ##   ####  ######
// ##        ##  ##  #### ##          ######### ##   ##   ##    ##        ##
// ##        ##  ##   ### ##          ##     ## ##    ##  ##    ##  ##    ##
// ######## #### ##    ## ########    ##     ## ##     ##  ######    ######
//--------------------------------------------

//Ask user to choose options at command line for script generation
void Choose_Arguments_From_CommandLine(char& include_systematics, char& include_statistical, TString& template_name, TString& region, TString& lumiName)
{
    cout<<endl<<FBLU("Choose the luminosity ('Run2'/'2016'/'2017'/'2018'/'201617'/'201618'/'201718') ")<<endl;
    cout<<DIM("('0' <--> 'Run2')")<<endl;
    cin>>lumiName;
    while(lumiName != "0" && lumiName != "Run2" && lumiName != "2016" && lumiName != "2017" && lumiName != "2018" && lumiName != "201617" && lumiName != "201618" && lumiName != "201718")
    {
    	cin.clear();
    	cin.ignore(1000, '\n');

    	cout<<" Wrong answer ! Retry :"<<endl;
    	cin>>lumiName;
    }
    if(lumiName == "0") {lumiName = "Run2";}

	//Choose whether to include shape syst or not
	cout<<endl<<FYEL("--- Do you want to include the *shape* systematics as nuisances in the datacards ? (y/n)")<<endl;
	cin>>include_systematics;
	while(include_systematics != 'y' && include_systematics != 'n')
	{
		cin.clear();
		cin.ignore(1000, '\n');

		cout<<" Wrong answer ! Need to type 'y' or 'n' ! Retry :"<<endl;
		cin>>include_systematics;
	}

    //Choose whether to include statistic uncert. or not (e.g. for pulls of NPs, need to remove them, else too long)
    cout<<endl<<FYEL("--- Do you want to include the *statistical* uncert. in the datacards ? (y/n)")<<endl;
    cin>>include_statistical;
    while(include_statistical != 'y' && include_statistical != 'n')
    {
        cin.clear();
        cin.ignore(1000, '\n');

        cout<<" Wrong answer ! Need to type 'y' or 'n' ! Retry :"<<endl;
        cin>>include_statistical;
    }

    return;
}







//--------------------------------------------
// ##     ##    ###    #### ##    ##
// ###   ###   ## ##    ##  ###   ##
// #### ####  ##   ##   ##  ####  ##
// ## ### ## ##     ##  ##  ## ## ##
// ##     ## #########  ##  ##  ####
// ##     ## ##     ##  ##  ##   ###
// ##     ## ##     ## #### ##    ##
//--------------------------------------------

int main()
{
//Can set options here (can be modified at command line)
//--------------------------------------------
    TString lumiName = "Run2"; //'2016','2017','2018','201617','201618','201718','Run2'

    vector<TString> v_templates; //'NN', 'BDT'
    v_templates.push_back("NN");
    // v_templates.push_back("NN0");

    vector<TString> v_channel; //'all', 'uuu', 'eeu', 'uue', 'eee'
    v_channel.push_back("all");
    // v_channel.push_back("uuu");
    // v_channel.push_back("eeu");
    // v_channel.push_back("uue");
    // v_channel.push_back("eee");

    vector<TString> v_regions;
	v_regions.push_back("SR");

    char include_systematics = 'n';
    char include_statistical = 'n';
//--------------------------------------------


//Automated
//--------------------------------------------
    TString template_name = "0", region = "0";
    Choose_Arguments_From_CommandLine(include_systematics, include_statistical, template_name, region, lumiName);

	Script_Datacards_TemplateFit(include_systematics, include_statistical, template_name, region, v_templates, v_channel, v_regions, lumiName);

	return 0;
}
