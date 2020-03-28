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

#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
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
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>

#include <cassert>     //Can be used to terminate program if argument is not true.
//Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> // to be able to use mkdir

using namespace std;

//--------------------------------------------
// ##     ## ######## ##       ########  ######## ########
// ##     ## ##       ##       ##     ## ##       ##     ##
// ##     ## ##       ##       ##     ## ##       ##     ##
// ######### ######   ##       ########  ######   ########
// ##     ## ##       ##       ##        ##       ##   ##
// ##     ## ##       ##       ##        ##       ##    ##
// ##     ## ######## ######## ##        ######## ##     ##
//--------------------------------------------

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




//--------------------------------------------
//  ######  ####  ######       ######  ######## ########  ######## ##    ##  ######   ######## ##     ##
// ##    ##  ##  ##    ##     ##    ##    ##    ##     ## ##       ###   ## ##    ##     ##    ##     ##
// ##        ##  ##           ##          ##    ##     ## ##       ####  ## ##           ##    ##     ##
//  ######   ##  ##   ####     ######     ##    ########  ######   ## ## ## ##   ####    ##    #########
//       ##  ##  ##    ##           ##    ##    ##   ##   ##       ##  #### ##    ##     ##    ##     ##
// ##    ##  ##  ##    ##     ##    ##    ##    ##    ##  ##       ##   ### ##    ##     ##    ##     ##
//  ######  ####  ######       ######     ##    ##     ## ######## ##    ##  ######      ##    ##     ##
//--------------------------------------------

void Draw_Compatibility_SigStrength(TString signal)
{
    const unsigned int n=5;   // number of mass points
    double x[n];
    double y[n];
    double x_errDown[n];
    double y_errDown[n];
    double x_errUp[n];
    double y_errUp[n];

    if(signal == "Hut")
    {
		x[4] = 0.108679;
        x[3] = -0.0272309;
        x[2] = -0.108448;
        x[1] = -0.0273704;
        x[0] = -0.0365462;
        y[4]	 = 5-0.5;
        y[3]	 = 4-0.5;
        y[2]	 = 3-0.5;
        y[1]	 = 2-0.5;
        y[0]	 = 1-0.5;
		x_errDown[4] = 0.0992573;
        x_errDown[3] = 0.0744194;
        x_errDown[2] = 0.0929672;
        x_errDown[1] = 0.163189;
        x_errDown[0] = 0.0486704;
        y_errDown[4] = 0;
        y_errDown[3] = 0;
        y_errDown[2] = 0;
        y_errDown[1] = 0;
        y_errDown[0] = 0;
		x_errUp[4] = 0.10787;
		x_errUp[3] = 0.0741001;
		x_errUp[2] = 0.083669;
		x_errUp[1] = 0.136481;
		x_errUp[0] = 0.0483706;
        y_errUp[4] = 0;
        y_errUp[3] = 0;
        y_errUp[2] = 0;
        y_errUp[1] = 0;
        y_errUp[0] = 0;
    }
    else if(signal == "Hct")
    {
		x[4] = 0.0745259;
        x[3] = -0.0320763;
        x[2] = -0.0770098;
        x[1] = -0.0166602;
        x[0] = -0.0430082;
        y[4]	 = 5-0.5;
        y[3]	 = 4-0.5;
        y[2]	 = 3-0.5;
        y[1]	 = 2-0.5;
        y[0]	 = 1-0.5;
        x_errDown[4] = 0.131673;
        x_errDown[3] = 0.0902493;
        x_errDown[2] = 0.0911303;
        x_errDown[1] = 0.168721;
        x_errDown[0] = 0.0567197;
        y_errDown[4] = 0;
        y_errDown[3] = 0;
        y_errDown[2] = 0;
        y_errDown[1] = 0;
        y_errDown[0] = 0;
		x_errUp[4] = 0.140377;
		x_errUp[3] = 0.0854555;
		x_errUp[2] = 0.0830211;
		x_errUp[1] = 0.154542;
		x_errUp[0] = 0.05636;
        y_errUp[4] = 0;
        y_errUp[3] = 0;
        y_errUp[2] = 0;
        y_errUp[1] = 0;
        y_errUp[0] = 0;
    }
	else {cout<<"Error ! "<<endl; return;}

    TCanvas* c1 = new TCanvas("", "", 1000, 800);
	c1->SetTopMargin(0.1);

	TH2F* h = new TH2F("", "", 1, -0.3, 0.3, 5, 0, 5);
	h->Draw("same");

	//Text labels
	const char *labels[n]  = {"Combined", "ee", "#mue", "#mu#mu", "3l"};
	for(int i=1;i<=n;i++) {h->GetYaxis()->SetBinLabel(i,labels[i-1]);}
	h->GetYaxis()->LabelsOption("v");
	h->GetXaxis()->SetTitle("#mu_{"+signal+"}");
	h->GetXaxis()->SetLabelSize(0.04);
	h->GetYaxis()->SetLabelSize(0.05);

    TGraphAsymmErrors* gr_obs = new TGraphAsymmErrors(n, x, y, x_errDown, x_errUp, y_errDown, y_errUp);
    // TGraphAsymmErrors(Int_t n, const Double_t* x, const Double_t* y, const Double_t* exl = 0, const Double_t* exh = 0, const Double_t* eyl = 0, const Double_t* eyh = 0)
	gr_obs->SetMarkerStyle(20);
	gr_obs->SetMarkerSize(1.5);
	gr_obs->SetLineColor(kBlue);
	gr_obs->SetLineWidth(2);
    gr_obs->Draw("E0 P same");

	TLine l;
	l.SetLineStyle(2);
	l.DrawLine(0, 0, 0, 5);

    TLatex latex;
    latex.SetNDC();
    TString lumi_13TeV = "41.5 fb^{-1} (13 TeV)";
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.95, 0.92,lumi_13TeV);
    // latex.DrawLatex(0.78, 0.92,lumi_13TeV);

    c1->SaveAs("ccc_"+signal+".png");

    delete c1;
    delete gr_obs;

    return;
}

//--------------------------------------------
// ##       #### ##     ## #### ########  ######
// ##        ##  ###   ###  ##     ##    ##    ##
// ##        ##  #### ####  ##     ##    ##
// ##        ##  ## ### ##  ##     ##     ######
// ##        ##  ##     ##  ##     ##          ##
// ##        ##  ##     ##  ##     ##    ##    ##
// ######## #### ##     ## ####    ##     ######
//--------------------------------------------
void Draw_Compatibility_Limits(TString signal)
{
    const unsigned int n=5;   // number of points
	double x_exp[n];
	double x_obs[n];
    double y[n];
    double y_errDown[n];
	double y_errUp[n];
	double x_errDown_1sigma[n];
    double x_errUp_1sigma[n];
	double x_errDown_2sigma[n];
    double x_errUp_2sigma[n];

    if(signal == "Hut")
    {
		x_exp[0] = 0.0864; //comb
		x_exp[4] = 0.1831; //3l
        x_exp[3] = 0.1387; //uu
        x_exp[2] = 0.1226; //ue
        x_exp[1] = 0.2373; //ee

		x_obs[0] = 0.0718;
		x_obs[4] = 0.3005;
		x_obs[3] = 0.1260;
		x_obs[2] = 0.0957;
		x_obs[1] = 0.2393;

		y[0]	 = 1-0.5;
        y[4]	 = 5-0.5;
        y[3]	 = 4-0.5;
        y[2]	 = 3-0.5;
        y[1]	 = 2-0.5;
		y_errDown[0] = 0.5;
        y_errDown[4] = 0.5;
        y_errDown[3] = 0.5;
        y_errDown[2] = 0.5;
        y_errDown[1] = 0.5;
		y_errUp[0] = 0.5;
		y_errUp[4] = 0.5;
		y_errUp[3] = 0.5;
		y_errUp[2] = 0.5;
		y_errUp[1] = 0.5;

		x_errDown_1sigma[0] = x_exp[0] - 0.0624;
		x_errDown_1sigma[4] = x_exp[4] - 0.1284;
		x_errDown_1sigma[3] = x_exp[3] - 0.1007;
		x_errDown_1sigma[2] = x_exp[2] - 0.0890;
		x_errDown_1sigma[1] = x_exp[1] - 0.1718;
		x_errUp_1sigma[0] = 0.1198 - x_exp[0];
		x_errUp_1sigma[4] = 0.2663 - x_exp[4];
		x_errUp_1sigma[3] = 0.1912 - x_exp[3];
		x_errUp_1sigma[2] = 0.1690 - x_exp[2];
		x_errUp_1sigma[1] = 0.3328 - x_exp[1];

		x_errDown_2sigma[0] = x_exp[0] - 0.0469;
		x_errDown_2sigma[4] = x_exp[4] - 0.0944;
		x_errDown_2sigma[3] = x_exp[3] - 0.0764;
		x_errDown_2sigma[2] = x_exp[2] - 0.0675;
		x_errDown_2sigma[1] = x_exp[1] - 0.1298;
		x_errUp_2sigma[0] = 0.1595 - x_exp[0];
		x_errUp_2sigma[4] = 0.3741 - x_exp[4];
		x_errUp_2sigma[3] = 0.2533 - x_exp[3];
		x_errUp_2sigma[2] = 0.2239 - x_exp[2];
		x_errUp_2sigma[1] = 0.4528 - x_exp[1];
    }
    else if(signal == "Hct")
    {
		x_exp[0] = 0.1050; //comb
		x_exp[4] = 0.2344; //3l
		x_exp[3] = 0.1694; //uu
		x_exp[2] = 0.1426; //ue
		x_exp[1] = 0.2812; //ee

		x_obs[0] = 0.0853;
		x_obs[4] = 0.3323;
		x_obs[3] = 0.1497;
		x_obs[2] = 0.1161;
		x_obs[1] = 0.2882;

		y[0]	 = 1-0.5;
		y[4]	 = 5-0.5;
		y[3]	 = 4-0.5;
		y[2]	 = 3-0.5;
		y[1]	 = 2-0.5;
		y_errDown[0] = 0.5;
		y_errDown[4] = 0.5;
		y_errDown[3] = 0.5;
		y_errDown[2] = 0.5;
		y_errDown[1] = 0.5;
		y_errUp[0] = 0.5;
		y_errUp[4] = 0.5;
		y_errUp[3] = 0.5;
		y_errUp[2] = 0.5;
		y_errUp[1] = 0.5;

		x_errDown_1sigma[0] = x_exp[0] - 0.0764;
		x_errDown_1sigma[4] = x_exp[4] - 0.1643;
		x_errDown_1sigma[3] = x_exp[3] - 0.1240;
		x_errDown_1sigma[2] = x_exp[2] - 0.1035;
		x_errDown_1sigma[1] = x_exp[1] - 0.2023;
		x_errUp_1sigma[0] = 0.1456 - x_exp[0];
		x_errUp_1sigma[4] = 0.3390 - x_exp[4];
		x_errUp_1sigma[3] = 0.2322 - x_exp[3];
		x_errUp_1sigma[2] = 0.1977 - x_exp[2];
		x_errUp_1sigma[1] = 0.3967 - x_exp[1];

		x_errDown_2sigma[0] = x_exp[0] - 0.0574;
		x_errDown_2sigma[4] = x_exp[4] - 0.1208;
		x_errDown_2sigma[3] = x_exp[3] - 0.0940;
		x_errDown_2sigma[2] = x_exp[2] - 0.0785;
		x_errDown_2sigma[1] = x_exp[1] - 0.1516;
		x_errUp_2sigma[0] = 0.1937 - x_exp[0];
		x_errUp_2sigma[4] = 0.4748 - x_exp[4];
		x_errUp_2sigma[3] = 0.3063 - x_exp[3];
		x_errUp_2sigma[2] = 0.2631 - x_exp[2];
		x_errUp_2sigma[1] = 0.5400 - x_exp[1];
    }
	else {cout<<"Error ! "<<endl; return;}

	//-- CHANGE manually scale (then write on plot E-3)
	for(int i=0; i<n; i++)
	{
		x_exp[i]*= 10;
		x_obs[i]*= 10;
		x_errDown_1sigma[i]*= 10;
		x_errUp_1sigma[i]*= 10;
		x_errDown_2sigma[i]*= 10;
		x_errUp_2sigma[i]*= 10;
	}


    TCanvas* c1 = new TCanvas("", "", 1000, 800);
	c1->SetTopMargin(0.1);
	c1->SetRightMargin(0.07);

	// TH2F* h = new TH2F("", "", 1, 0, 0.65, 5, 0, 5.2);
	TH2F* h = new TH2F("", "", 1, 0, 6.5, 5, 0, 5.2);
	h->SetMaximum(5.1);
	h->Draw("same");

	//Text labels
	const char *labels[n]  = {"Combined", "ee", "#mue", "#mu#mu", "3l"};
	for(int i=1;i<=n;i++) {h->GetYaxis()->SetBinLabel(i,labels[i-1]);}
	h->GetYaxis()->LabelsOption("v");
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetLabelSize(0.04);
	h->GetYaxis()->SetLabelSize(0.05);

	// h->GetXaxis()->SetTitle("95 % CL limit on \\mathscr{B}(t #rightarrow Hq)");

	TGraphAsymmErrors* gr_exp_2sigma = new TGraphAsymmErrors(n, x_exp, y, x_errDown_2sigma, x_errUp_2sigma, y_errDown, y_errUp);
	gr_exp_2sigma->SetFillColor(kYellow);
	gr_exp_2sigma->SetFillStyle(1001);
	gr_exp_2sigma->SetLineWidth(3);
	gr_exp_2sigma->SetLineStyle(2);

	gr_exp_2sigma->Draw("e2 same");

	TGraphAsymmErrors* gr_exp_1sigma = new TGraphAsymmErrors(n, x_exp, y, x_errDown_1sigma, x_errUp_1sigma, y_errDown, y_errUp);
	gr_exp_1sigma->SetFillColor(kGreen);
	gr_exp_1sigma->SetFillStyle(1001);
	gr_exp_1sigma->SetLineStyle(2);
	gr_exp_1sigma->SetLineWidth(3);
	gStyle->SetErrorX(0.0);
	gr_exp_1sigma->Draw("e2 same");

	// TLine l;
	// l.SetLineStyle(2);
	// l.DrawLine(0, 0, 0, 5);

	vector<TLine*> v_line_exp(n);
	for(int i=0; i<n; i++)
	{
		v_line_exp[i] = new TLine(x_exp[i], y[i]-0.5, x_exp[i], y[i]+0.5);
		v_line_exp[i]->SetLineStyle(2);
		v_line_exp[i]->SetLineWidth(3);
		v_line_exp[i]->Draw("same");
	}

	vector<TLine*> v_line_obs(n);
	for(int i=0; i<n; i++)
	{
		v_line_obs[i] = new TLine(x_obs[i], y[i]-0.5, x_obs[i], y[i]+0.5);
		v_line_obs[i]->SetLineWidth(3);
		v_line_obs[i]->Draw("same");
	}

	TLegend* qw = new TLegend(0.70,.50,0.92,0.70);
	qw->AddEntry(gr_exp_1sigma, "Expected #pm 1 #sigma", "LF");
	qw->AddEntry(gr_exp_2sigma, "Expected #pm 2 #sigma", "LF");
	qw->AddEntry(v_line_obs[0], "Observed", "L");
	// qw->AddEntry(v_line_exp[0], "Expected", "L");
	qw->SetBorderSize(0);
	qw->Draw("same");

	TLatex latex;
    latex.SetNDC();
    TString lumi_13TeV = "41.5 fb^{-1} (13 TeV)";
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.93, 0.92,lumi_13TeV);
    // latex.DrawLatex(0.78, 0.92,lumi_13TeV);

	TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextFont(42);
    latex2.SetTextAlign(31);
    latex2.SetTextSize(0.04);
	if(signal == "Hut") latex2.DrawLatex(0.94, 0.02,"95 %\ CL limit on Br(t #rightarrow Hu)");
	else if(signal == "Hct") latex2.DrawLatex(0.94, 0.02,"95 %\ CL limit on Br(t #rightarrow Hc)");

	TLatex power;
	power.SetNDC();
	power.SetTextSize(0.035);
	power.SetTextFont(42);
	power.DrawLatex(0.935,0.13,"#times 10^{-3}");

    c1->SaveAs("limitCompability_"+signal+".png");

    delete c1;
	delete gr_exp_1sigma;
	delete gr_exp_2sigma;
	delete qw;

	for(int i=0; i<n; i++)
	{
		delete v_line_obs[i];
	}
	for(int i=0; i<n; i++)
	{
		delete v_line_exp[i];
	}

    return;
}



//--------------------------------------------
//--------------------------------------------
int main()
{
	TString signal = "Hut";
	// TString signal = "Hct";

    Load_Canvas_Style();

    Draw_Compatibility_SigStrength(signal);
	// Draw_Compatibility_Limits(signal);
}
