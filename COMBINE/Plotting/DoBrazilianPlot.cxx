/**
 * NB : HARD-CODED SCRIPT TO PRODUCE BRAZILIAN PLOTS FOR THQ-FCNC ANALYSIS. ADAPT TO YOUR NEEDS/RESULTS !
*/

#include "TString.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include <iostream>
#include "TString.h"
#include "TMultiGraph.h"
#include "TLine.h"

using namespace std;

void Set_Style()
{
    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");

    // For the axis:

    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(1); // kTRUE
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);

    // Change for log plots:
    gStyle->SetOptLogx(0);
    gStyle->SetOptLogy(0);
    gStyle->SetOptLogz(0);
}

void Do_Brazilian_Plot(TString signal)
{
    const unsigned int n=2;   // number of mass points
    Double_t x[n] = {0.01, 0.1};

    double y_obs[n]    ;
    double y_obs_err[n];
    double y_exp2up[n] ;
    double y_expup[n]  ;
    double y_exp[n]    ;
    double y_expdn[n]  ;;
    double y_exp2dn[n] ;

    //-- Cross sections
    // SF_ST_hut_FCNC = 0.3311;
    // SF_ST_hct_FCNC = 0.0456;
    // SF_TT_FCNC = 1.093;

    if(signal == "hut")
    {
        TString name = "Yields";
        double cross_sect = 1.093+0.3311;
        y_obs[0]	 = 0.38956*cross_sect;
        y_obs[1]	 = 0.38956*cross_sect;
        y_obs_err[0]= 0.0091583*cross_sect;
        y_obs_err[1]= 0.0091583*cross_sect;
        y_exp2up[0] = 0.96284*cross_sect;
        y_exp2up[1] = 0.96284*cross_sect;
        y_expup[0]  = 0.63606*cross_sect;
        y_expup[1]  = 0.63606*cross_sect;
        y_exp[0]	 = 0.41744*cross_sect;
        y_exp[1]	 = 0.41744*cross_sect;
        y_expdn[0]  = 0.2842*cross_sect;
        y_expdn[1]  = 0.2842*cross_sect;
        y_exp2dn[0] = 0.20599*cross_sect;
        y_exp2dn[1] = 0.20599*cross_sect;
    }
    else if(signal == "hct")
    {
        TString name = "Yields";
        double cross_sect = 1.093+0.0456;
        y_obs[0]	 = 0.38672*cross_sect;
        y_obs[1]	 = 0.38672*cross_sect;
        y_obs_err[0]= 0.011028*cross_sect;
        y_obs_err[1]= 0.011028*cross_sect;
        y_exp2up[0] = 0.90769*cross_sect;
        y_exp2up[1] = 0.90769*cross_sect;
        y_expup[0]  = 0.61699*cross_sect;
        y_expup[1]	 = 0.61699*cross_sect;
        y_exp[0]	 = 0.41554*cross_sect;
        y_exp[1]	 = 0.41554*cross_sect;
        y_expdn[0]  = 0.2804*cross_sect;
        y_expdn[1]	 = 0.2804*cross_sect;
        y_exp2dn[0] = 0.20639*cross_sect;
        y_exp2dn[1] = 0.20639*cross_sect;
    }

    Double_t y_sigmaup[n], y_sigmadn[n], y_sigma2up[n], y_sigma2dn[n];

    for(int i = 0; i<n; i++)
    {
        y_sigmadn[i]  = y_exp[i] - y_expdn[i];
        y_sigmaup[i]  = y_expup[i] - y_exp[i];
        y_sigma2dn[i] = y_exp[i] - y_exp2dn[i];
        y_sigma2up[i] = y_exp2up[i] - y_exp[i];
    }

    TCanvas *c1 = new TCanvas("c1","c1", 1000, 800);
    c1->cd();

    TMultiGraph *mg = new TMultiGraph();
    //mg->SetName("");

    TGraphErrors* gr_obs = new TGraphErrors(n, x, y_obs, 0, y_obs_err);
    gr_obs->SetLineWidth(2);
    gr_obs->SetMarkerStyle(20);

    TGraphErrors* gr_exp = new TGraphErrors(n, x, y_exp, 0, 0);
    gr_exp->SetLineWidth(2);
    gr_exp->SetLineStyle(2);
    gr_exp->SetMarkerStyle(22);

    TGraphAsymmErrors* sigma1 = new TGraphAsymmErrors(n, x, y_exp, 0, 0, y_sigmadn, y_sigmaup);
    sigma1->SetFillColor(kGreen);
    sigma1->SetLineColor(kGreen);

    TGraphAsymmErrors* sigma2 = new TGraphAsymmErrors(n, x, y_exp, 0, 0, y_sigma2dn, y_sigma2up);
    sigma2->SetFillColor(kYellow);
    sigma2->SetLineColor(kYellow);

    TGraphErrors* gr_expup = new TGraphErrors(n, x, y_expup, 0, 0);
    gr_expup->SetLineColor(kRed);
    gr_expup->SetLineWidth(2);
    gr_expup->SetLineStyle(2);

    TGraphErrors* gr_expdn = new TGraphErrors(n, x, y_expdn, 0, 0);
    gr_expdn->SetLineColor(kRed);
    gr_expdn->SetLineWidth(2);
    gr_expdn->SetLineStyle(2);

    TGraphErrors* gr_exp2up = new TGraphErrors(n, x, y_exp2up, 0, 0);
    gr_exp2up->SetLineColor(kBlue);
    gr_exp2up->SetLineWidth(2);
    gr_exp2up->SetLineStyle(4);

    TGraphErrors* gr_exp2dn = new TGraphErrors(n, x, y_exp2dn, 0, 0);
    gr_exp2dn->SetLineColor(kBlue);
    gr_exp2dn->SetLineWidth(2);
    gr_exp2dn->SetLineStyle(4);

    mg->Add(sigma2);
    mg->Add(sigma1);
    //mg->Add(gr_exp);
    //mg->Add(gr_expup);
    //mg->Add(gr_expdn);
    //mg->Add(gr_exp2up);
    //mg->Add(gr_exp2dn);
    //mg->Add(gr_obs);
    mg->Draw("3AC");
    mg->GetXaxis()->SetTitleFont(62);
    mg->GetYaxis()->SetTitleFont(62);
    mg->GetXaxis()->SetTitleOffset(.8);
    mg->GetYaxis()->SetTitleOffset(.7);
    mg->GetXaxis()->SetTitleSize(0.055);
    mg->GetYaxis()->SetTitleSize(0.055);
    mg->GetXaxis()->SetLabelSize(0.04);
    mg->GetYaxis()->SetLabelSize(0.04);
    if(signal == "hut") {mg->GetXaxis()->SetTitle("Br(t#rightarrow Hu) %");}
    if(signal == "hct") {mg->GetXaxis()->SetTitle("Br(t#rightarrow Hc) %");}
    mg->GetYaxis()->SetTitle("95% CL limit on #sigma");
    mg->GetYaxis()->SetTitleOffset(1.1);
    mg->GetXaxis()->SetRangeUser(0.01,0.1);
    if(signal == "hct") mg->GetYaxis()->SetRangeUser(0.,0.15);
    if(signal == "hut") mg->GetYaxis()->SetRangeUser(0.,0.08);
    //sigma2  ->Draw("3A");
    //sigma1  ->Draw("same3A");
    gr_exp  ->Draw("samel");
    gr_obs  ->Draw("samel");
    //gr_expup->Draw("samec");
    //gr_expdn->Draw("samec");
    //gr_exp2up->Draw("samec");
    //gr_exp2dn->Draw("samec");

    TLine line (162.5, 1., 200., 1.);
    line.SetLineWidth(1);
    line.SetLineColor(1);
    line.SetLineStyle(2);
    line.Draw("same");

    TF1* thefunc = 0;
    if(signal == "hut") {thefunc = new TF1("funct2", "73.6424*x/100", 0.001,0.2);}
    else if(signal == "hct") {thefunc = new TF1("funct2", "50.034*x/100", 0.001,0.2);}
    thefunc->SetLineColor(2);
    thefunc->SetLineWidth(2);
    thefunc->Draw("same");

    TLegend *leg = new TLegend(0.3,0.65,0.55,0.85);
    leg->AddEntry(gr_obs  ," Observed","lpe");
    leg->AddEntry(gr_exp  ," Expected","lp");
    leg->AddEntry(sigma1  ," Expected #pm1#sigma","f");
    leg->AddEntry(sigma2  ," Expected #pm2#sigma","f");
    leg->AddEntry( thefunc ," theory","l");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw();

    // TLatex* text = new TLatex(0.37, 0.92, "CMS Preliminary, L = 20 fb^{-1}, #sqrt{s} = 8 TeV");
    // //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
    // text->SetNDC(true);
    // text->SetTextFont(62);
    // text->SetTextSize(0.045);
    // text->Draw();

//----------------
// CAPTIONS //
//----------------
// -- using https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines

	float l = c1->GetLeftMargin();
	float t = c1->GetTopMargin();

	TString cmsText = "CMS";
	TLatex latex;
	latex.SetNDC();
	latex.SetTextColor(kBlack);
	latex.SetTextFont(61);
	latex.SetTextAlign(11);
	latex.SetTextSize(0.06);
	// latex.DrawLatex(l + 0.01, 0.92, cmsText);

	TString extraText = "Preliminary";
	latex.SetTextFont(52);
	latex.SetTextSize(0.05);
    // latex.DrawLatex(l + 0.12, 0.92, extraText);

	TString lumi_13TeV = "41.5 fb^{-1} (13 TeV)";
	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04);
    latex.DrawLatex(0.95, 0.92,lumi_13TeV);
    // latex.DrawLatex(0.78, 0.92,lumi_13TeV);

	//------------------
	//-- Channel info
	// TLatex text2 ;
	// text2.SetNDC();
	// text2.SetTextAlign(13);
	// text2.SetTextSize(0.045);
	// text2.SetTextFont(42);

    c1->SaveAs( ("limits_"+signal+".png").Data() );

    delete c1;
    delete mg;
    delete gr_obs;
    delete gr_exp;
    delete gr_expdn;
    delete gr_expup;
    delete gr_exp2dn;
    delete gr_exp2up;
    // delete sigma1;
    // delete sigma2;
    delete thefunc;
    delete leg;
    delete text;

    return;
}

int main()
{
    TString signal = "hut";
    // TString signal = "hct";

    Set_Style();
    Do_Brazilian_Plot(signal);

    return 1;
}
