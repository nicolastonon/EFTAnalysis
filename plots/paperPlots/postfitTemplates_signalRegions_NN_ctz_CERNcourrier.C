void postfitTemplates_signalRegions_NN_ctz_CERNcourrier()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Jul  5 22:23:19 2021) by ROOT version 6.18/04
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,500,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.16);
   c1->SetRightMargin(0.03);
   c1->SetTopMargin(0.07);
   c1->SetBottomMargin(0.13);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: c1_1
   TPad *c1_1 = new TPad("c1_1", "c1_1",1e-11,9.999999e-12,1,1);
   c1_1->Draw();
   c1_1->cd();
   c1_1->Range(0.3907407,-3.205901,0.7611111,4.194764);
   c1_1->SetFillColor(0);
   c1_1->SetBorderMode(0);
   c1_1->SetBorderSize(2);
   c1_1->SetLogy();
   c1_1->SetTickx(1);
   c1_1->SetTicky(1);
   c1_1->SetLeftMargin(0.16);
   c1_1->SetRightMargin(0.03);
   c1_1->SetTopMargin(0.07);
   c1_1->SetBottomMargin(0.4);
   c1_1->SetFrameFillStyle(0);
   c1_1->SetFrameBorderMode(0);
   c1_1->SetFrameFillStyle(0);
   c1_1->SetFrameBorderMode(0);
   
   THStack * = new THStack();
   ->SetName("");
   ->SetTitle("");
   ->SetMinimum(1.5);
   ->SetMaximum(2868);
   
   TH1F *_stack_1 = new TH1F("_stack_1","",8,0.45,0.75);
   _stack_1->SetMinimum(0.5680218);
   _stack_1->SetMaximum(4750.261);
   _stack_1->SetDirectory(0);
   _stack_1->SetStats(0);
   _stack_1->SetLineStyle(0);
   _stack_1->GetXaxis()->SetLabelFont(42);
   _stack_1->GetXaxis()->SetLabelOffset(0.007);
   _stack_1->GetXaxis()->SetLabelSize(0);
   _stack_1->GetXaxis()->SetTitleSize(0.06);
   _stack_1->GetXaxis()->SetTickLength(0);
   _stack_1->GetXaxis()->SetTitleOffset(0.9);
   _stack_1->GetXaxis()->SetTitleFont(42);
   _stack_1->GetYaxis()->SetTitle("Events / bin");
   _stack_1->GetYaxis()->SetNdivisions(506);
   _stack_1->GetYaxis()->SetLabelFont(42);
   _stack_1->GetYaxis()->SetLabelOffset(0.007);
   _stack_1->GetYaxis()->SetLabelSize(0.048);
   _stack_1->GetYaxis()->SetTitleSize(0.06);
   _stack_1->GetYaxis()->SetTickLength(0.04);
   _stack_1->GetYaxis()->SetTitleOffset(1.2);
   _stack_1->GetYaxis()->SetTitleFont(42);
   _stack_1->GetZaxis()->SetLabelFont(42);
   _stack_1->GetZaxis()->SetLabelOffset(0.007);
   _stack_1->GetZaxis()->SetLabelSize(0.05);
   _stack_1->GetZaxis()->SetTitleSize(0.06);
   _stack_1->GetZaxis()->SetTitleOffset(1);
   _stack_1->GetZaxis()->SetTitleFont(42);
   ->SetHistogram(_stack_1);
   
   
   TH1F *_stack_1 = new TH1F("_stack_1","",8,0.45,0.75);
   _stack_1->SetBinContent(1,16.38008);
   _stack_1->SetBinContent(2,0.4185609);
   _stack_1->SetBinContent(3,0.3543069);
   _stack_1->SetBinContent(4,3.232592e-05);
   _stack_1->SetBinContent(5,0.3213319);
   _stack_1->SetBinContent(6,0.1734881);
   _stack_1->SetBinContent(7,0.07504268);
   _stack_1->SetBinContent(8,0.1393902);
   _stack_1->SetEntries(24);
   _stack_1->SetDirectory(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#cccccc");
   _stack_1->SetFillColor(ci);

   ci = TColor::GetColor("#cccccc");
   _stack_1->SetLineColor(ci);
   _stack_1->SetLineStyle(0);
   _stack_1->GetXaxis()->SetLabelFont(42);
   _stack_1->GetXaxis()->SetLabelOffset(0.007);
   _stack_1->GetXaxis()->SetLabelSize(0.05);
   _stack_1->GetXaxis()->SetTitleSize(0.06);
   _stack_1->GetXaxis()->SetTitleOffset(0.9);
   _stack_1->GetXaxis()->SetTitleFont(42);
   _stack_1->GetYaxis()->SetLabelFont(42);
   _stack_1->GetYaxis()->SetLabelOffset(0.007);
   _stack_1->GetYaxis()->SetLabelSize(0.05);
   _stack_1->GetYaxis()->SetTitleSize(0.06);
   _stack_1->GetYaxis()->SetTitleOffset(1.25);
   _stack_1->GetYaxis()->SetTitleFont(42);
   _stack_1->GetZaxis()->SetLabelFont(42);
   _stack_1->GetZaxis()->SetLabelOffset(0.007);
   _stack_1->GetZaxis()->SetLabelSize(0.05);
   _stack_1->GetZaxis()->SetTitleSize(0.06);
   _stack_1->GetZaxis()->SetTitleOffset(1);
   _stack_1->GetZaxis()->SetTitleFont(42);
   ->Add(_stack_1,"");
   
   TH1F *_stack_2 = new TH1F("_stack_2","",8,0.45,0.75);
   _stack_2->SetBinContent(1,3.67764);
   _stack_2->SetBinContent(2,0.1007296);
   _stack_2->SetBinContent(3,0.007791365);
   _stack_2->SetBinContent(4,0.03771733);
   _stack_2->SetBinContent(5,0.01258882);
   _stack_2->SetBinContent(6,2.984036e-05);
   _stack_2->SetBinContent(7,2.984036e-05);
   _stack_2->SetBinContent(8,2.984036e-05);
   _stack_2->SetEntries(24);
   _stack_2->SetDirectory(0);

   ci = TColor::GetColor("#cc6699");
   _stack_2->SetFillColor(ci);

   ci = TColor::GetColor("#cc6699");
   _stack_2->SetLineColor(ci);
   _stack_2->SetLineStyle(0);
   _stack_2->GetXaxis()->SetLabelFont(42);
   _stack_2->GetXaxis()->SetLabelOffset(0.007);
   _stack_2->GetXaxis()->SetLabelSize(0.05);
   _stack_2->GetXaxis()->SetTitleSize(0.06);
   _stack_2->GetXaxis()->SetTitleOffset(0.9);
   _stack_2->GetXaxis()->SetTitleFont(42);
   _stack_2->GetYaxis()->SetLabelFont(42);
   _stack_2->GetYaxis()->SetLabelOffset(0.007);
   _stack_2->GetYaxis()->SetLabelSize(0.05);
   _stack_2->GetYaxis()->SetTitleSize(0.06);
   _stack_2->GetYaxis()->SetTitleOffset(1.25);
   _stack_2->GetYaxis()->SetTitleFont(42);
   _stack_2->GetZaxis()->SetLabelFont(42);
   _stack_2->GetZaxis()->SetLabelOffset(0.007);
   _stack_2->GetZaxis()->SetLabelSize(0.05);
   _stack_2->GetZaxis()->SetTitleSize(0.06);
   _stack_2->GetZaxis()->SetTitleOffset(1);
   _stack_2->GetZaxis()->SetTitleFont(42);
   ->Add(_stack_2,"");
   
   TH1F *_stack_3 = new TH1F("_stack_3","",8,0.45,0.75);
   _stack_3->SetBinContent(1,11.16937);
   _stack_3->SetBinContent(2,2.20121);
   _stack_3->SetBinContent(3,0.6491132);
   _stack_3->SetBinContent(4,0.8741133);
   _stack_3->SetBinContent(5,0.8021057);
   _stack_3->SetBinContent(6,0.3217089);
   _stack_3->SetBinContent(7,0.1033402);
   _stack_3->SetBinContent(8,0.142947);
   _stack_3->SetEntries(24);
   _stack_3->SetDirectory(0);

   ci = TColor::GetColor("#35b863");
   _stack_3->SetFillColor(ci);

   ci = TColor::GetColor("#35b863");
   _stack_3->SetLineColor(ci);
   _stack_3->SetLineStyle(0);
   _stack_3->GetXaxis()->SetLabelFont(42);
   _stack_3->GetXaxis()->SetLabelOffset(0.007);
   _stack_3->GetXaxis()->SetLabelSize(0.05);
   _stack_3->GetXaxis()->SetTitleSize(0.06);
   _stack_3->GetXaxis()->SetTitleOffset(0.9);
   _stack_3->GetXaxis()->SetTitleFont(42);
   _stack_3->GetYaxis()->SetLabelFont(42);
   _stack_3->GetYaxis()->SetLabelOffset(0.007);
   _stack_3->GetYaxis()->SetLabelSize(0.05);
   _stack_3->GetYaxis()->SetTitleSize(0.06);
   _stack_3->GetYaxis()->SetTitleOffset(1.25);
   _stack_3->GetYaxis()->SetTitleFont(42);
   _stack_3->GetZaxis()->SetLabelFont(42);
   _stack_3->GetZaxis()->SetLabelOffset(0.007);
   _stack_3->GetZaxis()->SetLabelSize(0.05);
   _stack_3->GetZaxis()->SetTitleSize(0.06);
   _stack_3->GetZaxis()->SetTitleOffset(1);
   _stack_3->GetZaxis()->SetTitleFont(42);
   ->Add(_stack_3,"");
   
   TH1F *_stack_4 = new TH1F("_stack_4","",8,0.45,0.75);
   _stack_4->SetBinContent(1,36.31194);
   _stack_4->SetBinContent(2,5.365836);
   _stack_4->SetBinContent(3,2.686675);
   _stack_4->SetBinContent(4,2.689535);
   _stack_4->SetBinContent(5,2.655447);
   _stack_4->SetBinContent(6,1.260464);
   _stack_4->SetBinContent(7,0.612062);
   _stack_4->SetBinContent(8,0.484611);
   _stack_4->SetEntries(24);
   _stack_4->SetDirectory(0);

   ci = TColor::GetColor("#b2df8a");
   _stack_4->SetFillColor(ci);

   ci = TColor::GetColor("#b2df8a");
   _stack_4->SetLineColor(ci);
   _stack_4->SetLineStyle(0);
   _stack_4->GetXaxis()->SetLabelFont(42);
   _stack_4->GetXaxis()->SetLabelOffset(0.007);
   _stack_4->GetXaxis()->SetLabelSize(0.05);
   _stack_4->GetXaxis()->SetTitleSize(0.06);
   _stack_4->GetXaxis()->SetTitleOffset(0.9);
   _stack_4->GetXaxis()->SetTitleFont(42);
   _stack_4->GetYaxis()->SetLabelFont(42);
   _stack_4->GetYaxis()->SetLabelOffset(0.007);
   _stack_4->GetYaxis()->SetLabelSize(0.05);
   _stack_4->GetYaxis()->SetTitleSize(0.06);
   _stack_4->GetYaxis()->SetTitleOffset(1.25);
   _stack_4->GetYaxis()->SetTitleFont(42);
   _stack_4->GetZaxis()->SetLabelFont(42);
   _stack_4->GetZaxis()->SetLabelOffset(0.007);
   _stack_4->GetZaxis()->SetLabelSize(0.05);
   _stack_4->GetZaxis()->SetTitleSize(0.06);
   _stack_4->GetZaxis()->SetTitleOffset(1);
   _stack_4->GetZaxis()->SetTitleFont(42);
   ->Add(_stack_4,"");
   
   TH1F *_stack_5 = new TH1F("_stack_5","",8,0.45,0.75);
   _stack_5->SetBinContent(1,33.51057);
   _stack_5->SetBinContent(2,2.079117);
   _stack_5->SetBinContent(3,0.919536);
   _stack_5->SetBinContent(4,0.7591751);
   _stack_5->SetBinContent(5,0.4841586);
   _stack_5->SetBinContent(6,0.2385893);
   _stack_5->SetBinContent(7,0.1861717);
   _stack_5->SetBinContent(8,0.1284063);
   _stack_5->SetEntries(24);
   _stack_5->SetDirectory(0);

   ci = TColor::GetColor("#a6cee3");
   _stack_5->SetFillColor(ci);

   ci = TColor::GetColor("#a6cee3");
   _stack_5->SetLineColor(ci);
   _stack_5->SetLineStyle(0);
   _stack_5->GetXaxis()->SetLabelFont(42);
   _stack_5->GetXaxis()->SetLabelOffset(0.007);
   _stack_5->GetXaxis()->SetLabelSize(0.05);
   _stack_5->GetXaxis()->SetTitleSize(0.06);
   _stack_5->GetXaxis()->SetTitleOffset(0.9);
   _stack_5->GetXaxis()->SetTitleFont(42);
   _stack_5->GetYaxis()->SetLabelFont(42);
   _stack_5->GetYaxis()->SetLabelOffset(0.007);
   _stack_5->GetYaxis()->SetLabelSize(0.05);
   _stack_5->GetYaxis()->SetTitleSize(0.06);
   _stack_5->GetYaxis()->SetTitleOffset(1.25);
   _stack_5->GetYaxis()->SetTitleFont(42);
   _stack_5->GetZaxis()->SetLabelFont(42);
   _stack_5->GetZaxis()->SetLabelOffset(0.007);
   _stack_5->GetZaxis()->SetLabelSize(0.05);
   _stack_5->GetZaxis()->SetTitleSize(0.06);
   _stack_5->GetZaxis()->SetTitleOffset(1);
   _stack_5->GetZaxis()->SetTitleFont(42);
   ->Add(_stack_5,"");
   
   TH1F *_stack_6 = new TH1F("_stack_6","",8,0.45,0.75);
   _stack_6->SetBinContent(1,298.248);
   _stack_6->SetBinContent(2,48.51447);
   _stack_6->SetBinContent(3,21.11051);
   _stack_6->SetBinContent(4,22.98304);
   _stack_6->SetBinContent(5,20.31575);
   _stack_6->SetBinContent(6,6.561287);
   _stack_6->SetBinContent(7,4.050829);
   _stack_6->SetBinContent(8,2.232787);
   _stack_6->SetEntries(24);
   _stack_6->SetDirectory(0);

   ci = TColor::GetColor("#1f78b4");
   _stack_6->SetFillColor(ci);

   ci = TColor::GetColor("#1f78b4");
   _stack_6->SetLineColor(ci);
   _stack_6->SetLineStyle(0);
   _stack_6->GetXaxis()->SetLabelFont(42);
   _stack_6->GetXaxis()->SetLabelOffset(0.007);
   _stack_6->GetXaxis()->SetLabelSize(0.05);
   _stack_6->GetXaxis()->SetTitleSize(0.06);
   _stack_6->GetXaxis()->SetTitleOffset(0.9);
   _stack_6->GetXaxis()->SetTitleFont(42);
   _stack_6->GetYaxis()->SetLabelFont(42);
   _stack_6->GetYaxis()->SetLabelOffset(0.007);
   _stack_6->GetYaxis()->SetLabelSize(0.05);
   _stack_6->GetYaxis()->SetTitleSize(0.06);
   _stack_6->GetYaxis()->SetTitleOffset(1.25);
   _stack_6->GetYaxis()->SetTitleFont(42);
   _stack_6->GetZaxis()->SetLabelFont(42);
   _stack_6->GetZaxis()->SetLabelOffset(0.007);
   _stack_6->GetZaxis()->SetLabelSize(0.05);
   _stack_6->GetZaxis()->SetTitleSize(0.06);
   _stack_6->GetZaxis()->SetTitleOffset(1);
   _stack_6->GetZaxis()->SetTitleFont(42);
   ->Add(_stack_6,"");
   
   TH1F *_stack_7 = new TH1F("_stack_7","",8,0.45,0.75);
   _stack_7->SetBinContent(1,24.17974);
   _stack_7->SetBinContent(2,3.802644);
   _stack_7->SetBinContent(3,1.397188);
   _stack_7->SetBinContent(4,1.731064);
   _stack_7->SetBinContent(5,1.449534);
   _stack_7->SetBinContent(6,0.5996611);
   _stack_7->SetBinContent(7,0.385654);
   _stack_7->SetBinContent(8,0.1679598);
   _stack_7->SetEntries(24);
   _stack_7->SetDirectory(0);

   ci = TColor::GetColor("#fce181");
   _stack_7->SetFillColor(ci);

   ci = TColor::GetColor("#fce181");
   _stack_7->SetLineColor(ci);
   _stack_7->SetLineStyle(0);
   _stack_7->GetXaxis()->SetLabelFont(42);
   _stack_7->GetXaxis()->SetLabelOffset(0.007);
   _stack_7->GetXaxis()->SetLabelSize(0.05);
   _stack_7->GetXaxis()->SetTitleSize(0.06);
   _stack_7->GetXaxis()->SetTitleOffset(0.9);
   _stack_7->GetXaxis()->SetTitleFont(42);
   _stack_7->GetYaxis()->SetLabelFont(42);
   _stack_7->GetYaxis()->SetLabelOffset(0.007);
   _stack_7->GetYaxis()->SetLabelSize(0.05);
   _stack_7->GetYaxis()->SetTitleSize(0.06);
   _stack_7->GetYaxis()->SetTitleOffset(1.25);
   _stack_7->GetYaxis()->SetTitleFont(42);
   _stack_7->GetZaxis()->SetLabelFont(42);
   _stack_7->GetZaxis()->SetLabelOffset(0.007);
   _stack_7->GetZaxis()->SetLabelSize(0.05);
   _stack_7->GetZaxis()->SetTitleSize(0.06);
   _stack_7->GetZaxis()->SetTitleOffset(1);
   _stack_7->GetZaxis()->SetTitleFont(42);
   ->Add(_stack_7,"");
   
   TH1F *_stack_8 = new TH1F("_stack_8","",8,0.45,0.75);
   _stack_8->SetBinContent(1,32.28824);
   _stack_8->SetBinContent(2,4.08455);
   _stack_8->SetBinContent(3,1.643061);
   _stack_8->SetBinContent(4,1.773579);
   _stack_8->SetBinContent(5,1.771113);
   _stack_8->SetBinContent(6,0.6682047);
   _stack_8->SetBinContent(7,0.485924);
   _stack_8->SetBinContent(8,0.1825291);
   _stack_8->SetEntries(24);
   _stack_8->SetDirectory(0);

   ci = TColor::GetColor("#e3413d");
   _stack_8->SetFillColor(ci);

   ci = TColor::GetColor("#e3413d");
   _stack_8->SetLineColor(ci);
   _stack_8->SetLineStyle(0);
   _stack_8->GetXaxis()->SetLabelFont(42);
   _stack_8->GetXaxis()->SetLabelOffset(0.007);
   _stack_8->GetXaxis()->SetLabelSize(0.05);
   _stack_8->GetXaxis()->SetTitleSize(0.06);
   _stack_8->GetXaxis()->SetTitleOffset(0.9);
   _stack_8->GetXaxis()->SetTitleFont(42);
   _stack_8->GetYaxis()->SetLabelFont(42);
   _stack_8->GetYaxis()->SetLabelOffset(0.007);
   _stack_8->GetYaxis()->SetLabelSize(0.05);
   _stack_8->GetYaxis()->SetTitleSize(0.06);
   _stack_8->GetYaxis()->SetTitleOffset(1.25);
   _stack_8->GetYaxis()->SetTitleFont(42);
   _stack_8->GetZaxis()->SetLabelFont(42);
   _stack_8->GetZaxis()->SetLabelOffset(0.007);
   _stack_8->GetZaxis()->SetLabelSize(0.05);
   _stack_8->GetZaxis()->SetTitleSize(0.06);
   _stack_8->GetZaxis()->SetTitleOffset(1);
   _stack_8->GetZaxis()->SetTitleFont(42);
   ->Add(_stack_8,"");
   ->Draw("hist");
   
   TH1F *__1 = new TH1F("__1","",8,0.45,0.75);
   __1->SetBinContent(1,478);
   __1->SetBinContent(2,65);
   __1->SetBinContent(3,27);
   __1->SetBinContent(4,30);
   __1->SetBinContent(5,20);
   __1->SetBinContent(6,10);
   __1->SetBinContent(7,6);
   __1->SetBinContent(8,3);
   __1->SetBinError(1,21.86321);
   __1->SetBinError(2,8.062258);
   __1->SetBinError(3,5.196152);
   __1->SetBinError(4,5.477226);
   __1->SetBinError(5,4.472136);
   __1->SetBinError(6,3.162278);
   __1->SetBinError(7,2.44949);
   __1->SetBinError(8,1.732051);
   __1->SetEntries(23);
   __1->SetDirectory(0);
   __1->SetLineStyle(0);
   __1->SetMarkerStyle(20);
   __1->GetXaxis()->SetLabelFont(42);
   __1->GetXaxis()->SetLabelOffset(0.007);
   __1->GetXaxis()->SetLabelSize(0.05);
   __1->GetXaxis()->SetTitleSize(0.06);
   __1->GetXaxis()->SetTitleOffset(0.9);
   __1->GetXaxis()->SetTitleFont(42);
   __1->GetYaxis()->SetLabelFont(42);
   __1->GetYaxis()->SetLabelOffset(0.007);
   __1->GetYaxis()->SetLabelSize(0.05);
   __1->GetYaxis()->SetTitleSize(0.06);
   __1->GetYaxis()->SetTitleOffset(1.25);
   __1->GetYaxis()->SetTitleFont(42);
   __1->GetZaxis()->SetLabelFont(42);
   __1->GetZaxis()->SetLabelOffset(0.007);
   __1->GetZaxis()->SetLabelSize(0.05);
   __1->GetZaxis()->SetTitleSize(0.06);
   __1->GetZaxis()->SetTitleOffset(1);
   __1->GetZaxis()->SetTitleFont(42);
   __1->Draw("e0p same");
   
   Double_t Graph0_fx3001[8] = {
   0.46875,
   0.50625,
   0.54375,
   0.58125,
   0.61875,
   0.65625,
   0.69375,
   0.73125};
   Double_t Graph0_fy3001[8] = {
   455.7656,
   66.56712,
   28.76818,
   30.84826,
   27.81202,
   9.823434,
   5.899054,
   3.47866};
   Double_t Graph0_felx3001[8] = {
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875};
   Double_t Graph0_fely3001[8] = {
   18.39299,
   3.124905,
   1.572905,
   1.622252,
   1.586518,
   0.5954024,
   0.4066325,
   0.3290225};
   Double_t Graph0_fehx3001[8] = {
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875};
   Double_t Graph0_fehy3001[8] = {
   18.39299,
   3.124905,
   1.572905,
   1.622252,
   1.586518,
   0.5954024,
   0.4066325,
   0.3290225};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(8,Graph0_fx3001,Graph0_fy3001,Graph0_felx3001,Graph0_fehx3001,Graph0_fely3001,Graph0_fehy3001);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetFillStyle(3254);
   
   TH1F *Graph_Graph03001 = new TH1F("Graph_Graph03001","Graph",100,0.42,0.78);
   Graph_Graph03001->SetMinimum(2.834674);
   Graph_Graph03001->SetMaximum(521.2594);
   Graph_Graph03001->SetDirectory(0);
   Graph_Graph03001->SetStats(0);
   Graph_Graph03001->SetLineStyle(0);
   Graph_Graph03001->GetXaxis()->SetLabelFont(42);
   Graph_Graph03001->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph03001->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph03001->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph03001->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph03001->GetXaxis()->SetTitleFont(42);
   Graph_Graph03001->GetYaxis()->SetLabelFont(42);
   Graph_Graph03001->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph03001->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph03001->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph03001->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph03001->GetYaxis()->SetTitleFont(42);
   Graph_Graph03001->GetZaxis()->SetLabelFont(42);
   Graph_Graph03001->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph03001->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph03001->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph03001->GetZaxis()->SetTitleOffset(1);
   Graph_Graph03001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph03001);
   
   grae->Draw("e2 ");
  
// ------------>Primitives in pad: pad_ratio
   TPad *pad_ratio = new TPad("pad_ratio", "pad_ratio",0,0,1,1);
   pad_ratio->Draw();
   pad_ratio->cd();
   pad_ratio->Range(0.3907407,-2.815385,0.7611111,8.723077);
   pad_ratio->SetFillColor(0);
   pad_ratio->SetFillStyle(4000);
   pad_ratio->SetBorderMode(0);
   pad_ratio->SetBorderSize(2);
   pad_ratio->SetGridy();
   pad_ratio->SetTickx(1);
   pad_ratio->SetTicky(1);
   pad_ratio->SetLeftMargin(0.16);
   pad_ratio->SetRightMargin(0.03);
   pad_ratio->SetTopMargin(0.6);
   pad_ratio->SetBottomMargin(0.27);
   pad_ratio->SetFrameFillStyle(0);
   pad_ratio->SetFrameBorderMode(0);
   pad_ratio->SetFrameFillStyle(0);
   pad_ratio->SetFrameBorderMode(0);
   
   TH1F *__2 = new TH1F("__2","",8,0.45,0.75);
   __2->SetBinContent(1,1.048785);
   __2->SetBinContent(2,0.9764582);
   __2->SetBinContent(3,0.9385368);
   __2->SetBinContent(4,0.9725024);
   __2->SetBinContent(5,0.7191135);
   __2->SetBinContent(6,1.017974);
   __2->SetBinContent(7,1.017112);
   __2->SetBinContent(8,0.862401);
   __2->SetBinError(1,0.0479703);
   __2->SetBinError(2,0.1211147);
   __2->SetBinError(3,0.1806215);
   __2->SetBinError(4,0.1775538);
   __2->SetBinError(5,0.1607987);
   __2->SetBinError(6,0.3219117);
   __2->SetBinError(7,0.4152344);
   __2->SetBinError(8,0.4979075);
   __2->SetMinimum(0.3);
   __2->SetMaximum(1.8);
   __2->SetEntries(90.41521);
   __2->SetLineStyle(0);
   __2->SetMarkerStyle(20);
   __2->SetMarkerSize(1.2);
   __2->GetXaxis()->SetNdivisions(-8);
   __2->GetXaxis()->SetLabelFont(42);
   __2->GetXaxis()->SetLabelOffset(0.007);
   __2->GetXaxis()->SetLabelSize(0);
   __2->GetXaxis()->SetTitleSize(0);
   __2->GetXaxis()->SetTickLength(0.01);
   __2->GetXaxis()->SetTitleOffset(0.9);
   __2->GetXaxis()->SetTitleFont(42);
   __2->GetYaxis()->SetTitle("#frac{Data}{Pred.}");
   __2->GetYaxis()->CenterTitle(true);
   __2->GetYaxis()->SetNdivisions(303);
   __2->GetYaxis()->SetLabelFont(42);
   __2->GetYaxis()->SetLabelOffset(0.007);
   __2->GetYaxis()->SetTickLength(0.15);
   __2->GetYaxis()->SetTitleOffset(1.9);
   __2->GetYaxis()->SetTitleFont(42);
   __2->GetZaxis()->SetLabelFont(42);
   __2->GetZaxis()->SetLabelOffset(0.007);
   __2->GetZaxis()->SetLabelSize(0.05);
   __2->GetZaxis()->SetTitleSize(0.06);
   __2->GetZaxis()->SetTitleOffset(1);
   __2->GetZaxis()->SetTitleFont(42);
   __2->Draw("E1 X0 P");
   
   TH1F *__3 = new TH1F("__3","",8,0.45,0.75);
   __3->SetBinContent(1,-999);
   __3->SetBinContent(2,-999);
   __3->SetBinContent(3,-999);
   __3->SetBinContent(4,-999);
   __3->SetBinContent(5,-999);
   __3->SetBinContent(6,-999);
   __3->SetBinContent(7,-999);
   __3->SetBinContent(8,-999);
   __3->SetBinError(1,0.0479703);
   __3->SetBinError(2,0.1211147);
   __3->SetBinError(3,0.1806215);
   __3->SetBinError(4,0.1775538);
   __3->SetBinError(5,0.1607987);
   __3->SetBinError(6,0.3219117);
   __3->SetBinError(7,0.4152344);
   __3->SetBinError(8,0.4979075);
   __3->SetEntries(98.41521);
   __3->SetLineStyle(0);
   __3->SetMarkerStyle(26);
   __3->SetMarkerSize(1.5);
   __3->GetXaxis()->SetNdivisions(-8);
   __3->GetXaxis()->SetLabelFont(42);
   __3->GetXaxis()->SetLabelOffset(0.007);
   __3->GetXaxis()->SetLabelSize(0);
   __3->GetXaxis()->SetTitleSize(0);
   __3->GetXaxis()->SetTickLength(0.01);
   __3->GetXaxis()->SetTitleOffset(0.9);
   __3->GetXaxis()->SetTitleFont(42);
   __3->GetYaxis()->SetTitle("#frac{Data}{Pred.}");
   __3->GetYaxis()->CenterTitle(true);
   __3->GetYaxis()->SetNdivisions(303);
   __3->GetYaxis()->SetLabelFont(42);
   __3->GetYaxis()->SetLabelOffset(0.007);
   __3->GetYaxis()->SetTickLength(0.15);
   __3->GetYaxis()->SetTitleOffset(1.9);
   __3->GetYaxis()->SetTitleFont(42);
   __3->GetZaxis()->SetLabelFont(42);
   __3->GetZaxis()->SetLabelOffset(0.007);
   __3->GetZaxis()->SetLabelSize(0.05);
   __3->GetZaxis()->SetTitleSize(0.06);
   __3->GetZaxis()->SetTitleOffset(1);
   __3->GetZaxis()->SetTitleFont(42);
   __3->Draw("E1 X0 P same");
   
   TH1F *__4 = new TH1F("__4","",8,0.45,0.75);
   __4->SetBinContent(1,-999);
   __4->SetBinContent(2,-999);
   __4->SetBinContent(3,-999);
   __4->SetBinContent(4,-999);
   __4->SetBinContent(5,-999);
   __4->SetBinContent(6,-999);
   __4->SetBinContent(7,-999);
   __4->SetBinContent(8,-999);
   __4->SetBinError(1,0.0479703);
   __4->SetBinError(2,0.1211147);
   __4->SetBinError(3,0.1806215);
   __4->SetBinError(4,0.1775538);
   __4->SetBinError(5,0.1607987);
   __4->SetBinError(6,0.3219117);
   __4->SetBinError(7,0.4152344);
   __4->SetBinError(8,0.4979075);
   __4->SetEntries(98.41521);
   __4->SetLineStyle(0);
   __4->SetMarkerStyle(32);
   __4->SetMarkerSize(1.5);
   __4->GetXaxis()->SetNdivisions(-8);
   __4->GetXaxis()->SetLabelFont(42);
   __4->GetXaxis()->SetLabelOffset(0.007);
   __4->GetXaxis()->SetLabelSize(0);
   __4->GetXaxis()->SetTitleSize(0);
   __4->GetXaxis()->SetTickLength(0.01);
   __4->GetXaxis()->SetTitleOffset(0.9);
   __4->GetXaxis()->SetTitleFont(42);
   __4->GetYaxis()->SetTitle("#frac{Data}{Pred.}");
   __4->GetYaxis()->CenterTitle(true);
   __4->GetYaxis()->SetNdivisions(303);
   __4->GetYaxis()->SetLabelFont(42);
   __4->GetYaxis()->SetLabelOffset(0.007);
   __4->GetYaxis()->SetTickLength(0.15);
   __4->GetYaxis()->SetTitleOffset(1.9);
   __4->GetYaxis()->SetTitleFont(42);
   __4->GetZaxis()->SetLabelFont(42);
   __4->GetZaxis()->SetLabelOffset(0.007);
   __4->GetZaxis()->SetLabelSize(0.05);
   __4->GetZaxis()->SetTitleSize(0.06);
   __4->GetZaxis()->SetTitleOffset(1);
   __4->GetZaxis()->SetTitleFont(42);
   __4->Draw("E1 X0 P same");
   
   Double_t Graph0_fx3002[8] = {
   0.46875,
   0.50625,
   0.54375,
   0.58125,
   0.61875,
   0.65625,
   0.69375,
   0.73125};
   Double_t Graph0_fy3002[8] = {
   1,
   1,
   1,
   1,
   1,
   1,
   1,
   1};
   Double_t Graph0_felx3002[8] = {
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875};
   Double_t Graph0_fely3002[8] = {
   0.04035625,
   0.04694367,
   0.05467514,
   0.05258812,
   0.05704433,
   0.06061042,
   0.06893182,
   0.09458312};
   Double_t Graph0_fehx3002[8] = {
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875,
   0.01875};
   Double_t Graph0_fehy3002[8] = {
   0.04035625,
   0.04694367,
   0.05467514,
   0.05258812,
   0.05704433,
   0.06061042,
   0.06893182,
   0.09458312};
   grae = new TGraphAsymmErrors(8,Graph0_fx3002,Graph0_fy3002,Graph0_felx3002,Graph0_fehx3002,Graph0_fely3002,Graph0_fehy3002);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetFillStyle(3254);
   
   TH1F *Graph_Graph03002 = new TH1F("Graph_Graph03002","Graph",100,0.42,0.78);
   Graph_Graph03002->SetMinimum(0.8865003);
   Graph_Graph03002->SetMaximum(1.1135);
   Graph_Graph03002->SetDirectory(0);
   Graph_Graph03002->SetStats(0);
   Graph_Graph03002->SetLineStyle(0);
   Graph_Graph03002->GetXaxis()->SetLabelFont(42);
   Graph_Graph03002->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph03002->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph03002->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph03002->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph03002->GetXaxis()->SetTitleFont(42);
   Graph_Graph03002->GetYaxis()->SetLabelFont(42);
   Graph_Graph03002->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph03002->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph03002->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph03002->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph03002->GetYaxis()->SetTitleFont(42);
   Graph_Graph03002->GetZaxis()->SetLabelFont(42);
   Graph_Graph03002->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph03002->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph03002->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph03002->GetZaxis()->SetTitleOffset(1);
   Graph_Graph03002->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph03002);
   
   grae->Draw("e2 ");
  
// ------------>Primitives in pad: pad_ratio
   TPad *pad_ratio = new TPad("pad_ratio", "pad_ratio",0,0,1,1);
   pad_ratio->Draw();
   pad_ratio->cd();
   pad_ratio->Range(0.3907407,-6.37059,0.7611111,65.33531);
   pad_ratio->SetFillColor(0);
   pad_ratio->SetFillStyle(4000);
   pad_ratio->SetBorderMode(0);
   pad_ratio->SetBorderSize(2);
   pad_ratio->SetTickx(1);
   pad_ratio->SetTicky(1);
   pad_ratio->SetLeftMargin(0.16);
   pad_ratio->SetRightMargin(0.03);
   pad_ratio->SetTopMargin(0.73);
   pad_ratio->SetFrameFillStyle(0);
   pad_ratio->SetFrameBorderMode(0);
   pad_ratio->SetFrameFillStyle(0);
   pad_ratio->SetFrameBorderMode(0);
   
   TH1F *__5 = new TH1F("__5","",8,0.45,0.75);
   __5->SetBinContent(1,1.141862);
   __5->SetBinContent(2,1.290886);
   __5->SetBinContent(3,1.351332);
   __5->SetBinContent(4,1.504208);
   __5->SetBinContent(5,1.663718);
   __5->SetBinContent(6,2.000915);
   __5->SetBinContent(7,2.359776);
   __5->SetBinContent(8,3.494187);
   __5->SetBinError(1,0.01201137);
   __5->SetBinError(2,0.02496382);
   __5->SetBinError(3,0.04298839);
   __5->SetBinError(4,0.04619402);
   __5->SetBinError(5,0.06705618);
   __5->SetBinError(6,0.09360625);
   __5->SetBinError(7,0.1474305);
   __5->SetBinError(8,0.2929997);
   __5->SetMinimum(0.8);
   __5->SetMaximum(12.99);
   __5->SetEntries(1745.675);

   ci = TColor::GetColor("#990099");
   __5->SetLineColor(ci);
   __5->SetLineStyle(0);
   __5->SetLineWidth(3);
   __5->SetMarkerStyle(20);
   __5->SetMarkerSize(1.2);
   __5->GetXaxis()->SetTitle("NN-C_{tZ}-t#bar{t}Z output");
   __5->GetXaxis()->SetBinLabel(1,"1");
   __5->GetXaxis()->SetBinLabel(2,"2");
   __5->GetXaxis()->SetBinLabel(3,"3");
   __5->GetXaxis()->SetBinLabel(4,"4");
   __5->GetXaxis()->SetBinLabel(5,"5");
   __5->GetXaxis()->SetBinLabel(6,"6");
   __5->GetXaxis()->SetBinLabel(7,"7");
   __5->GetXaxis()->SetBinLabel(8,"8");
   __5->GetXaxis()->SetNdivisions(-8);
   __5->GetXaxis()->SetLabelFont(42);
   __5->GetXaxis()->SetLabelOffset(0.02);
   __5->GetXaxis()->SetLabelSize(0.07);
   __5->GetXaxis()->SetTitleSize(0.045);
   __5->GetXaxis()->SetTickLength(0.01);
   __5->GetXaxis()->SetTitleOffset(0.95);
   __5->GetXaxis()->SetTitleFont(42);
   __5->GetYaxis()->SetTitle("#frac{SM+EFT}{SM}");
   __5->GetYaxis()->CenterTitle(true);
   __5->GetYaxis()->SetNdivisions(505);
   __5->GetYaxis()->SetLabelFont(42);
   __5->GetYaxis()->SetLabelOffset(0.007);
   __5->GetYaxis()->SetTickLength(0.1);
   __5->GetYaxis()->SetTitleOffset(1.9);
   __5->GetYaxis()->SetTitleFont(42);
   __5->GetZaxis()->SetLabelFont(42);
   __5->GetZaxis()->SetLabelOffset(0.007);
   __5->GetZaxis()->SetLabelSize(0.05);
   __5->GetZaxis()->SetTitleSize(0.06);
   __5->GetZaxis()->SetTitleOffset(1);
   __5->GetZaxis()->SetTitleFont(42);
   __5->Draw("hist");
   
   TH1F *__6 = new TH1F("__6","",8,0.45,0.75);
   __6->SetBinContent(1,1.57717);
   __6->SetBinContent(2,2.156509);
   __6->SetBinContent(3,2.410824);
   __6->SetBinContent(4,3.022161);
   __6->SetBinContent(5,3.635533);
   __6->SetBinContent(6,5.021362);
   __6->SetBinContent(7,6.403821);
   __6->SetBinContent(8,10.86679);
   __6->SetBinError(1,0.0165602);
   __6->SetBinError(2,0.0404947);
   __6->SetBinError(3,0.07363709);
   __6->SetBinError(4,0.09454215);
   __6->SetBinError(5,0.1372023);
   __6->SetBinError(6,0.2350384);
   __6->SetBinError(7,0.4105678);
   __6->SetBinError(8,0.9957046);
   __6->SetEntries(985.0163);

   ci = TColor::GetColor("#ff99cc");
   __6->SetLineColor(ci);
   __6->SetLineStyle(0);
   __6->SetLineWidth(3);
   __6->GetXaxis()->SetLabelFont(42);
   __6->GetXaxis()->SetLabelOffset(0.007);
   __6->GetXaxis()->SetLabelSize(0.05);
   __6->GetXaxis()->SetTitleSize(0.06);
   __6->GetXaxis()->SetTitleOffset(0.9);
   __6->GetXaxis()->SetTitleFont(42);
   __6->GetYaxis()->SetLabelFont(42);
   __6->GetYaxis()->SetLabelOffset(0.007);
   __6->GetYaxis()->SetLabelSize(0.05);
   __6->GetYaxis()->SetTitleSize(0.06);
   __6->GetYaxis()->SetTitleOffset(1.25);
   __6->GetYaxis()->SetTitleFont(42);
   __6->GetZaxis()->SetLabelFont(42);
   __6->GetZaxis()->SetLabelOffset(0.007);
   __6->GetZaxis()->SetLabelSize(0.05);
   __6->GetZaxis()->SetTitleSize(0.06);
   __6->GetZaxis()->SetTitleOffset(1);
   __6->GetZaxis()->SetTitleFont(42);
   __6->Draw("hist same");
   
   TLegend *leg = new TLegend(0.2,0.17,0.82,0.26,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","C_{tZ} /#Lambda^{2} [TeV^{-2}]","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","1.5","L");

   ci = TColor::GetColor("#990099");
   entry->SetLineColor(ci);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","3","L");

   ci = TColor::GetColor("#ff99cc");
   entry->SetLineColor(ci);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","t#bar{t}Z","L");
   entry->SetLineColor(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","Total prediction","L");
   entry->SetLineColor(1);
   entry->SetLineStyle(2);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   TH1F *__7 = new TH1F("__7","",8,0.45,0.75);
   __7->SetBinContent(1,1.092524);
   __7->SetBinContent(2,1.212903);
   __7->SetBinContent(3,1.259505);
   __7->SetBinContent(4,1.378209);
   __7->SetBinContent(5,1.493965);
   __7->SetBinContent(6,1.672434);
   __7->SetBinContent(7,1.948337);
   __7->SetBinContent(8,2.604354);
   __7->SetEntries(12.66223);

   ci = TColor::GetColor("#990099");
   __7->SetLineColor(ci);
   __7->SetLineStyle(2);
   __7->SetLineWidth(3);
   __7->GetXaxis()->SetLabelFont(42);
   __7->GetXaxis()->SetLabelOffset(0.007);
   __7->GetXaxis()->SetLabelSize(0.05);
   __7->GetXaxis()->SetTitleSize(0.06);
   __7->GetXaxis()->SetTitleOffset(0.9);
   __7->GetXaxis()->SetTitleFont(42);
   __7->GetYaxis()->SetLabelFont(42);
   __7->GetYaxis()->SetLabelOffset(0.007);
   __7->GetYaxis()->SetLabelSize(0.05);
   __7->GetYaxis()->SetTitleSize(0.06);
   __7->GetYaxis()->SetTitleOffset(1.25);
   __7->GetYaxis()->SetTitleFont(42);
   __7->GetZaxis()->SetLabelFont(42);
   __7->GetZaxis()->SetLabelOffset(0.007);
   __7->GetZaxis()->SetLabelSize(0.05);
   __7->GetZaxis()->SetTitleSize(0.06);
   __7->GetZaxis()->SetTitleOffset(1);
   __7->GetZaxis()->SetTitleFont(42);
   __7->Draw("hist same");
   
   TH1F *__8 = new TH1F("__8","",8,0.45,0.75);
   __8->SetBinContent(1,1.376438);
   __8->SetBinContent(2,1.846465);
   __8->SetBinContent(3,2.04208);
   __8->SetBinContent(4,2.516835);
   __8->SetBinContent(5,2.961465);
   __8->SetBinContent(6,3.701627);
   __8->SetBinContent(7,4.768741);
   __8->SetBinContent(8,7.346688);
   __8->SetEntries(26.56034);

   ci = TColor::GetColor("#ff99cc");
   __8->SetLineColor(ci);
   __8->SetLineStyle(2);
   __8->SetLineWidth(3);
   __8->GetXaxis()->SetLabelFont(42);
   __8->GetXaxis()->SetLabelOffset(0.007);
   __8->GetXaxis()->SetLabelSize(0.05);
   __8->GetXaxis()->SetTitleSize(0.06);
   __8->GetXaxis()->SetTitleOffset(0.9);
   __8->GetXaxis()->SetTitleFont(42);
   __8->GetYaxis()->SetLabelFont(42);
   __8->GetYaxis()->SetLabelOffset(0.007);
   __8->GetYaxis()->SetLabelSize(0.05);
   __8->GetYaxis()->SetTitleSize(0.06);
   __8->GetYaxis()->SetTitleOffset(1.25);
   __8->GetYaxis()->SetTitleFont(42);
   __8->GetZaxis()->SetLabelFont(42);
   __8->GetZaxis()->SetLabelOffset(0.007);
   __8->GetZaxis()->SetLabelSize(0.05);
   __8->GetZaxis()->SetTitleSize(0.06);
   __8->GetZaxis()->SetTitleOffset(1);
   __8->GetZaxis()->SetTitleFont(42);
   __8->Draw("hist same");
   pad_ratio->Modified();
   pad_ratio->cd();
   TLatex *   tex = new TLatex(0.2,0.87,"CMS");
tex->SetNDC();
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.96,0.94,"138 fb^{-1} (13 TeV)");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.85,"SR-t#bar{t}Z");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   pad_ratio->Modified();
   c1_1->cd();
   
   leg = new TLegend(0.35,0.8,0.94,0.92,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("h_uncert","Unc.","F");
   entry->SetFillColor(1);
   entry->SetFillStyle(3254);
   entry->SetLineColor(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("","Data","ep");
   entry->SetLineColor(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("_stack_8","tZq","f");

   ci = TColor::GetColor("#e3413d");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#e3413d");
   entry->SetLineColor(ci);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("_stack_7","tWZ","f");

   ci = TColor::GetColor("#fce181");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#fce181");
   entry->SetLineColor(ci);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("_stack_6","t#bar{t}Z","f");

   ci = TColor::GetColor("#1f78b4");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#1f78b4");
   entry->SetLineColor(ci);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("_stack_5","t(#bar{t})X","f");

   ci = TColor::GetColor("#a6cee3");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#a6cee3");
   entry->SetLineColor(ci);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("_stack_4","WZ","f");

   ci = TColor::GetColor("#b2df8a");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#b2df8a");
   entry->SetLineColor(ci);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("_stack_3","VV(V)","f");

   ci = TColor::GetColor("#35b863");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#35b863");
   entry->SetLineColor(ci);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("_stack_2","X#gamma","f");

   ci = TColor::GetColor("#cc6699");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#cc6699");
   entry->SetLineColor(ci);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("_stack_1","NPL","f");

   ci = TColor::GetColor("#cccccc");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#cccccc");
   entry->SetLineColor(ci);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   c1_1->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
