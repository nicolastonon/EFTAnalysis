
#include <iostream>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"

using namespace std;

//Draw plots of nuisance params from file "output.root" in same directory
void Draw_Canvas()
{
	TString filename = "./nuisances.root"; //Change path here

	TFile* f = 0;
	f = TFile::Open(filename); if(!f) {cout<<"File not found !"<<endl; return;}

	TCanvas* c1 = 0;
	c1 = (TCanvas*) f->Get("nuisances"); if(!c1) {cout<<"c1 null!"<<endl; return;}
	c1->Draw();
	c1->SaveAs("./nuisances.png");

	TCanvas* c2 = 0;
	c2 = (TCanvas*) f->Get("post_fit_errs"); if(!c2) {cout<<"c2 null!"<<endl; return;}
	c2->Draw();
	c2->SaveAs("./post_fit_errs.png");
}
