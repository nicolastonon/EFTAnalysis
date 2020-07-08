// https://root.cern.ch/root/html534/guides/users-guide/AddingaClass.html

#include "TH1EFT.h"

using namespace std;

// ROOT needs this here:
ClassImp(TH1EFT) //-- Needed to include custom class within ROOT

TH1EFT::TH1EFT() {}
TH1EFT::~TH1EFT() {}

TH1EFT::TH1EFT(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup)
 : TH1D (name, title, nbinsx, xlow, xup)
{
    // Create/Initialize a fit function for each bin in the histogram
    WCFit new_fit;
    for (Int_t i = 0; i < nbinsx; i++) {
        this->hist_fits.push_back(new_fit);
    }
}

void TH1EFT::SetBins(Int_t nx, Double_t xmin, Double_t xmax)
{
    // Use this function with care! Non-over/underflow bins are simply
    // erased and hist_fits re-sized with empty fits.

    hist_fits.clear();
    WCFit new_fit;
    for (Int_t i = 0; i < nx; i++) {this->hist_fits.push_back(new_fit);}

    TH1::SetBins(nx, xmin, xmax);

    return;
}

Bool_t TH1EFT::Add(const TH1 *h1, Double_t c1)
{
    // check whether the object pointed to inherits from (or is a) TH1EFT:
    // if (true) { //Obtain class with IsA(), check if inherits from TH1EFT
    if (h1->IsA()->InheritsFrom(TH1EFT::Class())) { //Obtain class with IsA(), check if inherits from TH1EFT //should not work without ClassDef/ClassImp... ?
        if (this->hist_fits.size() == ((TH1EFT*)h1)->hist_fits.size()) {
            for (unsigned int i = 0; i < this->hist_fits.size(); i++) {
                // assumes this hist and the one whose fits we're adding have the same bins!
                this->hist_fits[i].addFit( ((TH1EFT*)h1)->hist_fits[i] );
            }
        } else {
            std::cout << "Attempt to add 2 TH1EFTs with different # of fits!" << std::endl;
            std::cout << this->hist_fits.size() << ", " << ((TH1EFT*)h1)->hist_fits.size() << std::endl;
        }
        this->overflow_fit.addFit( ((TH1EFT*)h1)->overflow_fit );
        this->underflow_fit.addFit( ((TH1EFT*)h1)->underflow_fit );
    }

    return TH1::Add(h1,c1);
}

// Custom merge function for using hadd
Long64_t TH1EFT::Merge(TCollection* list)
{
    TIter nexthist(list);
    TH1EFT *hist;
    while ((hist = (TH1EFT*)nexthist.Next())) {
        if (this->hist_fits.size() != hist->hist_fits.size()) {
            std::cout << "[WARNING] Skipping histogram with different # of fits" << std::endl;
            continue;
        }
        for (unsigned int i = 0; i < this->hist_fits.size(); i++) {
            this->hist_fits.at(i).addFit(hist->hist_fits.at(i));
        }
        this->overflow_fit.addFit(hist->overflow_fit);
        this->underflow_fit.addFit(hist->underflow_fit);
    }

    return TH1::Merge(list);
}

Int_t TH1EFT::Fill(Double_t x, Double_t w, WCFit& fit)
{
    Int_t bin_idx = this->FindFixBin(x) - 1;
    Int_t nhists  = this->hist_fits.size();

    if (bin_idx >= nhists) {
        // For now ignore events which enter overflow bin
        this->overflow_fit.addFit(fit);
        return Fill(x,w);
    } else if (bin_idx < 0) {
        // For now ignore events which enter underflow bin
        this->underflow_fit.addFit(fit);
        return Fill(x,w);
    }
    this->hist_fits.at(bin_idx).addFit(fit);
    return Fill(x,w); // the original TH1D member function
}

// Returns a fit function for a particular bin (no checks are made if the bin is an over/underflow bin)
WCFit TH1EFT::GetBinFit(Int_t bin)
{
    Int_t nhists = this->hist_fits.size();
    if (bin <= 0) {
        return this->underflow_fit;
    } else if (bin > nhists) {
        return this->overflow_fit;
    }
    return this->hist_fits.at(bin - 1);
}

// Returns a WCFit whose structure constants are determined by summing structure constants from all bins
WCFit TH1EFT::GetSumFit()
{
    WCFit summed_fit;
    for (unsigned int i = 0; i < this->hist_fits.size(); i++) {
        summed_fit.addFit(this->hist_fits.at(i));
    }
    return summed_fit;
}

// Returns a bin scaled by the the corresponding fit evaluated at a particular WC point
Double_t TH1EFT::GetBinContent(Int_t bin, WCPoint wc_pt)
{
    if(GetBinContent(bin) == 0) {return 0.;}
    if(this->GetBinFit(bin).getDim() <= 0) {return GetBinContent(bin);} // We don't have a fit for this bin, return regular bin contents

    return this->GetBinFit(bin).evalPoint(&wc_pt);
}

Double_t TH1EFT::GetBinContent(Int_t bin, std::string wc_pt_name)
{
    if(wc_pt_name=="") {return this->TH1D::GetBinContent(bin);}
    WCPoint wc_pt(wc_pt_name);
    return this->TH1EFT::GetBinContent(bin, wc_pt);
}

//NB : only name of WCPoint matters (not its weight)
void TH1EFT::Scale(WCPoint wc_pt)
{
    // Warning: calling GetEntries after a call to this function will return a
    // non-zero value, even if the histogram was never filled.

    if(!Check_WCPoint_Operators(wc_pt))
    {
        for (Int_t i = 1; i <= this->GetNbinsX(); i++)
        {
            this->SetBinContent(i,0);
        }
        return;
    }

    for (Int_t i = 1; i <= this->GetNbinsX(); i++)
    {
        Double_t new_content = this->GetBinContent(i,wc_pt);
        Double_t new_error = (GetBinFit(i)).evalPointError(&wc_pt);
        this->SetBinContent(i,new_content);
        this->SetBinError(i,new_error);
    }

    return;
}

//Overload regular TH1D::Scale() function, because not always found for some reason. Trick: had to use a different name to avoid ambiguity...
void TH1EFT::Scaler(Double_t sf)
{
    this->TH1D::Scale(sf);
    return;
}

// Uniformly scale all fits by amt
void TH1EFT::ScaleFits(double amt)
{
    for (uint i = 0; i < this->hist_fits.size(); i++) {
        this->hist_fits.at(i).scale(amt);
    }
}

// Display the fit parameters for all bins
void TH1EFT::DumpFits()
{
    cout<<"this->hist_fits.size() "<<this->hist_fits.size()<<endl;

    for (uint i = 0; i < this->hist_fits.size(); i++) {
        this->hist_fits.at(i).dump();
    }
}

bool TH1EFT::Check_WCPoint_Operators(WCPoint& pt)
{
    WCFit sumfit = this->GetSumFit();

    for(auto& kv: pt.inputs)
    {
        // cout<<"kv.first "<<kv.first<<endl;
        // cout<<"kv.second "<<kv.second<<endl;
        bool found = false;
        vector<std::string> names = sumfit.getNames();
        for(int i=0; i < names.size(); i++) //For each pair of coeffs
        {
            // cout<<"name "<<names[i]<<endl;
            if(kv.first == names[i]) {found = true;}
        }

        if(!found) {cout<<endl<<"ERROR ! Operator ["<<kv.first<<"] not found in the weight parameterization..."<<endl<<endl; return false;}
    }

    return true;
}

//Set the WCFit object associated with the given bin
void TH1EFT::Set_WCFit_Bin(int ibin, WCFit fit)
{
    this->hist_fits.at(ibin-1) = fit; //Bins start at 1; hist_fits indices start at 0

    return;
}
