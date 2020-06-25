// https://root.cern.ch/root/html534/guides/users-guide/AddingaClass.html

//TODO:
//- Ways to deal with under/overflow ?


#ifndef TH1EFT_H_
#define TH1EFT_H_

#include "TH1D.h"
#include <vector>
#include "TClass.h"
#include "WCFit.h"
#include "WCPoint.h"

using namespace std;

class TH1EFT : public TH1D
{
    public:

        // ROOT needs these
        TH1EFT();
        ~TH1EFT();

        // usual constructor:
        TH1EFT(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup);

        std::vector<WCFit> hist_fits;
        WCFit overflow_fit;
        WCFit underflow_fit;

        using TH1D::Fill;           // Bring the TH1D Fill fcts into scope
        using TH1D::GetBinContent;  // Bring the TH1D GetBinContent fcts into scope
        using TH1D::Scale;          // Bring the TH1D Scale fcts into scope (likely not needed)

        Int_t Fill(Double_t x, Double_t w, WCFit& fit);
        WCFit GetBinFit(Int_t bin);
        WCFit GetSumFit();
        Double_t GetBinContent(Int_t bin, WCPoint wc_pt);
        void Scale(WCPoint wc_pt);
        void ScaleFits(double amt);
        void DumpFits();
        bool Check_WCPoint_Operators(WCPoint&);

        void SetBins (Int_t nx, Double_t xmin, Double_t xmax);  // overriding virtual function from TH1
        Bool_t Add(const TH1 *h1, Double_t c1=1); // overriding virtual function from TH1
        Long64_t Merge(TCollection* list);
        void Set_WCFit_Bin(int, WCFit);

        ClassDef(TH1EFT,1) // Needed to include custom class within ROOT
};

/*
// ROOT needs this here:
// ClassImp(TH1EFT) //-- Needed to include custom class within ROOT

TH1EFT::TH1EFT() {}
TH1EFT::~TH1EFT() {}
*/

#endif
