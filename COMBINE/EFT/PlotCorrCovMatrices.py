import ROOT
from ROOT import *

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBlackBody)
ROOT.TColor.InvertPalette()


def Get_CovMatrix_TH2F(filepath, wcs):

    fIn=ROOT.TFile.Open(filepath)
    fit_mdf=fIn.Get("fit_mdf")
    poilist = []
    for wc in wcs:
	    poilist.append(ROOT.RooRealVar(wc,wc,-10,10))

    pois = ROOT.RooArgList(*poilist)

    covMat = fit_mdf.reducedCovarianceMatrix(pois)
    fIn.Close()

    hist_cov = ROOT.TH2F("cov", "cov", len(poilist), 0.5 , len(poilist)+0.5, len(poilist), 0.5, len(poilist)+0.5)
    for i in range(len(poilist)):
	    hist_cov.GetXaxis().SetBinLabel(i+1, wcs[i])
	    for j in range(len(poilist)):
	        #print("covMat({},{})".format(i+1,j+1), covMat(i,j))
	        hist_cov.GetYaxis().SetBinLabel(j+1, wcs[j])
	        hist_cov.SetBinContent(i+1, j+1, covMat(i,j))
	
    #print(covMat(0,0))
    #print(covMat(5,5))
    
    #hists_cov.append(hist_cov)
    #c = ROOT.TCanvas("", "")
    #hist_cov.Draw()
    #c.SaveAs("test.png")

    return hist_cov


def Get_CorrMatrix_TH2F(filepath, wcs):

    fIn=ROOT.TFile.Open(filepath)
    fit_mdf=fIn.Get("fit_mdf")
    fIn.Close()
    hist_corr = ROOT.TH2F("corr", "corr", len(wcs), 0.5 , len(wcs)+0.5, len(wcs), 0.5, len(wcs)+0.5)
    for i, p1 in enumerate(wcs):
        #print('p1', p1)
		hist_corr.GetXaxis().SetBinLabel(i+1, wcs[i])
		for j, p2 in enumerate(wcs):
			#print('p2', p2)
			#print('fit_mdf.correlation(p1,p2)', fit_mdf.correlation(p1,p2))
			hist_corr.SetBinContent(i+1,j+1, fit_mdf.correlation(p1,p2))
			hist_corr.GetYaxis().SetBinLabel(j+1, wcs[j])

    #print(hist_corr)
    return hist_corr


def Plot_TH2F(h2, outputname, roundN=2, in_percent=True, xTitle=r"", yTitle=r""):
    """
    h2: histogram to plot
    """

    #print('h2', h2)

    leftmargin = 0.125
    bottommargin = 0.15
    rightmargin = 0.2
    topmargin = 0.075

    canvas = ROOT.TCanvas("c", " ", 800, 600)
    pad1 = ROOT.TPad("pad1", "pad1", 0., 0., 1, 1.0)
    pad1.SetLeftMargin(leftmargin)
    pad1.SetRightMargin(rightmargin)
    pad1.SetTopMargin(topmargin)
    pad1.SetBottomMargin(bottommargin)
    pad1.Draw()
    pad1.cd()

    textsize = 32./(pad1.GetWh()*pad1.GetAbsHNDC())
    if in_percent: h2.Scale(100)

    h2.GetXaxis().SetTitle(xTitle)
    h2.GetXaxis().SetTitleFont(42)
    h2.GetXaxis().SetTitleSize(textsize)
    h2.GetYaxis().SetTitle(yTitle)
    h2.GetYaxis().SetTitleFont(42)
    h2.GetYaxis().SetTitleSize(textsize)
    h2.GetZaxis().SetTitleSize(textsize)
    h2.GetXaxis().SetLabelSize(textsize*1.5)
    h2.GetYaxis().SetLabelSize(textsize*1.5)
    h2.GetZaxis().SetLabelSize(textsize)
    h2.GetZaxis().SetTitleSize(textsize)
    h2.GetZaxis().SetTitleOffset(1.4)
    h2.GetZaxis().SetTitle("")
    h2.SetMarkerSize(2.0)

    h2.SetTitle("")
    if in_percent:
        ROOT.gStyle.SetPaintTextFormat("3.0f")
    else:
        ROOT.gStyle.SetPaintTextFormat("4.{0}f".format(roundN))

    h2.Draw("COLZ TEXT")

    canvas.SaveAs(outputname)
    canvas.Close()

    return


if __name__ == "__main__":

    filepath = "fit_results/multidimfit.5DObs.root"

    wcs = ["ctz", "ctw", "cpq3", "cpqm", "cpt"]

    h2 = Get_CovMatrix_TH2F(filepath, wcs)
    Plot_TH2F(h2, "cov_matrix.png", in_percent=False)
    
    h2 = Get_CorrMatrix_TH2F(filepath, wcs)
    Plot_TH2F(h2, "corr_matrix.png", in_percent=True)

