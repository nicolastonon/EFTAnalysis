import ROOT
from ROOT import *
from EFTPlotter import Load_Canvas_Style, Get_Parameter_LegName, Set_Custom_ColorGradient_Zaxis

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBlackBody)
ROOT.TColor.InvertPalette()


#============================
#============================


def Get_CovMatrix_TH2F(filepath, wcs=[]):

	nskip_par = 0 #Skip MCstat nuisances
	fIn=ROOT.TFile.Open(filepath)
	fit_mdf=fIn.Get("fit_mdf")
	poilist = []
	for wc in wcs:
		poilist.append(ROOT.RooRealVar(wc,wc,-10,10))

	if len(wcs)==0:
		covMat = fit_mdf.covarianceMatrix()
		floatParsFinal = fit_mdf.floatParsFinal() #RooArgList of parameters
		size = covMat.GetNrows()	
		for i in range(size):
			if 'prop_bin' in floatParsFinal[i].GetName(): nskip_par+= 1

	else:
		pois = ROOT.RooArgList(*poilist)
		covMat = fit_mdf.reducedCovarianceMatrix(pois)

	fIn.Close()

	size = covMat.GetNrows()
	npars = size-nskip_par
	#print('Size = ', size)

	hist_cov = ROOT.TH2F("cov", "cov", npars, 0.5 , npars+0.5, npars, 0.5, npars+0.5)
	i_tmp = 0; j_tmp = 0 #Without skip par
	for i in range(size):
		j_tmp = 0 #Reinit
		if len(wcs)>0: 
			hist_cov.GetXaxis().SetBinLabel(i+1, Get_Parameter_LegName(wcs[i], fullname=False))
		else: 
			if 'prop_bin' in floatParsFinal[i].GetName(): continue
			#print('i', floatParsFinal[i].GetName())
			hist_cov.GetXaxis().SetBinLabel(i_tmp+1, floatParsFinal[i].GetName())
			#print('i_tmp', i_tmp)

		for j in range(size):
			if len(wcs)>0: #WCs only
				hist_cov.GetYaxis().SetBinLabel(j+1, Get_Parameter_LegName(wcs[j], fullname=False))
				hist_cov.SetBinContent(i+1, j+1, covMat(i,j))
			else: #All params
				if 'prop_bin' in floatParsFinal[j].GetName(): continue
				#print('j', floatParsFinal[j].GetName())
				hist_cov.SetBinContent(i_tmp+1, j_tmp+1, covMat(i,j))
				hist_cov.GetYaxis().SetBinLabel(j_tmp+1, floatParsFinal[j].GetName())
				#print('j_tmp', j_tmp)
				j_tmp+= 1
			#print(covMat(i,j))

		if len(wcs)==0: i_tmp+= 1

	return hist_cov


def Get_CorrMatrix_TH2F(filepath, wcs=[]):

	fIn=ROOT.TFile.Open(filepath)
	fit_mdf=fIn.Get("fit_mdf")

	idx_skip_pars = [] #Skip MCstat nuisances (and more)

	size = len(wcs)
	if size==0:
		corrMat = fit_mdf.correlationMatrix()
		size = corrMat.GetNrows()
		floatParsFinal = fit_mdf.floatParsFinal() #RooArgList of parameters
		for i in range(size):
			if 'prop_bin' in floatParsFinal[i].GetName(): 
				idx_skip_pars.append(i)
				continue
			has_large_corr = False
			for j in range(size): 
				if abs(corrMat(i,j)) > 0.30 and corrMat(i,j) != 1: has_large_corr = True #Skip parameters that don't have any large correl
			if has_large_corr == False: idx_skip_pars.append(i)
		#print('len(idx_skip_pars)', len(idx_skip_pars))

	fIn.Close()

	npars = size-len(idx_skip_pars)

	hist_corr = ROOT.TH2F("corr", "corr", npars, 0.5 , npars+0.5, npars, 0.5, npars+0.5)
	#if npars > 15: 
	#	hist_corr.GetXaxis().SetLabelSize(0.01)
	#	hist_corr.GetYaxis().SetLabelSize(0.01)

	i_tmp = 0; j_tmp = 0 #Without skip par
	for i in range(size):
		j_tmp = 0 #Reinit
		if len(wcs)>0: hist_corr.GetXaxis().SetBinLabel(i+1, Get_Parameter_LegName(wcs[i], fullname=False))
		else: 
			#if 'prop_bin' in floatParsFinal[i].GetName(): continue
			if i in idx_skip_pars: continue
			#print('i', i)
			#print('i_tmp', i_tmp)
			hist_corr.GetXaxis().SetBinLabel(i_tmp+1, floatParsFinal[i].GetName())

		for j in range(size):
			if len(wcs)>0: #WCs only
				hist_corr.SetBinContent(i+1,j+1, fit_mdf.correlation(wcs[i],wcs[j]))
				hist_corr.GetYaxis().SetBinLabel(j+1, Get_Parameter_LegName(wcs[j], fullname=False))
				#print(fit_mdf.correlation(wcs[i],wcs[j])*100)
			else: #All params
				#if 'prop_bin' in floatParsFinal[j].GetName(): continue
				if j in idx_skip_pars: continue
				#print('j', j)
				#print('j_tmp', j_tmp)
				hist_corr.SetBinContent(i_tmp+1,j_tmp+1, corrMat(i,j))
				hist_corr.GetYaxis().SetBinLabel(j_tmp+1, floatParsFinal[j].GetName())
				j_tmp+= 1
				#print(corrMat(i,j))

		if len(wcs)==0: i_tmp+= 1

	#print(hist_corr)
	return hist_corr


def Plot_TH2F(h2, is_cov, outputname, roundN=1, xTitle=r"", yTitle=r"", zTitle=r"", textlabel=False):
	"""
	h2: histogram to plot
	"""

	if h2 is None: 
		print 'ERROR: h2 is None'
		return

	if is_cov: 
		covMatrix = True #Different settings
		roundN = 3

	#print('h2', h2)

	leftmargin = 0.125
	bottommargin = 0.15
	rightmargin = 0.2
	topmargin = 0.075

	canvas = ROOT.TCanvas("c", " ", 800, 600)
	pad1 = ROOT.TPad("pad1", "pad1", 0., 0., 1, 1.0)
	pad1.SetLeftMargin(leftmargin)
	if not is_cov: pad1.SetRightMargin(rightmargin)
	pad1.SetTopMargin(topmargin)
	pad1.SetBottomMargin(bottommargin)
	pad1.Draw()
	pad1.cd()

	textsize = 32./(pad1.GetWh()*pad1.GetAbsHNDC())
	if not is_cov: h2.Scale(100)

	h2.GetXaxis().SetTitle(xTitle)
	h2.GetXaxis().SetTitleFont(42)
	h2.GetXaxis().SetTitleSize(textsize)
	h2.GetYaxis().SetTitle(yTitle)
	h2.GetYaxis().SetTitleFont(42)
	h2.GetYaxis().SetTitleSize(textsize)
	h2.GetZaxis().SetTitleSize(textsize)
	h2.GetZaxis().SetLabelSize(textsize)
	h2.GetZaxis().SetTitleSize(0.07)
	h2.GetZaxis().SetTitleOffset(1)
	h2.GetZaxis().SetTitle(zTitle)
	h2.SetMarkerSize(2.0)
        h2.LabelsOption("v", "X")

	#h2.GetXaxis().SetLabelSize(textsize*1.)
	#h2.GetYaxis().SetLabelSize(textsize*1.)
	h2.GetXaxis().SetLabelSize(0.1)
	h2.GetYaxis().SetLabelSize(0.1)

	h2.SetTitle("")
	ROOT.gStyle.SetPaintTextFormat("4.{0}f".format(roundN)) #"4.{0}f".format(roundN) <-> up to 4 digits before point, and 'roundN' after

	if is_cov: h2.Draw("TEXT")
	elif textlabel: h2.Draw("COLZ TEXT")
	else: h2.Draw("COLZ")

	#-- Labels
        left = 0.13; top = 0.94
        CMS_text = ROOT.TLatex(left, top, "CMS")
        CMS_text.SetNDC()
        CMS_text.SetTextColor(ROOT.kBlack)
        CMS_text.SetTextFont(61)
        CMS_text.SetTextAlign(11)
        CMS_text.SetTextSize(0.07)
	CMS_text.Draw('same')

        extraText = ROOT.TLatex(left+0.12, top, "Supplementary")
        extraText.SetNDC()
        extraText.SetTextFont(52)
        extraText.SetTextSize(0.05)
	extraText.Draw('same')

        lumiText = ROOT.TLatex(0.80, top, "138 fb^{-1} (13 TeV)")
        if is_cov: lumiText = ROOT.TLatex(0.96, top, "138 fb^{-1} (13 TeV)")
        lumiText.SetNDC()
        lumiText.SetTextFont(42)
        lumiText.SetTextAlign(31)
        lumiText.SetTextSize(0.05)
	lumiText.Draw('same')

	canvas.SaveAs(outputname)
	canvas.Close()

	return


#============================
#============================


if __name__ == "__main__":

	filepath = "fit_results/multidimfit.5DExp.root"

	wcs = ["ctz", "ctw", "cpq3", "cpqm", "cpt"]

	#Styling
	Load_Canvas_Style() 
	#ROOT.gStyle.SetPalette(ROOT.kBird)
	#ROOT.gStyle.SetPalette(ROOT.kAquamarine)
	ROOT.gStyle.SetPalette(ROOT.kBlackBody)
	#ROOT.TColor.InvertPalette()

#-- WCs only

	h2 = Get_CovMatrix_TH2F(filepath, wcs)
	Plot_TH2F(h2, True, "cov_matrix.png", textlabel=True, zTitle="Covariance [TeV^{-4}]")

	h2 = Get_CorrMatrix_TH2F(filepath, wcs)
	Plot_TH2F(h2, False, "corr_matrix.png", textlabel=True, zTitle="Correlation [%]")

#-- All parameters

	#h2 = Get_CovMatrix_TH2F(filepath)
	#Plot_TH2F(h2, True, "cov_matrix_all.png")

	#h2 = Get_CorrMatrix_TH2F(filepath)
	#Plot_TH2F(h2, False, "corr_matrix_all.png")

