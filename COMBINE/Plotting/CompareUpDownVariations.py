#-- Plot up/down variations (relative or absolute) for all bins/processes/nuisances found in datacard/workspace

import ROOT
import os
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)


#-------------------------------------

plotSingleNuisance = '' #If non-empty, only plot this specific nuisance; else plot all nuisances found in WS
# plotSingleNuisance = 'bcTagSF'

excludeNuisanceGroup = '' #If non-empty, don't plot nuisances whose name include the defined string

relative = True #True <-> Plot relative bin-per-bin deviations; else, superimpose 'absolute' histograms

#-- List of input ROOT files containing the template histograms (use only 1 for now)
files = [
    #'/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/2016_llBg1.root',
    # '/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/TOP-18-010/cards/shapecard_2016_differential.root',
    #'/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/EFTAnalysis/COMBINE/templates/Templates_Zpt_EFT2_Run2.root',
    #'/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/TOP-18-010/cards/shapecard_2016_inclusive.root',
    #'/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/TOP-18-010/cards/shapecard_2017_inclusive.root',
    #'/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/datacards_200101_2018v39_ttbbccRateBjets/Hut/forCombine_FCNC_Hut_Discriminant_DNN_Hut_all.root',
    '~/Templates_NN_ctz_EFT2_Run2.root'
]

#-- List of input ROOT worskpaces containing the physics model and all parameters (use only 1 for now)
workspaces = [
    #'/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/LFVScalarU_2016.root',
    # '/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/TOP-18-010/cards/diff_2016.root',
    # '/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/EFTAnalysis/COMBINE/datacards/COMBINED_Datacard_TemplateFit_noStat_Run2.root',
    #'/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/TOP-18-010/cards/incl_2016.root',
    #'/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/TOP-18-010/cards/incl_2017.root',
    '/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/TOPCombine_contact/Datacards/datacards_200101_2018v39_ttbbccRateBjets/Hut/combine_workspace.root',
]

#-- Protection: consider only 1 file and WS for now, to avoid ambiguities #Otherwise, make sure the files/WS names match !
if len(files) != 1 or len(workspaces) != 1: print('ERROR: please consider only 1 input file and WS for now')

#-------------------------------------

nplots = 0 #Counter
ifile = 0 #File counter
for fname, wsname in zip(files,workspaces): #Open each file and corresponding WS
    f = ROOT.TFile.Open(fname)
    #f.cd('xxx') #Enter subdir if necessary
    ws = ROOT.TFile.Open(wsname)
    ifile = ifile + 1

    print('Open file:', fname)
    print('Open workspace:', wsname)

    #-- Get full list of nuisance parameters
    list_nuisances = []
    w = ws.Get("w")
    # w.Print()
    mc = w.genobj("ModelConfig")
    # mc->GetParametersOfInterest()->Print("V") #POI
    # mc->GetNuisanceParameters()->Print("V") #Nuisances (remove 'V' to access RooArgSet* list)
    set_nuisances = mc.GetNuisanceParameters() #Get list of nuisances
    iter = set_nuisances.createIterator()
    var = iter.Next()
    while var :
        # print var.GetName()
        list_nuisances.append(var.GetName())
        var = iter.Next()

    #-- Set user-defined list of bins (may be used in naming convention; else, insert 1 dummy element)
    list_bins = ['']
    #list_bins = ['Zpt_SRtZq_2016','Zpt_SRttZ_2016','mTW_CR_2016']
    #list_bins = ['dc_2016','dc_2017','dc_2018']    

    #-- Set user-defined list of processes
    #list_processes = ['LFVTtScalarU', 'LFVStScalarU', 'tt', 'Jets', 'tW', 'DY', 'Others']
    #list_processes = ['signal', 'WG', 'TT_misID', 'ZG', 'fakes', 'other', 'QCD', 'other_misID', 'DY_misID']
    #list_processes = ['PrivMC_tZq','PrivMC_ttZ','VVV','tX','NPL','WZ']
    list_processes = ['Hut125', 'ttlf', 'ttcc', 'ttbb', 'other', 'qcd']

    subdir_name = 'DNN_Hut_b2j3'

#-- Will loop on bins/processes/nuisances and compare nom/up/down histos (make sure to follow the relevant conventions)

    ibin = 0 #Bin counter
    for bin in list_bins:
        ibin = ibin + 1
        
        # print('Bin', bin)

        #-- histo naming conventions 
        include_bin_prefix = False #True <-> look for histogram names starting with bin name; else, if nbins>1, will store corresponding plots separately #NB: if false, make sure you read the proper RootFile (hence, may read histos corresponding to a different bin)
        separator = '_' #Define how the keys are separated in histo names (e.g. '_', '__', etc.)

        if include_bin_prefix == False: #If bin name not part of histo name, multiple bins must refer to multiple files
            if len(list_bins) > 1 and len(list_bins) != len(files):
                print('ERROR: wrong list of bins ?')
                exit(1)
            elif ibin != ifile: continue #1 bin per file

        for proc in list_processes:

            # print('Process', proc)

            #-- Loop on list of nuisances
            for nuisance in list_nuisances:
                if plotSingleNuisance != '' and nuisance != plotSingleNuisance: continue #Check individual syst
                if excludeNuisanceGroup != '' and excludeNuisanceGroup in nuisance: continue #Don't consider these nuisance group

                # print('Nuisance', nuisance)

                #-- Histo names #Make sure to use relevant naming conventions !
                hname = proc
                if include_bin_prefix: hname = bin + separator + hname
                if subdir_name != '': hname = subdir_name + '/' + hname
                hnameDown = hname+separator+nuisance+'Down'
                hnameUp = hname+separator+nuisance+'Up'

                if subdir_name=='' and hname not in f.GetListOfKeys():
                    #print(key.GetName(), key.GetClassName())
                    print('ERROR: ', hname, 'not in ', fname, '! Skip !')
                    continue
                if subdir_name=='' and hnameDown not in f.GetListOfKeys():
                    #print('ERROR: ', hnameDown, 'not in ', fname, '! Skip !')
                    continue
                if subdir_name=='' and hnameUp not in f.GetListOfKeys():
                    #print('ERROR: ', hnameUp, 'not in ', fname, '! Skip !')
                    continue

                nplots+= 1

                print(hname)
                n_tmp = f.Get(hname)
                #print(n_tmp.GetClassName())
                print('INTEGRAL', hname, ':', n_tmp.Integral())
                # print(hnameDown)
                d_tmp = f.Get(hnameDown)
                print('INTEGRAL', hnameDown, ':', d_tmp.Integral())
                # print(hnameUp)
                u_tmp = f.Get(hnameUp)
                print('INTEGRAL', hnameUp, ':', u_tmp.Integral())

                #-- Protections
                if n_tmp.Integral()==0 and d_tmp.Integral()==0 and u_tmp.Integral()==0: continue

                #-- Clone objects, so that we can modify without e.g. impacting the nominal histo (reused for all plots for given process)
                n = n_tmp.Clone(); u = u_tmp.Clone(); d = d_tmp.Clone()
                n.Sumw2(); u.Sumw2(); d.Sumw2()

                if relative: #Plot deviations relative to nominal --> Substract nominal

                    rel_mode = 1 #0 <-> substract nominal; 1 <-> divide by nominal #NB: if substract nominal, then errors of (nom+variations) are added in quadrature... misleading ? While if divide(var/nom), it is properly computed

                    # print('BEFORE [%s --- %s --- %s] / [ %s --- %s --- %s]'%(d.Integral(), n.Integral(), u.Integral(),d.GetEntries(), n.GetEntries(), u.GetEntries()))
                    if rel_mode == 0: u.Add(n,-1); d.Add(n,-1); n.Add(n,-1)
                    else: u.Divide(n); d.Divide(n); n.Divide(n)
                    # print('AFTER [%s --- %s --- %s]'%(d.Integral(), n.Integral(), u.Integral()))
                    #u.Scale(1./u.Integral())
                    #d.Scale(1./d.Integral())
                    #u.Print()
                    #d.Print()
                    themax = abs(u.GetMaximum()*1.1)
                    if abs(d.GetMaximum()*1.1) > themax:
                        themax=abs(d.GetMaximum()*1.1)
                    themax=1.2 #hardcode y range
                    if rel_mode == 0: u.GetYaxis().SetRangeUser(-themax, themax) #Relative diff: expect ~symmetric deviations around 0
                    else: u.GetYaxis().SetRangeUser(1-(themax-1), themax)

                c = ROOT.TCanvas('c', 'c', 800, 800)
                c.cd()
                u.SetTitle('File %s'%fname)
                u.GetXaxis().SetTitle(str('Bin '+bin+' / Proc '+proc+' / Syst '+nuisance))
                u.GetXaxis().SetTitleSize(0.04)
                n.SetLineColor(ROOT.kBlack); n.SetFillColorAlpha(ROOT.kBlack, 0)
                u.SetLineColor(ROOT.kRed); u.SetFillColorAlpha(ROOT.kRed, 0)
                d.SetLineColor(ROOT.kBlue); d.SetFillColorAlpha(ROOT.kBlue, 0)
                n.SetMarkerSize(0); u.SetMarkerSize(0); d.SetMarkerSize(0)

                u.Draw("hist e")
                d.Draw("hist same e")
                #u.Draw("hist")
                #d.Draw("hist same")
                n.Draw("hist same")

                l = ROOT.TLegend(0.6, 0.75, 0.9, 0.9)
                l.AddEntry(n, 'Nominal', 'l')
                l.AddEntry(u, nuisance+'Up', 'l')
                l.AddEntry(d, nuisance+'Down', 'l')
                l.Draw()

                outputdir = './plots/'
                prefix = ''
                if not os.path.exists(outputdir): os.makedirs(outputdir)
                if relative:
                    outputdir = './plots/rel/'
                    if not os.path.exists(outputdir): os.makedirs(outputdir)
                else:
                    outputdir = './plots/abs/'
                    if not os.path.exists(outputdir): os.makedirs(outputdir)
                if include_bin_prefix == False and len(list_bins) > 1: 
                    outputdir = outputdir + bin + '/'
                    if not os.path.exists(outputdir): os.makedirs(outputdir)
                else: prefix = (bin+'_') if bin != '' else ''
                suffix = '_rel' if relative else '_abs'
                c.Print(str(outputdir+prefix+proc+'_'+nuisance+suffix+'.png'))

if nplots == 0: print('ERROR: no histograms found to plot ! Check naming conventions, bin/proc names, etc. !')

#-------------------------------------
