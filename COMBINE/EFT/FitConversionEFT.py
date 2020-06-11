print "Importing libraries..."
import numpy as np
import itertools
import os
import ROOT
import sys
ROOT.gSystem.Load('/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/EFTAnalysis/myLib.so')

#Dict that will hold the parameterizations of the cross-sections
fits = {}

if len(sys.argv) != 2:
    hist_file = "../templates/Templates_NN_SR_2017.root"
else:
    hist_file = sys.argv[1]

#List of operators to extract parameterizations for
operators = ['sm']+['ctZ']

#Load file
print "Loading Root file..."
readfile = ROOT.TFile.Open(hist_file)

#Crawl through file
print "Extracting parameterizations..."
for key in readfile.GetListOfKeys():
    # print('Key', key)
    hist = readfile.Get(key.GetName())

    #Get categorical information
    histname = hist.GetName().split('__')
    category,systematic,process = '','',''
    if(len(histname)==3): [category,process,systematic] = histname
    if(len(histname)==2): [category,process] = histname
    #process = process.replace('tllq','tZq')

    #Skip systematic histograms
    if systematic != '': continue

    #Only use histograms from WC samples
    if 'TH1EFT' in category:
        #process = process.split('_',1)[0]

        #Loop through bins and extract parameterization
        for bin in range(1,10): # Doesn't include last bin
            fit = hist.GetBinFit(bin)

            category_njet=''
			#Debug
            #if bin==4:
            #    print category_njet,fit.getCoefficient('sm','sm'),fit.getCoefficient('ctW','sm'),fit.getCoefficient('ctW','ctW')
            #print category_njet,fit.getCoefficient('sm','sm')
            names = fit.getNames()
            if len(names)!=0 and fit.getCoefficient('sm','sm')==0:
                for op1 in operators:
                    for op2 in operators:
                        if fit.getCoefficient(op1,op2)!=0:
                            print "Error! SM yield is 0, but this bin has a nonzero contribution from EFT effects! The parameterization for this bin will be ignored."
                            print "    "+process,category_njet
                            print "    "+op1,op2," ",round(fit.getCoefficient(op1,op2),8)

            if len(names)==0 or fit.getCoefficient('sm','sm')==0: continue
            if (process,category_njet) not in fits.keys(): fits[(process,category_njet)]={}
#            for pair in itertools.combinations_with_replacement(operators,2):
#                fits[(process,category_njet)][pair] = round(fit.getCoefficient(pair[0],pair[1])/fit.getCoefficient('sm','sm'),8)
            for op1 in operators:
                for op2 in operators:
                    #Debug
                    #if category_njet=='C_2lss_m_ee_2b_ge7j' and op2 == 'sm' and op1 != 'sm': print process,[op1,op2],fit.getCoefficient(op1,op2)
                    #if category_njet=='C_2lss_m_ee_2b_ge7j' and op2 == 'cbW' and op1 != 'sm': print process,[op1,op2],fit.getCoefficient(op1,op2)
                    #if category_njet=='C_2lss_m_ee_2b_ge7j' and op2 == 'cbW' and op1 != 'cbW': print process,[op1,op2],fit.getCoefficient(op1,op2)
                    fits[(process,category_njet)][(op1,op2)] = round(fit.getCoefficient(op1,op2)/fit.getCoefficient('sm','sm'),8)

#print fits
#print "Categories:",[key[1] for key in fits.keys()]
#print "Processes:",[key[0] for key in fits.keys()]
#print "Keys:",fits.keys()
np.save('EFT_Parameterization.npy', fits)
print "\nSaved numpy file {}".format("EFT_Parameterization.npy !\n")
