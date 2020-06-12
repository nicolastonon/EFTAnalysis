#-- Extract EFT parameterization from TH1EFT objects
# Adapted from: https://github.com/cms-govner/EFTFit

from Utils.ColoredPrintout import colors
import numpy as np
import os
import ROOT
import sys
print(colors.fg.lightblue + "Importing libraries..." + colors.reset)
ROOT.gSystem.Load('/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/EFTAnalysis/myLib.so')

#USER OPTIONS
# //--------------------------------------------
SM_name = 'SM' #SM point naming convention
verbose = 1 #(Dis)activate printouts
operators = [SM_name]+['ctz'] #List of operators for which to extract parameterizations #Ex.: ['SM']+['ctZ','ctW']
# //--------------------------------------------

#Setup
# //--------------------------------------------
fits = {} #Dict that will hold the parameterizations of the cross-sections

if len(sys.argv) != 2: hist_file = "../templates/Templates_NN_SR_2017.root" #Default rootfile path
else: hist_file = sys.argv[1] #Rootfile path given in arg
# //--------------------------------------------

#Read file and extract parameterizations for each process / category
# //--------------------------------------------
#Load file
print(colors.fg.lightblue + "Loading Root file..." + colors.reset)
readfile = ROOT.TFile.Open(hist_file)
print(colors.fg.lightblue + "... Done !\n" + colors.reset)
print(colors.fg.lightblue + "Extracting parameterizations..." + colors.reset)
for key in readfile.GetListOfKeys():
    print('Key', key)
    print('key.GetName()', key.GetName())
    hist = readfile.Get(key.GetName())

    #Get categorical information
    histname = key.GetName().split('__') #Naming convention: 'categ__proc__syst'
    print('histname', histname)
    bin_name,process,systematic = '','',''
    if(len(histname)==3): [bin_name,process,systematic] = histname
    if(len(histname)==2): [bin_name,process] = histname
    #process = process.replace('tllq','tZq')

    #Skip systematic histograms
    if systematic != '': continue

    if verbose:
        print('bin_name', bin_name)
        print('process', process)
        # print('systematic', systematic)

    #Extract parameterization from TH1EFT objects only
    if 'TH1EFT' in bin_name:
        # print('process.split(_,1)', process.split('_',1))
        process = process.split('_')[1] #split(separator, maxsplit) #'PrivMC_tZq_xxx' --> 'tZq'
        if verbose: print('process', process)

        categ = bin_name.split('_')[2] #'TH1EFT_NN_SR_2017_xxx' --> 'SR'
        if verbose: print('categ', categ)

        #Loop through bins and extract parameterization
        if verbose: print('hist.GetNbinsX()', hist.GetNbinsX())
        for ibin in range(1, hist.GetNbinsX()+1):
            if verbose: print('ibin', ibin)
            fit = hist.GetBinFit(ibin)

            names = fit.getNames()
            if len(names)!=0 and fit.getCoefficient(SM_name,SM_name)==0:
                for op1 in operators:
                    for op2 in operators:
                        if fit.getCoefficient(op1,op2)!=0:
                            print "Error! SM yield is 0, but this bin has a nonzero contribution from EFT effects! The parameterization for this bin will be ignored."
                            print "    "+process,categ
                            print "    "+op1,op2," ",round(fit.getCoefficient(op1,op2),8)

            elif (process,categ,str(ibin)) not in fits.keys(): fits[(process,categ,str(ibin))]={}
            elif len(names)==0 or fit.getCoefficient(SM_name,SM_name)==0: continue

            #For given bin,
            for op1 in operators:
                for op2 in operators:
                    if verbose: print process, [op1,op2], fit.getCoefficient(op1,op2), round(fit.getCoefficient(op1,op2)/fit.getCoefficient(SM_name,SM_name), 8)
                    fits[(process,categ,str(ibin))][(op1,op2)] = round(fit.getCoefficient(op1,op2)/fit.getCoefficient(SM_name,SM_name), 8)

print(colors.fg.lightblue + "... Done !\n" + colors.reset)

#Summary printout and save results
# //--------------------------------------------
if verbose:
    print "Processes:",[key[0] for key in fits.keys()]
    print "Categories:",[key[1] for key in fits.keys()]
    print "Bins:",[key[2] for key in fits.keys()]
    print "Keys:",fits.keys()
    print "Fits:", fits

#Store fits
np.save('EFT_Parameterization.npy', fits)
print(colors.fg.lightblue + "\nStored fits in numpy file {}\n".format("EFT_Parameterization.npy") + colors.reset)
