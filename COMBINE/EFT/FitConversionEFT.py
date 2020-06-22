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
verbose = 0 #(Dis)activate printouts
operators = [SM_name]+['ctz'] #List of operators for which to extract parameterizations #Ex.: ['SM']+['ctZ','ctW']
# //--------------------------------------------

#Setup
# //--------------------------------------------
if len(sys.argv) != 2: hist_file = "../templates/Templates_NN_SR_2017.root" #Default rootfile path
else: hist_file = sys.argv[1] #Rootfile path given in arg

fits = {} #Dict that will hold the parameterizations of the cross-sections
# //--------------------------------------------

#Read file and extract parameterizations for each process / category
# //--------------------------------------------
print(colors.fg.lightblue + "Loading Root file..." + colors.reset)
readfile = ROOT.TFile.Open(hist_file)
print(colors.fg.lightblue + "... Done !\n" + colors.reset)
print(colors.fg.lightblue + "Extracting parameterizations..." + colors.reset)

for key in readfile.GetListOfKeys():

    name = key.GetName()

    if 'PrivMC' not in name or 'bin' in name: continue #Get parametrization from 'full' TH1EFT histograms only

    if verbose: print('\nkey.GetName()', name)

    hist = readfile.Get(name)

    #Get categorical information
    histname = name.split('__') #Naming convention: 'categ__proc__syst'
    if verbose: print('histname', histname)
    full_bin_name,process,systematic = '','',''
    if(len(histname)==3): [full_bin_name,process,systematic] = histname
    if(len(histname)==2): [full_bin_name,process] = histname
    #process = process.replace('tllq','tZq')

    #Skip systematic histograms (only use nominal histograms for parametrization)
    if systematic != '': continue

    if verbose:
        print('full_bin_name', full_bin_name)
        print('process', process)
        # print('systematic', systematic)

    #-- Extract parameterization
    if 'PrivMC' in process and 'bin' not in full_bin_name: #Get parametrization from 'full' TH1EFT histograms only
        if verbose: print('process', process)

        # bin_name_tmp = full_bin_name.split('_', 1)[1] #'TH1EFT_NN_SR_2017_xxx' --> 'NN_SR_2017' #split(separator, maxsplit)
        # if verbose: print('bin_name_tmp', bin_name_tmp)

        #Loop through bins and extract parameterization
        if verbose: print('hist.GetNbinsX()', hist.GetNbinsX())
        for ibin in range(1, hist.GetNbinsX()+1):
            if verbose: print('ibin', ibin)
            fit = hist.GetBinFit(ibin)
            # fit = hist.GetSumFit()

            # bin_name = bin_name_tmp
            # bin_name = bin_name_tmp + '_{0}'.format(str(ibin))
            bin_name = 'bin{0}_'.format(str(ibin)) + full_bin_name
            if verbose: print('bin_name', bin_name)

            names = fit.getNames()
            if len(names)==0 or fit.getCoefficient(SM_name,SM_name)==0: continue
            elif len(names)!=0 and fit.getCoefficient(SM_name,SM_name)==0:
                for op1 in operators:
                    for op2 in operators:
                        if fit.getCoefficient(op1,op2)!=0:
                            print "Error! SM yield is 0, but this bin has a nonzero contribution from EFT effects! The parameterization for this bin will be ignored."
                            print "    "+process,categ
                            print "    "+op1,op2," ",round(fit.getCoefficient(op1,op2), 8)
            elif (process,bin_name) not in fits.keys(): fits[(process,bin_name)]={}

            #For given bin,
            for op1 in operators:
                for op2 in operators:
                    if verbose: print process, [op1,op2], fit.getCoefficient(op1,op2), round(fit.getCoefficient(op1,op2)/fit.getCoefficient(SM_name,SM_name), 8)
                    fits[(process,bin_name)][(op1,op2)] = round(fit.getCoefficient(op1,op2)/fit.getCoefficient(SM_name,SM_name), 8)

print(colors.fg.lightblue + "... Done !\n" + colors.reset)

#Summary printout and save results
# //--------------------------------------------
if verbose:
    # print "Processes:",[key[0] for key in fits.keys()]
    # print "Bins:",[key[1] for key in fits.keys()]
    print "Keys:",fits.keys()
    # print "Fits:", fits

#Store fits
np.save('EFT_Parameterization.npy', fits)
print(colors.fg.lightblue + "\n---> Stored fits in numpy file {}\n".format("EFT_Parameterization.npy") + colors.reset)
print('(Inspect dictionnary content in file: EFT_Parameterization.txt)\n')
txt_file = open("EFT_Parameterization.txt", "w")
txt_file.write(str(fits) + '\n')
txt_file.close()
