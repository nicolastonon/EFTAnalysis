'''
=========
Store RooParametricHist (for EFT signal histos) and RooDataHist (for all other histograms) objects in a temporary workspace, which can then be called in a text datacard defining the full shape analysis.
=========

- See RooParametricHist example: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#rooparametrichist-gamman-for-shapes

- Naming conventions: rrv=RooRealVar, rph=RooParametricHist, rdh=RooDataHist, ...

- Example running command:
[[ python create_WS.py -d ../datacards/COMBINED_Datacard_TemplateFit_Run2.txt -f ./Parameterization_EFT.npy -t ../templates/Templates_Zpt_EFT2_Run2.root -o ../templates/Templates_otherRegions_Run2.root ]]

Written by Nicolas Tonon (DESY), 2021
'''

#- test 'useless' removal ?

import os
import stat
import sys
import logging
import subprocess as sp
from optparse import OptionParser
import numpy as np
from collections import defaultdict
from collections import OrderedDict
import getopt # command line parser
import argparse
import ROOT
from ROOT import *

from Utils.ColoredPrintout import colors
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from settings import opts #Custom dictionnary of settings
# ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")

#-- USER OPTIONS
# //--------------------------------------------

min_threshold_EFTcoeff = 0.0001 #Don't consider EFT terms whose coefficients are below this threshold (negligible)
min_threshold_MCstat = 0.01 #Don't include MCstat nuisance in RPH bin parameterization if relError is less than this relative threshold


# //--------------------------------------------
# //--------------------------------------------


##     ## ######## ##       ########  ######## ########
##     ## ##       ##       ##     ## ##       ##     ##
##     ## ##       ##       ##     ## ##       ##     ##
######### ######   ##       ########  ######   ########
##     ## ##       ##       ##        ##       ##   ##
##     ## ##       ##       ##        ##       ##    ##
##     ## ######## ######## ##        ######## ##     ##

def PrintBanner():
    '''
    Print some info
    '''

    print('\n' + colors.bg.orange + '                                           ' + colors.reset)
    print(colors.fg.orange + '-- EFTWorkspace --' + colors.reset + '\n')
    print('* NB: make sure you have extracted the relevant EFT parameterizations')
    print(colors.bg.orange + '                                           ' + colors.reset + '\n')

    return


# //--------------------------------------------
# //--------------------------------------------


 ######  ########  ########    ###    ######## ########    ##      ##  ######
##    ## ##     ## ##         ## ##      ##    ##          ##  ##  ## ##    ##
##       ##     ## ##        ##   ##     ##    ##          ##  ##  ## ##
##       ########  ######   ##     ##    ##    ######      ##  ##  ##  ######
##       ##   ##   ##       #########    ##    ##          ##  ##  ##       ##
##    ## ##    ##  ##       ##     ##    ##    ##          ##  ##  ## ##    ##
 ######  ##     ## ######## ##     ##    ##    ########     ###  ###   ######

class EFTWorkspace(object):
    '''
    xxx
    '''

 # #    # # #####
 # ##   # #   #
 # # #  # #   #
 # #  # # #   #
 # #   ## #   #
 # #    # #   #

    def __init__(self, opts):
        '''
        Initialize
        '''

        self.logger = logging.getLogger(__name__)

        self.wcs = opts['wcs']
        self.wc_ranges = opts['wc_ranges']
        self.allVars = [] #Trick -- need to save objects purely for python memory management #If get segfault when importing object, you may have forgotten to store some object in memory #See: https://root-forum.cern.ch/t/roofit-in-python-segmentation-faults-all-over-the-place/27955/3
        self.fit_file = fits_file # File containing WC parameterizations of each process+bin (if non-empty)!
        self.input_dc = '' #Path of input text datacard
        self.DC = None #Datacard object
        self.shape_systs = [] #List of shape systs
        self.dict_shapesyst_chan_proc = {} #Nested dict mapping syst/chan/proc/value
        self.dict_rrv_shapeSyst = {} #Dict whose keys are shapeSyst names, and values are the corresponding RooRealVars (nuisance params) --> So that we can always reuse the same RRV where needed

        return


# //--------------------------------------------
# //--------------------------------------------

    #-- Gives invalid fit results...
    def Get_PolyCoeffs_FromUpDownVariations(self, binc_nom, binc_up, binc_down):
        '''
        Combine can't interpolate from up/down variations automatically as for TH1s --> Need to parameterize the bins of the RPHs also on all shape nuisances.
        A good approximation to enforce this interpolation 'manually' is to parameterize the bins on some polynomial of the form (1+a*@0+b*@0*@0), where a and b can be determined from the +-1 sigma templates.

        This function aims to determine the values of the a and b coefficients based on the 'nom/up/down relative values' of the bins.
        In that notation, if e.g. bin_nom=2, bin_up=2.4, bin_down=1.6 --> nom=1, up=1.2, down=0.8 (<-> relative).
        The system of equations is then:
        (1 +a*1  +b+1*1)      = up    <=> a+b=1-up   ==U
        (1 +a*-1 +b*-1*-1)    = down  <=> a-b=1-down ==D
        which becomes (sum both equations):
        a = (U+D)/2
        b = (U - a)
        '''

        U = binc_up / binc_nom
        D = binc_down / binc_nom

        a = (U + D) / 2.
        b = (U - a)

        # print('binc_nom', binc_nom, 'binc_up', binc_up, 'binc_down', binc_down)
        # print('U', U, 'D', D)
        # print('--> a', a, ' /  b', b)

        return a, b


# //--------------------------------------------
# //--------------------------------------------

    def Load_Fits(self, verb=0):
        '''
        Load all the fit coefficients in memory (from numpy dict. previously created with 'DumEFTParameterization.py' code)
        '''

        if not os.path.exists(self.fit_file):
            print(colors.bg.red + "File {} does not exist!".format(hist_file) + colors.reset)
            return

        #-- Load fit coefficients for all bins/EFT signals
        self.fits = np.load(self.fit_file)[()] #Create dict object

        print(colors.fg.orange + 'Reading EFT parameterization dict: ' + colors.reset + str(self.fit_file))

        if verb>=2: print('-- Contents of [self.fits] :', self.fits)

        return


# //--------------------------------------------
# //--------------------------------------------


    def Parse_Input_Datacard(self, datacard_path, verb=0):
        '''
        Use combine methods to automatically parse the input datacard, to extract info about processes, channels, etc.
        '''

        #-- Need to defile_ine all the options expected by parseCard(), but don't care about their values
        parser = OptionParser()
        addDatacardParserOptions(parser)
        (options, args) = parser.parse_args(['dummy']) #Use dummy arg to avoid reading the command-line args

        file = open(datacard_path, "r")
        self.DC = parseCard(file, options) #Parse datacard textfile

        #-- Build list of names of all shape systematics (will need to loop over the corresponding shifted histos)
        #Ex: DC.systs = 	[('Lumi16', False, 'lnN', [], {'Zpt_SRtZq_2016': {'NPL': 0.0, 'tX': 1.022, 'PrivMC_tZq': 1.022, 'VVV': 1.022, 'XG': 1.022, 'PrivMC_tWZ': 1.022, 'WZ': 1.022, 'PrivMC_ttZ': 1.022}}), ('PU', False, 'shape', [], {'Zpt_SRtZq_2016': {'NPL': 1.0, 'tX': 1.0, 'PrivMC_tZq': 1.0, 'VVV': 1.0, 'XG': 1.0, 'PrivMC_tWZ': 1.0, 'WZ': 1.0, 'PrivMC_ttZ': 1.0}})] # <type 'list'>
        for isyst, syst in enumerate(self.DC.systs):
            if syst[2] == 'shape': #If shape systematic found
                self.shape_systs.append(syst[0]) #Append name of variation
                # self.shape_systs.append('{}Up'.format(syst[0])) #Append name of up variation
                # self.shape_systs.append('{}Down'.format(syst[0])) #Append name of down variation

        if verb>0:
            print(colors.fg.lightblue + '\nhasShapes --->' + colors.reset + str(self.DC.hasShapes))
            print(colors.fg.lightblue + 'bins --->' + colors.reset + str(self.DC.bins))
            print(colors.fg.lightblue + 'processes --->' + colors.reset + str(self.DC.processes))
            print(colors.fg.lightblue + 'isSignal --->' + colors.reset + str(self.DC.isSignal))
            print(colors.fg.lightblue + 'shapeMap --->' + colors.reset + str(self.DC.shapeMap))
            print(colors.fg.lightblue + 'self.shape_systs' + colors.reset + self.shape_systs)

            # print(colors.fg.lightblue + 'systs --->' + colors.reset + str(self.DC.systs)) #List of systs #For each element: [0]=name, [1]=bool(?), [2]=type, [3]=list(?), [4]=dict whose keys are channels, and values are dicts (whose keys are processes, and values are the corresponding syst values)
            # print(colors.fg.lightblue + 'obs --->' + colors.reset + str(self.DC.obs))
            # print(colors.fg.lightblue + 'signals --->' + colors.reset + str(self.DC.signals))
            # print(colors.fg.lightblue + 'keyline --->' + colors.reset + str(self.DC.keyline))
            # print(colors.fg.lightblue + 'exp --->' + colors.reset + str(self.DC.exp))
            # print(colors.fg.lightblue + 'flatParamNuisances --->' + colors.reset + str(self.DC.flatParamNuisances))
            # print(colors.fg.lightblue + 'rateParams --->' + colors.reset + str(self.DC.rateParams))
            # print(colors.fg.lightblue + 'systematicsShapeMap --->' + colors.reset + str(self.DC.systematicsShapeMap))
            # print(colors.fg.lightblue + 'binParFlags --->' + colors.reset + str(self.DC.binParFlags))
            # print(colors.fg.lightblue + 'groups --->' + colors.reset + str(self.DC.groups))

        return


# //--------------------------------------------
# //--------------------------------------------

    def Fill_Dict_ShapeSyst(self, verb=0):
        '''
        Create a simplified nest dictionnary as such: {syst1:{chan1:{proc1:val1,proc2:val2}, chan2:{proc1:val1,proc2:val2}}, syst2...}
        All shapesyst/chan/proc combinations are included. Can thus access the value of a shape syst for a given chan/proc as: self.dict_shapesyst_chan_proc[syst][chan][proc]; val=0 <-> the syst does not apply in that case
        '''

        #Ex: DC.systs = [('Lumi16', False, 'lnN', [], {'Zpt_SRtZq_2016': {'NPL': 0.0, 'tX': 1.022, 'PrivMC_tZq': 1.022, 'VVV': 1.022, 'XG': 1.022, 'PrivMC_tWZ': 1.022, 'WZ': 1.022, 'PrivMC_ttZ': 1.022}}), ('PU', False, 'shape', [], {'Zpt_SRtZq_2016': {'NPL': 1.0, 'tX': 1.0, 'PrivMC_tZq': 1.0, 'VVV': 1.0, 'XG': 1.0, 'PrivMC_tWZ': 1.0, 'WZ': 1.0, 'PrivMC_ttZ': 1.0}})] # <type 'list'>

        for isyst, syst in enumerate(self.shape_systs): #Loop on shape systematics

            self.dict_shapesyst_chan_proc[syst] = {}

            for ichan, chan in enumerate(self.DC.bins): #Loop on channels ('bins')

                if 'CR' in chan: continue

                self.dict_shapesyst_chan_proc[syst][chan] = {}

                for iproc, proc in enumerate(self.DC.processes): #Loop on MC physics processes (samples)

                    self.dict_shapesyst_chan_proc[syst][chan][proc] = 0. #Default = 0 <-> this shape syst does not apply to this chan/proc

                    syst_appliesToChanProc = False #True <-> apply this syst to this chan/proc

                    for systtmp in self.DC.systs:
                        # print('systtmp', systtmp)
                        if systtmp[0]==syst and systtmp[2]=='shape' and not syst_appliesToChanProc:
                            # print('systtmp[0]', systtmp[0])
                            for chantmp in systtmp[4]: #Keys in dict of dict (<-> channels)
                                # print('chantmp', chantmp)
                                if chan == chantmp and not syst_appliesToChanProc: #Channel found for this syst
                                    for procs_systtmp in systtmp[4][chantmp]: #Keys in dict (<-> processes)
                                        # print('procs_systtmp', procs_systtmp)
                                        # print('val systtmp[4][chantmp][procs_systtmp]', systtmp[4][chantmp][procs_systtmp])
                                        if proc == procs_systtmp and systtmp[4][chantmp][procs_systtmp] != 0. and not syst_appliesToChanProc:
                                            self.dict_shapesyst_chan_proc[syst][chan][proc] = systtmp[4][chantmp][procs_systtmp]
                                            syst_appliesToChanProc = True
                                            break
                                            # print(colors.fb.red + 'self.dict_shapesyst_chan_proc[' + syst + '][' + chan + '][' + proc + '] = ' + self.dict_shapesyst_chan_proc[syst][chan][proc] + colors.reset)

        #-- Order dict keys alphabetically, for easier reading
        self.dict_shapesyst_chan_proc = OrderedDict(sorted(self.dict_shapesyst_chan_proc.items()))

        # print('self.dict_shapesyst_chan_proc', self.dict_shapesyst_chan_proc)
        return


# //--------------------------------------------
# //--------------------------------------------


 # #    # #####   ####  #####  #####
 # ##  ## #    # #    # #    #   #
 # # ## # #    # #    # #    #   #
 # #    # #####  #    # #####    #
 # #    # #      #    # #   #    #
 # #    # #       ####  #    #   #

 #####  #####  #    #
 #    # #    # #    #
 #    # #    # ######
 #####  #####  #    #
 #   #  #      #    #
 #    # #      #    #

                 #     #  #####
 #####  ####     #  #  # #     #
   #   #    #    #  #  # #
   #   #    #    #  #  #  #####
   #   #    #    #  #  #       #
   #   #    #    #  #  # #     #
   #    ####      ## ##   #####

    def Import_RooParametricHists_ToWorkspace(self, datacard_path, template_file, template_file_otherRegions, read_UpDownShiftedHist_forRPH=False, verb=0):
        '''
        Create a rootfile containing a temporary workspace including all the necessary RooParametricHist objects (1 per EFT signal process and channel).

        Logic: for each channel and relevant process, read the corresponding EFT fit coefficients and define a RooFormulaVar describing the EFT parameterization of each histogram bin --> Create and store RooParametricHist
        '''

        #-- Open output file, create workspace object
        file_out = TFile("WS.root","RECREATE")
        wspace = RooWorkspace("w","w")

        #-- Create RooRealVar for each WC #Store the list of WCs contributing to the formula --> Will provide this same list of args for all EFT parameterizations formulas
        ral_wcs = RooArgList()
        for wc in self.wcs:
            rrv_wc = RooRealVar(wc, wc, 0., self.wc_ranges[wc][0], self.wc_ranges[wc][1])
            # print('rrv_wc', rrv_wc)
            getattr(wspace,'import')(rrv_wc, RooFit.RecycleConflictNodes(), RooFit.Silence())
            self.allVars.append(rrv_wc) #Trick (memory-management)
            ral_wcs.add(rrv_wc)
        self.allVars.append(ral_wcs) #Trick (memory-management)

        #-- For each shapeSyst, create a corresponding RRV (nuisance par) once, which we could reuse later as needed
        for isyst, syst in enumerate(self.shape_systs): #Loop on shape systematics
            rrv_shapeSyst = RooRealVar(syst, syst, 0., -4., 4.) #Create a dedicated RooRealVar for each shape syst (so that it will be floated coherently between regions/processes...)
            # rrv_shapeSyst = RooRealVar(syst, syst, 0.) #Create a dedicated RooRealVar for each shape syst (so that it will be floated coherently between regions/processes...)
            self.dict_rrv_shapeSyst[syst] = rrv_shapeSyst
            getattr(wspace,'import')(rrv_shapeSyst, RooFit.RecycleConflictNodes(), RooFit.Silence()) # Add 'RooFit.Silence()' arrg to remove printout

        include_MCstat_nuisances = False #Default: don't consider MC stat errors
        if self.DC.binParFlags: include_MCstat_nuisances = True #If this dict is not empty, it means we want to include autoMCStats


 #       ####   ####  #####
 #      #    # #    # #    #
 #      #    # #    # #    #
 #      #    # #    # #####
 #      #    # #    # #
 ######  ####   ####  #

#-- Main loop on chan / proc / nuisances / ...
# //--------------------------------------------

        for ichan, chan in enumerate(self.DC.bins): #Loop on channels ('bins')

            if 'CR' in chan: continue #Only need to create RooParametricHist/RooDataHist objects in SRs, i.e. in any region including RooParametricHists; in CRs, can use only regular TH1s as usual

            # print('-- chan', chan)

            #-- Get prefix for relative path based on datacard path #Obsolete -- now provide filepath directly as command-line arg
            # split_datacard_path = datacard_path.split('/')
            # prefix_path = ''
            # for s in split_datacard_path:
            #     if '.txt' not in s: prefix_path+= s + '/'
            # print(prefix_path)
            # rootfile_path = prefix_path + self.DC.shapeMap[chan]['*'][0] #Full template path = prefix + relative template path
            # rootfile_path = '../templates/Templates_Zpt_EFT2_Run2.root'
            # print(rootfile_path)

            rootfile_path = template_file
            if 'SRttZ4l' in chan: rootfile_path = template_file_otherRegions #SRttZ4l templates are stored in a separate file (which does not depend on NN-EFT trainings)

            file_in = TFile.Open(rootfile_path) #Open input template file

            rrv_xaxis = RooRealVar() #RRV corresponding to the current x-axis variable (shared by different processes/systs/... in same region)
            rrv_xaxis_alreadyDefined = False #Will define the RRV when reading first histo in this channel #Then recycle the same object, since all other histos should reuse the same

            for iproc, proc in enumerate(self.DC.processes + ['data_obs']): #Loop on physics processes (samples) #Must also consider data hist

                # print('-- proc', proc)

                #-- NOMINAL
                hname = '{}__{}'.format(chan,proc) #Hard-coded histo naming convention
                if hname not in file_in.GetListOfKeys():
                    print(colors.fg.red + 'Histo ' + hname + ' not found in rootfile ' + rootfile_path + ' ! Did you store the correct EFT parameterization (with DumEFTParameterization.py) ? Abort !' + colors.reset)
                    exit(1)

                h_tmp = file_in.Get(hname)
                # print('h_tmp.GetNbinsX()', h_tmp.GetNbinsX())

                #-- Define RRV for x-axis variable (common to different processes in same channel -- define only first time it's used)
                if not rrv_xaxis_alreadyDefined:
                    rrv_xaxis = RooRealVar("rrv_{}".format(chan), "rrv_{}".format(chan), h_tmp.GetXaxis().GetXmin(), h_tmp.GetXaxis().GetXmax())
                    self.allVars.append(rrv_xaxis) #Trick (memory-management)
                    rrv_xaxis_alreadyDefined = True

                list_bins = RooArgList() #Will store the complete RooFormulaVar describing the full bin parameterization for each histo bin

# EFT signals
# //--------------------------------------------
                if proc in self.DC.isSignal and self.DC.isSignal[proc]: #Only need to create RooParametricHist for 'EFT signals'

                    for ibin in range(1, h_tmp.GetNbinsX()+1): #Loop on histogram bins

                        bin_name = 'bin{}_{}'.format(ibin,chan) if 'countExp' not in chan else chan #Ex: bin8_mTW_SRother_2017, countExp_SRttZ4l_2016

                        bin_content = max(1e-5, h_tmp.GetBinContent(ibin)) #Make sure bin content is never negative
                        # print('ibin', ibin, 'bin_content', bin_content)


 ###### ###### #####    #####    ##   #####    ##   #    #
 #      #        #      #    #  #  #  #    #  #  #  ##  ##
 #####  #####    #      #    # #    # #    # #    # # ## #
 #      #        #      #####  ###### #####  ###### #    # ###
 #      #        #      #      #    # #   #  #    # #    # ###
 ###### #        #      #      #    # #    # #    # #    # ###

                        eft_analytical_expr = '1' #Default bin parameterization = '1' (<-> 1*SM)

                        ral_tmp = RooArgList() #By default, EFT expression will depend on 0 parameters (e.g. if parameterization not found in dict)

                        if (proc,bin_name) in self.fits: #EFT parameterization found

                            ral_tmp = ral_wcs #Parameterization found --> Make expression depend on all WCs

                            fitCoeffs_procchan = self.fits[proc,bin_name] #Read all fit coefficients for this chan/proc
                            # print('fitCoeffs_procchan', fitCoeffs_procchan)

                            for icoeff, coeff1 in enumerate(self.wcs): #First loop on WCs
                                # print('icoeff', icoeff, 'coeff1', coeff1)

                                #-- Get fit coefficient for interference term (SM-EFT)
                                coeff = 0.
                                if ('SM',coeff1) in fitCoeffs_procchan: coeff = fitCoeffs_procchan[('SM',coeff1)]
                                elif (coeff1,'SM') in fitCoeffs_procchan: coeff = fitCoeffs_procchan[(coeff1,'SM')]
                                else:
                                    print(colors.fg.red + 'Key ' + ('SM',coeff1) + ' not found in  fitCoeffs_procchan ! Abort !' + colors.reset)
                                    exit(1)

                                #-- Update full formula for the current histo bin
                                if abs(coeff) < min_threshold_EFTcoeff: coeff=0 #If coefficient is below min threshold, neglect it, but still write the corresponding (null) term because we always want to give the full list of WCs as args (so even if a WC is actually neglected, we include it in the rfv)
                                eft_analytical_expr+= '+{}*@{}'.format(coeff,icoeff) if coeff>=0 else '{}*@{}'.format(coeff,icoeff) #NB: Negative coefficient already contains the minus sign

                                for jcoeff, coeff2 in enumerate(self.wcs): #Second loop on WCs
                                    # print('jcoeff', jcoeff, 'coeff2', coeff2)

                                    if jcoeff < icoeff: continue #<-> if term (wc1,wc2) already included, don't include equivalent term (wc2,wc1)

                                    #-- Get fit coefficient for quadratic term (EFT-EFT)
                                    coeff = 0.
                                    if (coeff1,coeff2) in fitCoeffs_procchan: coeff = fitCoeffs_procchan[(coeff1,coeff2)]
                                    elif (coeff2,coeff1) in fitCoeffs_procchan: coeff = fitCoeffs_procchan[(coeff2,coeff1)]
                                    else:
                                        print(colors.fg.red + 'Key ' + (coeff1,coeff2) + ' not found in  fitCoeffs_procchan ! Abort !' + colors.reset)
                                        exit(1)

                                    #-- Update full formula for the current histo bin (if above given threshold)
                                    if abs(coeff) < min_threshold_EFTcoeff: coeff=0 #If coefficient is below min threshold, neglect it, but still write the corresponding (null) term because we always want to give the full list of WCs as args (so even if a WC is actually neglected, we include it in the rfv)
                                    eft_analytical_expr+= '+{}*@{}*@{}'.format(coeff,icoeff,jcoeff) if coeff>=0 else '{}*@{}*@{}'.format(coeff,icoeff,jcoeff) #NB: Negative coefficient already contains the minus sign


                        else: #EFT parameterization not found (this is not necessarily an error: it will happen e.g. if a given EFT signal does not populate a given bin at all) #Just define a constant expression '1'
                            print(colors.fg.red + 'Key ' + '({},{})'.format(proc,bin_name) + ' not found in parameterization dict ' + self.fit_file + ' ! Abort !' + colors.reset)

                        #-- Define the formulas defining the bin content
                        # print('eft_analytical_expr', eft_analytical_expr)
                        eft_expr_bin_name = 'eft_expr_{}_{}'.format(bin_name,proc)
                        eft_expr_bin = RooFormulaVar('eft_expr_{}_{}'.format(bin_name,proc), 'EFT parameterization in {}_{}'.format(bin_name,proc), eft_analytical_expr, ral_tmp)
                        self.allVars.append(eft_expr_bin) #Trick (memory-management)


  ####  #    #   ##   #####  ######     ####  #   #  ####  #####
 #      #    #  #  #  #    # #         #       # #  #        #
  ####  ###### #    # #    # #####      ####    #    ####    #
      # #    # ###### #####  #              #   #        #   #
 #    # #    # #    # #      #         #    #   #   #    #   #
  ####  #    # #    # #      ######     ####    #    ####    #

                        list_shapeSyst_expr_bin = [] #Contains the RooFormulaVars describing the shape syst parameterization in the bin

                        #-- Shape systs (also convert TH1s to RDHs)
                        for isyst, syst in enumerate(self.shape_systs): #Loop on shape systematics
                            # print('-- syst', syst)

                            if self.dict_shapesyst_chan_proc[syst][chan][proc] == 0: continue

                            rrv_shapeSyst = self.dict_rrv_shapeSyst[syst]

                            #-- Get shifted histos (hardcoded naming conventions)
                            hnameup = '{}__{}__{}Up'.format(chan,proc,syst)
                            hnamedown = '{}__{}__{}Down'.format(chan,proc,syst)
                            if hnameup not in file_in.GetListOfKeys() or hnamedown not in file_in.GetListOfKeys():
                                print(colors.fg.red + 'Histo [' + hnameup + '] or [' + hnamedown + '] not found in rootfile ' + rootfile_path + ' ! Abort !' + colors.reset)
                                exit(1)

                            hup = file_in.Get(hnameup)
                            hdown = file_in.Get(hnamedown)
                            # print('hup.GetNbinsX()', hup.GetNbinsX())
                            # print('hdown.GetNbinsX()', hdown.GetNbinsX())

                            #For each relevant syst, store the corresponding RooFormulaVar (and associated RooRealVar) for the current bin; 1 RFV <-> 1 RRV
                            if read_UpDownShiftedHist_forRPH == False:

                                # a, b = self.Get_PolyCoeffs_FromUpDownVariations(bin_content, hup.GetBinContent(ibin), hdown.GetBinContent(ibin))
                                # rfv_shapesyst_name = 'rfv_shapesyst_{}_{}_{}'.format(bin_name,proc,syst)
                                # rfv_shapesyst = RooFormulaVar(rfv_shapesyst_name, 'Parameterization of shapeSyst {} in bin {}'.format(syst,bin_name), "1+{}*@0+{}*@0*@0".format(a,b), RooArgList(rrv_shapeSyst)) #NB: having '+-xxx' in expression is not a problem ?
                                # list_shapeSyst_expr_bin.append(rfv_shapesyst)
                                # self.allVars.append(rfv_shapesyst) #Useless ?

                                #-- See: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part2/settinguptheanalysis/#rates-for-shape-analysis
                                #-- See: https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/102x/src/AsymQuad.cc#L94-L126
                                rfv_shapesyst_name = 'rfv_shapesyst_{}_{}_{}'.format(bin_name,proc,syst)
                                shapesyst_expr = "( @0<=-1 )? ( @0*({0}-{2})/{0} ):( ( @0>=1 )? ( @0*({1}-{0})/{0} ):( (@0/2.)*(({1}-{2})+(0.125*@0*(TMath::Power(@0,2)*(3.*TMath::Power(@0,2)-1.)+15.))*({1}+{2}-2.*{0}))/{0} ) )".format(bin_content,hup.GetBinContent(ibin),hdown.GetBinContent(ibin))
                                rfv_shapesyst = RooFormulaVar(rfv_shapesyst_name, 'Parameterization of shapeSyst {} in bin {}'.format(syst,bin_name), "1+({})".format(shapesyst_expr), RooArgList(rrv_shapeSyst))
                                # print('rfv_shapesyst.evaluate', rfv_shapesyst.evaluate()) #Evaluate RDV at initial RRV value
                                list_shapeSyst_expr_bin.append(rfv_shapesyst)
                                self.allVars.append(rfv_shapesyst) #Useless ?

                            elif ibin==1: #Convert up/down (shape syst) histograms for RPHs, just like for RDHs (uses PR #634) #Do that only once !

                                rdh = RooDataHist('rdh_{}'.format(hnameup), 'rdh_{}'.format(hnameup), RooArgList(rrv_xaxis), hup)
                                #print(colors.fg.lightblue + 'Import rdh_{}'.format(hnameup) + colors.reset)
                                getattr(wspace,'import')(rdh)
                                rdh = RooDataHist('rdh_{}'.format(hnamedown), 'rdh_{}'.format(hnamedown), RooArgList(rrv_xaxis), hdown)
                                #print(colors.fg.lightblue + 'Import rdh_{}'.format(hnamedown) + colors.reset)
                                getattr(wspace,'import')(rdh) #getattr(wspace,'import', [RooFit.RecycleConflictNodes(), RooFit.Silence()])(rdh)?

                        #-- End shapeSyst loop


 #    #  ####      ####  #####   ##   #####
 ##  ## #    #    #        #    #  #    #
 # ## # #          ####    #   #    #   #
 #    # #              #   #   ######   #
 #    # #    #    #    #   #   #    #   #
 #    #  ####      ####    #   #    #   #

                        include_MCstatError_rph = False
                        if include_MCstat_nuisances:

                            bin_error = h_tmp.GetBinError(ibin)
                            rel_bin_error = abs(bin_error/bin_content)
                            if rel_bin_error > min_threshold_MCstat:
                                include_MCstatError_rph = True

                                #-- Must create a dedicated RooRealVar (nuisance parameter) for each chan/proc/histbin, and tie it to the corresponding RooRealVar (derived from the corresponding MC stat uncertainty)
                                rrv_stat_name = 'MCstat_{}_{}'.format(bin_name,proc)
                                rrv_stat = RooRealVar(rrv_stat_name, rrv_stat_name, 0., -4., 4.)
                                # rrv_stat = RooRealVar(rrv_stat_name, rrv_stat_name, 0.)
                                self.allVars.append(rrv_stat) #Trick (memory-management)

                                # print(rrv_stat_name, ' --> include_MCstatError_rph=True !')

                                rfv_stat_name = 'rfv_{}'.format(rrv_stat_name)
                                rfv_stat = RooFormulaVar(rfv_stat_name, 'Parameterization of MC stat uncertainty in bin {}'.format(bin_name), "TMath::Power(1+{}, @0)".format(rel_bin_error), RooArgList(rrv_stat))
                                self.allVars.append(rfv_stat) #Trick (memory-management)

                                #-- Need to import these rrv manually, because in the end the RPH will take the rfv as argument (not the rrv itself)
                                getattr(wspace,'import')(rrv_stat, RooFit.RecycleConflictNodes(), RooFit.Silence())


 #####  # #    #    #####  ###### ###### # #    # # ##### #  ####  #    #
 #    # # ##   #    #    # #      #      # ##   # #   #   # #    # ##   #
 #####  # # #  #    #    # #####  #####  # # #  # #   #   # #    # # #  #
 #    # # #  # #    #    # #      #      # #  # # #   #   # #    # #  # #
 #    # # #   ##    #    # #      #      # #   ## #   #   # #    # #   ##
 #####  # #    #    #####  ###### #      # #    # #   #   #  ####  #    #

                        #-- Now that we have defined the EFT / shapeSyst / MCstat parameterizations, we can tie everything together to define the full parameterization of the current bin
                        counter_args = 0 #Must increment the args @XX coherently
                        total_list_args_bin = RooArgList()

                        #-- EFT
                        total_bin_expr = '{}*@{}'.format(bin_content,counter_args) #Default = bin_content_SM * EFT_param
                        total_list_args_bin.add(eft_expr_bin) #Minimal argument: rfv for EFT parameterization
                        counter_args+= 1

                        #-- Shape
                        if read_UpDownShiftedHist_forRPH == False:
                            for shape_expr in list_shapeSyst_expr_bin:
                                total_bin_expr+= '*@{}'.format(counter_args)
                                total_list_args_bin.add(shape_expr)
                                counter_args+= 1

                        #-- MC stat
                        if include_MCstatError_rph:
                            total_bin_expr+= '*@{}'.format(counter_args)
                            total_list_args_bin.add(rfv_stat)
                            counter_args+= 1

                        rfv_bin_name = '{}_{}'.format(bin_name,proc)
                        rfv_bin = RooFormulaVar(rfv_bin_name, 'Total parameterized yield in {}'.format(rfv_bin_name), total_bin_expr, RooArgList(total_list_args_bin))
                        list_bins.add(rfv_bin)
                        self.allVars.append(rfv_bin) #Trick (memory-management)

                    #-- End histo bins loop


 ######  ######  #     #
 #     # #     # #     #
 #     # #     # #     #
 ######  ######  #######
 #   #   #       #     #
 #    #  #       #     #
 #     # #       #     #

                    #-- Create RPH object
                    rph_name = "rph_{}__{}".format(chan,proc) #Hard-coded naming convention
                    rph = RooParametricHist(rph_name, rph_name, rrv_xaxis, list_bins, h_tmp)
                    self.allVars.append(rph) #Trick (memory-management)
                    self.allVars.append(list_bins) #Trick (memory-management)

                    #-- "Always include a _norm term which should be the sum of the yields (thats how combine likes to play with pdfs)"
                    rphnorm_name = "rph_{}__{}_norm".format(chan,proc) #Hard-coded naming convention
                    rph_norm = RooAddition(rphnorm_name, rphnorm_name, list_bins)
                    self.allVars.append(rph_norm) #Trick (memory-management)

                    #-- Import objects to workspace
                    # print(colors.fg.lightblue + 'Import {}'.format(rph_name_tmp) + colors.reset)
                    getattr(wspace,'import')(rph, RooFit.RecycleConflictNodes(), RooFit.Silence())
                    getattr(wspace,'import')(rph_norm, RooFit.RecycleConflictNodes(), RooFit.Silence())
                    # print('... done !')
# //--------------------------------------------

#Others (data, bkg)
# //--------------------------------------------
                else: #Not an EFT signal <-> still need to convert TH1 to RooDataHist (RDH) so that the same RRV is used to describe the x-axis of the different processes in the given region


 #####  #####  #    #
 #    # #    # #    #
 #    # #    # ######
 #####  #    # #    #
 #   #  #    # #    #
 #    # #####  #    #

                    #-- NOMINAL
                    rdh_name = "rdh_{}__{}".format(chan,proc) #Hard-coded naming convention
                    rdh = RooDataHist(rdh_name, rdh_name, RooArgList(rrv_xaxis), h_tmp) #Can ignore warnings of type 'INFO:DataHandling -- RooDataHist::adjustBinning(rdh_Zpt_SRtZq_2016__WZ): fit range of variable rrv_Zpt_SRtZq_2016 expanded to nearest bin boundaries: [0,350] --> [0,350]', as range does not actually change !
                    # print(colors.fg.lightblue + 'Import {}'.format(rdh_name) + colors.reset)
                    getattr(wspace,'import')(rdh) #commands 'RooFit.RecycleConflictNodes()' and 'RooFit.Silence()' only available for datasets starting from ROOT6.18 ?

                    #-- SHAPE SYSTS (convert corresponding TH1s to RDHs)
                    for isyst, syst in enumerate(self.shape_systs): #Loop on shape systematics
                        # print('-- syst', syst)

                        if proc=='data_obs' or ('NPL' in proc and not syst.startswith('FR')) or (syst.startswith('FR') and 'NPL' not in proc): continue
                        if self.dict_shapesyst_chan_proc[syst][chan][proc] == 0: continue

                        #-- For all shape systs which apply, need to store objects both for up/down variations
                        for ivariation, variation in enumerate(['Up', 'Down']):
                            systname_tmp = '{}{}'.format(syst,variation)

                            hname = '{}__{}__{}'.format(chan,proc,systname_tmp) #Hard-coded histo naming convention
                            if hname not in file_in.GetListOfKeys():
                                print(colors.fg.red + 'Histo ' + hname + ' not found in rootfile ' + rootfile_path + ' ! Abort !' + colors.reset)
                                exit(1)

                            h_tmp = file_in.Get(hname)
                            # print('h_tmp.GetNbinsX()', h_tmp.GetNbinsX())

                            rdh_name = "rdh_{}__{}__{}".format(chan,proc,systname_tmp) #Hard-coded naming convention
                            rdh = RooDataHist(rdh_name, rdh_name, RooArgList(rrv_xaxis), h_tmp)
                            # print(colors.fg.lightblue + 'Import {}'.format(rdh_name) + colors.reset)
                            getattr(wspace,'import')(rdh)

                    #-- End shape syst loop

                    #-- MC STAT (create new RDHs where a given histobin is shifted +-stat, and store it) --> MC stat will need to be treated as shape systematics for RDHs !
                    if include_MCstat_nuisances:
                        for ibin in range(1, h_tmp.GetNbinsX()+1): #Loop on histogram bins
                            if proc == 'data_obs': break #No MC stat uncert for data

                            bin_name = 'bin{}_{}'.format(ibin,chan) if 'countExp' not in chan else chan #Ex: bin8_mTW_SRother_2017, countExp_SRttZ4l_2016

                            #-- New RDHs with bin shifted with +- MC stat error (make sure the result is still >=0)
                            h_statUp_ibin = h_tmp.Clone(); h_statDown_ibin = h_tmp.Clone()
                            new_bin_content = max(1e-5, h_statUp_ibin.GetBinContent(ibin)+h_statUp_ibin.GetBinError(ibin))
                            h_statUp_ibin.SetBinContent(ibin, new_bin_content)
                            new_bin_content = max(1e-5, h_statDown_ibin.GetBinContent(ibin)-h_statDown_ibin.GetBinError(ibin))
                            h_statDown_ibin.SetBinContent(ibin, new_bin_content)

                            rrv_stat_name = 'MCstat_{}_{}'.format(bin_name,proc)
                            rdh_name = "rdh_{}__{}__{}Up".format(chan,proc,rrv_stat_name) #Hardcoded naming convention: 'rdh_binvar__proc__systUp/Down'
                            rdh_statUp = RooDataHist(rdh_name, rdh_name, RooArgList(rrv_xaxis), h_statUp_ibin)
                            rdh_name = "rdh_{}__{}__{}Down".format(chan,proc,rrv_stat_name) #Hardcoded naming convention: 'rdh_binvar__proc__systUp/Down'
                            rdh_statDown = RooDataHist(rdh_name, rdh_name, RooArgList(rrv_xaxis), h_statDown_ibin)
                            getattr(wspace,'import')(rdh_statUp)
                            getattr(wspace,'import')(rdh_statDown)

                #-- if/else EFT signal
# //--------------------------------------------
            #-- End proc loop
        #-- End chan loop

        #-- Write workspace object to output file
        file_out.cd()
        wspace.Write()

        #-- Close output file
        file_out.Close()

        return


# //--------------------------------------------
# //--------------------------------------------









##     ##    ###    #### ##    ##
###   ###   ## ##    ##  ###   ##
#### ####  ##   ##   ##  ####  ##
## ### ## ##     ##  ##  ## ## ##
##     ## #########  ##  ##  ####
##     ## ##     ##  ##  ##   ###
##     ## ##     ## #### ##    ##


# //--------------------------------------------
# //--------------------------------------------

if __name__ == "__main__":

    PrintBanner()

    read_UpDownShiftedHist_forRPH = True #True <-> convert up/down (shape syst) histograms for RPHs, just like for RDHs (uses PR #634); else: encode the shapeSyst effects directly into bin parameterizations

# User options -- Default values
# //--------------------------------------------
    datacard_path = './datacard.txt'
    verb = 0
    fits_file = './EFT_Parameterization.npy'
    template_file = '../templates/Templates_Zpt_EFT2_Run2.root'
    template_file_otherRegions = '../templates/Templates_otherRegions_Run2.root' #Also need path to 'otherRegions' rootfile containing SRttZ4l templates (and CRs)

# Set up the command line arguments
# //--------------------------------------------
    parser = argparse.ArgumentParser(description='Perform SM and EFT fits using custom Physics Model')
    parser.add_argument("-d", metavar="datacard path", help="Path to the txt datacard (to create RooWorkspace)", required=True)
    parser.add_argument("-f", metavar="parameterization file path", help="Path to the EFT parameterization dict", required=True)
    parser.add_argument("-t", metavar="template file path", help="Path to the rootfile containing templates", required=True)
    parser.add_argument("-o", metavar="otherRegions template file path", help="Path to the otherRegions rootfile containing SRttZ4l templates")
    parser.add_argument("-v", metavar="Combine verbosity level", help="Set combine output verbosity")

    args = parser.parse_args()

    if args.d: datacard_path = args.d
    if args.v: verb = int(args.v)
    if args.f: fits_file = args.f
    if args.t: template_file = args.t
    if args.o: template_file_otherRegions = args.o

# Create class object
# //--------------------------------------------
    eftws = EFTWorkspace(opts)

    eftws.Load_Fits(verb)
    eftws.Parse_Input_Datacard(datacard_path, verb)
    eftws.Fill_Dict_ShapeSyst(verb)
    eftws.Import_RooParametricHists_ToWorkspace(datacard_path, template_file, template_file_otherRegions, read_UpDownShiftedHist_forRPH, verb)

# //--------------------------------------------
