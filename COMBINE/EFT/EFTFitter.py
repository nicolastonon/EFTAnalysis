#-- Perform SM and EFT fits using custom Physics Model
# Adapted from: https://github.com/cms-govner/EFTFit
# Batch modes supported are: CRAB3 ('crab') and Condor ('condor')

'''
#'text2workspace will convert the datacard into a pdf which summaries the analysis'
- '--just-check-physics-model'
- '--dump-datacard'
- '--PO' Pass a given option to the physics model
- '-P' Pass Physics Model
- '-v 9' max verbose
'''

import os
import stat
import sys
import logging
import subprocess as sp
import ROOT
import itertools
import glob
import getpass
import array
from Utils.ColoredPrintout import colors
from collections import defaultdict
import getopt # command line parser
import argparse
from settings import opts #Custom dictionnary of settings
import shutil
import CombineHarvester.CombineTools.plotting as plot #Combine plotting utils
from EFTPlotter import BuildScan


def PrintBanner():
    #Print some info when starting this code

    print('\n' + colors.bg.orange + '                                           ' + colors.reset)
    print(colors.fg.orange + '-- EFTFitter --' + colors.reset + '\n')
    print('* NB1: make sure you are using the relevant template file')
    print('* NB2: make sure you have extracted the relevant EFT parameterizations')
    print(colors.bg.orange + '                                           ' + colors.reset + '\n')

    return


######## ######## ######## ######## #### ########
##       ##          ##    ##        ##     ##
##       ##          ##    ##        ##     ##
######   ######      ##    ######    ##     ##
##       ##          ##    ##        ##     ##
##       ##          ##    ##        ##     ##
######## ##          ##    ##       ####    ##

class EFTFit(object):

 # #    # # #####
 # ##   # #   #
 # # #  # #   #
 # #  # # #   #
 # #   ## #   #
 # #    # #   #

    def __init__(self, opts):
        self.logger = logging.getLogger(__name__)

        self.wcs = opts['wcs']
        self.wc = opts['wc']
        self.scan_wcs = opts['scan_wcs']
        self.wcs_tracked = opts['wcs_tracked']
        self.wc_ranges = opts['wc_ranges']
        self.wcs_pairs = opts['wcs_pairs']
        self.SM_mu = opts['SM_mu']
        self.SM_mus = opts['SM_mus']

        # Systematics names except for FR stats. Only used for debug
        self.systematics = []


 #       ####   ####   ####  ###### #####
 #      #    # #    # #    # #      #    #
 #      #    # #      #      #####  #    #
 #      #    # #  ### #  ### #      #####
 #      #    # #    # #    # #      #   #
 ######  ####   ####   ####  ###### #    #

        log_file = 'fitter.log'

        FORMAT1 = '%(message)s'
        FORMAT2 = '[%(levelname)s] %(message)s'
        FORMAT3 = '[%(levelname)s][%(name)s] %(message)s'

        frmt1 = logging.Formatter(FORMAT1)
        frmt2 = logging.Formatter(FORMAT2)
        frmt3 = logging.Formatter(FORMAT3)

        logging.basicConfig(
            level=logging.DEBUG,
            format=FORMAT2,
            filename=log_file,
            filemode='w'
        )

        # Configure logging to also output to stdout
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(frmt2)
        logging.getLogger('').addHandler(console)


 #    # ###### #      #####  ###### #####
 #    # #      #      #    # #      #    #
 ###### #####  #      #    # #####  #    #
 #    # #      #      #####  #      #####
 #    # #      #      #      #      #   #
 #    # ###### ###### #      ###### #    #

    def log_subprocess_output(self,pipe,level):
        ### Pipes Popen streams to logging class ###
        for line in iter(pipe.readline, ''):
            if level=='info': logging.info(line.rstrip('\n'))
            if level=='err': logging.error(line.rstrip('\n'))


 #    #  ####  #####  #    #  ####  #####    ##    ####  ######
 #    # #    # #    # #   #  #      #    #  #  #  #    # #
 #    # #    # #    # ####    ####  #    # #    # #      #####
 # ## # #    # #####  #  #        # #####  ###### #      #
 ##  ## #    # #   #  #   #  #    # #      #    # #    # #
 #    #  ####  #    # #    #  ####  #      #    #  ####  ######

    def makeWorkspaceSM(self, datacard='datacard.txt', verbosity=0):
        ### Generates a workspace from a datacard ###
        logging.info(colors.fg.lightblue + "Creating workspace" + colors.reset)
        if not os.path.isfile(datacard):
            logging.error("Datacard does not exist!")
            return
        # CMSSW_BASE = os.getenv('CMSSW_BASE')

        args = ['text2workspace.py',datacard,'-P','HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel', '-o','SMWorkspace.root']

        # Map signal strengths to processes in all bins
        for iproc,proc in enumerate(opts["processes"]):
            args.extend(['--PO', 'map=.*/'+proc+':'+opts["SM_mus"][iproc]+'[1,'+str(opts["SMmu_ranges"][opts["SM_mus"][iproc]][0])+','+str(opts["SMmu_ranges"][opts["SM_mus"][iproc]][1])+']'])
        if verbosity>0:
            args.extend(['-v', str(verbosity)])
        args.extend(['--channel-masks']) #Creates additional parameters allowing to later mask specific channels

        logging.info(colors.fg.purple + " ".join(args) + colors.reset)
        process = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
        with process.stdout,process.stderr:
            self.log_subprocess_output(process.stdout,'info')
            self.log_subprocess_output(process.stderr,'err')
        process.wait()


    def makeWorkspaceEFT(self, datacard='datacard.txt', verbosity=0):
        ### Generates a workspace from a datacard and fit parameterization file ###
        logging.info(colors.fg.lightblue + "Creating workspace" + colors.reset)
        if not os.path.isfile(datacard):
            logging.error("Datacard does not exist!")
            sys.exit()
        CMSSW_BASE = os.getenv('CMSSW_BASE')
        args = ['text2workspace.py',datacard,'-P','EFTModel:eftmodel','--PO','fits=./Parameterization_EFT.npy','-o','EFTWorkspace.root']
        if verbosity>0: args.extend(['-v', str(verbosity)])
        args.extend(['--channel-masks']) #Creates additional parameters allowing to later mask specific channels

        logging.info(colors.fg.purple + ' '.join(args) + colors.reset)
        process = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
        with process.stdout,process.stderr:
            self.log_subprocess_output(process.stdout,'info')
            self.log_subprocess_output(process.stderr,'err')
        process.wait()


 #####  ######  ####  #####    ###### # #####
 #    # #      #        #      #      #   #
 #####  #####   ####    #      #####  #   #
 #    # #           #   #      #      #   #
 #    # #      #    #   #      #      #   #
 #####  ######  ####    #      #      #   #

    def bestFitSM(self, name='.SM', params_POI=[], freeze=[], autoMaxPOIs=True, other=[], exp=False, verbosity=0, mask=[], antimask=[]):
        '''
        Perform a MLF to find the best fit value of the POI(s).
        '''

        ### Multidimensional fit ###
        logging.info(colors.fg.lightblue + "Enter function bestFitSM()\n" + colors.reset)

        if len(params_POI) == 0: params = self.SM_mus

        #-- #Define channel masking regexp pattern, if any
        maskPattern = []; antimaskPattern = []
        if len(mask)>0: maskPattern=[','.join('rgx{{mask.*_{}_.*}}=1'.format(chan) for chan in mask)] #Use '{{' to insert a litteral bracket, not a replacement field #More info on regexp meaning: https://regex101.com/
        #if len(antimask)>0: antimaskPattern=[','.join('rgx{{^mask_(?!.*{}).*$}}=1'.format(chan) for chan in antimask)]
        if len(antimask)>0: antimaskPattern=['rgx{^mask_(?!.*('+'|'.join('{}'.format(chan) for chan in antimask)+')).*$}=1'] #Opposite: mask all channels NOT matching ANY 'chan'

        args=['combine', '-d','./SMWorkspace.root', '-v','2', '-M','MultiDimFit', '--saveFitResult','--cminPoiOnlyFit','--do95','1','--robustFit','1']

        for mu in params_POI: args.extend(['-P', '{}'.format(mu)]) #Define signal strengths as POIs
        #args.extend(['--setParameters',','.join('{}=1'.format(mu) for mu in self.SM_mus)]) 
        args.extend(['--setParameters',','.join([','.join('{}=1'.format(mu) for mu in self.SM_mus)]+maskPattern+antimaskPattern)]) #Set default values to 1
        if freeze:
            frozen_pois = [wc for wc in self.wcs if wc not in params_POI] #Define which WCs are frozen
            args.extend(['--freezeParameters',','.join('{}'.format(mu) for mu in self.SM_mus if mu not in params_POI and len(frozen_pois)>0)]) #Freeze other parameters
        else: args.extend(['--floatOtherPOIs','1']) #Float other parameters defined in the physics model
        if name: args.extend(['-n','{}'.format(name)])
        if verbosity>0: args.extend(['-v', str(verbosity)])
        if exp: args.extend(['-t', '-1'])
        if other:args.extend(other)

        if debug: print('args --> ', args)
        logging.info(colors.fg.purple + " ".join(args) + colors.reset)
        process = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
        with process.stdout,process.stderr:
            self.log_subprocess_output(process.stdout,'info')
            self.log_subprocess_output(process.stderr,'err')
        process.wait()
        logging.info(colors.fg.lightblue + "Done with SM best fit !" + colors.reset)
        self.printBestFit(name=name, params=params_POI)
    
        return


    def bestFitEFT(self, name='.EFT', params_POI=[], freeze=False, startValue='', autoBounds=True, other=[], exp=False, verbosity=0, fixedPointNLL=False, mask=[], antimask=[]):
        '''
        Perform a (multi-dim.) MLF to find the best fit value of the POI(s).

        NB: the error 'Looks like the last fit did not float this parameter' could arise e.g. if the bin names differ between file/datacard/EFT parametrization file.
        '''

        logging.info(colors.fg.lightblue + "Enter function bestFitEFT()\n" + colors.reset)

        if params_POI == []: params_POI = self.wcs
        if name == '': name = '.EFT'

        #-- #Define channel masking regexp pattern, if any
        maskPattern = []; antimaskPattern = []
        if len(mask)>0: maskPattern=[','.join('rgx{{mask.*_{}_.*}}=1'.format(chan) for chan in mask)] #Use '{{' to insert a litteral bracket, not a replacement field #More info on regexp meaning: https://regex101.com/
        #if len(antimask)>0: antimaskPattern=[','.join('rgx{{^mask_(?!.*{}).*$}}=1'.format(chan) for chan in antimask)]
        if len(antimask)>0: antimaskPattern=['rgx{^mask_(?!.*('+'|'.join('{}'.format(chan) for chan in antimask)+')).*$}=1'] #Opposite: mask all channels NOT matching ANY 'chan'

        # CMSSW_BASE = os.getenv('CMSSW_BASE')
        # args=['combine','-d','./EFTWorkspace.root','-M','MultiDimFit','--saveNLL','--saveFitResult','-H','AsymptoticLimits','--cminPoiOnlyFit']
        args=['combine','-d','./EFTWorkspace.root','-M','MultiDimFit','--saveNLL','--saveFitResult','--do95','1','--robustFit','1']

        args.extend(['-n','{}'.format(name)])
        if fixedPointNLL:
            args.extend(['--X-rtd','REMOVE_CONSTANT_ZERO_POINT=1']) #Access absolute NLL #Necessary/useful ?
            args.extend(['--algo','fixed','--fixedPointPOIs','{}={}'.format(opts['wc'],startValue)])
        if params_POI:
            for wc in params_POI: args.extend(['-P','{}'.format(wc)])
        args.extend(['--setParameters',','.join([','.join('{}=0'.format(poi) for poi in self.wcs)]+maskPattern+antimaskPattern)]) #Set default values to 0 #Mask channels if needed
        if freeze:
            frozen_pois = [wc for wc in self.wcs if wc not in params_POI] #Define which WCs are frozen
            args.extend(['--freezeParameters',','.join('{}'.format(poi) for poi in opts['wcs'] if poi not in params_POI and len(frozen_pois)>0)]) #Freeze other parameters
        else: args.extend(['--floatOtherPOIs','1']) #Float other parameters defined in the physics model
        if autoBounds:          args.extend(['--autoBoundsPOIs=*']) #Auto adjust POI bounds if found close to boundary
        if exp:                 args.extend(['-t','-1']) #Assume MC expected (Asimov?)
        if verbosity>0:           args.extend(['-v', str(verbosity)])
        if other:               args.extend(other)
        check = True in (wc not in params_POI for wc in self.wcs)
        if check: args.extend(['--trackParameters',','.join(wc for wc in self.wcs_tracked if wc not in params_POI)]) #Save values of additional parameters (e.g. profiled nuisances)
        args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(wc,self.wc_ranges[wc][0],self.wc_ranges[wc][1]) for wc in self.wcs)]) #in params_POI

        if debug: print('args --> ', args)
        logging.info(colors.fg.purple + " ".join(args) + colors.reset)
        process = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
        with process.stdout,process.stderr:
            self.log_subprocess_output(process.stdout,'info')
            self.log_subprocess_output(process.stderr,'err')
        process.wait()
        logging.info(colors.fg.lightblue + "Done with bestFitEFT." + colors.reset)

        if not fixedPointNLL: self.printBestFit(name=name, params=params_POI)


  ####  #####  # #####      ####   ####    ##   #    #
 #    # #    # # #    #    #      #    #  #  #  ##   #
 #      #    # # #    #     ####  #      #    # # #  #
 #  ### #####  # #    #         # #      ###### #  # #
 #    # #   #  # #    #    #    # #    # #    # #   ##
  ####  #    # # #####      ####   ####  #    # #    #

    #FIXME -- could merge with EFT func ? (also harmonize batch mode, ...) #Check when will rerun it
    def gridScanSM(self, name='.SM', batch='', scan_params=['r_tzq'], params_tracked=[], points=300, freeze=False, other=[], exp=False, verbosity=0):

        ### Runs deltaNLL Scan in a parameter using CRAB ###
        ### Can be used to do 2D scans as well ###
        logging.info(colors.fg.lightblue + "Enter function gridScanSM()\n" + colors.reset)

        if(len(scan_params) == 2):
            # points*= points #2D --> square npoints
            name = name + '.2D'

        # CMSSW_BASE = os.getenv('CMSSW_BASE')
        args = ['combineTool.py','-d','./SMWorkspace.root','-M','MultiDimFit','--algo','grid','--cminPreScan','--cminDefaultMinimizerStrategy=0']

        for mu in scan_params: args.extend(['-P', '{}'.format(mu)]) #Define signal strengths as POIs
        args.extend(['--setParameters',','.join('{}=1'.format(mu) for mu in scan_params)]) #Set default values to 1
        args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(mu,opts["SMmu_ranges"][mu][0],opts["SMmu_ranges"][mu][1]) for mu in opts["SM_mus"])])

        args.extend(['--points','{}'.format(points)])
        if name:                args.extend(['-n','{}'.format(name)])
        if params_tracked:      args.extend(['--trackParameters',','.join(params_tracked)])
        if verbosity>0:           args.extend(['-v', str(verbosity)])
        if not freeze:          args.extend(['--floatOtherPOIs','1']) #Float other parameters defined in the physics model
        if exp: args.extend(['-t', '-1'])
        if other:             args.extend(other)

        if batch=='crab':      args.extend(['--job-mode','crab3','--task-name',name.replace('.',''),'--custom-crab','Utils/custom_crab.py','--split-points','2000'])
        if batch=='condor':    args.extend(['--job-mode','condor','--task-name',name.replace('.',''),'--split-points','2000','--sub-opts=','\'getenv = true\''])

        if debug: print('args --> ', args)
        logging.info(colors.fg.purple + ' '.join(args) + colors.reset)

        process = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
        with process.stdout,process.stderr:
            self.log_subprocess_output(process.stdout,'info')
            self.log_subprocess_output(process.stderr,'err')
        process.wait()
        logging.info(colors.fg.lightblue + "Done with gridScan batch submission." + colors.reset)

        if not batch:
            logging.info(colors.fg.lightblue + "Done with gridScan." + colors.reset)

        return


    def gridScanEFT(self, name='.EFT', batch='', freeze=False, scan_params=['ctz'], startValuesString='', params_tracked=[], points=1000, exp=False, other=[], verbosity=0):

        ### Runs deltaNLL Scan in two parameters using CRAB or Condor ###
        logging.info(colors.fg.lightblue + "Enter function gridScanEFT()\n" + colors.reset)

        # CMSSW_BASE = os.getenv('CMSSW_BASE')
        # Other options: '--cminFallbackAlgo Minuit2,Combined,2:0.3'
        args = ['combineTool.py','-d','./EFTWorkspace.root','-M','MultiDimFit','--algo','grid','--cminPreScan','--cminDefaultMinimizerStrategy=0']

        for wc in scan_params: args.extend(['-P', '{}'.format(wc)]) #Define signal strengths as POIs
        args.extend(['--setParameters',','.join('{}=0'.format(wc) for wc in opts["wcs"])]) #Set default values to 1
        args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(wc,self.wc_ranges[wc][0],self.wc_ranges[wc][1]) for wc in self.wcs)])

        args.extend(['--points','{}'.format(points)])
        if name: args.extend(['-n','{}'.format(name)])
        check = True in (wc not in self.wcs for wc in self.wcs_tracked)
        if check: args.extend(['--trackParameters',','.join([wc for wc in self.wcs_tracked if wc not in self.wcs])]) #Save values of additional parameters (e.g. profiled nuisances)
        # if startValuesString:   args.extend(['--setParameters',startValuesString])
        if freeze: args.extend(['--freezeParameters',','.join('{}'.format(poi) for poi in opts['wcs'] if poi not in scan_params)]) #Freeze other parameters
        else: args.extend(['--floatOtherPOIs','1']) #Float other parameters defined in the physics model
        if exp:               args.extend(['-t -1'])
        if verbosity>0:           args.extend(['-v', str(verbosity)])
        if other:             args.extend(other)
        # args.extend(['--fastScan']) #No profiling (speed up) of nuisances, kept to best fit value

        if batch=='crab': args.extend(['--job-mode','crab3','--task-name',name.replace('.',''),'--custom-crab','Utils/custom_crab.py','--split-points','100'])
        if batch=='condor': args.extend(['--job-mode','condor','--task-name',name.replace('.',''),'--split-points','100','--dry-run'])

        if debug: print('args --> ', args)
        logging.info(colors.fg.purple + ' '.join(args) + colors.reset)

        # Run the combineTool.py command
        process = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
        with process.stdout,process.stderr:
            self.log_subprocess_output(process.stdout,'info')
            self.log_subprocess_output(process.stderr,'err')
        process.wait()

        # Condor needs executab permissions on the .sh file, so we used --dry-run
        # Add the permission and complete the submission.
        if batch=='condor':
            if os.path.exists('condor{}'.format(name)):
                logging.error("Directory condor{} already exists!".format(name))
                logging.error("OVERWRITING !")
                shutil.rmtree('condor{}'.format(name), ignore_errors=True)
                #return

            sp.call(['mkdir','condor{}'.format(name)])
            sp.call(['chmod','a+x','condor_{}.sh'.format(name.replace('.',''))])
            logging.info('Now submitting condor jobs.')
            condorsub = sp.Popen(['condor_submit','-append','initialdir=condor{}'.format(name),'condor_{}.sub'.format(name.replace('.',''))], stdout=sp.PIPE, stderr=sp.PIPE)
            with condorsub.stdout,condorsub.stderr:
                self.log_subprocess_output(condorsub.stdout,'info')
                self.log_subprocess_output(condorsub.stderr,'err')
            condorsub.wait()

        if batch: logging.info(colors.fg.lightblue + "Done with gridScan batch submission." + colors.reset)
        else: logging.info(colors.fg.lightblue + "Done with gridScan." + colors.reset)

        if batch=='': 
            fitter.printIntervalFitsEFT(basename='.EFT', scan_params=param_tmp) #Print exclusion range #Obsolete

        return


 #####    ##   #####  ####  #    #
 #    #  #  #    #   #    # #    #
 #####  #    #   #   #      ######
 #    # ######   #   #      #    #
 #    # #    #   #   #    # #    #
 #####  #    #   #    ####  #    #

    def batch1DScanSM(self, basename='.test', batch='', scan_params=[], points=300, freeze=False, other=[]):
        ### For each SM signal strength, run a 1D deltaNLL Scan.
        if not scan_params:
            scan_params = ['r_tzq']

        for param in scan_params:
            scanmax = 3
            if param=='r_ttH': scanmax = 6
            if param=='r_tllq': scanmax = 4
            self.gridScanSM('{}.{}'.format(basename,param), batch, [param], self.systematics+[params for params in scan_params if params != param], points, freeze, ['--setParameterRanges','{}=0,{}'.format(param,scanmax)]+other)

    def batchRetrieve1DScansSM(self, basename='.test', batch='crab'):
        ### For each wc, retrieves finished 1D deltaNLL crab jobs, extracts, and hadd's into a single file ###
        for param in ['r_ttll','r_ttlnu','r_ttH','r_tllq']:
            self.retrieveGridScan('{}.{}'.format(basename,param),batch)

    def batch1DScanEFT(self, basename='.EFT', batch='crab', freeze=False, scan_wcs=[], points=300, other=[]):
        ### For each wc, run a 1D deltaNLL Scan.
        if not scan_wcs:
            scan_wcs = self.wcs

        for wc in scan_wcs:
            self.gridScanEFT('{}.{}'.format(basename,wc), batch, freeze, [wc], [wcs for wcs in self.wcs if wcs != wc], points, other)

    def batch2DScanEFT(self, basename='.EFT.gridScan', batch='crab', freeze=False, points=90000, allPairs=False, other=[]):
        ### For pairs of wcs, runs deltaNLL Scan in two wcs using CRAB or Condor ###

        # Use EVERY combination of wcs
        if allPairs:
            scan_wcs = self.wcs

            for wcs in itertools.combinations(scan_wcs,2):
                wcs_tracked = [wc for wc in self.wcs if wc not in wcs]
                #print pois, wcs_tracked
                self.gridScanEFT(name='{}.{}{}'.format(basename,wcs[0],wcs[1]), batch=batch, freeze=freeze, scan_params=list(wcs), params_tracked=wcs_tracked, points=points, other=other)

        # Use each wc only once
        if not allPairs:
            scan_wcs = [('ctw','ctG'),('ctz','ctG'),('ctp','ctG'),('cpQM','ctG'),('cbW','ctG'),('cpQ3','ctG'),('cptb','ctG'),('cpt','ctG'),('cQl3i','ctG'),('cQlMi','ctG'),('cQei','ctG'),('ctli','ctG'),('ctei','ctG'),('ctlSi','ctG'),('ctlTi','ctG')]
            #pairs from AN
            scan_wcs = [('cQlMi','cQei'),('cpQ3','cbW'),('cptb','cQl3i'),('ctG','cpQM'),('ctz','ctw'),('ctei','ctlTi'),('ctlSi','ctli'),('ctp','cpt')]

            for wcs in scan_wcs:
                wcs_tracked = [wc for wc in self.wcs if wc not in wcs]
                #print pois, wcs_tracked
                self.gridScanEFT(name='{}.{}{}'.format(basename,wcs[0],wcs[1]), batch=batch, freeze=freeze, scan_params=list(wcs), params_tracked=wcs_tracked, points=points, other=other)

    def batch3DScanEFT(self, basename='.EFT.gridScan', batch='crab', freeze=False, points=27000000, allPairs=False, other=[], wc_triplet=[]):
        ### For pairs of wcs, runs deltaNLL Scan in two wcs using CRAB or Condor ###

        # Use EVERY combination of wcs
        if allPairs:
            scan_wcs = self.wcs

            for wcs in itertools.combinations(scan_wcs,2):
                wcs_tracked = [wc for wc in self.wcs if wc not in wcs]
                #print pois, wcs_tracked
                self.gridScanEFT(name='{}.{}{}{}'.format(basename,wcs[0],wcs[1],wcs[2]), batch=batch, freeze=freeze, scan_params=list(wcs), params_tracked=wcs_tracked, points=points, other=other)

        # Use each wc only once
        if not allPairs:
            scan_wcs = [('ctz','ctp','cpt')]
            if len(wc_triplet)>0: scan_wcs = wc_triplet
            #scan_wcs = [('ctz','ctw'),('ctp','cpt'),('ctlSi','ctli'),('cptb','cQl3i'),('ctG','cpQM'),('ctei','ctlTi'),('cQlMi','cQei'),('cpQ3','cbW')]
            #scan_wcs = [('ctw','ctG'),('ctz','ctG'),('ctp','ctG'),('cpQM','ctG'),('cbW','ctG'),('cpQ3','ctG'),('cptb','ctG'),('cpt','ctG'),('cQl3i','ctG'),('cQlMi','ctG'),('cQei','ctG'),('ctli','ctG'),('ctei','ctG'),('ctlSi','ctG'),('ctlTi','ctG')]

            for wcs in scan_wcs:
                wcs_tracked = [wc for wc in self.wcs if wc not in wcs]
                #print pois, wcs_tracked
                self.gridScanEFT(name='{}.{}{}{}'.format(basename,wcs[0],wcs[1],wcs[2]), batch=batch, freeze=freeze, scan_params=list(wcs), params_tracked=wcs_tracked, points=points, other=other)

    def batchResubmit1DScansEFT(self, basename='.EFT.gridScan', scan_wcs=[]):
        ### For each wc, attempt to resubmit failed CRAB jobs ###
        if not scan_wcs:
            scan_wcs = self.wcs

        for wc in scan_wcs:
            process = sp.Popen(['crab','resubmit','crab_'+basename.replace('.','')+wc], stdout=sp.PIPE, stderr=sp.PIPE)
            with process.stdout,process.stderr:
                self.log_subprocess_output(process.stdout,'info')
                self.log_subprocess_output(process.stderr,'err')
            process.wait()

    def batchResubmit2DScansEFT(self, basename='.EFT.gridScan', allPairs=False):
        ### For pairs of wcs, attempt to resubmit failed CRAB jobs ###

        # Use EVERY combination of wcs
        if allPairs:
            scan_wcs = self.wcs

            for wcs in itertools.combinations(scan_wcs,2):
                process = sp.Popen(['crab','resubmit','crab_'+basename.replace('.','')+wcs[0]+wcs[1]], stdout=sp.PIPE, stderr=sp.PIPE)
                with process.stdout,process.stderr:
                    self.log_subprocess_output(process.stdout,'info')
                    self.log_subprocess_output(process.stderr,'err')
                process.wait()

        # Use each wc only once
        if not allPairs:
            scan_wcs = [('ctz','ctw'),('ctp','cpt'),('ctlSi','ctli'),('cptb','cQl3i'),('ctG','cpQM'),('ctei','ctlTi'),('cQlMi','cQei'),('cpQ3','cbW')]
            scan_wcs = [('cQlMi','cQei'),('cpQ3','cbW'),('cptb','cQl3i'),('ctG','cpQM'),('ctz','ctw'),('ctei','ctlTi'),('ctlSi','ctli'),('ctp','cpt')]

            for wcs in scan_wcs:
                process = sp.Popen(['crab','resubmit','crab_'+basename.replace('.','')+wcs[0]+wcs[1]], stdout=sp.PIPE, stderr=sp.PIPE)
                with process.stdout,process.stderr:
                    self.log_subprocess_output(process.stdout,'info')
                    self.log_subprocess_output(process.stderr,'err')
                process.wait()


  ####  ###### #####    #####  ######  ####  #####    #    #   ##   #
 #    # #        #      #    # #      #        #      #    #  #  #  #
 #      #####    #      #####  #####   ####    #      #    # #    # #
 #  ### #        #      #    # #           #   #      #    # ###### #
 #    # #        #      #    # #      #    #   #       #  #  #    # #
  ####  ######   #      #####  ######  ####    #        ##   #    # ######

    def getBestValues2D(self, name, scan_params=[], params_tracked=[]):
        ### Gets values of parameters for grid scan point with best deltaNLL ###

        bestDeltaNLL=1000000;
        bestEntry=-1;

        fitFile = './higgsCombine'+name+'.MultiDimFit.root'
        print fitFile

        if not os.path.isfile(fitFile):
            logging.error("fitFile does not exist!")
            sys.exit()
        rootFile = ROOT.TFile.Open(fitFile);
        limitTree = rootFile.Get("limit");

        for entry in range(limitTree.GetEntries()):
            limitTree.GetEntry(entry)
            if(bestDeltaNLL > limitTree.GetLeaf("deltaNLL").GetValue(0)):
              bestDeltaNLL = limitTree.GetLeaf("deltaNLL").GetValue(0)
              bestEntry=entry
              #cout << Form("Entry %i deltaNLL=%f, ctw=%f ctz=%f",bestEntry,bestDeltaNLL,limitTree.GetLeaf("ctw").GetValue(0),limitTree.GetLeaf("ctz").GetValue(0)) << endl;

        limitTree.GetEntry(bestEntry)
        startValues = []
        for param in scan_params:
            value = limitTree.GetLeaf(param).GetValue(0)
            startValues.append('{}={}'.format(param,value))
        for param in params_tracked:
            value = limitTree.GetLeaf('trackedParam_'+param).GetValue(0)
            startValues.append('{}={}'.format(param,value))
        return ','.join(startValues)

    def getBestValues1DEFT(self, basename=".EFT", wcs=[]):
        ### Gets values of WCs for grid scan point with best deltaNLL ###
        if not wcs:
            wcs = self.wcs

        startValues = []

        for wc in wcs:

            bestDeltaNLL=1000000;
            bestEntry=-1;

            fitFile = './higgsCombine{}.{}.MultiDimFit.root'.format(basename,wc)
            logging.info("Obtaining best value from {}".format(fitFile))

            if not os.path.isfile(fitFile):
                logging.error("fitFile does not exist!")
                sys.exit()
            rootFile = ROOT.TFile.Open(fitFile);
            limitTree = rootFile.Get("limit");

            for entry in range(limitTree.GetEntries()):
                limitTree.GetEntry(entry)
                if(bestDeltaNLL > limitTree.GetLeaf("deltaNLL").GetValue(0)):
                  bestDeltaNLL = limitTree.GetLeaf("deltaNLL").GetValue(0)
                  bestEntry=entry
                  #cout << Form("Entry %i deltaNLL=%f, ctw=%f ctz=%f",bestEntry,bestDeltaNLL,limitTree.GetLeaf("ctw").GetValue(0),limitTree.GetLeaf("ctz").GetValue(0)) << endl;

            limitTree.GetEntry(bestEntry)

            value = limitTree.GetLeaf(wc).GetValue(0)
            startValues.append('{}={}'.format(wc,value))

        return ','.join(startValues)


 #####  ###### ##### #####  # ###### #    # ######     ####  #####  # #####
 #    # #        #   #    # # #      #    # #         #    # #    # # #    #
 #    # #####    #   #    # # #####  #    # #####     #      #    # # #    #
 #####  #        #   #####  # #      #    # #         #  ### #####  # #    #
 #   #  #        #   #   #  # #       #  #  #         #    # #   #  # #    #
 #    # ######   #   #    # # ######   ##   ######     ####  #    # # #####

    def retrieveGridScan(self, name='.test', batch='crab', user='ntonon'):#getpass.getuser()):

        ### Retrieve grid jobs outputs, extracts/hadd them into a single file ###
        taskname = name.replace('.','')
        logging.info("Retrieving gridScan files. Task name: "+taskname)

        if batch=='crab':
            # Find crab output files (defaults to user directory)
            outputspath = '/store/user/{}/EFT/Combine/{}'.format(user, taskname)
            (tarpath,tardirs,tarfiles) = os.walk(outputspath)
            if not tarfiles[2]:
                logging.error("No files found in store!")
                sys.exit()

            # Make a temporary folder to hold the extracted root files
            if not os.path.isdir(taskname+'tmp'):
                sp.call(['mkdir',taskname+'tmp'])
            else:
                logging.error("Directory {}tmp/ already exists! Please rename this directory.".format(taskname))
                return

            # Extract the root files
            for tarfile in tarfiles[2]:
                if tarfile.endswith('.tar'):
                    print tarfiles[0]+'/'+tarfile
                    sp.call(['tar', '-xf', tarfiles[0]+'/'+tarfile,'-C', taskname+'tmp'])
            haddargs = ['hadd','-f','./higgsCombine'+name+'.MultiDimFit.root']+['{}tmp/{}'.format(taskname,rootfile) for rootfile in os.listdir(taskname+'tmp') if rootfile.endswith('.root')]
            process = sp.Popen(haddargs, stdout=sp.PIPE, stderr=sp.PIPE)
            with process.stdout,process.stderr:
                self.log_subprocess_output(process.stdout,'info')
                self.log_subprocess_output(process.stderr,'err')
            process.wait()

            # Remove the temporary directory and split root files
            sp.call(['rm','-r',taskname+'tmp'])

        elif batch=='condor':
            if not glob.glob('higgsCombine{}.POINTS*.root'.format(name)): #glob: find matching patterns
                logging.info("No files with names higgsCombine{}.POINTS*.root to hadd. Returning.".format(name))
                return
            #haddargs = ['hadd','-f','higgsCombine'+name+'.MultiDimFit.root']+sorted(glob.glob('higgsCombine{}.POINTS*.root'.format(name)))
            haddargs = ['hadd','-f','./higgsCombine'+name+'.MultiDimFit.mH120.root']+sorted(glob.glob('higgsCombine{}.POINTS*.root'.format(name)))
            process = sp.Popen(haddargs, stdout=sp.PIPE, stderr=sp.PIPE)
            with process.stdout,process.stderr:
                self.log_subprocess_output(process.stdout,'info')
                self.log_subprocess_output(process.stderr,'err')
            process.wait()
            for rootfile in glob.glob('higgsCombine{}.POINTS*.root'.format(name)):
                os.remove(rootfile)
            if os.path.isfile('condor_{}.sh'.format(name.replace('.',''))):
                os.rename('condor_{}.sh'.format(name.replace('.','')),'condor{0}/condor_{0}.sh'.format(name))
            if os.path.isfile('condor_{}.sub'.format(name.replace('.',''))):
                os.rename('condor_{}.sub'.format(name.replace('.','')),'condor{0}/condor_{0}.sub'.format(name))

    def batchRetrieve1DScansEFT(self, basename='.EFT', batch='crab', scan_wcs=[]):
        ### For each wc, retrieves finished 1D deltaNLL grid jobs, extracts, and hadd's into a single file ###
        if not scan_wcs:
            scan_wcs = self.wcs

        for wc in scan_wcs:
            self.retrieveGridScan('{}.{}'.format(basename,wc),batch)

    def batchRetrieve2DScansEFT(self, wc_pair=[], basename='.EFT', batch='crab', allPairs=False):
        ### For pairs of wcs, retrieves finished grid jobs, extracts, and hadd's into a single file ###

        # Use EVERY combination of wcs
        if allPairs:
            wc_pair = self.wcs
            for wcs in itertools.combinations(wc_pair,2):
                self.retrieveGridScan(name='{}.{}{}'.format(basename,wcs[0],wcs[1]), batch=batch)

        # Consider a single pair of WCs
        else:
            if wc_pair == []: wc_pair = self.wcs_pairs
            print wc_pair
            print '{}.{}'.format(batch, basename)
            self.retrieveGridScan(name=basename, batch=batch)


 #####  ###### #####  #    #  ####  ##### #  ####  #    #
 #    # #      #    # #    # #    #   #   # #    # ##   #
 #    # #####  #    # #    # #        #   # #    # # #  #
 #####  #      #    # #    # #        #   # #    # #  # #
 #   #  #      #    # #    # #    #   #   # #    # #   ##
 #    # ###### #####   ####   ####    #   #  ####  #    #

    def reductionFitEFT(self, name='.EFT', wc='ctz', final=True, from_wcs=[], alreadyRun=True):
        ### Extract a 1D scan from a higher-dimension scan to avoid discontinuities ###
        if not wc:
            logging.error("No WC specified!")
            return
        if final and not alreadyRun:
            os.system('hadd -f higgsCombine{}.MultiDimFit.mH120.root higgsCombine{}.POINTS*.{}reduced.MultiDimFit.root '.format(name,name,''.join(from_wcs)))
        if alreadyRun and not os.path.exists('./higgsCombine{}.MultiDimFit.root'.format(name)):
            logging.error("File higgsCombine{}.MultiDimFit.root does not exist!".format(name))
            return
        elif not alreadyRun and not os.path.exists('higgsCombine{}.MultiDimFit.mH120.root'.format(name)):
        #if not os.path.exists('./higgsCombine{}.MultiDimFit.root'.format(name)):
            logging.error("File higgsCombine{}.MultiDimFit.root does not exist!".format(name))
            return

        rootFile = []
        if alreadyRun:
            rootFile = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.root'.format(name))
        else: rootFile = ROOT.TFile.Open('higgsCombine{}.MultiDimFit.mH120.root'.format(name))
        limitTree = rootFile.Get('limit')

        # First loop through entries and get deltaNLL list for each value of the WC
        wc_dict = defaultdict(list)
        for entry in range(limitTree.GetEntries()):
            limitTree.GetEntry(entry)
            wc_dict[limitTree.GetLeaf(wc).GetValue(0)].append(limitTree.GetLeaf('deltaNLL').GetValue(0))
        rootFile.Close()

        # Next pick the best deltaNLL for each WC value
        wc_dict_reduced = {}
        for key in wc_dict:
            wc_dict_reduced[key] = min(wc_dict[key])

        # Now make a new .root file with the new TTree
        # Only the WC and deltaNLL will be branches
        # These can be directly used by EFTPlotter
        outFile = []
        if final:
            outFile = ROOT.TFile.Open('./higgsCombine{}.{}reduced.MultiDimFit.root'.format(name,wc),'RECREATE')
        else:
            outFile = ROOT.TFile.Open('higgsCombine{}.{}reduced.MultiDimFit.root'.format(name,wc),'RECREATE')
        outTree = ROOT.TTree('limit','limit')

        wc_branch = array.array('f',[0.])
        deltaNLL_branch = array.array('f',[0.])
        outTree.Branch(wc,wc_branch,wc+'/F')
        outTree.Branch('deltaNLL',deltaNLL_branch,'deltaNLL/F')

        # Fill the branches
        for event in range(len(wc_dict_reduced.keys())):
            wc_branch[0] = wc_dict_reduced.keys()[event]
            deltaNLL_branch[0] = wc_dict_reduced.values()[event]
            outTree.Fill()

        # Write the file
        outFile.Write()

    def reduction2DFitEFT(self, name='.EFT.Private.Unblinded.Nov16.28redo.Float.cptcpQM', wcs=['cpt','ctp'], final=True):
        ### Extract a 2D scan from a higher-dimension scan to avoid discontinuities ###
        if not wcs:
            logging.error("No WC specified!")
            return
        if final:
            os.system('hadd -f higgsCombine{}.MultiDimFit.mH120.root higgsCombine{}.POINTS*.{}reduced.MultiDimFit.root '.format(name,name,''.join(wcs)))
        if not os.path.exists('higgsCombine{}.MultiDimFit.mH120.root'.format(name)):
        #if not os.path.exists('./higgsCombine{}.MultiDimFit.root'.format(name)):
            logging.error("File higgsCombine{}.MultiDimFit.root does not exist!".format(name))
            return

        rootFile = []
        rootFile = ROOT.TFile.Open('higgsCombine{}.MultiDimFit.mH120.root'.format(name))
        #rootFile = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.root'.format(name))
        limitTree = rootFile.Get('limit')

        # First loop through entries and get deltaNLL list for each value of the WC
        wc_dict = defaultdict(list)
        for entry in range(limitTree.GetEntries()):
            limitTree.GetEntry(entry)
            wc_dict[limitTree.GetLeaf(wcs[0]).GetValue(0),limitTree.GetLeaf(wcs[1]).GetValue(0)].append(limitTree.GetLeaf('deltaNLL').GetValue(0))
        rootFile.Close()

        # Next pick the best deltaNLL for each WC value
        wc_dict_reduced = {}
        for key in wc_dict:
            wc_dict_reduced[key] = min(wc_dict[key])

        # Now make a new .root file with the new TTree
        # Only the WC and deltaNLL will be branches
        # These can be directly used by EFTPlotter
        outFile = []
        if final:
            outFile = ROOT.TFile.Open('./higgsCombine{}.{}reduced.MultiDimFit.root'.format(name,''.join(wcs)),'RECREATE')
        else:
            outFile = ROOT.TFile.Open('higgsCombine{}.{}{}reduced.MultiDimFit.root'.format(name,wcs[0],wcs[1]),'RECREATE')
        outTree = ROOT.TTree('limit','limit')

        wc_branch1 = array.array('f',[0.])
        wc_branch2 = array.array('f',[0.])
        deltaNLL_branch = array.array('f',[0.])
        outTree.Branch(wcs[0],wc_branch1,wcs[0]+'/F')
        outTree.Branch(wcs[1],wc_branch2,wcs[1]+'/F')
        outTree.Branch('deltaNLL',deltaNLL_branch,'deltaNLL/F')

        # Fill the branches
        for wc1,wc2 in wc_dict_reduced:
            wc_branch1[0] = wc1
            wc_branch2[0] = wc2
            deltaNLL_branch[0] = wc_dict_reduced[(wc1,wc2)]
            outTree.Fill()

        # Write the file
        outFile.Write()


  ####   ####  #    # #####    ##   #####  ######    ###### # #####  ####
 #    # #    # ##  ## #    #  #  #  #    # #         #      #   #   #
 #      #    # # ## # #    # #    # #    # #####     #####  #   #    ####
 #      #    # #    # #####  ###### #####  #         #      #   #        #
 #    # #    # #    # #      #    # #   #  #         #      #   #   #    #
  ####   ####  #    # #      #    # #    # ######    #      #   #    ####

    def compareFitsEFT(self,basename='.EFT.SM.Float'):
        ### Compare results of different 1D EFT scans ###
        tfiles = {}
        limits = {}
        bestFits = {} # Nested dict; bestFit of key1 according to key2
        # First get all scan files
        for wc in self.wcs:
            tfiles[wc] = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.root'.format(basename+'.'+wc))
            limits[wc] = tfiles[wc].Get('limit')
            bestFits[wc] = {}
        # Get best fits
        for poiwc in self.wcs:
            limit = limits[poiwc]
            # First get POI best fit
            bestNLL = (-1,1000000)
            for entry in range(limit.GetEntries()):
                limit.GetEntry(entry)
                currentNLL = limit.GetLeaf('deltaNLL').GetValue(0)
                if bestNLL[1] > currentNLL: bestNLL = (entry,currentNLL)
            print "Best entry for {} is {}.".format(poiwc,bestNLL[0])
            limit.GetEntry(bestNLL[0])
            bestFits[poiwc][poiwc] = limit.GetLeaf(poiwc).GetValue(0)
            # Second get corresponding fits for the other wcs
            trackedwcs = list(self.wcs)
            trackedwcs.remove(poiwc)
            for trackedwc in trackedwcs:
                bestFits[trackedwc][poiwc] = limit.GetLeaf('trackedParam_'+wc).GetValue(0)

        # Print full set of results
        for poiwc in self.wcs:
            trackedwcs = list(self.wcs)
            trackedwcs.remove(poiwc)
            print("Best value of {}: {}".format(poiwc,bestFits[poiwc][poiwc]))
            for trackedwc in trackedwcs:
                print("Value according to {}: {}".format(trackedwc,bestFits[poiwc][trackedwc]))


 #####  #####  # #    # #####
 #    # #    # # ##   #   #
 #    # #    # # # #  #   #
 #####  #####  # #  # #   #
 #      #   #  # #   ##   #
 #      #    # # #    #   #

    def printBestFit(self, name='.EFT', params = []):
        '''
        Print a table of SM signal strengths, their best fits, and their uncertainties ###

        Example to inspect file via command line:
        root multidimfit.EFT.root
        a = fit_mdf->floatParsFinal().find("ctz")
        b = (RooAbsReal*) a
        b->Print()
        rf = dynamic_cast<RooRealVar*>(a)
        rf->getMin("err68") #...

        Example in python:
        import ROOT
        fit_file = ROOT.TFile.Open('multidimfit.EFT.root')
        fit = fit_file.Get('fit_mdf')
        roorealvar = fit.floatParsFinal().find('ctz')
        value = round(roorealvar.getVal(),2) #Best fit value
        value = round(roorealvar.getMin('err68'),2) #Lower 68% limit; idem with getMax, err95
        '''

        logging.info(colors.fg.lightblue + "\nEnter function printBestFit()" + colors.reset)

        if len(params)==0: 
            print('params is empty... Return !')

        fit_array = []

        logging.info("Obtaining result of fit: multidimfit{}.root".format(name))
        fit_file = ROOT.TFile.Open('./multidimfit{}.root'.format(name))
        fit = fit_file.Get('fit_mdf')

        for param in params:
            roorealvar = fit.floatParsFinal().find(param)
            if not roorealvar: continue

            value = round(roorealvar.getVal(),3)
            #err_sym =  round(roorealvar.getError(),3)
            err_low = round(roorealvar.getErrorLo(),3)
            err_high = round(roorealvar.getErrorHi(),3)
            err_low_95 = -9; err_high_95 = -9
            if roorealvar.hasRange('err95'): #If 95% CL errors available (using --do95 --robustFit 1 options)
                err_low_95 = round(roorealvar.getMin('err95'),3)
                err_high_95 = round(roorealvar.getMax('err95'),3)

            fit_array.append((param,value,err_low,err_high,err_low_95,err_high_95))

        logging.info('\n' + colors.fg.orange + "Fit result:" + colors.reset)
        logging.info(colors.fg.orange + "Param | Best Fit Value | [68% interval] | [95% interval]" + colors.reset)
        for row in fit_array:
            logging.info(colors.fg.orange + row[0] + ' | ' + str(row[1]) + ' | ' + "[" + str(row[2]) + ";" + str(row[3]) + "] | [" + str(row[4]) + ';' + str(row[5]) + ']' + colors.reset + '\n')
            #logging.debug("{} {} +/- {}".format(row[0],row[1],row[2]))
                        
        return


    def printIntervalFitsEFT(self, basename='.EFT', scan_params=[]):
        ### Print a table of wcs, their best fits, and their uncertainties ###
        ### Use 1D scans instead of regular MultiDimFit ###
        logging.info(colors.fg.lightblue + "Enter function printIntervalFitsEFT()\n" + colors.reset)

        if not scan_params: scan_params = [self.wc]
        ROOT.gROOT.SetBatch(True)
        fit_array = []
        canvas = ROOT.TCanvas()
        for wc in scan_params:

            canvas.Clear()

            logging.info("Obtaining result of scan: higgsCombine{}.MultiDimFit.mH120.root for WC {}".format(basename,wc))

            #-- Get scan TTree
            rootFile = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.mH120.root'.format(basename))
            limitTree = rootFile.Get('limit')

            #-- Use CombineTool utils (see: https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools/python/plotting.py)
            graph = plot.TGraphFromTree(limitTree, wc, '2*deltaNLL', 'quantileExpected > -1.5')

            yvals = [1., 3.84] #1sigma, 95%CL intervals
            #func, crossings, val, val_2sig, cross_1sig, cross_2sig, other_1sig, other_2sig = BuildScan(graph, ROOT.kBlack, yvals)
            main_scan = BuildScan(graph, ROOT.kBlack, yvals)

            crossings = main_scan['crossings'][yvals[1]][0]
            if crossings['valid_lo'] == True and crossings['valid_hi'] == True:
                print(colors.fg.orange + wc + ": 95% interval: [" + str(crossings['lo']) + ", " + str(crossings['hi']) + "]" + colors.reset)
            else: print('Error: invalid crossing X-values...')

        return


    '''
    def printIntervalFitsEFT(self, basename='.EFT', scan_params=[]):
        ### Print a table of wcs, their best fits, and their uncertainties ###
        ### Use 1D scans instead of regular MultiDimFit ###
        #-- NB: this method is actually very imprecise (takes crossings as (x1+x2)/2), relies on very fine scans. Better to use Combine function 'BuildScan' like in EFTPlotter. See new function above.
        logging.info(colors.fg.lightblue + "Enter function printIntervalFitsEFT()\n" + colors.reset)

        NLL_threshold = 3.84 #Define NN threshold to determine exclusion boundaries #3.84 <-> 95%

        if not scan_params: scan_params = self.wc
        ROOT.gROOT.SetBatch(True)
        fit_array = []
        canvas = ROOT.TCanvas()
        for wc in scan_params:

            canvas.Clear()

            logging.info("Obtaining result of scan: higgsCombine{}.MultiDimFit.mH120.root".format(basename))
            fit_file = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.mH120.root'.format(basename))
            limit_tree = fit_file.Get('limit')

            limit_tree.Draw('2*deltaNLL:{}>>{}1DNLL(50,{},{})'.format(wc,wc,self.wc_ranges[wc][0],self.wc_ranges[wc][1]),'2*deltaNLL>-1','same')
            graph = canvas.GetPrimitive('Graph')
            # print(graph.GetN())

            graph.Sort() #If error here, problem with limit tree

            lowedges=[]
            highedges=[]
            minimums=[]
            true_minimums=[]
            best = [-1000,1000]
            prev = 1000
            for idx in range(1, graph.GetN()):
                # print('idx', idx)
                y_val = graph.GetY()[idx]
                if prev>NLL_threshold and NLL_threshold>y_val:
                    lowedges.append((graph.GetX()[idx-1]+graph.GetX()[idx+1])/2)
                if prev<NLL_threshold and NLL_threshold<y_val:
                    highedges.append((graph.GetX()[idx-1]+graph.GetX()[idx+1])/2)
                if y_val < best[1]:
                    best = [graph.GetX()[idx],y_val]
                if y_val<prev and y_val<graph.GetY()[idx+1]:
                    minimums.append((graph.GetX()[idx],y_val))
                prev = y_val
            if not len(lowedges) == len(highedges):
                logging.error("Something is strange! Interval is missing endpoint!")

            def sortkey(elem):
                return elem[1]

            for interval in zip(lowedges,highedges):
                true_min = [-1000,1000]
                for minimum in minimums:
                    if minimum[1]<true_min[1] and interval[0]<minimum[0] and minimum[0]<interval[1]:
                        true_min = minimum
                true_minimums.append(true_min[0])

            fit_array.append([wc,[list(l) for l in zip(true_minimums,lowedges,highedges)]])

        for line in fit_array:
            print line
        
        return
    ''' 







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

# User options -- Default values
# //--------------------------------------------
    SM = False #True <-> consider SM scenario (rather than SMEFT)
    datacard_path = './datacard.txt'
    exp = False #True <-> Asimov a-priori expected; False <-> observed
    scan_dim = '1D' #'1d' <-> 1D scan (default); '2d' <-> 2D scan
    verb=0
    name = '' #Suffix added to output filenames
    startValue = '' #Starting value for the POI
    fixedPointNLL = False #True <-> perform the NLL scan at a fixed point (different Combine options)
    debug=False
    freeze=False
    createWS = 0 #0 <-> create WS and proceed ; 1 <-> create WS and exit ; 2 <-> don't create WS and proceed
    POI=[]
    mode = '' #Can choose to run only specific functions (not all) #'','grid','bestfit'
    batch = '' #Can choose to run jobs on 'crab' or 'condor'
    dryrun = '' #Perform dry run (don't submit jobs)
    points = -1 #Choose npoints for grid scans #-1 <-> use default values set below (different for 1D/2D)
    mask = [] #Can choose to mask specific channels from the likelihood #Mask channels matching pattern(s)
    antimask = [] #Mask channels NOT matching any pattern

# Set up the command line arguments
# //--------------------------------------------
    parser = argparse.ArgumentParser(description='Perform SM and EFT fits using custom Physics Model')
    parser.add_argument("-d", metavar="datacard path", help="Path to the txt datacard (to create RooWorkspace)")
    parser.add_argument("-v", metavar="Combine verbosity level", help="Set combine output verbosity")
    parser.add_argument("--sm", metavar="SM", help="Consider SM scenario (rather than SMEFT)")
    parser.add_argument("-dim", metavar="dim", help="1D or 2D scan")
    parser.add_argument('--exp', help='Use MC predictions only (no data)', nargs='?', const=1) #nargs='?' <-> 0 or 1 arg (default value is const=1)
    parser.add_argument("--fixed", metavar="fixed", help="Get NLL for fixed point", nargs='?', const=1)
    parser.add_argument("-name", metavar="name", help="add suffix to output filename")
    parser.add_argument("-m", metavar="mode", help="Can choose to run only specific functions (grid, bestfit, etc.)")
    parser.add_argument("-batch", metavar="batchmode", help="crab or condor")
    parser.add_argument("-val", metavar="val", help="Starting value for the POI")
    parser.add_argument("--debug", metavar="debug", help="Activate code debug printouts", nargs='?', const=1)
    parser.add_argument("--freeze", metavar="freeze", help="Freeze other POIs", nargs='?', const=1)
    parser.add_argument('-P','--POI', metavar="POI", nargs='+', help='Define POI(s)', required=False) #Takes >=0 args
    parser.add_argument("--noworkspace", metavar="noworkspace", help="Don't recreate workspace", nargs='?', const=1)
    parser.add_argument("--onlyworkspace", metavar="onlyworkspace", help="Only create workspace", nargs='?', const=1)
    parser.add_argument("--dryrun", metavar="dryrun", help="Perform dry run (don't submit jobs)", nargs='?', const=1)
    parser.add_argument("-points", metavar="points", help="Number of points for grid scans")
    parser.add_argument('--mask', metavar="mask", nargs='+', help='Mask channels matching pattern', required=False) #Takes >=0 args
    parser.add_argument('--antimask', metavar="antimask", nargs='+', help='Mask channels NOT matching any pattern', required=False) #Takes >=0 args

    args = parser.parse_args()
    if args.sm: sm = True
    if args.d: datacard_path = args.d
    if args.v: verb = int(args.v)
    if args.exp: exp = True
    if args.name: name = args.name
    if args.val: startValue = args.val
    if args.dim == '2D' or args.dim == '2D': scan_dim = '2D'
    if args.fixed: fixedPointNLL = True
    if args.debug: debug = True
    if args.freeze: freeze = True
    if args.POI: POI = args.POI
    if args.noworkspace: createWS = 2
    if args.onlyworkspace: createWS = 1
    if args.m: mode = args.m
    if args.batch: batch = args.batch
    if args.dryrun: dryrun = '--dry-run'
    if args.points: points = args.points
    if args.mask: mask = args.mask
    if args.antimask: antimask = args.antimask #FIXME also SM

    fitter = EFTFit(opts) #Create EFTFit object

    if batch=='crab':
        print('ERROR: crab mode is not supported yet ! Use condor to submit jobs !')
        exit(1)

    PrintBanner() #Init printouts


# SM fit
# //--------------------------------------------
    if SM:
        if '.root' not in datacard_path and (createWS<2 or '.txt' in datacard_path): fitter.makeWorkspaceSM(datacard_path, verbosity=verb)
        if createWS==1: exit(1)

        if mode in ['','bestfit']: fitter.bestFitSM(params_POI=POI, exp=exp, verbosity=verb, mask=mask, antimask=antimask)
        elif mode is 'printbestfit': #Only print best fit results
            fitter.printBestFit(name='.EFT', params=POI)
            exit(1)

        if mode in ['','grid']:
            if scan_dim=='1D': fitter.gridScanSM(scan_params=[opts["SM_mu"]], points=100, exp=exp, verbosity=verb, batch=batch) #1D
            elif scan_dim=='2D': fitter.gridScanSM(scan_params=opts["SM_mus"], points=1000, exp=exp, verbosity=verb, batch=batch) #2D

# SMEFT fit
# //--------------------------------------------
    else:
        #-- Create Combine Workspace
        if '.root' not in datacard_path and (createWS<2 or '.txt' in datacard_path): fitter.makeWorkspaceEFT(datacard_path, verbosity=verb)
        if createWS==1: exit(1)
        if name == '': name = '.EFT' #Default

        #-- Maximum Likelihood Fit
        if mode in ['','bestfit']: fitter.bestFitEFT(params_POI=POI, exp=exp, verbosity=verb, name=name, startValue=startValue, fixedPointNLL=fixedPointNLL, freeze=freeze, mask=mask, antimask=antimask)
        elif mode is 'printbestfit': #Only print best fit results
            fitter.printBestFit(name='.EFT', params=POI)
            exit(1)

        #-- Grid Scan
        if not fixedPointNLL and mode in ['','grid']:
            if scan_dim=='1D':
                param_tmp = POI if len(POI) == 1 else [opts['wc']]
                points = points if points != -1 else 50
                fitter.gridScanEFT(name=name, scan_params=param_tmp, exp=exp, points=points, verbosity=verb, freeze=freeze, batch=batch, other=[dryrun,])
                # fitter.batchRetrieve1DScansEFT(basename=name, batch=batch, scan_wcs=param_tmp)
            elif scan_dim=='2D':
                param_tmp = POI if len(POI) == 2 else [opts['wcs_pairs']]
                points = points if points != -1 else 500
                fitter.gridScanEFT(scan_params=param_tmp, exp=exp, points=points, verbosity=verb, freeze=freeze, batch=batch, name=name, other=[dryrun,])
                # fitter.batchRetrieve2DScansEFT(basename=name, batch=batch, wc_pair=param_tmp, allPairs=False)

        #-- OTHERS
        # fitter.batch1DScanEFT() # Freeze other WCs
        # fitter.batchRetrieve1DScansEFT()
        # fitter.getBestValues1DEFT()
# //--------------------------------------------
