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
        self.wc_ranges_scan1D = opts['wc_ranges_scan1D']
        self.wc_ranges_scan2D = opts['wc_ranges_scan2D']
        self.wcs_pairs = opts['wcs_pairs']
        self.SM_mu = opts['SM_mu']
        self.SM_mus = opts['SM_mus']
        self.SMmu_ranges = opts['SMmu_ranges']

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

    def makeWorkspace(self, ws_output, SM=False, datacard='datacard.txt', verbosity=0, wcs=[]):
        ### Generates a workspace from a datacard and fit parameterization file ###
        logging.info(colors.fg.lightblue + "Creating workspace" + colors.reset)
        if not os.path.isfile(datacard):
            logging.error("Datacard does not exist!")
            sys.exit()
        CMSSW_BASE = os.getenv('CMSSW_BASE')

        if SM:
            args = ['text2workspace.py',datacard,'-P','HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel', '-o',ws]
            # Map signal strengths to processes in all bins
            for iproc,proc in enumerate(opts["processes"]):
                args.extend(['--PO', 'map=.*/'+proc+':'+opts["SM_mus"][iproc]+'[1,'+str(self.SMmu_ranges[opts["SM_mus"][iproc]][0])+','+str(self.SMmu_ranges[opts["SM_mus"][iproc]][1])+']'])
        else:
            args = ['text2workspace.py',datacard, '-P','EFTModel:eftmodel', '-o',ws,'--PO','fits=./Parameterization_EFT.npy']
            if len(wcs) != 0: #Useless ?
                for wc in wcs: args.extend(['--PO','wcs=',' '.join('{}'.format(wc) for wc in wcs)])

        if verbosity>0: args.extend(['-v', str(verbosity)])
        args.extend(['--channel-masks']) #Creates additional parameters allowing to later mask specific channels

        # Remove pre-existing WS
        sp.call(['rm',ws])

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

    def bestFit(self, datacard='./EFTWorkspace.root', SM=False, name='.EFT', params_POI=[], freeze=False, startValues=[], autoBounds=True, other=[], exp=False, verbosity=0, fixedPointNLL=False, mask=[], antimask=[], ntoys=-1):
        '''
        Perform a (multi-dim.) MLF to find the best fit value of the POI(s).

        NB: the error 'Looks like the last fit did not float this parameter' could arise e.g. if the bin names differ between file/datacard/EFT parametrization file. For some reason the parameter is not freely floating,  or not associated to any process, etc.
        '''

        logging.info(colors.fg.lightblue + "Enter function bestFit()\n" + colors.reset)

        if params_POI == []:
            if SM: params_POI = self.SM_mus
            else: params_POI = self.wcs
        if name == '':
            if SM: name = '.SM'
            else: name = '.EFT'

        #-- #Define channel masking regexp pattern, if any
        maskPattern = []; antimaskPattern = []
        if len(mask)>0: maskPattern=[','.join('rgx{{mask_.*{}.*}}=1'.format(chan) for chan in mask)] #Use '{{' to insert a litteral bracket, not a replacement field #More info on regexp meaning: https://regex101.com/
        if len(antimask)>0: antimaskPattern=['rgx{^mask_(?!.*('+'|'.join('{}'.format(chan) for chan in antimask)+')).*$}=1'] #Opposite: mask all channels NOT matching ANY 'chan'
        # if len(mask)>0: maskPattern=[','.join('rgx{{mask.*_{}_.*}}=1'.format(chan) for chan in mask)] #Use '{{' to insert a litteral bracket, not a replacement field #More info on regexp meaning: https://regex101.com/
        # if len(antimask)>0: antimaskPattern=['rgx{^mask_(?!.*('+'|'.join('{}'.format(chan) for chan in antimask)+')).*$}=1'] #Opposite: mask all channels NOT matching ANY 'chan'

        # CMSSW_BASE = os.getenv('CMSSW_BASE')
        # args=['combine','-d',datacard,'-M','MultiDimFit','--saveNLL','--saveFitResult','-H','AsymptoticLimits','--cminPoiOnlyFit']
        args=['combine','-d',datacard,'-M','MultiDimFit','--saveNLL','--saveFitResult','--do95','1','--robustFit','1']

        args.extend(['-n','{}'.format(name)])
        if fixedPointNLL:
            #args.extend(['--X-rtd','REMOVE_CONSTANT_ZERO_POINT=1']) #Access absolute NLL #Necessary/useful ?
            # if SM: args.extend(['--algo','fixed','--fixedPointPOIs','{}={}'.format(opts['SM_mu'],startValues)])
            # else: args.extend(['--algo','fixed','--fixedPointPOIs','{}={}'.format(opts['wc'],startValues)])
            args.extend(['--algo','fixed'])
            if SM:
                args.extend(['--fixedPointPOIs',','.join([','.join('{}={}'.format(mu,startValues[imu]) for imu,mu in enumerate(self.SM_mus))])])
            else:
                args.extend(['--fixedPointPOIs',','.join([','.join('{}={}'.format(wc,startValues[iwc]) for iwc,wc in enumerate(self.wcs))])])
        if params_POI:
            for param in params_POI: args.extend(['-P','{}'.format(param)])
        if SM: args.extend(['--setParameters',','.join([','.join('{}=1'.format(mu) for mu in self.SM_mus)]+maskPattern+antimaskPattern)]) #Set default values to 1
        else: args.extend(['--setParameters',','.join([','.join('{}=0'.format(poi) for poi in self.wcs)]+maskPattern+antimaskPattern)]) #Set default values to 0 #Mask channels if needed
        if freeze: #NB: looks like using option 'floatOtherPOIs' overrides 'freezeParameters'
            frozen_pois = []
            if SM:
	            frozen_pois = [par for par in self.SM_mus if par not in params_POI] #Define which WCs are frozen
                #print(colors.bg.red + "WARNING: preventing to use --freeze for SM... if needed, remove this protection !" + colors.reset)
            else: frozen_pois = [wc for wc in self.wcs if wc not in params_POI]
            if len(frozen_pois)>0: args.extend(['--freezeParameters',','.join('{}'.format(poi) for poi in frozen_pois)]) #Freeze other parameters
        else: args.extend(['--floatOtherPOIs','1']) #Float other parameters defined in the physics model
        #if autoBounds: args.extend(['--autoBoundsPOIs', '*']) #Auto adjust POI bounds if found close to boundary
        # if autoBounds: args.extend(['--autoBoundsPOIs=',','.join([','.join('{}'.format(poi) for poi in self.wcs)])]) #Auto adjust POI bounds if found close to boundary
        if exp: args.extend(['-t','-1']) #Assume MC expected (Asimov?)
        if verbosity>0: args.extend(['-v', str(verbosity)])
        if other: args.extend(other)
        check = True in (wc not in params_POI for wc in self.wcs_tracked)
        if check and not SM: args.extend(['--trackParameters',','.join(wc for wc in self.wcs_tracked if wc not in params_POI)]) #Save values of additional parameters (e.g. profiled nuisances)
        if SM: args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(mu,self.SMmu_ranges[mu][0],self.SMmu_ranges[mu][1]) for mu in self.SM_mus)]) #in params_POI
        else: args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(wc,self.wc_ranges[wc][0],self.wc_ranges[wc][1]) for wc in self.wcs)]) #in params_POI
        if ntoys>0: args.extend(['-t',str(ntoys)])

        if debug: print('args --> ', args)
        logging.info(colors.fg.purple + " ".join(args) + colors.reset)
        process = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE)
        with process.stdout,process.stderr:
            self.log_subprocess_output(process.stdout,'info')
            self.log_subprocess_output(process.stderr,'err')
        process.wait()
        logging.info(colors.fg.lightblue + "Done with bestFit !" + colors.reset)

        if not fixedPointNLL: self.printBestFit(name=name, params=params_POI, SM=SM)


  ####  #####  # #####      ####   ####    ##   #    #
 #    # #    # # #    #    #      #    #  #  #  ##   #
 #      #    # # #    #     ####  #      #    # # #  #
 #  ### #####  # #    #         # #      ###### #  # #
 #    # #   #  # #    #    #    # #    # #    # #   ##
  ####  #    # # #####      ####   ####  #    # #    #

    def gridScan(self, datacard='./EFTWorkspace.root', SM=False, name='', batch='', freeze=False, scan_params=['ctz'], startValuesString='', params_tracked=[], points=1000, exp=False, other=[], verbosity=0, collect=False, mask=[], antimask=[], ntoys=-1):

        ### Runs deltaNLL Scan in two parameters using CRAB or Condor ###
        logging.info(colors.fg.lightblue + "Enter function gridScan()\n" + colors.reset) #'You need to supply shell=True to execute the command through a shell interpreter. If you do that however, you can no longer supply a list as the first argument, because the arguments will get quoted then. Instead, specify the raw commandline as you want it to be passed to the shell'

        if collect==False: #Run grid scan

            if batch != '': #Delete pre-existing scan files with same name, to avoid unvoluntary mixing
                print('rm higgsCombine{}.POINTS*root'.format(name))
                sp.call('rm higgsCombine{}.POINTS*root'.format(name), shell=True, stdout=sp.PIPE, stderr=sp.PIPE)

            args = ['combineTool.py','-d',datacard,'-M','MultiDimFit','--algo','grid','--cminPreScan','--cminDefaultMinimizerStrategy=0']

            for wc in scan_params: args.extend(['-P', '{}'.format(wc)]) #Define signal strengths as POIs

            #-- #Define channel masking regexp pattern, if any
            maskPattern = []; antimaskPattern = []
            if len(mask)>0: maskPattern=[','.join('rgx{{mask_.*{}.*}}=1'.format(chan) for chan in mask)] #Use '{{' to insert a litteral bracket, not a replacement field #More info on regexp meaning: https://regex101.com/
            if len(antimask)>0: antimaskPattern=['rgx{^mask_(?!.*('+'|'.join('{}'.format(chan) for chan in antimask)+')).*$}=1'] #Opposite: mask all channels NOT matching ANY 'chan'
            # if len(mask)>0: maskPattern=[','.join('rgx{{mask.*_{}_.*}}=1'.format(chan) for chan in mask)] #Use '{{' to insert a litteral bracket, not a replacement field #More info on regexp meaning: https://regex101.com/
            # if len(antimask)>0: antimaskPattern=['rgx{^mask_(?!.*('+'|'.join('{}'.format(chan) for chan in antimask)+')).*$}=1'] #Opposite: mask all channels NOT matching ANY 'chan'

            if SM:
                args.extend(['--setParameters',','.join(['{}=1'.format(mu) for mu in scan_params]+maskPattern+antimaskPattern)]) #Set default values to 1
                args.extend(['--setParameterRanges', ':'.join(['{}={},{}'.format(mu,self.SMmu_ranges[mu][0],self.SMmu_ranges[mu][1]) for mu in self.SM_mus]+maskPattern+antimaskPattern)])
            else:
                args.extend(['--setParameters',','.join(['{}=0'.format(wc) for wc in self.wcs]+maskPattern+antimaskPattern)]) #Set default values to 1

                if freeze:
                    if len(scan_params)==1: #1D
                        args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(wc,self.wc_ranges_scan1D[wc][0],self.wc_ranges_scan1D[wc][1]) for wc in self.wcs)])
                    else: #n-D
                        args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(wc,self.wc_ranges_scan2D[wc][0],self.wc_ranges_scan2D[wc][1]) for wc in self.wcs)])
                else:
                    args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(wc,self.wc_ranges_scan2D[wc][0],self.wc_ranges_scan2D[wc][1]) for wc in self.wcs)])

                #Obsolete? -- want to adjust 1D ranges for grid scans only (not best fits), else loosing too many points !
                '''
                if freeze:
                    if len(scan_params)==1: #1D
                        min = -1; max = 1 #ctw
                        if scan_params[0]  == 'ctz': min = -1.5; max = 1
                        elif scan_params[0]  == 'cpq3': min = -2; max = 2
                        elif scan_params[0]  == 'cpqm': min = -4; max = 4
                        elif scan_params[0]  == 'cpt': min = -7; max = 7
                        args.extend([':{}={},{}'.format(scan_params[0],min,max)])
                        # args.extend(['--setParameterRanges', ':{}={},{}'.format(scan_params[0],min,max)])
                    else: #n-D
                        for iparam,param in enumerate(scan_params):
                            min = -2; max = 2 #For ctz, ctw, cpq3
                            if scan_params[iparam]  == 'cpqm': min = -20; max = 20
                            elif scan_params[iparam]  == 'cpt': min = -25; max = 20
                            args.extend([':{}={},{}'.format(scan_params[iparam],min,max)])
                            #args.extend(['--setParameterRanges', ':{}={},{}'.format(scan_params[iparam],min,max)])
                            # args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(wc,self.wc_ranges[wc][0],self.wc_ranges[wc][1]) for wc in self.wcs)]) #in params_POI
                '''

                # else:
                #     for iparam,param in enumerate(scan_params):
                #         min = -2; max = 2 #For ctz, ctw, cpq3
                #         if scan_params[iparam] in ['cpqm','cpt']: min = -15; max = 15 #For cpqm,cpt
                #         args.extend(['--setParameterRanges', ':{}={},{}'.format(scan_params[iparam],min,max)])
                #     #args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(wc,self.wc_ranges[wc][0],self.wc_ranges[wc][1]) for wc in self.wcs)])

                # args.extend(['--setParameterRanges', ':'.join('{}={},{}'.format(wc,self.wc_ranges_scan[wc][0],self.wc_ranges_scan[wc][1]) for wc in self.wcs)]) #in params_POI

            args.extend(['--points','{}'.format(points)])
            if name: args.extend(['-n','{}'.format(name)])
            check = True in (wc not in scan_params for wc in self.wcs_tracked)
            if check:
                tracked_pois = []
                if SM: tracked_pois = [par for par in self.SM_mus if par not in scan_params] #Define which params are frozen
                else: tracked_pois = [par for par in self.wcs if par not in scan_params]
                if len(tracked_pois)>0: args.extend(['--trackParameters',','.join([par for par in tracked_pois if par not in scan_params])]) #Save values of additional parameters (e.g. profiled nuisances) --> Can investigate their behjaviour while scanning a different POI
            # if startValuesString:   args.extend(['--setParameters',startValuesString])
            if freeze:
                frozen_pois = []
                if SM: frozen_pois = [par for par in self.SM_mus if par not in scan_params] #Define which params are frozen
                else: frozen_pois = [par for par in self.wcs if par not in scan_params]
                if len(frozen_pois)>0: args.extend(['--freezeParameters',','.join('{}'.format(par) for par in frozen_pois if par not in scan_params)]) #Freeze other parameters
            else: args.extend(['--floatOtherPOIs','1']) #Float other parameters defined in the physics model
            if exp:               args.extend(['-t -1'])
            if verbosity>0:           args.extend(['-v', str(verbosity)])
            if other:             args.extend(other)
            # args.extend(['--fastScan']) #No profiling (speed up) of nuisances, kept to best fit value
            if ntoys>0: args.extend(['-t',str(ntoys)])

            if batch=='crab': args.extend(['--job-mode','crab3','--task-name',name.replace('.',''),'--custom-crab','Utils/custom_crab.py','--split-points','50'])
            if batch=='condor': args.extend(['--job-mode','condor','--task-name',name.replace('.',''),'--split-points','50','--dry-run'])

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
                fitter.printInterval_fromScan(name=name, scan_params=scan_params, SM=SM) #Print exclusion range #Obsolete

        else: #Only collect/hadd grid scan outputs (produced via batch)
            outfilename = 'scan_'
            if len(scan_params)==1: outfilename+= scan_params[0]
            elif len(scan_params)==2: outfilename+= scan_params[0]+'_'+scan_params[1]
            else: outfilename+= '5D'
            outfilename+= '.root'
            print(colors.fg.orange + '-- Collecting grid scan output files from batch --> [{}] ...'.format(outfilename) + colors.reset)
            hadd_cmd = 'hadd -f {} higgsCombine{}.POINTS*.MultiDimFit.mH120.root'.format(outfilename,name)
            proc = sp.Popen(hadd_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
            print(colors.fg.orange + '... Done !' + colors.reset)

        return


 #####    ##   #####  ####  #    #
 #    #  #  #    #   #    # #    #
 #####  #    #   #   #      ######
 #    # ######   #   #      #    #
 #    # #    #   #   #    # #    #
 #####  #    #   #    ####  #    #

    '''
    def batch1DScanSM(self, name='.test', batch='', scan_params=[], points=300, freeze=False, other=[]):
        ### For each SM signal strength, run a 1D deltaNLL Scan.
        if not scan_params:
            scan_params = ['r_tzq']

        for param in scan_params:
            scanmax = 3
            if param=='r_ttH': scanmax = 6
            if param=='r_tllq': scanmax = 4
            self.gridScanSM('{}.{}'.format(name,param), batch, [param], self.systematics+[params for params in scan_params if params != param], points, freeze, ['--setParameterRanges','{}=0,{}'.format(param,scanmax)]+other)

    def batchRetrieve1DScansSM(self, name='.test', batch='crab'):
        ### For each wc, retrieves finished 1D deltaNLL crab jobs, extracts, and hadd's into a single file ###
        for param in ['r_ttll','r_ttlnu','r_ttH','r_tllq']:
            self.retrieveGridScan('{}.{}'.format(name,param),batch)

    def batch1DScanEFT(self, name='', batch='crab', freeze=False, scan_wcs=[], points=300, other=[]):
        ### For each wc, run a 1D deltaNLL Scan.
        if not scan_wcs:
            scan_wcs = self.wcs

        for wc in scan_wcs:
            self.gridScanEFT('{}.{}'.format(name,wc), batch, freeze, [wc], [wcs for wcs in self.wcs if wcs != wc], points, other)

    def batch2DScanEFT(self, name='EFT.gridScan', batch='crab', freeze=False, points=90000, allPairs=False, other=[]):
        ### For pairs of wcs, runs deltaNLL Scan in two wcs using CRAB or Condor ###

        # Use EVERY combination of wcs
        if allPairs:
            scan_wcs = self.wcs

            for wcs in itertools.combinations(scan_wcs,2):
                wcs_tracked = [wc for wc in self.wcs if wc not in wcs]
                #print pois, wcs_tracked
                self.gridScanEFT(name='{}.{}{}'.format(name,wcs[0],wcs[1]), batch=batch, freeze=freeze, scan_params=list(wcs), params_tracked=wcs_tracked, points=points, other=other)

        # Use each wc only once
        if not allPairs:
            scan_wcs = [('ctw','ctG'),('ctz','ctG'),('ctp','ctG'),('cpQM','ctG'),('cbW','ctG'),('cpQ3','ctG'),('cptb','ctG'),('cpt','ctG'),('cQl3i','ctG'),('cQlMi','ctG'),('cQei','ctG'),('ctli','ctG'),('ctei','ctG'),('ctlSi','ctG'),('ctlTi','ctG')]
            #pairs from AN
            scan_wcs = [('cQlMi','cQei'),('cpQ3','cbW'),('cptb','cQl3i'),('ctG','cpQM'),('ctz','ctw'),('ctei','ctlTi'),('ctlSi','ctli'),('ctp','cpt')]

            for wcs in scan_wcs:
                wcs_tracked = [wc for wc in self.wcs if wc not in wcs]
                #print pois, wcs_tracked
                self.gridScanEFT(name='{}.{}{}'.format(name,wcs[0],wcs[1]), batch=batch, freeze=freeze, scan_params=list(wcs), params_tracked=wcs_tracked, points=points, other=other)

    def batch3DScanEFT(self, name='EFT.gridScan', batch='crab', freeze=False, points=27000000, allPairs=False, other=[], wc_triplet=[]):
        ### For pairs of wcs, runs deltaNLL Scan in two wcs using CRAB or Condor ###

        # Use EVERY combination of wcs
        if allPairs:
            scan_wcs = self.wcs

            for wcs in itertools.combinations(scan_wcs,2):
                wcs_tracked = [wc for wc in self.wcs if wc not in wcs]
                #print pois, wcs_tracked
                self.gridScanEFT(name='{}.{}{}{}'.format(name,wcs[0],wcs[1],wcs[2]), batch=batch, freeze=freeze, scan_params=list(wcs), params_tracked=wcs_tracked, points=points, other=other)

        # Use each wc only once
        if not allPairs:
            scan_wcs = [('ctz','ctp','cpt')]
            if len(wc_triplet)>0: scan_wcs = wc_triplet
            #scan_wcs = [('ctz','ctw'),('ctp','cpt'),('ctlSi','ctli'),('cptb','cQl3i'),('ctG','cpQM'),('ctei','ctlTi'),('cQlMi','cQei'),('cpQ3','cbW')]
            #scan_wcs = [('ctw','ctG'),('ctz','ctG'),('ctp','ctG'),('cpQM','ctG'),('cbW','ctG'),('cpQ3','ctG'),('cptb','ctG'),('cpt','ctG'),('cQl3i','ctG'),('cQlMi','ctG'),('cQei','ctG'),('ctli','ctG'),('ctei','ctG'),('ctlSi','ctG'),('ctlTi','ctG')]

            for wcs in scan_wcs:
                wcs_tracked = [wc for wc in self.wcs if wc not in wcs]
                #print pois, wcs_tracked
                self.gridScanEFT(name='{}.{}{}{}'.format(name,wcs[0],wcs[1],wcs[2]), batch=batch, freeze=freeze, scan_params=list(wcs), params_tracked=wcs_tracked, points=points, other=other)

    def batchResubmit1DScansEFT(self, name='EFT.gridScan', scan_wcs=[]):
        ### For each wc, attempt to resubmit failed CRAB jobs ###
        if not scan_wcs:
            scan_wcs = self.wcs

        for wc in scan_wcs:
            process = sp.Popen(['crab','resubmit','crab_'+name.replace('.','')+wc], stdout=sp.PIPE, stderr=sp.PIPE)
            with process.stdout,process.stderr:
                self.log_subprocess_output(process.stdout,'info')
                self.log_subprocess_output(process.stderr,'err')
            process.wait()

    def batchResubmit2DScansEFT(self, name='EFT.gridScan', allPairs=False):
        ### For pairs of wcs, attempt to resubmit failed CRAB jobs ###

        # Use EVERY combination of wcs
        if allPairs:
            scan_wcs = self.wcs

            for wcs in itertools.combinations(scan_wcs,2):
                process = sp.Popen(['crab','resubmit','crab_'+name.replace('.','')+wcs[0]+wcs[1]], stdout=sp.PIPE, stderr=sp.PIPE)
                with process.stdout,process.stderr:
                    self.log_subprocess_output(process.stdout,'info')
                    self.log_subprocess_output(process.stderr,'err')
                process.wait()

        # Use each wc only once
        if not allPairs:
            scan_wcs = [('ctz','ctw'),('ctp','cpt'),('ctlSi','ctli'),('cptb','cQl3i'),('ctG','cpQM'),('ctei','ctlTi'),('cQlMi','cQei'),('cpQ3','cbW')]
            scan_wcs = [('cQlMi','cQei'),('cpQ3','cbW'),('cptb','cQl3i'),('ctG','cpQM'),('ctz','ctw'),('ctei','ctlTi'),('ctlSi','ctli'),('ctp','cpt')]

            for wcs in scan_wcs:
                process = sp.Popen(['crab','resubmit','crab_'+name.replace('.','')+wcs[0]+wcs[1]], stdout=sp.PIPE, stderr=sp.PIPE)
                with process.stdout,process.stderr:
                    self.log_subprocess_output(process.stdout,'info')
                    self.log_subprocess_output(process.stderr,'err')
                process.wait()
    '''


  ####  ###### #####    #####  ######  ####  #####    #    #   ##   #
 #    # #        #      #    # #      #        #      #    #  #  #  #
 #      #####    #      #####  #####   ####    #      #    # #    # #
 #  ### #        #      #    # #           #   #      #    # ###### #
 #    # #        #      #    # #      #    #   #       #  #  #    # #
  ####  ######   #      #####  ######  ####    #        ##   #    # ######

    '''
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

    def getBestValues1DEFT(self, name="EFT", wcs=[]):
        ### Gets values of WCs for grid scan point with best deltaNLL ###
        if not wcs:
            wcs = self.wcs

        startValues = []

        for wc in wcs:

            bestDeltaNLL=1000000;
            bestEntry=-1;

            fitFile = './higgsCombine{}.{}.MultiDimFit.root'.format(name,wc)
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
    '''


 #####  ###### ##### #####  # ###### #    # ######     ####  #####  # #####
 #    # #        #   #    # # #      #    # #         #    # #    # # #    #
 #    # #####    #   #    # # #####  #    # #####     #      #    # # #    #
 #####  #        #   #####  # #      #    # #         #  ### #####  # #    #
 #   #  #        #   #   #  # #       #  #  #         #    # #   #  # #    #
 #    # ######   #   #    # # ######   ##   ######     ####  #    # # #####

    '''
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


    def batchRetrieve1DScansEFT(self, name='', batch='crab', scan_wcs=[]):
        ### For each wc, retrieves finished 1D deltaNLL grid jobs, extracts, and hadd's into a single file ###
        if not scan_wcs:
            scan_wcs = self.wcs

        for wc in scan_wcs:
            self.retrieveGridScan('{}.{}'.format(name,wc),batch)

    def batchRetrieve2DScansEFT(self, wc_pair=[], name='', batch='crab', allPairs=False):
        ### For pairs of wcs, retrieves finished grid jobs, extracts, and hadd's into a single file ###

        # Use EVERY combination of wcs
        if allPairs:
            wc_pair = self.wcs
            for wcs in itertools.combinations(wc_pair,2):
                self.retrieveGridScan(name='{}.{}{}'.format(name,wcs[0],wcs[1]), batch=batch)

        # Consider a single pair of WCs
        else:
            if wc_pair == []: wc_pair = self.wcs_pairs
            print wc_pair
            print '{}.{}'.format(batch, name)
            self.retrieveGridScan(name=name, batch=batch)
    '''


 #####  ###### #####  #    #  ####  ##### #  ####  #    #
 #    # #      #    # #    # #    #   #   # #    # ##   #
 #    # #####  #    # #    # #        #   # #    # # #  #
 #####  #      #    # #    # #        #   # #    # #  # #
 #   #  #      #    # #    # #    #   #   # #    # #   ##
 #    # ###### #####   ####   ####    #   #  ####  #    #

#-- FIXME SHOULD USE THAT TO PRODUCE CORRECT SCAN ? EVEN 1D ?

    '''
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
            logging.error("File higgsCombine{}.MultiDimFit.root does not exist!".format(name))
            return

        rootFile = []
        rootFile = ROOT.TFile.Open('higgsCombine{}.MultiDimFit.mH120.root'.format(name))
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
    '''


  ####   ####  #    # #####    ##   #####  ######    ###### # #####  ####
 #    # #    # ##  ## #    #  #  #  #    # #         #      #   #   #
 #      #    # # ## # #    # #    # #    # #####     #####  #   #    ####
 #      #    # #    # #####  ###### #####  #         #      #   #        #
 #    # #    # #    # #      #    # #   #  #         #      #   #   #    #
  ####   ####  #    # #      #    # #    # ######    #      #   #    ####

    '''
    def compareFitsEFT(self,name='EFT.SM.Float'):
        ### Compare results of different 1D EFT scans ###
        tfiles = {}
        limits = {}
        bestFits = {} # Nested dict; bestFit of key1 according to key2
        # First get all scan files
        for wc in self.wcs:
            tfiles[wc] = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.root'.format(name+'.'+wc))
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
    '''


 #####  #####  # #    # #####
 #    # #    # # ##   #   #
 #    # #    # # # #  #   #
 #####  #####  # #  # #   #
 #      #   #  # #   ##   #
 #      #    # # #    #   #

    def printBestFit(self, name='', params=[], SM=False):
        '''
        Print best fit values and uncertainties of parameters.
        Read multidimfit.root file (output from max. likelihood fit)

        Example to inspect file via command line:
        root multidimfit.EFT.root
        a = fit_mdf->floatParsFinal().find("ctz") #POI
        b = (RooAbsReal*) a; b->Print() #Printout
        rf = dynamic_cast<RooRealVar*>(a);
        rf->getErrorLo(); f->getErrorHi(); #Get +-1sigma uncertainties
        rf->getMin("err68")-roorealvar.getVal(); rf->getMax("err68")-roorealvar.getVal(); #Idem
        rf->getMin("err68"); rf->getMax("err68"); #Get 68%CL INTERVAL
        rf->getMin("err95"); rf->getMax("err95"); #Get 95%CL INTERVAL (if used option '--do95')

        Example in python:
        import ROOT
        fit_file = ROOT.TFile.Open('multidimfit.EFT.root')
        fit = fit_file.Get('fit_mdf')
        roorealvar = fit.floatParsFinal().find('ctz')
        value = round(roorealvar.getVal(),2) #Best fit value
        value = round(roorealvar.getMin('err68'),2) #Lower 68% interval bound
        '''

        logging.info(colors.fg.lightblue + "\nEnter function printBestFit()" + colors.reset)

        if len(params)==0: print('params is empty... Return !')
        fit_array = []

        for param in params:

            filepath = "multidimfit{}.root".format(name) #Default naming convention
            if name=='': #Convention: look for WC-specific file (hardcoded naming convention)
                suffix = 'Obs' #'Exp' / 'Obs'
                filepath = "multidimfit.{}{}.root".format(param,suffix)

            if not os.path.exists(filepath):
                logging.error(colors.fg.red + "File {} does not exist!".format(filepath) + colors.reset)
                continue

            logging.info("Obtaining result of fit: {}".format(filepath)) #NB: no '.' before 'name', since this file is produced automatically by combine (not by this code)
            fit_file = ROOT.TFile.Open(filepath)
            fit = fit_file.Get('fit_mdf')

            roorealvar = fit.floatParsFinal().find(param)
            if not roorealvar:
                print('ERROR: RooRealVar not found for param:', param)
                continue

            value = round(roorealvar.getVal(),3)

            #-- Get 68% CL interval
            err_low = round(roorealvar.getMin('err68'),3); err_high = round(roorealvar.getMax('err68'),3)
            # err_low = round(roorealvar.getErrorLo(),3) #Returns error directly, no need to substract bestFitVal
            # err_high = round(roorealvar.getErrorHi(),3) #Returns error directly, no need to substract bestFitVal

            #-- Get 95% CL interval
            err_low_95 = -9; err_high_95 = -9
            if roorealvar.hasRange('err95'): #If 95% CL errors available (using --do95 1 --robustFit 1 options)
                err_low_95 = round(roorealvar.getMin('err95'),3)
                err_high_95 = round(roorealvar.getMax('err95'),3)

            if SM == True: #DIFFERENT CONVENTION FOR SM RESULTS --> Instead of quoting INTERVALS (like for EFT results), quote UNCERTAINTIES directly !

                #-- To convert interval->uncertainty, simply substract best fit value
                err_low-= value
                err_high-= value
                err_low_95-= value
                err_high_95-= value

            fit_array.append((param,value,err_low,err_high,err_low_95,err_high_95))

        logging.info('\n' + colors.fg.orange + "Fit result:" + colors.reset)
        if SM==True: logging.info(colors.fg.orange + "Param | Best Fit Value | [68% ERROR] | [95% ERROR]" + colors.reset)
        else: logging.info(colors.fg.orange + "Param | Best Fit Value | [68% interval] | [95% interval]" + colors.reset)

        for row in fit_array:
            logging.info(colors.fg.orange + row[0] + ' | ' + str(row[1]) + ' | ' + "[" + str(row[2]) + ";" + str(row[3]) + "] | [" + str(row[4]) + ';' + str(row[5]) + ']' + colors.reset + '\n')
            #logging.debug("{} {} +/- {}".format(row[0],row[1],row[2]))

        return


    def printInterval_fromScan(self, name='.EFT', scan_params=[], SM=False):
        '''
        Print best fit values and uncertainties of parameters.
        Read higgsCombine.root file (output from grid scan)
        '''

        logging.info(colors.fg.lightblue + "Enter function printInterval_fromScan()\n" + colors.reset)

        if not scan_params:
            if SM: scan_params = [self.wc]
            else: scan_params = [self.SM_mu]

        ROOT.gROOT.SetBatch(True)
        fit_array = []
        canvas = ROOT.TCanvas()
        for param in scan_params:

            canvas.Clear()
            logging.info("Obtaining result of scan: higgsCombine{}.MultiDimFit.mH120.root for POI {}".format(name,param))

            #-- Get scan TTree
            rootFile = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.mH120.root'.format(name))
            limitTree = rootFile.Get('limit')

            #-- Use CombineTool utils (see: https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools/python/plotting.py)
            graph = plot.TGraphFromTree(limitTree, param, '2*deltaNLL', 'quantileExpected > -1.5')

            yvals = [1., 3.84] #1sigma, 95%CL intervals
            #func, crossings, val, val_2sig, cross_1sig, cross_2sig, other_1sig, other_2sig = BuildScan(graph, ROOT.kBlack, yvals)
            main_scan = BuildScan(graph, ROOT.kBlack, yvals)

            crossings = main_scan['crossings'][yvals[1]][0]
            if crossings['valid_lo'] == True and crossings['valid_hi'] == True:
                print(colors.fg.orange + ": 95% CL interval: [" + str(crossings['lo']) + ", " + str(crossings['hi']) + "]" + colors.reset)
            else: print('Error: invalid crossing X-values...')

        return


######## ########  ######  ######## #### ##    ##  ######
   ##    ##       ##    ##    ##     ##  ###   ## ##    ##
   ##    ##       ##          ##     ##  ####  ## ##
   ##    ######    ######     ##     ##  ## ## ## ##   ####
   ##    ##             ##    ##     ##  ##  #### ##    ##
   ##    ##       ##    ##    ##     ##  ##   ### ##    ##
   ##    ########  ######     ##    #### ##    ##  ######

    def Testing(self, name):

        logging.info("Obtaining result of scan: higgsCombine.{}.MultiDimFit.mH120.root".format(name))

        #-- Get scan TTree
        rootFile = ROOT.TFile.Open('./higgsCombine.{}.MultiDimFit.mH120.root'.format(name))
        limitTree = rootFile.Get('limit')
        list_br = limitTree.GetListOfBranches()
        # print(len(list_br))

        # for ientry in range(limitTree.GetEntries()): #Loop on entries
        for ientry in range(1): #Loop just once, and compare entries 0 and 1
            # limitTree.GetEntry(ientry)
            # print('=== Entry', ientry, '===')
            for ibr in range(len(list_br)):
                br_name = list_br.At(ibr).GetName()
                # print(br_name, ' = ', limitTree.GetLeaf(br_name).GetValue(ientry))
                limitTree.GetEntry(0)
                val0 = limitTree.GetLeaf(br_name).GetValue(0)
                limitTree.GetEntry(1)
                val1 = limitTree.GetLeaf(br_name).GetValue(0)
                if val0 != 0:
                    if (val1/val0>1.5 or val1/val0<0.5):
                        print(colors.fg.red + br_name + ' = ' + str(val0) + ' / ' + str(val1) + colors.reset)
                    elif (val1/val0>1.1 or val1/val0<0.9):
                        print(colors.fg.orange + br_name + ' = ' + str(val0) + ' / ' + str(val1) + colors.reset)
                    else: print(br_name, ' = ', val0, ' / ', val1)
                else:
                    print(br_name, ' = ', val0, ' / ', val1)
                # print(br_name, ' = ', limitTree.GetLeaf(br_name).GetValue(0), ' / ', limitTree.GetLeaf(br_name).GetValue(1))

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

# User options -- Default values
# //--------------------------------------------
    SM = False #True <-> consider SM scenario (rather than SMEFT)
    datacard_path = './datacard.txt'
    exp = False #True <-> Asimov a-priori expected; False <-> observed
    scan_dim = '1D' #'1d' <-> 1D scan (default); '2d' <-> 2D scan
    verb=0
    name = '' #Suffix added to output filenames
    startValues = [] #Starting value(s) for the POI(s)
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
    only_collect_outputs = False #True <-> only collect/hadd output files previously produced e.g. on the grid
    wcs = [] #For EFTModel... useless ?
    ws = 'EFTWorkspace.root' #Name for output workspace
    ntoys = -1

# Set up the command line arguments
# //--------------------------------------------
    parser = argparse.ArgumentParser(description='Perform SM and EFT fits using custom Physics Model')
    parser.add_argument("-d", metavar="datacard path", help="Path to the txt datacard (to create RooWorkspace)")
    parser.add_argument("-v", metavar="Combine verbosity level", help="Set combine output verbosity")
    parser.add_argument("--sm", metavar="SM", help="Consider SM scenario (rather than SMEFT)", nargs='?', const=1)
    parser.add_argument("-dim", metavar="dim", help="1D or 2D scan")
    parser.add_argument('--exp', help='Use MC predictions only (no data)', nargs='?', const=1) #nargs='?' <-> 0 or 1 arg (default value is const=1)
    parser.add_argument("--fixed", metavar="fixed", help="Get NLL for fixed point", nargs='?', const=1)
    parser.add_argument('-val', metavar="val", nargs='+', help='Starting value(s) for the POI(s)', required=False) #Takes >=0 args
    parser.add_argument("-name", metavar="name", help="add suffix to output filename")
    parser.add_argument("-m", metavar="mode", help="Can choose to run only specific functions (grid, bestfit, etc.)")
    parser.add_argument("--batch", metavar="batchmode", help="(crab or condor -- only condor for now)", nargs='?', const=1)
    parser.add_argument("--debug", metavar="debug", help="Activate code debug printouts", nargs='?', const=1)
    parser.add_argument("--freeze", metavar="freeze", help="Freeze other POIs", nargs='?', const=1)
    parser.add_argument('-P','--POI', metavar="POI", nargs='+', help='Define POI(s)', required=False) #Takes >=0 args
    parser.add_argument("--noworkspace", metavar="noworkspace", help="Don't recreate workspace", nargs='?', const=1)
    parser.add_argument("--onlyworkspace", metavar="onlyworkspace", help="Only create workspace", nargs='?', const=1)
    parser.add_argument("--dryrun", metavar="dryrun", help="Perform dry run (don't submit jobs)", nargs='?', const=1)
    parser.add_argument("-points", metavar="points", help="Number of points for grid scans")
    parser.add_argument('--mask', metavar="mask", nargs='+', help='Mask channels matching pattern', required=False) #Takes >=0 args
    parser.add_argument('--antimask', metavar="antimask", nargs='+', help='Mask channels NOT matching any pattern', required=False) #Takes >=0 args
    parser.add_argument("--collect", metavar="collect", help="Impact plot: only collect output files previously produced e.g. on the grid", nargs='?', const=1)
    parser.add_argument('--wcs', metavar="wcs", nargs='+', help='testing', required=False)
    parser.add_argument("-o", metavar="output workspace name", help="Name for output workspace")
    parser.add_argument("-t", metavar="number of MC toys", help="Number of MC toys")

    args = parser.parse_args()
    if args.sm: SM = True
    if args.d: datacard_path = args.d
    if args.v: verb = int(args.v)
    if args.exp: exp = True
    if args.name: name = args.name
    if args.val: startValues = args.val
    if args.dim == '2D': scan_dim = '2D'
    if args.fixed: fixedPointNLL = True
    if args.debug: debug = True
    if args.freeze: freeze = True
    if args.POI:
        POI = args.POI
        if len(POI)==2: scan_dim = '2D'
    if args.noworkspace: createWS = 2
    if args.onlyworkspace: createWS = 1
    if args.m: mode = args.m
    if args.batch: batch = 'condor' #Only use condor
    if args.dryrun: dryrun = '--dry-run'
    if args.points: points = args.points
    if args.mask: mask = args.mask
    if args.antimask: antimask = args.antimask #-- implement this also for SM fits
    if args.collect: only_collect_outputs = True
    if args.wcs: wcs = args.wcs
    if args.o: ws = args.o
    if args.t and not exp: ntoys = args.t

    fitter = EFTFit(opts) #Create EFTFit object

    if batch=='crab':
        print('ERROR: crab mode is not supported yet ! Use condor to submit jobs !')
        exit(1)

    if mode not in ['', 'bestfit', 'printbestfit', 'scan', 'grid', 'other']:
        print('ERROR: mode not recognized !')
        exit(1)

    if name != '' and name[0] != '.': name = '.{}'.format(name) #Convention: add a '.' before custom name

    PrintBanner() #Init printouts


# SM fit
# //--------------------------------------------
    if SM:
        if ws=='EFTWorkspace.root': ws = 'SMWorkspace.root' #SM name
        if '.root' not in datacard_path and (createWS<2 or '.txt' in datacard_path) and only_collect_outputs == False and mode != 'other': fitter.makeWorkspace(SM=SM, datacard=datacard_path, verbosity=verb, ws_output=ws)
        if createWS==1: exit(1)
        if name == '': name = '.SM' #Default
        if '.txt' in datacard_path: datacard_path = ws #If WS was created, make sure to update path

        if mode in ['','bestfit']: fitter.bestFit(datacard=datacard_path, SM=SM, params_POI=POI, exp=exp, verbosity=verb, name=name, startValues=startValues, fixedPointNLL=fixedPointNLL, freeze=freeze, mask=mask, antimask=antimask, ntoys=ntoys)
        elif mode is 'printbestfit': #Only print best fit results
            fitter.printBestFit(name=name, params=POI, SM=SM)
            exit(1)

        if mode in ['','grid','scan']:
            if scan_dim=='1D': fitter.gridScan(datacard=datacard_path, SM=SM, name=name, scan_params=[opts["SM_mu"]], points=50, exp=exp, verbosity=verb, batch=batch, collect=only_collect_outputs, mask=mask, antimask=antimask, ntoys=ntoys) #1D
            elif scan_dim=='2D': fitter.gridScan(datacard=datacard_path, SM=SM, name=name, scan_params=opts["SM_mus"], points=50*50, exp=exp, verbosity=verb, batch=batch, collect=only_collect_outputs, mask=mask, antimask=antimask, ntoys=ntoys) #2D

# SMEFT fit
# //--------------------------------------------
    else:
        #-- Create Combine Workspace
        if '.root' not in datacard_path and (createWS<2 or '.txt' in datacard_path) and only_collect_outputs == False and mode != 'other': fitter.makeWorkspace(SM=SM, datacard=datacard_path, verbosity=verb, wcs=wcs, ws_output=ws)
        if createWS==1: exit(1)
        if name == '': name = '.EFT' #Default
        if '.txt' in datacard_path: datacard_path = ws #If WS was created, make sure to update path

        #-- Maximum Likelihood Fit
        if mode in ['','bestfit']: fitter.bestFit(datacard=datacard_path, SM=SM, params_POI=POI, exp=exp, verbosity=verb, name=name, startValues=startValues, fixedPointNLL=fixedPointNLL, freeze=freeze, mask=mask, antimask=antimask, ntoys=ntoys)
        elif mode == 'printbestfit': #Only print best fit results
            fitter.printBestFit(name=name, params=POI, SM=SM)
            exit(1)

        #-- Grid Scan
        if not fixedPointNLL and mode in ['','grid','scan']:
            if scan_dim=='1D':
                param_tmp = POI if len(POI) == 1 else [opts['wc']]
                points = points if points != -1 else 50
                fitter.gridScan(datacard=datacard_path, SM=SM, name=name, scan_params=param_tmp, exp=exp, points=points, verbosity=verb, freeze=freeze, batch=batch, other=[dryrun,], collect=only_collect_outputs, mask=mask, antimask=antimask, ntoys=ntoys)

            elif scan_dim=='2D':
                param_tmp = POI if len(POI) == 2 else [opts['wcs_pairs']]
                points = points if points != -1 else 50*50 # 2500 pts
                fitter.gridScan(datacard=datacard_path, SM=SM, name=name, scan_params=param_tmp, exp=exp, points=points, verbosity=verb, freeze=freeze, batch=batch, other=[dryrun,], collect=only_collect_outputs, mask=mask, antimask=antimask, ntoys=ntoys)

        #-- OTHERS
        if mode == 'other':
            fitter.printBestFit(name='', params=opts['wcs'], SM=SM) #Print intervals for ALL wcs (if files found) #Keep name='' for naming conventions
            fitter.Testing('test')

# //--------------------------------------------
