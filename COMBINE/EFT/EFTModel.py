#-- Define custom Physics Model
# Adapted from: https://github.com/cms-govner/EFTFit
# See: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part2/physicsmodels/

import numpy as np
import ROOT
import pprint
from Utils.ColoredPrintout import colors
from settings import opts #Custom dictionnary of settings
ROOT.gSystem.Load('/afs/cern.ch/work/n/ntonon/private/Combine/CMSSW_10_2_13/src/EFTAnalysis/myLib.so') #Library for custom classes WCPoint, WCFit, TH1EFT

#Based on 'Quadratic' model: HiggsAnalysis.CombinedLimit.QuadraticScaling
from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel

verbose = 0 #(Dis)activate printouts for this code

# //--------------------------------------------
# //--------------------------------------------


######## ######## ######## ##     ##  #######  ########  ######## ##
##       ##          ##    ###   ### ##     ## ##     ## ##       ##
##       ##          ##    #### #### ##     ## ##     ## ##       ##
######   ######      ##    ## ### ## ##     ## ##     ## ######   ##
##       ##          ##    ##     ## ##     ## ##     ## ##       ##
##       ##          ##    ##     ## ##     ## ##     ## ##       ##
######## ##          ##    ##     ##  #######  ########  ######## ########

class EFTModel(PhysicsModel):
    """
    Apply process scaling due to EFT operators.
    This class takes a dictionary of quadratic fits describing how processes are
    scaled as a function of an EFT operator's Wilson coefficient and adds it to
    the workspace. For an example coefficient x, dictionary values should have
    the form `(a, b, c)` where `xsec_NP(x) / xsec_SM = a + bx + cx^2`.

    To produce an example dictionary, for coefficient `cuW`:
    >>> import numpy as np
    >>> scales = {'cuW': {'ttZ': (1, 0.322778, 653.371), 'ttW': (1, 1.20998, 205.528)}}
    >>> np.save('scales.npy', scales)

    Example for running:
    text2workspace.py ttV.txt -P EFTModel:eftmodel --PO scaling=scales.npy --PO process=ttZ --PO process=ttW --PO coefficient=cuW -o cuW.root
    combine -M MultiDimFit cuW.root --setParameterRanges=-4,4
    """


  ####  #####  ##### #  ####  #    #  ####
 #    # #    #   #   # #    # ##   # #
 #    # #    #   #   # #    # # #  #  ####
 #    # #####    #   # #    # #  # #      #
 #    # #        #   # #    # #   ## #    #
  ####  #        #   #  ####  #    #  ####

    def setPhysicsOptions(self, options):

        self.fits = None # File containing WC parameterizations of each process+bin *with events*!
        wcs_override = [] # WCs specified by arguments
        self.procbins = [] # Process+bin combinations (tuple) that we have events for
        procbin_override = [] # Process+bin combinations (tuple) specified by arguments

        # Option specified when creating workspace, e.g.: '... EFTModel:eftmodel --PO fits=./xxx.npy'
        for option, value in [x.split('=') for x in options]:
            if option == 'fits': # .npy fit file created with DumpEFTParametrization.py
                self.fits = value
            elif option == 'wcs': # Override to fit only a subset of WCs
                wcs_override = value.split(',')
            elif option == 'procbins': # Override to fit only a subset of proc+category combinations
                procbin_override = value.split(',')
            else:
                print "Unknown option", option

        #If procbins are specified, only use subset that we have fits for.
        #Otherwise, use all of the process+bin combinations that we have fits for.
        fits = np.load(self.fits)[()]
        self.procbins.extend(fits.keys())
        if len(wcs_override)>0: self.wcs = np.intersect1d(self.wcs,wcs_override)
        if len(procbin_override)>0: self.procbins = np.intersect1d(self.procbins,procbins_override)

        if verbose: print(colors.fg.orange + 'procbins: ' + str(self.procbins) + colors.reset)


  ####  ###### ##### #    # #####
 #      #        #   #    # #    #
  ####  #####    #   #    # #    #
      # #        #   #    # #####
 #    # #        #   #    # #
  ####  ######   #    ####  #

    def setup(self):

        print(colors.fg.lightblue + "Setting up fits..." + colors.reset)
        fits = np.load(self.fits)[()]
        for procbin in self.procbins:
            name = 'r_{0}_{1}'.format(procbin[0],procbin[1])
            procbin_name = '_'.join(procbin)
            if verbose:
                print('name', name)
                print('procbin_name', procbin_name)

            if not self.modelBuilder.out.function(name): #If r_proc_cat_bin not yet setup
                # Initialize function pieces
                constant = '{}'.format(fits[procbin][(SM_name,SM_name)]) # constant term (should be 1)
                lin_name = procbin_name+"_L" # Name of linear function
                lin_term = [] # Linear term
                lin_args = [] # List of wcs in linear term
                quartic_names = [procbin_name+"_Q"+str(idx) for idx,wc in enumerate(self.wcs)] # Names of quadratic functions
                quartic_terms = [[] for wc in self.wcs] # Quartic terms, but split into chunks
                quartic_args = [[] for wc in self.wcs] # List of wcs in quartic terms
                fit_terms = [constant] # List of fit terms

                # Fill function pieces
                for idx,wc1 in enumerate(self.wcs):
                    if abs(fits[procbin][(SM_name,wc1)]) >= 0.001:
                        lin_term.append('{0}*{1}'.format(round(fits[procbin][(SM_name,wc1)],4),wc1))
                        lin_args.append(wc1)
                    for idy,wc2 in enumerate(self.wcs):
                        # print('wc1', wc1, 'wc2', wc2)
                        if (wc1,wc2) not in fits[procbin]:
                            print(colors.fg.red + "ERROR: WCs " + str((wc1,wc2)) + " not found in EFT parameterization... Did you run [DumpEFTParametrization.py] on the proper template file ?" + colors.reset); exit(1)
                        elif (idy >= idx) and (abs(fits[procbin][(wc1,wc2)]) >= 0.001):
                            quartic_terms[idx].append('{0}*{1}*{2}'.format(round(fits[procbin][(wc1,wc2)],4),wc1,wc2))
                            quartic_args[idx].extend([wc1,wc2])

                if verbose: print('lin_term', lin_term)
                if verbose: print('quartic_terms', quartic_terms)

                # Compile linear function for combine
                if lin_term:
                    lin_expr = "expr::{lin_name}('{lin_term}',{lin_args})".format(lin_name=lin_name,lin_term="+".join(lin_term),lin_args=",".join(lin_args))
                    lin_func = self.modelBuilder.factory_(lin_expr)
                    self.modelBuilder.out._import(lin_func)
                    fit_terms.append(lin_name)

                # Compile quartic functions separately first
                for idx,fn in enumerate(quartic_terms):
                    if not fn: continue # Skip empty quartic functions
                    quartic_expr = "expr::{quartic_names}('{quartic_terms}',{quartic_args})".format(quartic_names=quartic_names[idx], quartic_terms="+".join(fn), quartic_args=",".join(list(set(quartic_args[idx]))))
                    quartic_func = self.modelBuilder.factory_(quartic_expr)
                    self.modelBuilder.out._import(quartic_func)
                    fit_terms.append(quartic_names[idx])

                # Compile the full function
                fit_function = "expr::{name}('{fit_terms}',{fit_args})".format(name=name,fit_terms="+".join(fit_terms),fit_args=",".join(fit_terms[1:]))
                if verbose: print('fit_function', str(fit_function))
                quadratic = self.modelBuilder.factory_(fit_function) #Create RooFormulaVar object (from string)

                # Export fit function
                self.modelBuilder.out._import(quadratic)

        print(colors.fg.lightblue + "... Done !" + colors.reset)


 #####   ####  #
 #    # #    # #
 #    # #    # #
 #####  #    # #
 #      #    # #
 #       ####  #

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""

        # user can call combine with `--setPhysicsModelParameterRanges` to set to sensible ranges
        for wc in self.wcs:
            self.modelBuilder.doVar('{0}[0, {1}, {2}]'.format(wc,self.wc_ranges[wc][0],self.wc_ranges[wc][1]))
            print(colors.fg.lightblue + ' WCs to fit: {0}[0, {1}, {2}]'.format(wc,self.wc_ranges[wc][0],self.wc_ranges[wc][1]))
        # print(colors.fg.lightblue + "WCs to fit for: "+colors.reset+",".join(self.wcs))
        self.modelBuilder.doSet('POI', ','.join(self.wcs))
        self.setup()


 #   # # ###### #      #####     ####   ####    ##   #      ######
  # #  # #      #      #    #   #      #    #  #  #  #      #
   #   # #####  #      #    #    ####  #      #    # #      #####
   #   # #      #      #    #        # #      ###### #      #
   #   # #      #      #    #   #    # #    # #    # #      #
   #   # ###### ###### #####     ####   ####  #    # ###### ######

    def getYieldScale(self, bin, process):
        "Return the name of a RooAbsReal to scale this yield by or the two special values 1 and 0 (don't scale, and set to zero)"

        if (process,bin) not in self.procbins:
            if verbose: print(colors.fg.orange + '* ({0},{1}) not in self.procbins ! => Do not scale'.format(process, bin) + colors.reset)
            return 1
        else:
            if verbose: print(colors.fg.orange + '* Scaling ({0},{1}) by r_{0}_{1}'.format(process, bin) + colors.reset)
            return 'r_{0}_{1}'.format(process,bin)

# //--------------------------------------------
# //--------------------------------------------

print(colors.fg.lightblue + 'Creating custom Physics model...' + colors.reset)

#-- Create PhysicsModel object & read user options
eftmodel = EFTModel()
SM_name = opts["SM_name"]
eftmodel.wcs = opts["wcs"]
eftmodel.wc_ranges = opts["wc_ranges"]
