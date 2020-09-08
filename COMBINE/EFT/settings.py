opts = {

#=== GENERAL ===#
'SM_name': 'SM', #SM point naming convention
'verbosity': 0, #Verbosity level

# "wcs": ['ctz', 'ctw', 'cpq3', 'cpqm', 'cpt'], #Full list of EFT operators
# "wcs": ['ctz', 'ctw'],
"wcs": ['ctz'],
"wc": 'ctz', #If want to float a single operator
"wc_ranges": { #Ranges for scans, plots, etc.
                'ctz':  [-4,4],
                'ctw':  [-6,6],
                'cpq3': [-6,6],
                'cpqm': [-6,6],
                'cpt':  [-6,6],
              },
# "wcs_tracked": ['ctz','ctw'], #Save values of these profiled nuisances, if not considered as POIs #FIXME -- should instead define as POIs and use '--redefineSignalPOIs' (and float/freeze/...)?
"wcs_tracked": ['ctz'],
"scan_wcs": ['ctz','ctw'], # Default pair of wcs for 2D scans
"wcs_pairs": ['ctz','ctw'], #Pair of operators for 2D plots

"processes": ['PrivMC_tZq','PrivMC_ttZ'], #Names of signal processes
# "processes": ['PrivMC_tZq'],
#"SM_mus": ['r_tzq','r_ttz'], #Names of SM signal strengths for processes of interest
"SM_mus": ['r_tzq'],
"SM_mu": 'r_tzq', #If want to float a single process
# "processes": ['tzq'], #SM signal processes
# "SM_mus": ['r_tzq'], #Names of SM signal strengths for processes of interest
"SMmu_ranges": { #Ranges for scans, plots, etc.
                'r_tzq':  [0,5],
                'r_ttz':  [0,5],
              },
}
