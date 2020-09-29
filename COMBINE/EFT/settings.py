opts = {

#=== GENERAL ===#
'SM_name': 'SM', #SM point naming convention
'verbosity': 0, #Verbosity level

"wcs": ['ctz', 'ctw', 'cpq3', 'cpqm', 'cpt'], #Full list of EFT operators #by default, include all operators in the physics model
# "wcs": ['ctz', 'ctw'],
# "wcs": ['ctz'],
"wc": 'ctz', #If want to float a single operator
"wc_ranges": { #Ranges for scans, plots, etc.
                'ctz':  [-4,4],
                'ctw':  [-4,4],
                'cpq3': [-8,8],
                'cpqm': [-8,8],
                'cpt':  [-8,8],
              },
"wcs_tracked": ['ctz', 'ctw', 'cpq3', 'cpqm', 'cpt'], #WCs to track *when not considered as POIs* (all by default)
# "wcs_tracked": ['ctz','ctw'], 
# "wcs_tracked": ['ctz'],
"scan_wcs": ['ctz','ctw'], #Default pair of wcs for 2D scans
"wcs_pairs": ['ctz','ctw'], #Default pair of operators for 2D plots

"processes": ['PrivMC_tZq','PrivMC_ttZ'], #Names of signal processes
# "processes": ['PrivMC_tZq'],
# "processes": ['tzq'], #SM signal processes

"SM_mus": ['r_tzq'],
#"SM_mus": ['r_tzq','r_ttz'], #Names of SM signal strengths for processes of interest
# "SM_mus": ['r_tzq'], #Names of SM signal strengths for processes of interest
"SM_mu": 'r_tzq', #If want to float a single process

"SMmu_ranges": { #Ranges for scans, plots, etc.
                'r_tzq':  [0,3],
                'r_ttz':  [0,3],
              },
}
