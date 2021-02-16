opts = {

#=== GENERAL ===#
'SM_name': 'SM', #SM point naming convention
'verbosity': 0, #Verbosity level

#-- Different lists / associations of WCs
"wcs": ['ctz', 'ctw', 'cpq3', 'cpqm', 'cpt'], #Full list of EFT operators #by default, include all operators in the physics model
# "wcs": ['ctz', 'ctw'],
#"wcs": ['ctz'],
"wc": 'ctz', #If want to float a single operator
"wc_ranges": { #Ranges for fits
                'ctz':  [-4,4],
                'ctw':  [-4,4],
                'cpq3': [-8,8],
                'cpqm': [-25,25], #Must increase range for 5D fits
                'cpt':  [-25,25], #Must increase range for 5D fits
                # 'cpqm': [-10,15],
                # 'cpt':  [-20,15],
              },
"wc_ranges_scan1D": { #Ranges for 1D grid scans
                'ctz':  [-1.5,1.5],
                'ctw':  [-1,1],
                'cpq3': [-2,2],
                'cpqm': [-4,4],
                'cpt':  [-7,7],
              },
"wc_ranges_scan2D": { #Ranges for 2D grid scans
                'ctz':  [-2,2],
                'ctw':  [-2,2],
                'cpq3': [-2,2],
                'cpqm': [-20,20],
                'cpt':  [-25,20],
              },
#"wcs_tracked": ['ctz', 'ctw', 'cpq3', 'cpqm', 'cpt'], #WCs to track *when not considered as POIs* (all by default)
# "wcs_tracked": ['ctz','ctw'],
"wcs_tracked": [],
"scan_wcs": ['ctz','ctw'], #Default pair of wcs for 2D scans
"wcs_pairs": ['ctz','ctw'], #Default pair of operators for 2D plots

#-- List of SM processes
"processes": ['tZq'],
#"processes": ['PrivMC_tWZ'],
"processes": ['PrivMC_tZq','PrivMC_ttZ','PrivMC_tWZ'], #Names of signal processes
#"processes": ['tZq','ttZ','tWZ'], #SM signal processes

#-- Names of SM signal strengths for processes of interest
#"SM_mus": ['r_tzq'], #Names of SM signal strengths for processes of interest
"SM_mus": ['r_tzq','r_ttz','r_twz'], #Names of SM signal strengths for processes of interest

#-- Name of SM signal strengths for single process
"SM_mu": 'r_tzq',

#-- Ranges of SM signal strengths
"SMmu_ranges": { #Ranges for scans, plots, etc.
                'r_tzq':  [0,3],
                'r_ttz':  [0,3],
                'r_twz':  [-5,15],
              },
}
