# //--------------------------------------------
# SCANS 1D OBS

#python EFTPlotter.py -P ctz -n ctzObs -m 1D
#python EFTPlotter.py -P ctz -n ctzObs -m 1D --paper

#python EFTPlotter.py -P ctw -n ctwObs -m 1D
#python EFTPlotter.py -P ctw -n ctwObs -m 1D --paper

#python EFTPlotter.py -P cpq3 -n cpq3Obs -m 1D
#python EFTPlotter.py -P cpq3 -n cpq3Obs -m 1D --paper

#python EFTPlotter.py -P cpqm -n cpqmObs -m 1D
#python EFTPlotter.py -P cpqm -n cpqmObs -m 1D --paper

#python EFTPlotter.py -P cpt -n cptObs -m 1D
#python EFTPlotter.py -P cpt -n cptObs -m 1D --paper

# //--------------------------------------------
# SCANS 1D EXP

#python EFTPlotter.py -P ctz -n ctzExp -m 1D
#python EFTPlotter.py -P ctz -n ctzExp -m 1D --paper

#python EFTPlotter.py -P ctw -n ctwExp -m 1D
#python EFTPlotter.py -P ctw -n ctwExp -m 1D --paper

#python EFTPlotter.py -P cpq3 -n cpq3Exp -m 1D
#python EFTPlotter.py -P cpq3 -n cpq3Exp -m 1D --paper

#python EFTPlotter.py -P cpqm -n cpqmExp -m 1D
#python EFTPlotter.py -P cpqm -n cpqmExp -m 1D --paper

#python EFTPlotter.py -P cpt -n cptExp -m 1D
#python EFTPlotter.py -P cpt -n cptExp -m 1D --paper

# //--------------------------------------------
# SCANS 2D OBS

python EFTPlotter.py -P ctz ctw -n ctzctwObs -f scan_ctz_ctw.root
python EFTPlotter.py -P ctz ctw -n ctzctwObs -f scan_ctz_ctw.root --paper

python EFTPlotter.py -P cpqm cpt -n cpqmcptObs -f scan_cpqm_cpt.root
python EFTPlotter.py -P cpqm cpt -n cpqmcptObs -m 1D -f scan_cpqm_cpt.root --paper

# //--------------------------------------------
