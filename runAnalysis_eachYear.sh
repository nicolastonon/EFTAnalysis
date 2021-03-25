make

#---------------------------------------

#-- All years for region specified in main ('signal' by default)

# ./analysis_main.exe 2016
# ./analysis_main.exe 2017
# ./analysis_main.exe 2018
# ./analysis_main.exe Run2

#---------------------------------------

#-- All regions

#./analysis_main.exe 2016 ttZ4l
#./analysis_main.exe 2016 xg
#./analysis_main.exe 2016 zz
#./analysis_main.exe 2016 tt
#./analysis_main.exe 2016 wz
#./analysis_main.exe 2016 dy
#./analysis_main.exe 2017 ttZ4l
#./analysis_main.exe 2017 xg
#./analysis_main.exe 2017 zz
#./analysis_main.exe 2017 tt
#./analysis_main.exe 2017 wz
#./analysis_main.exe 2017 dy
#./analysis_main.exe 2018 ttZ4l
#./analysis_main.exe 2018 xg
#./analysis_main.exe 2018 zz
#./analysis_main.exe 2018 tt
#./analysis_main.exe 2018 wz
#./analysis_main.exe 2018 dy

# ./analysis_main.exe Run2 ttZ4l
# ./analysis_main.exe Run2 xg
# ./analysis_main.exe Run2 zz
# ./analysis_main.exe Run2 tt
# ./analysis_main.exe Run2 wz
# ./analysis_main.exe Run2 dy

#---------------------------------------
#-- Make yield tables

# ./Yield_Table.exe 2016 signal
# ./Yield_Table.exe 2016 ttZ4l
# ./Yield_Table.exe 2016 xg
# ./Yield_Table.exe 2016 zz
# ./Yield_Table.exe 2016 tt
# ./Yield_Table.exe 2016 wz
# ./Yield_Table.exe 2016 dy
# ./Yield_Table.exe 2017 signal
# ./Yield_Table.exe 2017 ttZ4l
# ./Yield_Table.exe 2017 xg
# ./Yield_Table.exe 2017 zz
# ./Yield_Table.exe 2017 tt
# ./Yield_Table.exe 2017 wz
# ./Yield_Table.exe 2017 dy
# ./Yield_Table.exe 2018 signal
# ./Yield_Table.exe 2018 ttZ4l
# ./Yield_Table.exe 2018 xg
# ./Yield_Table.exe 2018 zz
# ./Yield_Table.exe 2018 tt
# ./Yield_Table.exe 2018 wz
# ./Yield_Table.exe 2018 dy

# ./Yield_Table.exe Run2 signal
# ./Yield_Table.exe Run2 ttZ4l
# ./Yield_Table.exe Run2 xg
# ./Yield_Table.exe Run2 zz
# ./Yield_Table.exe Run2 tt
# ./Yield_Table.exe Run2 wz
# ./Yield_Table.exe Run2 dy

#---------------------------------------
#-- All templates

# ./analysis_main.exe NN_SM
# ./analysis_main.exe NN_ctz
# ./analysis_main.exe NN_ctw
# ./analysis_main.exe NN_cpq3
# ./analysis_main.exe NN_5D
# ./analysis_main.exe Zpt

#-- Submit all templates

# ./submit analyze --template NN_SM
# ./submit analyze --template NN_ctz
# ./submit analyze --template NN_ctw
# ./submit analyze --template NN_cpq3
# ./submit analyze --template NN_5D
# ./submit analyze --template Zpt

#---------------------------------------
