#Usage : 'python Parser_Datacard_Template.py [channel] [variable] [root file] ...'

import fileinput
import os
import sys, getopt
import string
import re
import os

# //--------------------------------------------
total = len(sys.argv)
cmdargs = str(sys.argv)

theVar = str(sys.argv[1])
channel = str(sys.argv[2])
year = str(sys.argv[3])
theFiletoRead= sys.argv[4]
systChoice = str(sys.argv[5])
statChoice = str(sys.argv[6])
datacard_dir = str(sys.argv[7])
# region = str(sys.argv[8]) #'SR', 'ttZ' (CR), etc.

print('\n * Creating datacard for year : '+year+' / channel : '+channel+' / variable : '+theVar)

fileToSearch = "Template_Datacard.txt" #TEMPLATE to parse

# //--------------------------------------------
# if(channel!="" and channel!="all" and channel!="uuu" and channel!="uue" and channel!="eeu" and channel!="eee" and channel!="ee" and channel!="uu" and channel!="ue"):
#     print("wrong channel")
#     print("channel should be '', 'all', 'uuu', 'uue', 'eeu' or 'uuu' or 'ee' or 'uu' or 'ue'")
#     exit()

if channel=="all": #Use the channel='all' keyword because parser needs to read some arg ! But don't want to appear in datacards. Also remove the "_" in front
    channel=""
else:
    channel = "_" + channel #Keep the "_" in front
    # varchan="varchan" #if no subcategorization, want to remove the "_" between 'var' and 'chan' !

#If don't want shape systematics, will comment them out
if systChoice=="withShape":
    shape = ""
elif systChoice == "noShape":
    shape = "#"
else:
    print("Wrong systChoice value ! should be 'withShape' or 'noShape'")
    exit()

#If don't want statistical uncertainties, will comment them out
if statChoice=="withStat":
    stat = ""
elif statChoice == "noStat":
    stat = "#"
else:
    print("Wrong statChoice value ! should be 'withStat' or 'noStat'")
    exit()

#Deal with all possible year correlations
is2016="#"
is2017="#"
is2018="#"
is201617="#"
is201718="#"
if year=="2016":
    is2016=""
    is201617=""
elif year=="2017":
    is2017=""
    is201617=""
    is201718=""
elif year=="2018":
    is2018=""
    is201718=""

#Hard-coded special cases: e.g. if a lN syst. is correlated between years with different values, use a marker replaced with year-specific values by parsing code
Lumi1617 = "-"
Lumi1718 = "-"
LumiXY = "-"
if year=="2016":
    Lumi1617 = "1.008"
    LumiXY = "1.009"
elif year=="2017":
    Lumi1617 = "1.006"
    Lumi1718 = "1.004"
    LumiXY = "1.008"
elif year=="2018":
    Lumi1718 = "1.003"
    LumiXY = "1.02"

#--------------------------------------------
# ele_sys = "" #Ele systematics only in ele channels

# if channel=="uuu" or "mmm" in channel:
#     ele_sys = "#"

#--------------------------------------------
#Can add a rateParam line to control normalization of processes from datacard
ratePar = "#" #empty to activate, or '#' to disactivate
sigPar = "tZq" #e.g. "tZq"
rateVal = "1" #or '2' to double the rate of the process ? (verify)


#--------------------------------------------
s = open(fileToSearch).read()

#-- REPLACE KEYWORDS
s = s.replace("[YEAR]", year)
s = s.replace("[2016]", is2016)
s = s.replace("[2017]", is2017)
s = s.replace("[2018]", is2018)
s = s.replace("[201617]", is201617)
s = s.replace("[201718]", is201718)
s = s.replace("[Lumi1617]", Lumi1617)
s = s.replace("[Lumi1718]", Lumi1718)
s = s.replace("[LumiXY]", LumiXY)

s = s.replace("[SHAPE]", shape)
s = s.replace("[STAT]", stat)
s = s.replace("filetoread", theFiletoRead)
s = s.replace("[VAR]",theVar)
s = s.replace("_[CHAN]", channel)

# s = s.replace("[ele_sys]",ele_sys)
# s = s.replace("sigPar", sigPar)
# s = s.replace("[ratePar]", ratePar)
# s = s.replace("rateVal", rateVal)




# //--------------------------------------------
#Replace some predefined markers with relevant values

#-- QFlip markers
# if nLep == "3l" or (nLep == "2l" and (channel == "uu" or "mm" in channel) ):
#     s = re.sub(r'\[qflip\].*?\[qflip\]', r'', s) #Erase stuff signaled by markers => Remove QFlip
    # print(s)
# s = s.replace("[qflip]", "") #Remove the remaining markers

#--------------------------------------------
print('==> Datacard created...')

outname = datacard_dir+"/datacard_"+theVar;
if channel != "":
    outname=outname+channel;
outname=outname+".txt"

f = open(outname, 'w')
f.write(s)
f.close()
