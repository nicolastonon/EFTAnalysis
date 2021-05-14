#Usage : 'python Parser_Datacard_Template.py [channel] [variable] [root file] ...'

import fileinput
import os
import sys, getopt
import string
import re
import os

# //--------------------------------------------
total = len(sys.argv)
#cmdargs = str(sys.argv)

theVar = str(sys.argv[1])
channel = str(sys.argv[2])
year = str(sys.argv[3])
theFiletoRead = sys.argv[4]
systChoice = str(sys.argv[5])
statChoice = str(sys.argv[6])
datacard_dir = str(sys.argv[7])
theWorkspaceToRead = str(sys.argv[8])
fileToSearch = str(sys.argv[9])
nbins = str(sys.argv[10])

# region = str(sys.argv[8]) #'SR', 'ttZ' (CR), etc.

print('\n * Creating datacard for year : '+year+' / channel : '+channel+' / variable : '+theVar)

# fileToSearch = "Template_Datacard.txt" #TEMPLATE card to parse
file = open(fileToSearch).read()
# print('fileToSearch', fileToSearch)

use_rph = False
if "#RPH" in file: use_rph = True #Detect keyword at b.o.f.
# print('use_rph', use_rph)

#if any(x in theVar for x in ['_CR']): fileToSearch = "Template_Datacard_noSig.txt" #Special case: don't consider any signal in CRs (approximation)

# //--------------------------------------------
# if(channel!="" and channel!="all" and channel!="uuu" and channel!="uue" and channel!="eeu" and channel!="eee" and channel!="ee" and channel!="uu" and channel!="ue"):
#     print("wrong channel")
#     print("channel should be '', 'all', 'uuu', 'uue', 'eeu' or 'uuu' or 'ee' or 'uu' or 'ue'")
#     exit()

#-- Channels
if channel=="all": #Use the channel='all' keyword because parser needs to read some arg ! But don't want to appear in datacards. Also remove the "_" in front
    channel=""
else:
    channel = "_" + channel #Keep the "_" in front
    # varchan="varchan" #if no subcategorization, want to remove the "_" between 'var' and 'chan' !

#-- Shape systs
if systChoice=="withShape":
    shape = ""
elif systChoice == "noShape":
    shape = "#"
else:
    print("Wrong systChoice value ! should be 'withShape' or 'noShape'")
    exit()

#-- MC stat uncertainty
if statChoice=="withStat":
    stat = ""
elif statChoice == "noStat":
    stat = "#"
else:
    print("Wrong statChoice value ! should be 'withStat' or 'noStat'")
    exit()

#-- SR/CR differences (when using RooParametricHist)
SR = '' #Default: consider region as SR
CR = '#'
PrivMC_CR = '1' #Default: activate systs also for PrivMC EFT samples
rate_sig = '-1'
if 'SR' in theVar:
    SR = ''
    rate_sig = '1' #Applied only to EFT signals in SR (must multiply PDF_norm by 1)
elif 'CR' in theVar:
    CR = ''
    SR = '#'
    PrivMC_CR = '-' #Choose to ignore shape systs for PrivMC samples in CRs (since nominals are dummies anyway)

#-- Year-to-year correlations
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

#-- Hard-coded special cases: e.g. if a lN syst. is correlated between years with different values, use a marker replaced with year-specific values by parsing code
#Lumi1617 = "-"
Lumi1718 = "-"
Lumi161718 = "-" #LumiXY
if year=="2016":
    #Lumi1617 = "1.008"
    Lumi161718 = "1.006"
elif year=="2017":
    #Lumi1617 = "1.006"
    Lumi1718 = "1.006"
    Lumi161718 = "1.009"
elif year=="2018":
    Lumi1718 = "1.002"
    Lumi161718 = "1.02"

#-- Region-specific flags
isSRtZq = '#'; isWZextrap = '#' #Disactive corresponding systematics by default (unless in these specific regions)
if 'SRtZq' in theVar: isSRtZq = ''
if 'SR' in theVar: isWZextrap = '' #Apply WZ extrapolation uncertainties in all SR, since they have b jet requirements not present in the WZ CR

theorygroup = ''
if 'CR' in theVar: theorygroup = '#' #'theory' systematics only apply to signals; but signals' nom/variation are ignored in CRs -> can't define theory nuisance group in CRs

#-- Obsolete, to remove
# isCRWZ = '#'; isCRZZ = '#'; isCRDY = '#'
# if 'CRWZ' in theVar: isCRWZ = ''
# elif 'CRZZ' in theVar: isCRZZ = ''
# elif 'CRDY' in theVar: isCRDY = ''

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

#-- REPLACE KEYWORDS #Order may be important (if subkeywords are included in larger keywords, etc.)
file = file.replace("[YEAR]", year)
file = file.replace("[2016]", is2016)
file = file.replace("[2017]", is2017)
file = file.replace("[2018]", is2018)
file = file.replace("[201617]", is201617)
file = file.replace("[201718]", is201718)
file = file.replace("[Lumi161718]", Lumi161718) #Was: LumiXY
file = file.replace("[Lumi1718]", Lumi1718)
file = file.replace("[SHAPE]", shape)
file = file.replace("[STAT]", stat)
file = file.replace("filetoread", theFiletoRead)
file = file.replace("workspacetoread", theWorkspaceToRead)
file = file.replace("[VAR]",theVar)
file = file.replace("_[CHAN]", channel)
file = file.replace("[SR]", SR)
file = file.replace("[CR]", CR)
file = file.replace("[RATE_SIG]", rate_sig)
file = file.replace("[PrivMC_CR]", PrivMC_CR)

#For region-specific uncertainties (e.g. CR extrapolation uncertainties)
file = file.replace("[SRtZq]", isSRtZq)
file = file.replace("[WZext]", isWZextrap)

file = file.replace("[THEORY]", theorygroup)

#Testing: make some nuisance unique in each bin
# if 'SRttZ' in theVar:
#     file = file.replace("PrivMC_ttZ_rate", "PrivMC_ttZ_rate_{}".format(theVar))

#Obsolete
# file = file.replace("[CRWZ]", isCRWZ)
# file = file.replace("[CRZZ]", isCRZZ)
# file = file.replace("[CRDY]", isCRDY)
#file = file.replace("[Lumi1617]", Lumi1617)

#--------------------------------------------

outname = datacard_dir+"/datacard_"+theVar;
if channel != "": outname+= channel;
if use_rph: outname+= '_rph';
outname+= ".txt"

print('==> Datacard ', outname, ' created...')

f = open(outname, 'w')
f.write(file)
f.close()


# //--------------------------------------------
#-- If using RPH/RDH (in SR regions), can't use autoMCstats
# --> Effect of MC stat error incorporated directly into RooParametricHists (for EFT signals)
# --> For RooDataHists (RDHs), must also create RDH corresponding to up/down MCstat variations in each bin for each chan/proc (treated as shape systematics)
# But the template datacard does not know about histo bins, etc. Below we duplicate/rename all the corresponding shape nuisance parameters for all histo bins, etc.

if use_rph==True and 'SR' in theVar: #Only need this in SRs, if using RPH/RDH objects

    #-- Read contents of the card we just produced
    f = open(outname, "r")
    contents = f.readlines()
    f.close()

    #-- Reopen card in write mode
    f = open(outname, "w")

    #-- Look for the 'MCstat' shape nuisances (attached to RDHs)
    idx_line_orig = 0 #Keep track of the number of the line we read in the original card
    counter_newlines = 0 #Keep track of the positions of the new lines (in the modified card)
    skip_newline = False
    tmp = {}
    for line in contents:

        if 'MCstat' in line:
            # print('idx_line_orig', idx_line_orig, 'line', line)

            idx_tmp = idx_line_orig+counter_newlines

            if 'countExp' in theVar: #Special case: single bin
                line_tmp = line #Duplicate template line
                line_tmp = line_tmp.replace("_[BIN]", "") #No bin number
                line_tmp = line_tmp.replace("#MCstat", "MCstat") #Original template line is commented out and will be ignored ; but modified line must be activated
                tmp[idx_tmp+1] = line_tmp #Will insert new line

            else: #Default: 1 MCstat nuisance per process --> Duplicate lines to have 1 nuisance per histo bin
                for ibin in range(1, int(nbins)+1): #Loop on histobins
                    line_tmp = line #Duplicate line
                    line_tmp = line_tmp.replace("[BIN]", "bin{}".format(ibin)) #Bin number info
                    line_tmp = line_tmp.replace("#MCstat", "MCstat") #Original template line is commented out and will be ignored ; but modified line must be activated
                    tmp[idx_tmp+ibin] = line_tmp #Will insert new line
                    # print('Inserted key', idx_tmp+ibin, ' / value', line_tmp, ' in dict...')
                    counter_newlines+= 1
        idx_line_orig+= 1

    #-- Insert all new lines in card
    for idx_line in tmp:
        contents.insert(idx_line, tmp[idx_line]) #Duplicate line, for each bin
        # print('... Inserted new line:', tmp[idx_line])

    f.writelines(contents) #Write modified card
    f.close()
