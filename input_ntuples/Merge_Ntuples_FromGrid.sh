#Example MC filename : "Analyzer3l-V05_X-Summer16-TTZToLLNuNu_M_10-job331-d20200116-t161919.root"
#Example DATA filename : "Analyzer3l-V05_X-Run2018A-DoubleMuon-job13-d20200116-t154927.root"

#Choose either to process Data, MC, or both
processData=true
processMC=true
force="" #"-f" <-> overwrite existing target files; "" <-> do not
tolerance="" #'If the option -k is used, hadd will not exit on corrupt or non-existant input files but skip the offending files instead'; "" <-> do not; Warning: you may want missing files to trigger a failure !

dataStr="SinglePhoton|EGamma|SingleElectron|SingleMuon|DoubleEG|DoubleMuon|MuonEG"
rootSuffix=".root"

#This command : lists all jobs outputs -> keeps only rootfiles -> looks at the 4th keyword (samplename) -> removes any duplicate in list -> check if it is a data/MC file
liDATA=$(ls | egrep -e "$rootSuffix" | cut -d"-" -f4 | awk '!a[$0]++' | egrep -e "$dataStr")
liMC=$(ls | egrep -e "$rootSuffix" | cut -d"-" -f4 | awk '!a[$0]++' | egrep -v "$dataStr")
# echo "$liDATA"
# echo "$liMC"

outDir=$PWD/merged_ntuples
echo "mkdir $outDir"
mkdir $outDir

# --- MC ---
if [ "$processMC" = true ]; then

    # rm /tmp/tempMC.txt
    for yearname in "Summer16" "Fall17" "Autumn18"
    do
        echo $yearname

        if [ $yearname == "Summer16" ]; then
            dirname="2016"
        elif [ $yearname == "Fall17" ]; then
            dirname="2017"
        elif [ $yearname == "Autumn18" ]; then
            dirname="2018"
        else
            echo "ERROR ! Wrong yearname value :" $yearname
        fi
        echo "mkdir $outDir/$dirname"
        mkdir $outDir/$dirname

        for line in $liMC
        do
          echo $line

          echo "hadd $force $outDir/$dirname/merged_$line.root *$yearname-$line-*" > /tmp/ScriptMerging_$yearname_$line.sh #Merge all rootfiles matching exact 'year' and 'samplename' patterns
          # cat /tmp/tempMC.txt >> /tmp/tmp_ScriptMerging_$yearname_$line.sh #Concatenate rootfiles
          # tr '\n' ' ' < /tmp/tmp_ScriptMerging_$yearname_$line.sh > /tmp/ScriptMerging_$yearname_$line.sh #Replace newlines with spaces, copy to new script
          # rm /tmp/tmp_ScriptMerging_$yearname_$line.sh #Remove tmp script
          chmod 755 /tmp/ScriptMerging_$yearname_$line.sh #Make script executable
          /tmp/ScriptMerging_$yearname_$line.sh #Run merging script
          rm /tmp/ScriptMerging_$yearname_$line.sh #Remove script
          # rm /tmp/tempMC.txt #Remove list
        done
    done
fi

# --- DATA ---
if [ "$processData" = true ]; then

    # rm /tmp/tempDATA.txt
    for yearname in "Run2016*" "Run2017*" "Run2018*"
    do
        echo $yearname

        if [ $yearname == "Run2016*" ]; then
            dirname="2016"
        elif [ $yearname == "Run2017*" ]; then
            dirname="2017"
        elif [ $yearname == "Run2018*" ]; then
            dirname="2018"
        else
            echo "ERROR ! Wrong yearname value :" $yearname
        fi
        echo "Directory name :" $dirname
	echo "mkdir $outDir/$dirname"
        mkdir $outDir/$dirname

        for line in $liDATA
        do
          echo $line

          echo "hadd $force $outDir/$dirname/merged_$line.root *$yearname-$line-*" > /tmp/ScriptMerging_$yearname_$line.sh #Merge all rootfiles matching exact 'yearname' and 'samplename' patterns
          chmod 755 /tmp/ScriptMerging_$yearname_$line.sh #Make script executable
          /tmp/ScriptMerging_$yearname_$line.sh #Run merging script
          rm /tmp/ScriptMerging_$yearname_$line.sh #Remove script
          # rm /tmp/tempDATA.txt #Remove list
        done
    done
fi

# //--------------------------------------------

#--- FURTHER HADDS & RENAME NTUPLES ---
echo ''
echo ''
echo ''
echo '=== FURTHER MERGING AND RENAMING ==='

# --- MC ---
if [ "$processMC" = true ]; then

    for yearname in "Summer16" "Fall17" "Autumn18"
    do
        echo $yearname

        if [ $yearname == "Summer16" ]; then
            dirname="2016"
        elif [ $yearname == "Fall17" ]; then
            dirname="2017"
        elif [ $yearname == "Autumn18" ]; then
            dirname="2018"
        else
            echo "ERROR ! Wrong yearname value :" $yearname
        fi
        echo "Directory name :" $dirname

# //--------------------------------------------
#Merge by subgroups (not needed anymore?)

        #-- DY
        # hadd $force $outDir/$dirname/DY.root $outDir/$dirname/merged_DY*.root
        # rm $outDir/$dirname/merged_DY*.root

        #-- ST
        # hadd $force $outDir/$dirname/ST.root $outDir/$dirname/merged_ST_tW_antitop*.root $outDir/$dirname/merged_ST_tW_top*.root
        # rm $outDir/$dirname/merged_ST_tW_antitop*.root $outDir/$dirname/merged_ST_tW_top*.root

        #-- TTZ -- don't merge anymore (1to10 has no PDfs, etc.)
        #hadd $force $outDir/$dirname/ttZ.root $outDir/$dirname/merged_TTZToLLNuNu_M_10*.root $outDir/$dirname/merged_TTZToLL_M_1to10*.root
        #rm $outDir/$dirname/merged_TTZToLLNuNu_M_10*.root $outDir/$dirname/merged_TTZToLL_M_1to10*.root

        #-- ZZTo4l (all prod modes) -- now merge directly as VV(V), see below
        # hadd $force $outDir/$dirname/ZZ4l.root $outDir/$dirname/merged_ZZTo4L*.root $outDir/$dirname/merged_GluGluToContinToZZTo*.root $outDir/$dirname/merged_GluGluHToZZTo4L*.root $outDir/$dirname/merged_VBF_HToZZTo4L*.root $outDir/$dirname/merged_VHToNonbb*.root
        # rm $outDir/$dirname/merged_ZZTo4L*.root $outDir/$dirname/merged_GluGluToContinToZZTo*.root $outDir/$dirname/merged_GluGluHToZZTo4L*.root $outDir/$dirname/merged_VBF_HToZZTo4L*.root $outDir/$dirname/merged_VHToNonbb*.root

        #-- ZGToLLG_01J (sums ZGToLLG_01J+ZGToLLG_01J_LoosePtlPtg if available, else only 1) -- now merge directly as XG, see below
        # hadd $force $outDir/$dirname/ZGToLLG_01J.root $outDir/$dirname/merged_ZGToLLG_01J*.root
        # rm $outDir/$dirname/merged_ZGToLLG_01J*.root

# //--------------------------------------------
#- Merge directly ntuple groups

        #-- t(t)X
    	hadd $force $outDir/$dirname/tX.root \
        $outDir/$dirname/merged_TTZToLL_M_1to10*.root \
        $outDir/$dirname/merged_THQ*.root \
        $outDir/$dirname/merged_THW*.root \
        $outDir/$dirname/merged_ttHToNonbb_M125*.root \
        $outDir/$dirname/merged_TTWJetsToLNu*.root \
        $outDir/$dirname/merged_TTZZ*.root \
        $outDir/$dirname/merged_TTWW*.root \
        $outDir/$dirname/merged_TTWZ*.root \
        $outDir/$dirname/merged_TTZH*.root \
        $outDir/$dirname/merged_TTWH*.root \
        $outDir/$dirname/merged_TTTT*.root \
        $outDir/$dirname/merged_TTHH*.root

        rm $outDir/$dirname/merged_TTZToLL_M_1to10*.root \
        $outDir/$dirname/merged_THQ*.root \
        $outDir/$dirname/merged_THW*.root \
        $outDir/$dirname/merged_ttHToNonbb_M125*.root \
        $outDir/$dirname/merged_TTWJetsToLNu*.root \
        $outDir/$dirname/merged_TTZZ*.root \
        $outDir/$dirname/merged_TTWW*.root \
        $outDir/$dirname/merged_TTWZ*.root \
        $outDir/$dirname/merged_TTZH*.root \
        $outDir/$dirname/merged_TTWH*.root \
        $outDir/$dirname/merged_TTTT*.root \
        $outDir/$dirname/merged_TTHH*.root

        #-- VV(V) #Very slow #Command line: [hadd VVV_18.root Analyzer3l-V12-*18-*ZZT* Analyzer3l-V12-*18-*VHT* Analyzer3l-V12-*18-*ZZZ* Analyzer3l-V12-*18-*WZZ* Analyzer3l-V12-*18-*WWW* Analyzer3l-V12-*18-*WZZ*]
        hadd $force $outDir/$dirname/VVV.root \
        $outDir/$dirname/merged_ZZTo4L*.root \
        $outDir/$dirname/merged_GluGluToContinToZZTo*.root \
        $outDir/$dirname/merged_GluGluHToZZTo4L*.root \
        $outDir/$dirname/merged_VBF_HToZZTo4L*.root \
        $outDir/$dirname/merged_VHToNonbb*.root \
        $outDir/$dirname/merged_ZZZ*.root \
        $outDir/$dirname/merged_WZZ*.root \
        $outDir/$dirname/merged_WWW*.root \
        $outDir/$dirname/merged_WWZ*.root

        rm $outDir/$dirname/merged_ZZTo4L*.root \
        $outDir/$dirname/merged_GluGluToContinToZZTo*.root \
        $outDir/$dirname/merged_GluGluHToZZTo4L*.root \
        $outDir/$dirname/merged_VBF_HToZZTo4L*.root \
        $outDir/$dirname/merged_VHToNonbb*.root \
        $outDir/$dirname/merged_ZZZ*.root \
        $outDir/$dirname/merged_WZZ*.root \
        $outDir/$dirname/merged_WWW*.root \
        $outDir/$dirname/merged_WWZ*.root

        #-- X+G #NB: will sum ZGToLLG_01J+ZGToLLG_01J_LoosePtlPtg+ZGToLLG_01J_lowMLL_Fall17 if available... if it is the case, make sure to combine samples properly
        hadd $force $outDir/$dirname/XG.root \
        $outDir/$dirname/merged_ZGToLLG_01J*.root \
        $outDir/$dirname/merged_TTGamma_Dilept*.root

        rm $outDir/$dirname/merged_ZGToLLG_01J*.root \
        $outDir/$dirname/merged_TTGamma_Dilept*.root

        #-- TTbar
        #hadd $force $outDir/$dirname/TTbar.root \
        #$outDir/$dirname/merged_TTTo2L2Nu*.root \
        #$outDir/$dirname/merged_TTToSemiLeptonic*.root

        #rm $outDir/$dirname/merged_TTTo2L2Nu*.root \
        #$outDir/$dirname/merged_TTToSemiLeptonic*.root

        #-- DY
        #hadd $force $outDir/$dirname/DY.root $outDir/$dirname/merged_DY*.root

        #rm $outDir/$dirname/merged_DY*.root

#//--------------------------------------------
#-- Merge files for single samples

        #-- t(t)X
    	# mv $outDir/$dirname/merged_TTZToLL_M_1to10*.root $outDir/$dirname/ttZ_M1to10.root
        # mv $outDir/$dirname/merged_THQ*.root $outDir/$dirname/tHq.root
        # mv $outDir/$dirname/merged_THW*.root $outDir/$dirname/tHW.root
        # mv $outDir/$dirname/merged_ttHToNonbb_M125*.root $outDir/$dirname/ttH.root
        # mv $outDir/$dirname/merged_TTWJetsToLNu*.root $outDir/$dirname/ttW.root
        # mv $outDir/$dirname/merged_TTZZ*.root $outDir/$dirname/ttZZ.root
        # mv $outDir/$dirname/merged_TTWW*.root $outDir/$dirname/ttWW.root
        # mv $outDir/$dirname/merged_TTWZ*.root $outDir/$dirname/ttWZ.root
        # mv $outDir/$dirname/merged_TTZH*.root $outDir/$dirname/ttZH.root
        # mv $outDir/$dirname/merged_TTWH*.root $outDir/$dirname/ttWH.root
        # mv $outDir/$dirname/merged_TTTT*.root $outDir/$dirname/tttt.root
        # mv $outDir/$dirname/merged_TTHH*.root $outDir/$dirname/ttHH.root

        #-- VV(V)
        # mv $outDir/$dirname/merged_WWW*.root $outDir/$dirname/WWW.root
        # mv $outDir/$dirname/merged_WWZ*.root $outDir/$dirname/WWZ.root
        # mv $outDir/$dirname/merged_WZZ*.root $outDir/$dirname/WZZ.root
        # mv $outDir/$dirname/merged_ZZZ*.root $outDir/$dirname/ZZZ.root

        #-- X+G
        # mv $outDir/$dirname/merged_TTGamma_Dilept*.root $outDir/$dirname/TTGamma_Dilep.root

        #-- TTbar
        # mv $outDir/$dirname/merged_TTTo2L2Nu*.root $outDir/$dirname/TTbar_DiLep.root
        # mv $outDir/$dirname/merged_TTToSemiLeptonic*.root $outDir/$dirname/TTbar_SemiLep.root

        #-- WZ
        mv $outDir/$dirname/merged_WZTo3LNu*.root $outDir/$dirname/WZ.root

        #-- Signals
        mv $outDir/$dirname/merged_TTZToLLNuNu_M_10*.root $outDir/$dirname/ttZ.root
        mv $outDir/$dirname/merged_tZq_ll*.root $outDir/$dirname/tZq.root
        mv $outDir/$dirname/merged_ST_tWll*.root $outDir/$dirname/tWZ.root
        mv $outDir/$dirname/merged_PrivMC_ttll*.root $outDir/$dirname/PrivMC_ttZ.root
    	mv $outDir/$dirname/merged_PrivMC_tllq*.root $outDir/$dirname/PrivMC_tZq.root
    	mv $outDir/$dirname/merged_PrivMC_twll*.root $outDir/$dirname/PrivMC_tWZ.root

        #-- OBSOLETE
        #mv $outDir/$dirname/merged_TGJets*.root $outDir/$dirname/tGJets.root
        #mv $outDir/$dirname/merged_ZZTo2L2Nu*.root $outDir/$dirname/ZZTo2L2Nu.root
        #mv $outDir/$dirname/merged_ZGTo2LG*.root $outDir/$dirname/ZG2l2g.root
        #mv $outDir/$dirname/merged_ZZTo2L2Q*.root $outDir/$dirname/ZZ2l2q.root
        #mv $outDir/$dirname/merged_WWTo2L2Nu*.root $outDir/$dirname/WWTo2L2Nu.root
        #mv $outDir/$dirname/merged_WZZTo2L2Nu*.root $outDir/$dirname/ZZTo2L2Nu.root
        #mv $outDir/$dirname/merged_WGToLNuG*.root $outDir/$dirname/WGToLNuG.root
        #mv $outDir/$dirname/merged_WZTo2L2Q*.root $outDir/$dirname/WZ2l2q.root
        #mv $outDir/$dirname/merged_ZGToLLG_01J*.root $outDir/$dirname/ZGToLLG_01J.root
    	#mv $outDir/$dirname/merged_WJetsToLNu*.root $outDir/$dirname/WJetsToLNu.root
        #mv $outDir/$dirname/merged_ZZTo4L*.root $outDir/$dirname/ZZ4l.root

    done
fi
#//--------------------------------------------

# --- DATA ---
if [ "$processData" = true ]; then

    for yearname in "Run2016*" "Run2017*" "Run2018*"
    do
        echo $yearname

        if [ $yearname == "Run2016*" ]; then
            dirname="2016"
        elif [ $yearname == "Run2017*" ]; then
            dirname="2017"
        elif [ $yearname == "Run2018*" ]; then
            dirname="2018"
        else
            echo "ERROR ! Wrong yearname value :" $yearname
        fi
        echo "Directory name :" $dirname

        hadd $tolerance $force $outDir/$dirname/DATA.root $outDir/$dirname/*EGamma* $outDir/$dirname/*SinglePhoton* $outDir/$dirname/*SingleElectron* $outDir/$dirname/*SingleMuon* $outDir/$dirname/*MuonEG* $outDir/$dirname/*DoubleEG* $outDir/$dirname/*DoubleMuon* #If the option -k is used, hadd will not exit on corrupt or non-existant input files but skip the offending files instead
        rm $outDir/$dirname/*EGamma*.root
        rm $outDir/$dirname/*SinglePhoton*.root
        rm $outDir/$dirname/*SingleElectron*.root
        rm $outDir/$dirname/*SingleMuon*.root
        rm $outDir/$dirname/*DoubleEG*.root
        rm $outDir/$dirname/*DoubleMuon*.root
        rm $outDir/$dirname/*MuonEG*.root
    done
fi
