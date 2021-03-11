# //--------------------------------------------
#FIXME -- remove mbjmax ? cf. different correl data/mc

#-- CARL tZq

#-- ctz
features_CARL_tZq_ctz = []
features_CARL_tZq_ctz.append("recoZ_Pt")
features_CARL_tZq_ctz.append("recoZ_Eta")
features_CARL_tZq_ctz.append("recoZ_dPhill")
features_CARL_tZq_ctz.append("dR_blW")
features_CARL_tZq_ctz.append("dR_tZ")
features_CARL_tZq_ctz.append("recoLepTop_Eta")
features_CARL_tZq_ctz.append("recoLepTop_Pt")
features_CARL_tZq_ctz.append("dEta_Zjprime")
features_CARL_tZq_ctz.append("dR_tClosestLep")
features_CARL_tZq_ctz.append("mTW")
features_CARL_tZq_ctz.append("jprime_Pt") #Added

#-- ctw
features_CARL_tZq_ctw = []
features_CARL_tZq_ctw.append("recoZ_Pt")
features_CARL_tZq_ctw.append("recoZ_Eta")
features_CARL_tZq_ctw.append("recoZ_dPhill")
features_CARL_tZq_ctw.append("dR_tClosestLep")
features_CARL_tZq_ctw.append("dR_tZ")
features_CARL_tZq_ctw.append("lAsymmetry")
features_CARL_tZq_ctw.append("cosThetaStarPolTop")
features_CARL_tZq_ctw.append("jprime_Pt")
features_CARL_tZq_ctw.append("recoLepTop_Eta")
features_CARL_tZq_ctw.append("dR_blW")

#features_CARL_tZq_ctw.append("dEta_bjprime") #7mar21
#features_CARL_tZq_ctw.append("jPrimeAbsEta") #7mar21
#features_CARL_tZq_ctw.append("mHT")
#features_CARL_tZq_ctw.append("TopZsystem_M")
#features_CARL_tZq_ctw.append("maxDiJet_Pt")
#features_CARL_tZq_ctw.append("Mass_3l")
#features_CARL_tZq_ctw.append("dEta_lWjprime")

#-- cpqm
features_CARL_tZq_cpqm = []
features_CARL_tZq_cpqm.append("recoZ_Pt")
features_CARL_tZq_cpqm.append("recoZ_Eta")
features_CARL_tZq_cpqm.append("dR_jprimeClosestLep")
features_CARL_tZq_cpqm.append("dR_tZ")
features_CARL_tZq_cpqm.append("recoLepTop_Pt")
features_CARL_tZq_cpqm.append("dEta_Zjprime")
features_CARL_tZq_cpqm.append("dR_bjprime")
features_CARL_tZq_cpqm.append("dR_lWjprime")

#-- cpq3
features_CARL_tZq_cpq3 = []
features_CARL_tZq_cpq3.append("recoZ_Pt")
features_CARL_tZq_cpq3.append("recoZ_Eta")
features_CARL_tZq_cpq3.append("recoZ_dPhill")
features_CARL_tZq_cpq3.append("cosThetaStarPolZ")
features_CARL_tZq_cpq3.append("dR_blW")
features_CARL_tZq_cpq3.append("dR_tClosestLep")
features_CARL_tZq_cpq3.append("recoLepTop_Eta")
features_CARL_tZq_cpq3.append("dR_tZ")
features_CARL_tZq_cpq3.append("recoLepTop_Pt")
features_CARL_tZq_cpq3.append("jprime_Pt")
#features_CARL_tZq_cpq3.append("Mass_3l")

#-- cpt
features_CARL_tZq_cpt = []
features_CARL_tZq_cpt.append("recoZ_Eta")
features_CARL_tZq_cpt.append("dR_tZ")
features_CARL_tZq_cpt.append("recoZ_dPhill")
features_CARL_tZq_cpt.append("recoZ_Pt")
features_CARL_tZq_cpt.append("dR_Zjprime")
features_CARL_tZq_cpt.append("TopZsystem_M")

#-- Multiple operators
features_CARL_tZq_5D = []
features_CARL_tZq_5D.append("recoZ_Pt")
features_CARL_tZq_5D.append("recoZ_Eta")
features_CARL_tZq_5D.append("recoZ_dPhill")
features_CARL_tZq_5D.append("dR_blW") #data/mc
features_CARL_tZq_5D.append("dR_tZ")
features_CARL_tZq_5D.append("recoLepTop_Pt")
features_CARL_tZq_5D.append("dEta_Zjprime")
features_CARL_tZq_5D.append("jPrimeAbsEta") #FIXME -- remove ?
#features_CARL_tZq_5D.append("dR_tClosestLep")
#features_CARL_tZq_5D.append("cosThetaStarPolZ") #diff behaviours for diff operators...?

#CHANGED 6Mar21
features_CARL_tZq_5D.append("dR_jprimeClosestLep")
features_CARL_tZq_5D.append("jprime_Pt")
features_CARL_tZq_5D.append("recoLepTop_Eta")
features_CARL_tZq_5D.append("lAsymmetry")
#features_CARL_tZq_5D.append("mHT") #correl Zpt
#features_CARL_tZq_5D.append("dEta_bjprime") #poor discr
#features_CARL_tZq_5D.append("Mass_3l") #~ data/mc
#features_CARL_tZq_5D.append("cosThetaStarPolTop") #poor discr
#features_CARL_tZq_5D.append("maxDiJet_Pt") #~data/mc
#features_CARL_tZq_5D.append("TopZsystem_M") #~data/mc



# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------



#-- CARL ttZ

#-- ctz
features_CARL_ttZ_ctz = []
features_CARL_ttZ_ctz.append("recoZ_Pt")
features_CARL_ttZ_ctz.append("recoZ_Eta")
features_CARL_ttZ_ctz.append("recoZ_dPhill")
features_CARL_ttZ_ctz.append("dR_tZ")
features_CARL_ttZ_ctz.append("mTW")
features_CARL_ttZ_ctz.append("cosThetaStarPolZ")
features_CARL_ttZ_ctz.append("recoLepTop_Eta")
features_CARL_ttZ_ctz.append("recoLepTop_Pt")
#features_CARL_ttZ_ctz.append("Mass_3l")
#features_CARL_ttZ_ctz.append("mHT")

#-- ctw
features_CARL_ttZ_ctw = []
features_CARL_ttZ_ctw.append("recoZ_Pt")
features_CARL_ttZ_ctw.append("recoZ_dPhill")
features_CARL_ttZ_ctw.append("recoLepTop_Eta")
features_CARL_ttZ_ctw.append("lAsymmetry")
features_CARL_ttZ_ctw.append("Mass_3l")
features_CARL_ttZ_ctw.append("recoLepTop_Pt") #~ data/mc

#-- Could add these, but be very careful about data/MC agreement in Run2 template...
#features_CARL_ttZ_ctw.append("mTW") #~ data/mc SRttZ (1 bin)
#features_CARL_ttZ_ctw.append("maxDiJet_Pt") #~ data/mc #correl mHT
#features_CARL_ttZ_ctw.append("dR_tClosestLep") #~ data/mc
#features_CARL_ttZ_ctw.append("mHT") #~ data/mc

#-- cpqm
features_CARL_ttZ_cpqm = []
features_CARL_ttZ_cpqm.append("recoZ_Pt")
features_CARL_ttZ_cpqm.append("recoZ_Eta")
features_CARL_ttZ_cpqm.append("recoZ_dPhill")
features_CARL_ttZ_cpqm.append("dR_jprimeClosestLep")
features_CARL_ttZ_cpqm.append("dR_lWjprime")
features_CARL_ttZ_cpqm.append("recoLepTop_Pt")
features_CARL_ttZ_cpqm.append("dR_ZlW")
features_CARL_ttZ_cpqm.append("cosThetaStarPolZ")

#-- cpq3
features_CARL_ttZ_cpq3 = []
features_CARL_ttZ_cpq3.append("recoZ_Pt")
features_CARL_ttZ_cpq3.append("recoZ_Eta")
features_CARL_ttZ_cpq3.append("recoZ_dPhill")
features_CARL_ttZ_cpq3.append("mHT")
features_CARL_ttZ_cpq3.append("Mass_3l")

#-- cpt
features_CARL_ttZ_cpt = []
features_CARL_ttZ_cpt.append("recoZ_Pt")
features_CARL_ttZ_cpt.append("dEta_tjprime")
features_CARL_ttZ_cpt.append("dR_bjprime")
features_CARL_ttZ_cpt.append("dR_jprimeClosestLep")
features_CARL_ttZ_cpt.append("TopZsystem_M")
features_CARL_ttZ_cpt.append("dR_Zjprime")
features_CARL_ttZ_cpt.append("recoZ_dPhill")

#-- Multiple operators
features_CARL_ttZ_5D = []
features_CARL_ttZ_5D.append("recoZ_Pt")
features_CARL_ttZ_5D.append("recoZ_Eta")
features_CARL_ttZ_5D.append("recoZ_dPhill")
features_CARL_ttZ_5D.append("dR_tZ")
features_CARL_ttZ_5D.append("mTW")
features_CARL_ttZ_5D.append("recoLepTop_Eta")
features_CARL_ttZ_5D.append("maxDiJet_Pt")
#Added 6Mar21
features_CARL_ttZ_5D.append("recoLepTop_Pt")
features_CARL_ttZ_5D.append("Mass_3l")
features_CARL_ttZ_5D.append("cosThetaStarPolZ")


# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------


# NN-SM: tZq vs ttZ vs Others
#TMP check: '#OK #NO!' <-> good agreement between priv/central for tZq, but not for ttZ
features_SM = []
features_SM.append("Mass_3l") #OK #OK
features_SM.append("lAsymmetry") #OK #OK
features_SM.append("maxDeepJet") #OK #OK
features_SM.append("nbjets") #OK #OK
features_SM.append("recoZ_Eta") #~ #OK
features_SM.append("recoZ_dPhill") #~ #OK
features_SM.append("dR_lWjprime") #~ #OK
features_SM.append("mbjMax") #~ #OK
features_SM.append("jprime_Pt") #OK #OK
features_SM.append("metEt") #OK #OK

#-- removed to improve data/mc in SRtZq
#features_SM.append("recoZ_Pt") #OK #OK #Imperfect data/MC, correlated with recoZ_dPhill
#features_SM.append("maxDiJet_M") #? imperfect data/mc
#features_SM.append("mTW") #OK #OK #imperfect data/mc
#features_SM.append("jPrimeAbsEta") #OK #OK #imperfect data/mc

#-- removing improves modeling ?
#features_SM.append("maxEtaJet") #Poor data/MC
# features_SM.append("dR_blW") #NO! #OK

#-- Testing
# features_SM.append("maxDelRbL")
# features_SM.append("maxDiJet_Pt")

# //--------------------------------------------




# //--------------------------------------------
# //--------------------------------------------
# //--------------------------------------------
#-- Full list of available high-level features

'''
list_features.append("recoZ_Pt")
list_features.append("recoZ_Eta")
list_features.append("mHT")

list_features.append("mTW")
list_features.append("Mass_3l")
list_features.append("recoZ_dPhill")

list_features.append("lAsymmetry")
list_features.append("jPrimeAbsEta")
list_features.append("maxEtaJet")
list_features.append("maxDeepJet")
list_features.append("njets")
list_features.append("nbjets")

list_features.append("recoLepTop_Pt")
list_features.append("recoLepTop_Eta")
list_features.append("TopZsystem_M")
list_features.append("recoLepTopLep_Pt")
list_features.append("mbjMax")
list_features.append("maxDiJet_Pt")
list_features.append("maxDelRbL")
list_features.append("minDelRbL")
list_features.append("dR_ZlW")
list_features.append("dR_blW")
list_features.append("dR_tClosestJet")
list_features.append("dR_bW")
list_features.append("dEta_jprimeClosestLep")

list_features.append("cosThetaStarPolTop")
list_features.append("cosThetaStarPolZ")
list_features.append("dR_tjprime")
list_features.append("dEta_bjprime")
list_features.append("dR_bjprime")
list_features.append("dR_lWjprime")
list_features.append("dR_Zjprime")
list_features.append("dEta_lWjprime")

list_features.append("dR_tClosestLep")
list_features.append("recoLepTopLep_Eta")
list_features.append("maxDiJet_dPhi")
list_features.append("dR_jprimeClosestLep")
list_features.append("jprime_Pt")
list_features.append("DeepJet_2nd")
list_features.append("maxDeepCSV")
list_features.append("deepCSV_2nd")
list_features.append("recoZ_Phi")
list_features.append("recoZ_M")
list_features.append("maxDiJet_M")
list_features.append("recoLepTop_M")
list_features.append("TopZsystem_Pt")
list_features.append("dR_tZ")
list_features.append("dEta_Zjprime")
list_features.append("dEta_tjprime")
list_features.append("maxDiJet_dEta")
list_features.append("maxDiJet_dR")
list_features.append("recoZLepMinus_Pt")
list_features.append("recoZLepMinus_Eta")
list_features.append("recoZLepMinus_Phi")
list_features.append("recoZLepPlus_Pt")
list_features.append("recoZLepPlus_Eta")
list_features.append("recoZLepPlus_Phi")
list_features.append("recoLepTopLep_Phi")
list_features.append("jprime_Eta")
list_features.append("jprime_Phi")
list_features.append("TopZsystem_Eta")
list_features.append("TopZsystem_Phi")
list_features.append("recoLepTopB_Pt")
list_features.append("recoLepTopB_Eta")
list_features.append("recoLepTopB_Phi")
list_features.append("recoLepTop_Phi")
'''
