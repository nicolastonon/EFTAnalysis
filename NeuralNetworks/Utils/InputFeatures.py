# //--------------------------------------------

#-- CARL_singlePoint tZq

#-- ctz
features_CARL_singlePoint_tZq_ctz = []
features_CARL_singlePoint_tZq_ctz.append("recoZ_Pt")
features_CARL_singlePoint_tZq_ctz.append("recoZ_Eta")
features_CARL_singlePoint_tZq_ctz.append("recoZ_dPhill")
features_CARL_singlePoint_tZq_ctz.append("jPrimeAbsEta")
features_CARL_singlePoint_tZq_ctz.append("cosThetaStarPolZ")

features_CARL_singlePoint_tZq_ctz.append("dEta_Zjprime")
features_CARL_singlePoint_tZq_ctz.append("dR_blW")
features_CARL_singlePoint_tZq_ctz.append("dR_jprimeClosestLep")
features_CARL_singlePoint_tZq_ctz.append("dR_tClosestLep")
features_CARL_singlePoint_tZq_ctz.append("dR_tZ")
features_CARL_singlePoint_tZq_ctz.append("lAsymmetry")
features_CARL_singlePoint_tZq_ctz.append("recoLepTop_Eta")
features_CARL_singlePoint_tZq_ctz.append("recoLepTop_Pt")
features_CARL_singlePoint_tZq_ctz.append("maxDiJet_Pt")
features_CARL_singlePoint_tZq_ctz.append("mTW")

# features_CARL_singlePoint_tZq_ctz.append("mHT")
# features_CARL_singlePoint_tZq_ctz.append("Mass_3l")
# features_CARL_singlePoint_tZq_ctz.append("recoLepTop_Pt")
# features_CARL_singlePoint_tZq_ctz.append("dR_tClosestJet")
# features_CARL_singlePoint_tZq_ctz.append("recoLepTop_Eta")
# features_CARL_singlePoint_tZq_ctz.append("lAsymmetry")

#-- ctw
features_CARL_singlePoint_tZq_ctw = []
features_CARL_singlePoint_tZq_ctw.append("recoZ_Pt")
features_CARL_singlePoint_tZq_ctw.append("recoZ_Eta")
features_CARL_singlePoint_tZq_ctw.append("recoZ_dPhill")
features_CARL_singlePoint_tZq_ctw.append("jPrimeAbsEta")
features_CARL_singlePoint_tZq_ctw.append("cosThetaStarPolTop")
features_CARL_singlePoint_tZq_ctw.append("lAsymmetry")
#
features_CARL_singlePoint_tZq_ctw.append("dEta_bjprime")
features_CARL_singlePoint_tZq_ctw.append("dR_blW")
features_CARL_singlePoint_tZq_ctw.append("dR_tClosestLep")
features_CARL_singlePoint_tZq_ctw.append("dR_tZ")
features_CARL_singlePoint_tZq_ctw.append("recoLepTop_Eta")
features_CARL_singlePoint_tZq_ctw.append("Mass_3l")
features_CARL_singlePoint_tZq_ctw.append("maxDiJet_Pt")
features_CARL_singlePoint_tZq_ctw.append("jprime_Pt")

#-- cpqm
features_CARL_singlePoint_tZq_cpqm = []
features_CARL_singlePoint_tZq_cpqm.append("recoZ_Pt")
features_CARL_singlePoint_tZq_cpqm.append("recoZ_Eta")
features_CARL_singlePoint_tZq_cpqm.append("recoZ_dPhill")
features_CARL_singlePoint_tZq_cpqm.append("dEta_tjprime")
features_CARL_singlePoint_tZq_cpqm.append("dR_bjprime")

features_CARL_singlePoint_tZq_cpqm.append("dR_jprimeClosestLep")
features_CARL_singlePoint_tZq_cpqm.append("dR_lWjprime")
features_CARL_singlePoint_tZq_cpqm.append("dR_tZ")
features_CARL_singlePoint_tZq_cpqm.append("lAsymmetry")
features_CARL_singlePoint_tZq_cpqm.append("jprime_Eta")
features_CARL_singlePoint_tZq_cpqm.append("recoLepTop_Pt")

#-- cpq3
features_CARL_singlePoint_tZq_cpq3 = []
# features_CARL_singlePoint_tZq_cpq3.append("recoZ_Pt")
features_CARL_singlePoint_tZq_cpq3.append("recoZ_Eta")
features_CARL_singlePoint_tZq_cpq3.append("cosThetaStarPolZ")
features_CARL_singlePoint_tZq_cpq3.append("dEta_bjprime")
features_CARL_singlePoint_tZq_cpq3.append("dR_Zjprime")

features_CARL_singlePoint_tZq_cpq3.append("mTW")
features_CARL_singlePoint_tZq_cpq3.append("dR_tZ")
features_CARL_singlePoint_tZq_cpq3.append("dR_blW")
features_CARL_singlePoint_tZq_cpq3.append("dR_jprimeClosestLep")
features_CARL_singlePoint_tZq_cpq3.append("lAsymmetry")
features_CARL_singlePoint_tZq_cpq3.append("recoLepTop_Eta")
features_CARL_singlePoint_tZq_cpq3.append("maxDiJet_Pt")

# features_CARL_singlePoint_tZq_cpq3.append("jPrimeAbsEta")
# features_CARL_singlePoint_tZq_cpq3.append("recoLepTop_Pt")
# features_CARL_singlePoint_tZq_cpq3.append("Mass_3l")

#-- cpt
features_CARL_singlePoint_tZq_cpt = []
features_CARL_singlePoint_tZq_cpt.append("recoZ_Pt")
features_CARL_singlePoint_tZq_cpt.append("recoZ_Eta")
features_CARL_singlePoint_tZq_cpt.append("recoZ_dPhill")
features_CARL_singlePoint_tZq_cpt.append("dEta_tjprime")
features_CARL_singlePoint_tZq_cpt.append("dR_bjprime")

# //--------------------------------------------

#-- CARL_singlePoint ttZ
features_CARL_singlePoint_ttZ_ctz = features_CARL_singlePoint_tZq_ctz
features_CARL_singlePoint_ttZ_ctw = features_CARL_singlePoint_tZq_ctw

# //--------------------------------------------

# NN-SM: tZq vs ttZ vs Others
features_SM = []
features_SM.append("recoZ_Pt")
features_SM.append("recoZ_Eta")
features_SM.append("mHT")
features_SM.append("mTW")
features_SM.append("Mass_3l")
features_SM.append("recoZ_dPhill")
features_SM.append("lAsymmetry")
features_SM.append("jPrimeAbsEta")
features_SM.append("maxEtaJet")
features_SM.append("maxDeepJet")
features_SM.append("nbjets")
features_SM.append("maxDiJet_M")
features_SM.append("dR_blW")
features_SM.append("dR_lWjprime")
# features_SM.append("njets") #FIXME -- remove because of PrivMC_tZq discrepancy ?

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
