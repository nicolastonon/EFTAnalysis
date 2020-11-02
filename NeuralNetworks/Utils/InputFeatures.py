# //--------------------------------------------

#-- CARL_singlePoint tZq
features_CARL_singlePoint_tZq = []
features_CARL_singlePoint_tZq.append("recoZ_Pt")
features_CARL_singlePoint_tZq.append("recoZ_Eta")
features_CARL_singlePoint_tZq.append("jprime_Pt")
features_CARL_singlePoint_tZq.append("jPrimeAbsEta")
features_CARL_singlePoint_tZq.append("maxDelPhiLL")
features_CARL_singlePoint_tZq.append("cosThetaStarPolZ")

# features_CARL_singlePoint_tZq.append("dR_bW")
# features_CARL_singlePoint_tZq.append("mHT")
# features_CARL_singlePoint_tZq.append("Mass_3l")
# features_CARL_singlePoint_tZq.append("TopZsystem_M")
# features_CARL_singlePoint_tZq.append("recoLepTop_Pt")
# features_CARL_singlePoint_tZq.append("dR_tClosestJet")

# //--------------------------------------------

#-- CARL_singlePoint ttZ
features_CARL_singlePoint_ttZ = features_CARL_singlePoint_tZq

# //--------------------------------------------

# NN-SM: tZq vs ttZ vs Others
features_SM = []
features_SM.append("recoZ_Pt")
features_SM.append("recoZ_Eta")
features_SM.append("mHT")
features_SM.append("mTW")
features_SM.append("Mass_3l")
features_SM.append("maxDelPhiLL")
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
list_features.append("maxDelPhiLL")

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
