#_______________________________________________________#
#_______________________________________________________#
# !!! AENMD ANALYSIS OF CLINVAR AND VEP ANALYSIS OF !!! #
# !!!        CLINVAR MUST BE COMPLETED FIRST        !!! #
# !      AENMD via: ./external_dataset/clinvar      !!! #
# !                 VEP via: ./vep/                 !!! #
#_______________________________________________________#
#_______________________________________________________#
#  Description Code for comparing the results of        #
# annotating the VEP dataset,                           #
# clinvar_20221015_Clnsig.txt, with AENMD and VEP       #
# NMD plugin.                                           #
#_______________________________________________________#
#   There are 2 parts, 					                #
# - Supp_Data_VEPxAENMD_OverlapAnalysis, which gives    #
# basic stats as to the overlap and differential calls  #
# between the 2 softwares.                              #
# - Supp_Data_VEPxAENMD_DiffNMDescCalls, which goes in  #
# and analyzes the differential calls between the two   #
# softwares to find reasons for the differential calls  #
# and check for problematic calls in either program.    #
#_______________________________________________________#
#_______________________________________________________#

##Start with the overlap analysis of VEP and AENMD NMD escape annotation
##results that returns basic stats:
#source the function:
source("./functions/Supp_Data_VEPxAENMD_OverlapAnalysis.R")

#import the VEP results, the AENMD results, and transcript information required for both functions called here:
clinvar_vep_out <- readr::read_table("./functions/accessory_files/SGandFS_resetID_clinvar_20221211_vep_out.txt",comment = '##')
aenmdVEPsettings_20221211clinvar105_results <- readRDS("./functions/accessory_files/aenmdVEPsettings_20221211clinvar105_results.rds")

#Get the transcript information required for both functions called here:
AENMD_base_transcriptset <- future::value(aenmd:::._EA_txs_grl)
AENMD_base_exonset  <- future::value(aenmd:::._EA_exn_grl)

#Run the Supp_Data_VEPxAENMD_OverlapAnalysis() command, this automatically writes the tables as well:
VEPxAENMD_OverlapAnalysis_res <- Supp_Data_VEPxAENMD_OverlapAnalysis(clinvar_vep_out, aenmdVEPsettings_20221211clinvar105_results, AENMD_base_transcriptset)

#Write results csv
readr::write_csv(VEPxAENMD_OverlapAnalysis_res$Supp_DT_11_AENMDxVEP_AllVarTxPairs, file="./sup_data/DT_11_AENMDxVEP_AllVarTxPairs.csv")
readr::write_csv(VEPxAENMD_OverlapAnalysis_res$Supp_DT_12_AENMDxVEPmerged_StopGainVarTxPairs, file="./sup_data/DT_12_AENMDxVEPmerged_StopGainVarTxPairs.csv")

##Next run the differential NMDesc call analysis 
##that checks all the differential called variant-transcript pairs:
#source:
source("./functions/Supp_Data_VEPxAENMD_DiffNMDescCallsAnalysis.R")

#Use output from Supp_Data_VEPxAENMD_OverlapAnalysis() as well as the above imported AENMD_base_transcriptset
#To run the differential NMDesc analysis, this automaticall writes the analysis results tables as well:
VEPxAENMD_DiffNMDescCalls_res <- Supp_Data_VEPxAENMD_DiffNMDescCalls(VEPxAENMD_OverlapAnalysis_res$VEP_AENMD_clinvar_res_merged, AENMD_base_transcriptset, AENMD_base_exonset)

#Write cvs of esults of manual check
readr::write_csv(VEPxAENMD_DiffNMDescCalls_res$AnalyzedByHand_Var_AENMD_NMDesc_NotVEP, file="./sup_data/DT_13_AnalyzedByHand_Var_AENMD_NMDesc_NotVEP.csv")
readr::write_csv(VEPxAENMD_DiffNMDescCalls_res$AnalyzedByHand_VEP_NMDesc_notAENMD_SingleExonCheck, file="./sup_data/DT_14_AnalyzedByHand_VEP_NMDesc_notAENMD_SingleExonCheck.csv")