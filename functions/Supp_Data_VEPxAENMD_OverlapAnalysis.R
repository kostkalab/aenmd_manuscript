############____________________________________#########
############____________________________________#########
############____________________________________#########
#   For comparing AENMD vs VEP results of analysis of   #
# clinvar. Both AENMD and VEP annotated clinvar results #
# is REQUIRED for running this, also transcript         #
# information for AENMD's base transcript set is        # 
# required, see /functions/Retrieve_AENMD_BaseTxSet.R   #
############____________________________________#########
############____________________________________#########
############____________________________________#########
#' @param clinvar_vep_out Tab delimited file with vcf like header. Results from VEP annoatation of Clinvar, 
#'                  including the NMD plugin (see /external_dataset/vep).
#' @param aenmdVEPsettings_clinvar_res GRangesList object returned by AENMD, with VEP settings (100bp start proximal rule),
#'                  from annotation of Clinvar (see /external_dataset/clinvar).
#' @param AENMD_base_transcriptset GRanges object. Granges for high-quality ENSTs analyzed by base AENMD with additional descriptive
#'                  information (see /functions/Retrieve_AENMD_BaseTxSetR).
#' @return List of 3 data.frames.
#'                  1. VEP_AENMD_clinvar_res_merged ; Data frame containing results from VEPxAENMD NMD annotation comparison of clinvar
#'              stop gain variant results. This is  also used by Supp_Data_VEPxAENMD_DiffNMDescCalls() .
#'                  2. Supp_DT_11_AENMDxVEP_AllVarTxPairs ; Data frame containing results of how many total variant-transcript pairs 
#'              were analyzed by AENMD and basic stats of the variant-transcript pairs analyzed by VEP NMD plugin.
#'                  3. Supp_DT_12_AENMDxVEPmerged_StopGainVarTxPairs: Data frame containing basic stats of the overlapping 
#'              variant-transcript pairs that were annotated by both AENMD and VEP.
#' @output 5 tab delimited data tables output to ./functions/accessory_files (1) or ./sup_data/raw_datatables (2,3,4,5)
#'                  1. VEPxAENMD_Clinvar_ResultsMerged.txt - All the overlapping variant-transcript pairs from clinvar with annotations 
#'                  from VEP and AENMD. $VEP_NMDesc - NMDesc calls by VEP ; $AENMD_NMDesc_VEPrules - NMD calls made by AENMD.
#'                  THE ONLY DATA TABLE WRITTEN TO ./functions/accessory_files 
#'                  2. VEPxAENMD_Clinvar_ResultsMerged_BothNMDesc.txt - Subset of VEPxAENMD_Clinvar_ResultsMerged.txt ; 
#'              All the variant-transcript pairs from clinvar called as NMD escape by both VEP and AENMD.
#'                  3. VEPxAENMD_Clinvar_ResultsMerged_BothNMDdeg.txt - Subset of VEPxAENMD_Clinvar_ResultsMerged.txt ; 
#'              All the variant-transcript pairs from clinvar called as not escaping NMD by both VEP and AENMD.
#'                  4. VEPxAENMD_Clinvar_ResultsMerged_AENMD_NMDesc_notVEP.txt - Subset of VEPxAENMD_Clinvar_ResultsMerged.txt ; 
#'              All the variant-transcript pairs from clinvar called as NMD escape by AENMD but not called as NMD escaping by VEP.
#'                  5. VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD.txt - Subset of VEPxAENMD_Clinvar_ResultsMerged.txt ;
#'              All the variant-transcript pairs from clinvar called as NMD escape by VEP but not called as NMD escaping by AENMD.
#'              
Supp_Data_VEPxAENMD_OverlapAnalysis <- function(clinvar_vep_out, aenmdVEPsettings_clinvar_res, AENMD_base_transcriptset) {
 
  ##___________________
  ##___________________
  ##___________________
  ###Starting with the AENMD results:
  #unlist and sort:
  aenmdVEPsettings_clinvar_res <-sort(GenomeInfoDb::sortSeqlevels(unlist(aenmdVEPsettings_clinvar_res)))
  #isolate just the PTC calls
  ind_PTC <- S4Vectors::mcols(aenmdVEPsettings_clinvar_res)$res_aenmd[,"is_ptc"]
  aenmdVEPsettings_clinvar_res_PTCs <- aenmdVEPsettings_clinvar_res[ind_PTC]
  #reformat to DF to allow for merging with VEP results:
  aenmdVEPsettings_clinvar_res <- as.data.frame(aenmdVEPsettings_clinvar_res, row.names = 1:length(aenmdVEPsettings_clinvar_res))
  aenmdVEPsettings_clinvar_res_PTCs <- as.data.frame(aenmdVEPsettings_clinvar_res_PTCs, row.names = 1:length(aenmdVEPsettings_clinvar_res_PTCs))

  #overall number of NMDesc annotated PTCs
  df_num_aenmdVEPsettings_clinvar_res_PTCs <- data.frame( functional_class = "AENMD_VEPsettings_Clinvar_AllPTCcausing",
                                                          value = nrow(aenmdVEPsettings_clinvar_res_PTCs), percent = NA )

  ##___________________
  ##___________________
  ##___________________
  ###Next, the Clinvar results:
  #filter for only Tx we we use for base AENMD analysis
  clinvar_vep_out_TSLtxOnly <- clinvar_vep_out %>% dplyr::filter(Feature %in% AENMD_base_transcriptset$tx_id)

  ##-filter out splice variants
  clinvar_vep_out_TSLtxOnly_noSplice <- clinvar_vep_out_TSLtxOnly %>% dplyr::filter(!stringr::str_detect(Consequence, "splice"))

  #expand info
  library(tidyr)
  colnames(clinvar_vep_out_TSLtxOnly_noSplice)[1] <- stringr::str_remove(colnames(clinvar_vep_out_TSLtxOnly_noSplice)[1],'#')
  clinvar_vep_out_TSLtxOnly_noSplice              <- janitor::clean_names(clinvar_vep_out_TSLtxOnly_noSplice)
  clinvar_vep_out_TSLtxOnly_noSplice <- clinvar_vep_out_TSLtxOnly_noSplice %>% dplyr::distinct() %>% separate_rows(extra, sep=";") %>%
    tidyr::separate(extra, into = c("key", "value"), sep="=") %>%
    tidyr::spread(key,value) %>%
    janitor::clean_names()

  #-Change the name of VEP's NMD escape annotation column to something more informative:
  clinvar_vep_out_TSLtxOnly_noSplice <- clinvar_vep_out_TSLtxOnly_noSplice %>% dplyr::mutate("VEP_NMDesc" = nmd) %>% dplyr::select(-nmd)
  
  ##___________________
  ##___________________
  ##___________________
  ##-Now Isolate Stop Gains, the only annotated variant type.
  clinvar_vep_out_stopgain_TSLtxOnly_noSplice <- clinvar_vep_out_TSLtxOnly_noSplice %>% 
    dplyr::filter(stringr::str_detect(consequence,"stop_gained"))

  ##___________________
  ##___________________
  ##___________________
  #Collecting stats for the unfiltered Clinvar Results.

  #Function that will be used twice, once for the entire VEP SG, and another time on the VEP-AENMD Overlap.
  AENMD_VEP_comparison_DT_fx <- function(VEP_annotation_containing_DF, VEP_or_AENMD = "VEP") {
    
    VEP_or_AENMD <- tolower(VEP_or_AENMD)

    if ( VEP_or_AENMD == "vep") {
      ###count all 
      num_clinvar_SG <- nrow(VEP_annotation_containing_DF)
      ##Strand overall:
      num_clinvar_posStrand <- nrow( VEP_annotation_containing_DF %>% dplyr::filter(strand == "1") )
      num_clinvar_negStrand <- nrow( VEP_annotation_containing_DF %>% dplyr::filter(strand == "-1") )

      prcnt_clinvar_posStrand <- ( nrow(VEP_annotation_containing_DF %>% dplyr::filter(strand == "1") ) / 
                                  nrow(VEP_annotation_containing_DF) * 100 )
      prcnt_clinvar_negStrand <- ( nrow(VEP_annotation_containing_DF %>% dplyr::filter(strand == "-1") ) / 
                                  nrow(VEP_annotation_containing_DF) * 100 )
    } else {print("AENMD data analysis, not analyzing the entire merge data set again.")}
    
    #set param based on user input
    NMDesc_fltr_param <- switch(VEP_or_AENMD, vep = "VEP_NMDesc == \'NMD_escaping_variant\'", aenmd = "AENMD_NMDesc_VEPrules == T")
    NMDdeg_fltr_param <- switch(VEP_or_AENMD, vep = "is.na(VEP_NMDesc) == T", aenmd = "AENMD_NMDesc_VEPrules == F")

    ##___________________
    ###NMDesc
    clinvar_NMDesc <- VEP_annotation_containing_DF %>% dplyr::filter( eval(rlang::parse_expr(NMDesc_fltr_param)) )
    num_clinvar_NMDesc <- nrow( clinvar_NMDesc )

    prcnt_clinvar_NMDesc <- ( nrow( clinvar_NMDesc ) / 
                                  nrow(VEP_annotation_containing_DF) * 100 )
    ##NMDesc strand:
    num_clinvar_NMDesc_posStrand <- nrow( clinvar_NMDesc %>% dplyr::filter(strand == "1") )
    num_clinvar_NMDesc_negStrand <- nrow( clinvar_NMDesc %>% dplyr::filter(strand == "-1") )

    prcnt_clinvar_NMDesc_posStrand <- ( nrow(clinvar_NMDesc %>% dplyr::filter(strand == "1") ) / 
                                  nrow(clinvar_NMDesc) * 100 )
    prcnt_clinvar_NMDesc_negStrand <- ( nrow(clinvar_NMDesc %>% dplyr::filter(strand == "-1") ) / 
                                  nrow(clinvar_NMDesc) * 100 )
    
    ##___________________
    #NMDdeg
    clinvar_NMDdeg <- VEP_annotation_containing_DF %>% dplyr::filter( eval(rlang::parse_expr(NMDdeg_fltr_param)) )
    num_clinvar_NMDdeg <- nrow( clinvar_NMDdeg )

    prcnt_clinvar_NMDdeg <-  ( nrow( clinvar_NMDdeg ) / 
                                    nrow(VEP_annotation_containing_DF) * 100 )
    #NMDdeg strand:
    num_clinvar_NMDdeg_posStrand <- nrow( clinvar_NMDdeg %>% dplyr::filter(strand == "1") )
    num_clinvar_NMDdeg_negStrand <- nrow( clinvar_NMDdeg %>% dplyr::filter(strand == "-1") )

    prcnt_clinvar_NMDdeg_posStrand <- ( nrow(clinvar_NMDdeg %>% dplyr::filter(strand == "1") ) / 
                                  nrow(clinvar_NMDdeg) * 100 )
    prcnt_clinvar_NMDdeg_negStrand <- ( nrow(clinvar_NMDdeg %>% dplyr::filter(strand == "-1") ) / 
                                  nrow(clinvar_NMDdeg) * 100 )

    #output:
    if ( VEP_or_AENMD == "vep") {
      df <- data.frame( functional_class = c( "vep_clinvar_stopgains_total", 
        "vep_clinvar_posStrand_total",  
        "vep_clinvar_negStrand_total", 
        "vep_clinvar_NMDesc_total", 
        "vep_clinvar_NMDesc_posStrand_total", 
        "vep_clinvar_NMDesc_negStrand_total",
        "vep_clinvar_NMDdeg_total", 
        "vep_clinvar_NMDdeg_posStrand_total", 
        "vep_clinvar_NMDdeg_negStrand_total"),
        value          = format( c( num_clinvar_SG, 
          num_clinvar_posStrand, 
          num_clinvar_negStrand,
          num_clinvar_NMDesc,
          num_clinvar_NMDesc_posStrand,
          num_clinvar_NMDesc_negStrand, 
          num_clinvar_NMDdeg,
          num_clinvar_NMDdeg_posStrand,
          num_clinvar_NMDdeg_negStrand ),
                                       scientific = F, digits = 2),  
        percent       = format( c( NA, 
          (prcnt_clinvar_posStrand),
          (prcnt_clinvar_negStrand),
          (prcnt_clinvar_NMDesc), 
          (prcnt_clinvar_NMDesc_posStrand), 
          (prcnt_clinvar_NMDesc_negStrand),
          (prcnt_clinvar_NMDdeg),
          (prcnt_clinvar_NMDdeg_posStrand),
          (prcnt_clinvar_NMDdeg_negStrand) ),
                                scientific = F, digits = 2) )
      } else if (VEP_or_AENMD == "aenmd") {
        df <- data.frame( functional_class = c( 
        "AENMD_overlapStopgains_NMDesc_total", 
        "AENMD_overlapStopgains_NMDesc_posStrand_total",  
        "AENMD_overlapStopgains_NMDesc_negStrand_total",
        "AENMD_overlapStopgains_NMDdeg_total",
        "AENMD_overlapStopgains_NMDdeg_posStrand_total", 
        "AENMD_overlapStopgains_NMDdeg_negStrand_total" ),
        value          = format( c(
          num_clinvar_NMDesc,
          num_clinvar_NMDesc_posStrand,
          num_clinvar_NMDesc_negStrand, 
          num_clinvar_NMDdeg,
          num_clinvar_NMDdeg_posStrand, 
          num_clinvar_NMDdeg_negStrand ),
                                scientific = F, digits = 2),
        percent          = format( c(
          (prcnt_clinvar_NMDesc), 
          (prcnt_clinvar_NMDesc_posStrand), 
          (prcnt_clinvar_NMDesc_negStrand),
          (prcnt_clinvar_NMDdeg),
          (prcnt_clinvar_NMDdeg_posStrand),
          (prcnt_clinvar_NMDdeg_negStrand) ),
                                scientific = F, digits = 2) )
      } else { stop("Please enter VEP or AENMD, case insensitive.") }
    return(df)
  }
  #Run the function to create the DT for all the variants (stop gains) that VEP analyzes:
  clinvar_vep_out_stopgain_all_DT <- AENMD_VEP_comparison_DT_fx(clinvar_vep_out_stopgain_TSLtxOnly_noSplice)

  ##___________________
  ##___________________
  ##___________________
  ##Stop gain analysis
  ##Now isolate all variant-transcript pairs where the variant is also in the aenmd dataset
  #1. intersection between VEP stop gains and all of the variants AENMD output (PTC causing SNVs and all indels + substiutions).
  SG_inVEP_and_Var_inAENMD  <-  clinvar_vep_out_stopgain_TSLtxOnly_noSplice[which(clinvar_vep_out_stopgain_TSLtxOnly_noSplice$uploaded_variation %in% aenmdVEPsettings_clinvar_res$id),]
  #2. intersection  between VEP stop gains and all of the variants AENMD output at PTC causing.
  SG_inVEP_and_AENMD  <-  clinvar_vep_out_stopgain_TSLtxOnly_noSplice[which(clinvar_vep_out_stopgain_TSLtxOnly_noSplice$uploaded_variation %in% aenmdVEPsettings_clinvar_res_PTCs$id),]

  ##1 and 2 should be equal:
  #If so, create the merged VEP and AENMD results
  if (identical(SG_inVEP_and_AENMD, SG_inVEP_and_Var_inAENMD)) {
    #-Shorten the dataframes to make them more readable:
    SG_inVEP_and_AENMD_shrtnd <- SG_inVEP_and_AENMD %>% dplyr::select(uploaded_variation,feature,consequence,c_dna_position,cds_position,protein_position, amino_acids,hgv_sp,VEP_NMDesc,strand)
    res_clinvar_PTC_shrtnd <- aenmdVEPsettings_clinvar_res_PTCs %>% dplyr::select(id, tx_id, res_aenmd.is_ptc, res_aenmd.is_last, res_aenmd.is_penultimate, res_aenmd.is_cssProximal, res_aenmd.is_single, res_aenmd.is_407plus, clinsig, type)
    
    #-merge them
    VEP_AENMD_clinvar_res_merged <- merge(SG_inVEP_and_AENMD_shrtnd, res_clinvar_PTC_shrtnd, by.x = c("uploaded_variation", "feature"), by.y = c("id", "tx_id"), all.x = TRUE)
    
    #-Get another column for the joining of the AENMD rules that VEP uses:
    VEP_AENMD_clinvar_res_merged$AENMD_NMDesc_VEPrules <- ifelse(VEP_AENMD_clinvar_res_merged$res_aenmd.is_last == T | 
                                                         VEP_AENMD_clinvar_res_merged$res_aenmd.is_penultimate ==T | 
                                                         VEP_AENMD_clinvar_res_merged$res_aenmd.is_cssProximal == T | 
                                                         VEP_AENMD_clinvar_res_merged$res_aenmd.is_single == T, 
                                                       TRUE, FALSE)
  } else {stop("There seems to be variants that VEP calls PTC that AENMD does not, this shouldnt be the case")}


  ##___________________
  ##___________________
  ##___________________
  ###Collecting stats.

  ##Individual dataset calls:

  ##___________________
  ##___________________
  #VEP and AENMD overlap, as well as just VEP:
  clinvar_VEP_AENMD_merged_DT_VEP <- AENMD_VEP_comparison_DT_fx(VEP_AENMD_clinvar_res_merged)
  #rename first column to notify this is the total overlap:
  clinvar_VEP_AENMD_merged_DT_VEP[1,1] <- "AENMDandVEP_overlap_PTC_causing"
  #rename the other columns to make it clear we are looking at the overlapping variants:
  clinvar_VEP_AENMD_merged_DT_VEP$functional_class <- stringr::str_replace(clinvar_VEP_AENMD_merged_DT_VEP$functional_class, "vep_clinvar_NMDesc", "VEP_overlapStopgains_NMDesc")
  clinvar_VEP_AENMD_merged_DT_VEP$functional_class <- stringr::str_replace(clinvar_VEP_AENMD_merged_DT_VEP$functional_class, "vep_clinvar_NMDdeg", "VEP_overlapStopgains_NMDdeg")
  clinvar_VEP_AENMD_merged_DT_VEP$functional_class <- stringr::str_replace(clinvar_VEP_AENMD_merged_DT_VEP$functional_class, "vep_clinvar", "AENMDandVEP_overlapStopgains")



  ##___________________
  ##___________________
  ##AENMD DT values:
  clinvar_VEP_AENMD_merged_DT_AENMD <- AENMD_VEP_comparison_DT_fx(VEP_AENMD_clinvar_res_merged, "aenmd")

  ##___________________
  ##___________________
    
  ##___________________
  ##Same calls:
   #both esc
  VEP_AENMD_clinvar_res_merged_BothNMDesc <- VEP_AENMD_clinvar_res_merged %>% dplyr::filter(is.na(VEP_NMDesc) == F, AENMD_NMDesc_VEPrules == T)
  num_VEPaenmdIntrsct_Both_NMDesc <- nrow(VEP_AENMD_clinvar_res_merged_BothNMDesc)
   #both deg
  VEP_AENMD_clinvar_res_merged_BothNMDdeg <- VEP_AENMD_clinvar_res_merged %>% dplyr::filter(is.na(VEP_NMDesc) == T, AENMD_NMDesc_VEPrules == F)
  num_VEPaenmdIntrsct_Both_NMDdeg <- nrow(VEP_AENMD_clinvar_res_merged_BothNMDdeg)

  ##___________________
  ##different calls:
    #AENMD_NMDesc_notVEP
  VEP_AENMD_clinvar_res_merged_AENMD_NMDesc_notVEP <- VEP_AENMD_clinvar_res_merged %>% dplyr::filter(is.na(VEP_NMDesc) == T, AENMD_NMDesc_VEPrules == T)
  num_VEPaenmdIntrsct_AENMD_NMDesc_notVEP <- nrow(VEP_AENMD_clinvar_res_merged_AENMD_NMDesc_notVEP)
   #VEP_NMDesc_notAENMD
  VEP_AENMD_clinvar_res_merged_VEP_NMDesc_notAENMD <- VEP_AENMD_clinvar_res_merged %>% dplyr::filter(is.na(VEP_NMDesc) == F, AENMD_NMDesc_VEPrules == F)
  num_VEPaenmdIntrsct_VEP_NMDesc_notAENMD <- nrow(VEP_AENMD_clinvar_res_merged_VEP_NMDesc_notAENMD)

  #Join into a DF
  df_VEPxAENMD_clinvar_merged_res <- data.frame( functional_class = c("VEP_AENMD_Intrsct_Both_NMDesc",  "VEP_AENMD_Intrsct_Both_NMDdeg",
                          "VEP_AENMD_Intrsct_AENMD_NMDesc_notVEP", "VEP_AENMD_Intrsct_VEP_NMDesc_notAENMD"),
                        value =   format( c(num_VEPaenmdIntrsct_Both_NMDesc, num_VEPaenmdIntrsct_Both_NMDdeg,
                                            num_VEPaenmdIntrsct_AENMD_NMDesc_notVEP, num_VEPaenmdIntrsct_VEP_NMDesc_notAENMD) ,
                               scientific = F, digits = 2),
                        percent = format( c(num_VEPaenmdIntrsct_Both_NMDesc / as.integer(clinvar_VEP_AENMD_merged_DT_VEP[1,2]) * 100 , num_VEPaenmdIntrsct_Both_NMDdeg / as.integer(clinvar_VEP_AENMD_merged_DT_VEP[1,2]) * 100,
                                            num_VEPaenmdIntrsct_AENMD_NMDesc_notVEP / as.integer(clinvar_VEP_AENMD_merged_DT_VEP[1,2]) * 100, num_VEPaenmdIntrsct_VEP_NMDesc_notAENMD / as.integer(clinvar_VEP_AENMD_merged_DT_VEP[1,2]) * 100) ,
                                scientific = F, digits = 2) ) 



  ##___________________
  ##___________________
  ##___________________
  ##formulate, save, and output the output:
  ##Save the raw data used to generate data as tab separated tables:
  write.table(VEP_AENMD_clinvar_res_merged, file="./functions/accessory_files/VEPxAENMD_Clinvar_ResultsMerged.txt", sep="\t", row.names=F, col.names = T,  quote = F)

  write.table(VEP_AENMD_clinvar_res_merged_BothNMDesc, file="./sup_data/raw_datatables/VEPxAENMD_Clinvar_ResultsMerged_BothNMDesc.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  write.table(VEP_AENMD_clinvar_res_merged_BothNMDdeg, file="./sup_data/raw_datatables/VEPxAENMD_Clinvar_ResultsMerged_BothNMDdeg.txt", sep="\t", row.names=F, col.names = T,  quote = F)

  write.table(VEP_AENMD_clinvar_res_merged_AENMD_NMDesc_notVEP, file="./sup_data/raw_datatables/VEPxAENMD_Clinvar_ResultsMerged_AENMD_NMDesc_notVEP.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  write.table(VEP_AENMD_clinvar_res_merged_VEP_NMDesc_notAENMD, file="./sup_data/raw_datatables/VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD.txt", sep="\t", row.names=F, col.names = T,  quote = F)

  ##Join data tables to create final output:
  Supp_DT_11_AENMDxVEP_AllVarTxPairs <- rbind(df_num_aenmdVEPsettings_clinvar_res_PTCs, clinvar_vep_out_stopgain_all_DT)
  Supp_DT_12_AENMDxVEPmerged_StopGainVarTxPairs <- rbind(clinvar_VEP_AENMD_merged_DT_VEP, clinvar_VEP_AENMD_merged_DT_AENMD, df_VEPxAENMD_clinvar_merged_res)

  return( list(VEP_AENMD_clinvar_res_merged = VEP_AENMD_clinvar_res_merged, Supp_DT_11_AENMDxVEP_AllVarTxPairs = Supp_DT_11_AENMDxVEP_AllVarTxPairs, Supp_DT_12_AENMDxVEPmerged_StopGainVarTxPairs = Supp_DT_12_AENMDxVEPmerged_StopGainVarTxPairs) )
  
  ##_____________________________________________
  ##_______End Output Generating code____________
  ##_____________________________________________
  ##_____________________________________________
  ##_______End Output Generating code____________
  ##_____________________________________________
  ##_____________________________________________
  ##_______End Output Generating code____________
  ##_____________________________________________


  ##_____________________________Supplemental code:____________________________________________________
  ##____frameshift - checking to see if all the frameshifts that standard VEP annotation_______________
  ##_________________finds is present in AENMD output. The below code does not yield any_______________
  ##_________________output and is present only for reference. In summary, all frameshifts_____________
  ##_________________that we would expect to be there, are there.______________________________________
  ##___________________________________________________________________________________________________

  ##### frameshifts:
  #get frameshifts without stop gains
  clinvar_vep_out_frameshift_TSLtxOnly_noSplice <- clinvar_vep_out_TSLtxOnly_noSplice %>% 
      dplyr::filter(stringr::str_detect(consequence,"frameshift_variant")) %>% 
      dplyr::filter(!stringr::str_detect(consequence,"stop_gained"))
  #check that they were isolated properly
  table(clinvar_vep_out_frameshift_TSLtxOnly_noSplice$consequence)

  #count all variant transcript pairs
  nrow(clinvar_vep_out_frameshift_TSLtxOnly_noSplice)

  ##________Now isolate all variant-transcript pairs where the variant is also in the aenmd dataset
  ##intersections between VEP frameshifts and all of the variants AENMD output (PTC causing variants and all indels).
  FS_inVEP_and_Var_inAENMD  <-  clinvar_vep_out_frameshift_TSLtxOnly_noSplice[which(clinvar_vep_out_frameshift_TSLtxOnly_noSplice$uploaded_variation %in% aenmdVEPsettings_clinvar_res$id),]

  FS_inVEP_and_Var_inAENMD_shrtnd <- FS_inVEP_and_Var_inAENMD %>% dplyr::select(uploaded_variation,feature,consequence,c_dna_position,cds_position,protein_position, amino_acids,hgv_sp,VEP_NMDesc,strand)
  #-Shorten the dataframes to make them more readable:
  res_clinvar_shrtnd <- aenmdVEPsettings_clinvar_res %>% dplyr::select(id, tx_id, res_aenmd.is_ptc, res_aenmd.is_last, res_aenmd.is_penultimate, res_aenmd.is_cssProximal, res_aenmd.is_single, res_aenmd.is_407plus, clinsig, type)
  #-merge them
  FS_VEP_AENMD_clinvar_res_merged <- merge(FS_inVEP_and_Var_inAENMD, res_clinvar_shrtnd, by.x = c("uploaded_variation", "feature"), by.y = c("id", "tx_id"), all.x = TRUE)

  #number of unique frameshifts found by VEP:
  length(unique(clinvar_vep_out_frameshift_TSLtxOnly_noSplice$uploaded_variation))
  #number of unique frameshifts overlapping:
  length(unique(FS_VEP_AENMD_clinvar_res_merged$uploaded_variation))

  #make protein length change a integrer:
  FS_VEP_AENMD_clinvar_res_merged$protein_length_change <- as.integer(FS_VEP_AENMD_clinvar_res_merged$protein_length_change)

  #Identify those that SHOULD be called PTCs but AENMD does not (none are found):
  aenmd_PTC_false <- FS_VEP_AENMD_clinvar_res_merged %>% dplyr::filter( res_aenmd.is_ptc == F )
  aenmd_PTC_false_shortenedProtein <- aenmd_PTC_false %>% dplyr::filter( protein_length_change < 0 ) 
  aenmd_PTC_false_shortenedProtein_notStartDisrupting <- aenmd_PTC_false_shortenedProtein %>% dplyr::filter( !stringr::str_detect(hgv_sp, "Met1") )
  aenmd_PTC_false_shortenedProtein_notStartDisrupting_notNonstop <- aenmd_PTC_false_shortenedProtein_notStartDisrupting %>% dplyr::filter(!stringr::str_detect(hgv_sp, "\\?") )

}