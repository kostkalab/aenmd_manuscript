###Retrieve basic stats for AENMD output.
###Used for at least the analysis of gnomAD or Clinvar in R.
#' @param data GRangesList object returned by annotate_nmd() funciton from package aenmd.
#'            Each object corrresponds to an ENST and within this objcet, each "sequence" corresponds to a SNP that affects the parent ENST
#' @param canonical_tx_set List. List of ENSTs that are from the ensembl canonical transcript set, code to create is stored in /functions/Retrieve_CanonTxSet.R
#' @param clinvar_extended_analysis Logical. Standard data frames or the data frames for the additional clinsig segmented analysis of the clinvar dataset?
#' @return List of 2 to 4 data frames containing basic stats from aenmd analysis. This list changes based on whether clinvar_extended_analysis = F (default) or
#'            clinvar_extended_analysis = T. 
#'        For clinvar_extended_analysis = F :
#'            1. Var_Tx_DF -Stats for variant transcript pairs, 
#'            2. Uni_Var_DF - Stats by unique variant
#'        For clinvar_extended_analysis = T :
#'            1. ext_clinvar_Uni_Var_DF - basic stats after filtering the results to remove any clinsig annotations that are not benign, likely benign,
#'               pathogenic, likely pathogenic, conflicting or uncertain. 
#'            2. Clinsig_benign_DF - analysis of clinsig benign variants, 
#'            3. Clinsig_confl_DF - analysis of clinsig conflicting or uncertain variants, 
#'            4. Clinsig_patho_DF - analysis of clinsig pathogenic variants.
Basic_Data_Table <- function(data, canonical_tx_set, var_in_coding, clinvar_extended_analysis = F) {
	##___________
	##___________
	##___________
	##___________
	#initialize not in:
	`%ni%` <- Negate(`%in%`)

	#Unlist and sort:
	aenmd_results <-sort(GenomeInfoDb::sortSeqlevels(unlist(data)))

	##___________
	##___________
	##___________
	##___________
	if (clinvar_extended_analysis == F) {
		print("Creating data tables for var-tx pairs and unique vairants for AENMD annotation.")
	} else if (clinvar_extended_analysis == T) {
		print("Creating a data table for the additional analysis of Clinvar AENMD annotated unique vairants." )
		#filter clinvar, removing variants with clinsig annotations that are not of interest.
		aenmd_results <- aenmd_results[(aenmd_results$clinsig != ".") & (aenmd_results$clinsig != "confers_sensitivity") & 
                (aenmd_results$clinsig != "association") & (aenmd_results$clinsig != "protective") &  
                (aenmd_results$clinsig != "Affects")   &   (aenmd_results$clinsig != "drug_response")   &   
                (aenmd_results$clinsig != "not_provided")   &  (aenmd_results$clinsig != "risk_factor") & 
                (aenmd_results$clinsig != "other")]
	} else {
		stop("Invalid entry. Options are: TRUE or FALSE")
	}

	##___________
	##___________
	##___________
	##___________
	#Create indexes for PTCs, NMDesc, NMDdeg based on all rules
	ind_PTC <- S4Vectors::mcols(aenmd_results)$res_aenmd[,"is_ptc"]
	ind_NMDesc <- S4Vectors::mcols(aenmd_results)$res_aenmd[,"is_ptc"] & (rowSums(S4Vectors::mcols(aenmd_results)$res_aenmd[1:6] |> as.matrix()) >1 )
	ind_NMDdeg <- S4Vectors::mcols(aenmd_results)$res_aenmd[,"is_ptc"] & (rowSums(S4Vectors::mcols(aenmd_results)$res_aenmd[1:6] |> as.matrix()) == 1 )

	#Create indexes for NMDesc, NMDdeg based canonical NMDesc rules:
	ind_NMDesc_canRules <- S4Vectors::mcols(aenmd_results)$res_aenmd[,"is_ptc"] & (rowSums(S4Vectors::mcols(aenmd_results)$res_aenmd[2:3] |> as.matrix()) >=1 )
	ind_NMDdeg_canRules <- S4Vectors::mcols(aenmd_results)$res_aenmd[,"is_ptc"] & 
	                        (rowSums(S4Vectors::mcols(aenmd_results)$res_aenmd[2:3] |> as.matrix()) == 0 )

	#Create indexes for NMDesc, NMDdeg based noncanonical NMDesc rules:
	ind_NMDesc_nonCanRules <- S4Vectors::mcols(aenmd_results)$res_aenmd[,"is_ptc"] & (rowSums(S4Vectors::mcols(aenmd_results)$res_aenmd[4:6] |> as.matrix()) >=1 )
	ind_NMDdeg_nonCanRules <- S4Vectors::mcols(aenmd_results)$res_aenmd[,"is_ptc"] & 
	                        (rowSums(S4Vectors::mcols(aenmd_results)$res_aenmd[4:6] |> as.matrix()) == 0 )

	#Create index for the BP proximal rules:
	ind_BPproxrules <- S4Vectors::mcols(aenmd_results)$res_aenmd[,"is_ptc"]  & (rowSums(S4Vectors::mcols(aenmd_results)$res_aenmd[3:4] |> as.matrix()) >=1 )

	#Isolate the PTCs and nonPTCs:
	PTCs_gr                  <- aenmd_results[ind_PTC]
	NonPTC_gr                <- aenmd_results[!ind_PTC]

	#Isolate the NMDesc and NMDdeg using all rules:
	PTCs_NMDesc_gr           <- aenmd_results[ind_NMDesc]
	PTCs_NMDdeg_gr           <- aenmd_results[ind_NMDdeg]

	#Isolate the NMDesc and NMDdeg using canonical rules:
	PTCs_NMDesc_CanRules_gr  <- aenmd_results[ind_NMDesc_canRules]
	PTCs_NMDdeg_CanRules_gr  <- aenmd_results[ind_NMDdeg_canRules]

	#Isolate the NMDesc and NMDdeg using noncanonical rules:
	PTCs_NMDesc_NonCanRules_gr  <- aenmd_results[ind_NMDesc_nonCanRules]
	PTCs_NMDdeg_NonCanRules_gr  <- aenmd_results[ind_NMDdeg_nonCanRules]


	#Total PTC-transcript pairs:
	num_PtcTxPair <- length(PTCs_gr)
	#Total non-PTC-transcript pairs:
	num_nonPtcTxPair <- length(NonPTC_gr)

	#function to get the adjusted number of unique PTCs and the number of Tx dependent PTCs
	PTC_analysis <- function(ptcs_gr, nonptcs_gr) {
		##Find the unique var numbers: 
		#1. Transcript dependent calls for PTC vs non PTC
		num_TxDepndnt_PTC_1  <- length(which(unique(ptcs_gr$id) %in% unique(nonptcs_gr$id)))
		num_TxDepndnt_PTC_2  <- length(which(unique(nonptcs_gr$id) %in% unique(ptcs_gr$id)))
		num_TxDepndnt_PTC_1 == num_TxDepndnt_PTC_2 #should be the same

		#2. Variants that are called enitrely one way or another regardless of tx
		if (num_TxDepndnt_PTC_1 == num_TxDepndnt_PTC_2) { #should be the same 
			#Total unique PTCs = # unique PTCs - # Variants whose PTC status is Tx Dependent.
			num_uni_PTC <- length(unique(ptcs_gr$id)) - num_TxDepndnt_PTC_1
			#Total unique non-PTCs = # unique PTCs - # Variants whose PTC status is Tx Dependent.
			num_uni_nonPTC  <- length(unique(nonptcs_gr$id)) - num_TxDepndnt_PTC_1
		} else {stop("Something is wrong with the tx dependent PTC numbers")}

		df <- data.frame( num_uni_PTC = num_uni_PTC, num_TxDepndnt_PTC = num_TxDepndnt_PTC_1)

		return(df)

	}
	#running PTC analysis 
	PTC_analysis_res <- PTC_analysis(PTCs_gr, NonPTC_gr)

	NMDesc_NMDdeg_analysis  <- function (PTCs_gr, PTCs_NMDesc_gr, PTCs_NMDdeg_gr) {

		#NMD escape predicted variant-transcript pairs:
		num_NMDesc_PtcTxPair <- length(PTCs_NMDesc_gr)
		#NMDdegraded variant-transcript pairs:
		num_NMDdeg_PtcTxPair_1 <- length(PTCs_gr) - length(PTCs_NMDesc_gr)
		num_NMDdeg_PtcTxPair_2 <- length(PTCs_NMDdeg_gr)
		num_NMDdeg_PtcTxPair_1 == num_NMDdeg_PtcTxPair_2 #should be the same

		#Find transcript dependent calls for NMDesc vs NMDdeg
		num_TxDepndnt_NMDesc_1  <- length(which(unique(PTCs_NMDesc_gr$id) %in% unique(PTCs_NMDdeg_gr$id)))
		num_TxDepndnt_NMDesc_2  <- length(which(unique(PTCs_NMDdeg_gr$id) %in% unique(PTCs_NMDesc_gr$id)))
		num_TxDepndnt_NMDesc_1 == num_TxDepndnt_NMDesc_2 #should be the same

		if (num_TxDepndnt_NMDesc_1 == num_TxDepndnt_NMDesc_2) { #should be the same 
			#Total NMD escape variant count:
			num_uni_NMDescPTC <- length(unique(PTCs_NMDesc_gr$id)) - num_TxDepndnt_NMDesc_1
			#NMDdegraded unique variant count:
			num_uni_NMDdegPTC <- length(unique(PTCs_NMDdeg_gr$id)) - num_TxDepndnt_NMDesc_1
		} 

		df <- data.frame(  num_NMDesc_PtcTxPair = num_NMDesc_PtcTxPair, num_NMDdeg_PtcTxPair = num_NMDdeg_PtcTxPair_1, 
			num_uni_NMDescPTC = num_uni_NMDescPTC, num_uni_NMDdegPTC = num_uni_NMDdegPTC,
			num_TxDepndnt_NMDesc = num_TxDepndnt_NMDesc_1)

		return(df)
	}
	#Running NMD call analysis
	NMDesc_NMDdeg_analysis_res <- NMDesc_NMDdeg_analysis(PTCs_gr, PTCs_NMDesc_gr, PTCs_NMDdeg_gr)

	#take an unlisted and sorted granges list and split it by the id. then get a distribution and mean of transcripts per PTC
	PTCs_gr_split <- split(PTCs_gr, f=PTCs_gr$id)
	table(lengths(PTCs_gr_split))
	av_Tx_per_Var <- mean(lengths(PTCs_gr_split))

	##-Getting data for all var-tx pairs that are insertion deletions >=3bp called as NMDesc by the BP  -##
	##- proximal rules to get an idea for how imapctful our incorperation of indel size is to NMDesc  	-##
	##- predicitons. In the future, annotate indels for their variant AA position and see how many   	-##
	##- would be called as NMDesc incorrectly if indel size was not considered.							-##
	#isolate all indels:
	PTCs_InDel_gr <- PTCs_gr[which(PTCs_gr$type == "ins" | PTCs_gr$type == "del")]
	#isolate all indels >=3bp
	PTCs_Indel3bpPlus_gr <- PTCs_gr[which(stringr::str_count(as.character(PTCs_gr$ref)) >=3 | stringr::str_count(as.character(PTCs_gr$alt)) >=3)]
	#ind for bp proximal rules
	ind_Indel3bpPlus_BPproxrules <- (rowSums(S4Vectors::mcols(PTCs_Indel3bpPlus_gr)$res_aenmd[3:4] |> as.matrix()) >=1 )
	#get indel >=3bp that NMDesc by BP proximal rules:
	PTCs_BPproxrules_Indel3bpPlus_gr <- (PTCs_Indel3bpPlus_gr[ind_Indel3bpPlus_BPproxrules])

	num_PTCs_BPproxrules_Indel3bpPlus <- length(PTCs_BPproxrules_Indel3bpPlus_gr)

	##___________
	##___________
	##___________
	##___________
	#########
	#####Numbers by rule:#########
	#########

	##___________
	#########
	##-Repeat analysis for canonical rules
	#########
	
	#Running NMD call analysis for canonical rules
	NMDesc_NMDdeg_CanonRules_analysis_res <- NMDesc_NMDdeg_analysis(PTCs_gr, PTCs_NMDesc_CanRules_gr, PTCs_NMDdeg_CanRules_gr)
	
	##-Look at overall numbers for both caonical rules in isloate
	#########
	######- Penult50:
	penult50 <- PTCs_gr[PTCs_gr$res_aenmd$is_penultimate == T]
	#Var-tx pairs
	num_penult50_PtcTxPair <- length(penult50)
	#Unique var - THIS HAS OVERLAP WITH ALL THE OTHER RULES
	num_uniq_penult50 <- length(unique(penult50$id))

	#########
	######- FinalExon:
	finalExon   <- PTCs_gr[PTCs_gr$res_aenmd$is_last ==T]
	#Var-tx pairs
	num_finalExon_PtcTxPair <- length(finalExon)
	#Unique var - THIS HAS OVERLAP WITH ALL THE OTHER RULES
	num_uniq_finalExon <- length(unique(finalExon$id))

	##___________
	#########
	##-Repeat analysis for noncanonical rules
	#########

	#Running NMD call analysis for noncanonical rules
	NMDesc_NMDdeg_NonCanonRules_analysis_res <- NMDesc_NMDdeg_analysis(PTCs_gr, PTCs_NMDesc_NonCanRules_gr, PTCs_NMDdeg_NonCanRules_gr)
	
	##-Look at overall numbers for both caonical rules in isloate
	#########
	######- >407bp:
	bp407       <- PTCs_gr[PTCs_gr$res_aenmd$is_407plus == T]
	#Var-tx pairs
	num_bp407_PtcTxPair <- length(bp407)
	#Unique var - THIS HAS OVERLAP WITH ALL THE OTHER RULES
	num_uniq_bp407 <- length(unique(bp407$id))

	#########
	######- Start proximal:
	CSSprox     <- PTCs_gr[PTCs_gr$res_aenmd$is_cssProximal == T]
	#Var-tx pairs
	num_CSSprox_PtcTxPair <- length(CSSprox)
	#Unique var - THIS HAS OVERLAP WITH ALL THE OTHER RULES
	num_uniq_CSSprox <- length(unique(CSSprox$id))

	#########
	######- Single exon:
	singleExon  <- PTCs_gr[PTCs_gr$res_aenmd$is_single == T]
	#Var-tx pairs
	num_singleExon_PtcTxPair <- length(singleExon)
	#Unique var - THIS HAS OVERLAP WITH ALL THE OTHER RULES
	num_uniq_singleExon <- length(unique(singleExon$id))

	##___________
	#########
	######- Finally, get NMDesc annotations per var-tx pair and unique variant by dividing the total from all rules and dividing by total unique NMDesc
	#########
	#Total unique NMDesc
	num_total_uniq_NMDescPTC <- length(unique(PTCs_NMDesc_gr$id))

	#Total NMDesc annotations:
	num_total_NMDescVarTxPair_eachRule <- (num_penult50_PtcTxPair + num_finalExon_PtcTxPair + num_bp407_PtcTxPair + num_CSSprox_PtcTxPair + num_singleExon_PtcTxPair)
	num_total_uniq_NMDescPTC_eachRule  <- (num_uniq_penult50 + num_uniq_finalExon + num_uniq_bp407 + num_uniq_CSSprox + num_uniq_singleExon) 

	#Do the math:
	num_NMDescAnnot_perNMD_VarTX_pair <- num_total_NMDescVarTxPair_eachRule / NMDesc_NMDdeg_analysis_res$num_NMDesc_PtcTxPair 
	num_NMDescAnnot_perUni_NMDesc_Var <- num_total_uniq_NMDescPTC_eachRule / num_total_uniq_NMDescPTC

	##___________
	##___________
	##___________
	##___________
	#########
	##repeat analysis for just the high quality transcripts that are also canonical transcripts
	#########

	#Filter all the PTCs and nonPTCs
	PTCs_gr_canontxFiltered <- PTCs_gr[PTCs_gr$tx_id %in% canonical_tx_set$tx_id]
	nonPTCs_gr_canontxFiltered <- NonPTC_gr[NonPTC_gr$tx_id %in% canonical_tx_set$tx_id]

	#Get all the NMDdeg and NMDesc var-tx pairs
	ind_NMDescTxPair_canontx <- (rowSums(S4Vectors::mcols(PTCs_gr_canontxFiltered)$res_aenmd[1:6] |> as.matrix()) >1 )
	ind_NMDdegTxPair_canontx <- (rowSums(S4Vectors::mcols(PTCs_gr_canontxFiltered)$res_aenmd[1:6] |> as.matrix()) ==1 )
	PTCs_NMDescTxPair_canontx <- PTCs_gr_canontxFiltered[ind_NMDescTxPair_canontx]
	PTCs_NMDdegTxPair_canontx <- PTCs_gr_canontxFiltered[ind_NMDdegTxPair_canontx]

	#Total PTC-transcript pairs:
	num_PtcTxPair_canontx <- length(PTCs_gr_canontxFiltered)
	#Total non-PTC-transcript pairs:
	num_nonPtcTxPair_canontx <- length(nonPTCs_gr_canontxFiltered)

	#NMD escape predicted variant-transcript pairs:s
	num_NMDesc_PtcTxPair_canontx <- length(PTCs_NMDescTxPair_canontx)
	#NMD degraded predicted variant-transcript pairs:
	num_NMDdeg_PtcTxPair_canontx_1 <- length(PTCs_gr_canontxFiltered) - length(PTCs_NMDescTxPair_canontx)
	num_NMDdeg_PtcTxPair_canontx_2 <- length(PTCs_NMDdegTxPair_canontx)
	num_NMDdeg_PtcTxPair_canontx_1 == num_NMDdeg_PtcTxPair_canontx_2 #should be the same

	#take an unlisted and sorted granges list and split it by the id. then get a distribution and mean of transcripts per PTC
	#Distribution of the number of tx per variant
	PTCs_gr_canontxFiltered_split <- split(PTCs_gr_canontxFiltered, f = PTCs_gr_canontxFiltered$id)
	table(lengths(PTCs_gr_canontxFiltered_split))
	av_Tx_per_Var_canontx <- mean(lengths(PTCs_gr_canontxFiltered_split))


	#unique var dt numbers, we want this regardless of whether we are just doing standard analysis, or clinvar extended analysis:
	Uni_Var_DF <- data.frame( functional_class = c("AENMD_annotated_PTC_causing", "Transcript_dependent_PTC", "NMD_escaping", "NMD_triggering", 
													            "Transcript_dependent",
													            "NMD_escaping_canonical_rules", "NMD_triggering_canonical_rules", "Transcript_dependent_canonical_rules",
													            "NMD_escaping_noncanonical_rules", "NMD_triggering_noncanonical_rules", "Transcript_dependent_noncanonical_rules",
													            "NMD_escaping_>407bp_rule", "NMD_escaping_start_proximal_rule", "NMD_escaping_single_exon_rule", 
													            "NMD_escaping_annotation_per_NMD_escape_variant"),
                  		          value        = format(c(PTC_analysis_res$num_uni_PTC, PTC_analysis_res$num_TxDepndnt_PTC, NMDesc_NMDdeg_analysis_res$num_uni_NMDescPTC, NMDesc_NMDdeg_analysis_res$num_uni_NMDdegPTC,
                  		          			NMDesc_NMDdeg_analysis_res$num_TxDepndnt_NMDesc,
                  		          			NMDesc_NMDdeg_CanonRules_analysis_res$num_uni_NMDescPTC, NMDesc_NMDdeg_CanonRules_analysis_res$num_uni_NMDdegPTC, NMDesc_NMDdeg_CanonRules_analysis_res$num_TxDepndnt_NMDesc,
                  		          		  NMDesc_NMDdeg_NonCanonRules_analysis_res$num_uni_NMDescPTC, NMDesc_NMDdeg_NonCanonRules_analysis_res$num_uni_NMDdegPTC, NMDesc_NMDdeg_NonCanonRules_analysis_res$num_TxDepndnt_NMDesc,
                  		          			num_uniq_bp407, num_uniq_CSSprox, num_uniq_singleExon,
                  		          			num_NMDescAnnot_perUni_NMDesc_Var), 
                  		          		scientific = F, digits = 2),
													      percent      = format(c( PTC_analysis_res$num_uni_PTC / (var_in_coding) *100, PTC_analysis_res$num_TxDepndnt_PTC / (var_in_coding) *100, NMDesc_NMDdeg_analysis_res$num_uni_NMDescPTC / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100, NMDesc_NMDdeg_analysis_res$num_uni_NMDdegPTC / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100,
													            NMDesc_NMDdeg_analysis_res$num_TxDepndnt_NMDesc / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100,
													            NMDesc_NMDdeg_CanonRules_analysis_res$num_uni_NMDescPTC / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) * 100, NMDesc_NMDdeg_CanonRules_analysis_res$num_uni_NMDdegPTC / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100,  NMDesc_NMDdeg_CanonRules_analysis_res$num_TxDepndnt_NMDesc / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100,
													            NMDesc_NMDdeg_NonCanonRules_analysis_res$num_uni_NMDescPTC / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100, NMDesc_NMDdeg_NonCanonRules_analysis_res$num_uni_NMDdegPTC  / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100, NMDesc_NMDdeg_NonCanonRules_analysis_res$num_TxDepndnt_NMDesc / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100,
													            num_uniq_bp407 / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100, num_uniq_CSSprox / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100,  num_uniq_singleExon / (PTC_analysis_res$num_uni_PTC + PTC_analysis_res$num_TxDepndnt_PTC) *100,
													            NA),
													        scientific = F, digits = 2)
													      )


	#now we branch according to whether we do basic analysis or extended analysis:
	if (clinvar_extended_analysis == F) {
		Var_Tx_DF <- data.frame(  functional_class = c("AENMD_annotated_PTC_causing", "NMD_escaping", "NMD_triggering", "Transcripts_per_variant", 
													          "NMD_escaping_canonical_rules", "NMD_triggering_canonical_rules",
													          "NMD_escaping_noncanonical_rules", "NMD_triggering_noncanonical_rules",
												          	"NMD_escaping_>407bp_rule", "NMD_escaping_start_proximal_rule", "NMD_escaping_single_exon_rule",
												          	"NMD_escaping_annotation_per_NMD_escape_variant_tx_pair",
												          	"PTC_causing_canonical_transcripts", "NMD_escaping_canon_tx", "NMD_triggering_canon_tx", "Transcripts_per_variant_canon_tx"),
                  		          value            = format(c(num_PtcTxPair, NMDesc_NMDdeg_analysis_res$num_NMDesc_PtcTxPair, NMDesc_NMDdeg_analysis_res$num_NMDdeg_PtcTxPair, av_Tx_per_Var, 
                  		    					 NMDesc_NMDdeg_CanonRules_analysis_res$num_NMDesc_PtcTxPair, NMDesc_NMDdeg_CanonRules_analysis_res$num_NMDdeg_PtcTxPair, 
                  		    					 NMDesc_NMDdeg_NonCanonRules_analysis_res$num_NMDesc_PtcTxPair, NMDesc_NMDdeg_NonCanonRules_analysis_res$num_NMDdeg_PtcTxPair,
                  		    					 num_bp407_PtcTxPair, num_CSSprox_PtcTxPair, num_singleExon_PtcTxPair,
                  		    					 num_NMDescAnnot_perNMD_VarTX_pair,
                  		    					 num_PtcTxPair_canontx, num_NMDesc_PtcTxPair_canontx, num_NMDdeg_PtcTxPair_canontx_1, av_Tx_per_Var_canontx), 
                  		          		scientific = F, digits = 2),
													      percent            = format(c(NA, NMDesc_NMDdeg_analysis_res$num_NMDesc_PtcTxPair / num_PtcTxPair *100, NMDesc_NMDdeg_analysis_res$num_NMDdeg_PtcTxPair / num_PtcTxPair *100, NA, 
													           NMDesc_NMDdeg_CanonRules_analysis_res$num_NMDesc_PtcTxPair / num_PtcTxPair *100, NMDesc_NMDdeg_CanonRules_analysis_res$num_NMDdeg_PtcTxPair / num_PtcTxPair *100, 
													           NMDesc_NMDdeg_NonCanonRules_analysis_res$num_NMDesc_PtcTxPair / num_PtcTxPair *100, NMDesc_NMDdeg_NonCanonRules_analysis_res$num_NMDdeg_PtcTxPair / num_PtcTxPair *100,
													           num_bp407_PtcTxPair / num_PtcTxPair *100, num_CSSprox_PtcTxPair / num_PtcTxPair *100, num_singleExon_PtcTxPair / num_PtcTxPair *100,
													           NA,
													           NA, num_NMDesc_PtcTxPair_canontx / num_PtcTxPair_canontx * 100, num_NMDdeg_PtcTxPair_canontx_1 / num_PtcTxPair_canontx * 100, NA), 
													         scientific = F, digits = 2)													           
													    )

		
		print("Done.")

		return(list(Var_Tx_DF = Var_Tx_DF, Uni_Var_DF = Uni_Var_DF))

	} else if (clinvar_extended_analysis == T) {
		
		print("Continuing to extended analysis...")

		#Creating function for per clinsig annotation class  analysis....
		Clinvar_Clinsig_Analysis <- function(clinsig_annotation, PTCs = PTCs_gr, NonPTC = NonPTC_gr, PTCs_NMDesc = PTCs_NMDesc_gr, PTCs_NMDdeg = PTCs_NMDdeg_gr, 
											PTCs_NMDesc_CanRules = PTCs_NMDesc_CanRules_gr, PTCs_NMDdeg_CanRules = PTCs_NMDdeg_CanRules_gr, 
											PTCs_NMDesc_NonCanRules = PTCs_NMDesc_NonCanRules_gr, PTCs_NMDdeg_NonCanRules = PTCs_NMDdeg_NonCanRules_gr) {
			if (tolower(clinsig_annotation) == "benign") {
				###Benign analysis:
				PTCs_clinsig <- ( PTCs_gr[stringr::str_detect(PTCs_gr$clinsig, "(?i)benign")] )

				#############################
				##All Rules NMDesc and NMDdeg
				#NMDesc
				PTCs_NMDesc_clinsig <- ( PTCs_NMDesc_gr[stringr::str_detect(PTCs_NMDesc_gr$clinsig, "(?i)benign")] )
				#NMDdeg
				PTCs_NMDdeg_clinsig <- ( PTCs_NMDdeg_gr[stringr::str_detect(PTCs_NMDdeg_gr$clinsig, "(?i)benign")] )

				#############################
				##Canonical NMDesc and NMDdeg rules
				#NMDesc
				PTCs_NMDesc_CanRules_clinsig <- ( PTCs_NMDesc_CanRules_gr[stringr::str_detect(PTCs_NMDesc_CanRules_gr$clinsig, "(?i)benign")] ) 
				#NMDdeg
				PTCs_NMDdeg_CanRules_clinsig <- ( PTCs_NMDdeg_CanRules_gr[stringr::str_detect(PTCs_NMDdeg_CanRules_gr$clinsig, "(?i)benign")] ) 

				#############################
				##Nonanonical NMDesc and NMDdeg rules
				#NMDesc
				PTCs_NMDesc_NonCanRules_clinsig <- ( PTCs_NMDesc_NonCanRules_gr[stringr::str_detect(PTCs_NMDesc_NonCanRules_gr$clinsig, "(?i)benign")] ) 
				#NMDdeg
				PTCs_NMDdeg_NonCanRules_clinsig <- ( PTCs_NMDdeg_NonCanRules_gr[stringr::str_detect(PTCs_NMDdeg_NonCanRules_gr$clinsig, "(?i)benign")] ) 

			} else if (tolower(clinsig_annotation) == "pathogenic") {
				###Pathogenic and likely pathogenic analysis:
				PTCs_clinsig <- ( PTCs_gr[stringr::str_detect(PTCs_gr$clinsig, "(?i)pathogenic") & !stringr::str_detect(PTCs_gr$clinsig, "(?i)conflicting")] )

				#############################
				##All Rules NMDesc and NMDdeg
				#NMDesc
				PTCs_NMDesc_clinsig <- ( PTCs_NMDesc_gr[stringr::str_detect(PTCs_NMDesc_gr$clinsig, "(?i)pathogenic") & !stringr::str_detect(PTCs_NMDesc_gr$clinsig, "(?i)conflicting")] )
				#NMDdeg
				PTCs_NMDdeg_clinsig <- ( PTCs_NMDdeg_gr[stringr::str_detect(PTCs_NMDdeg_gr$clinsig, "(?i)pathogenic") & !stringr::str_detect(PTCs_NMDdeg_gr$clinsig, "(?i)conflicting")] )

				#############################
				##Canonical NMDesc and NMDdeg rules
				#NMDesc
				PTCs_NMDesc_CanRules_clinsig <- ( PTCs_NMDesc_CanRules_gr[stringr::str_detect(PTCs_NMDesc_CanRules_gr$clinsig, "(?i)pathogenic") & !stringr::str_detect(PTCs_NMDesc_CanRules_gr$clinsig, "(?i)conflicting")] ) 
				#NMDdeg
				PTCs_NMDdeg_CanRules_clinsig <- ( PTCs_NMDdeg_CanRules_gr[stringr::str_detect(PTCs_NMDdeg_CanRules_gr$clinsig, "(?i)pathogenic") & !stringr::str_detect(PTCs_NMDdeg_CanRules_gr$clinsig, "(?i)conflicting")] ) 

				#############################
				##Nonanonical NMDesc and NMDdeg rules
				#NMDesc
				PTCs_NMDesc_NonCanRules_clinsig <- ( PTCs_NMDesc_NonCanRules_gr[stringr::str_detect(PTCs_NMDesc_NonCanRules_gr$clinsig, "(?i)pathogenic") & !stringr::str_detect(PTCs_NMDesc_NonCanRules_gr$clinsig, "(?i)conflicting")] ) 
				#NMDdeg
				PTCs_NMDdeg_NonCanRules_clinsig <- ( PTCs_NMDdeg_NonCanRules_gr[stringr::str_detect(PTCs_NMDdeg_NonCanRules_gr$clinsig, "(?i)pathogenic") & !stringr::str_detect(PTCs_NMDdeg_NonCanRules_gr$clinsig, "(?i)conflicting")] ) 

			} else if (tolower(clinsig_annotation) == "uncertain" | tolower(clinsig_annotation) == "conflicting") {
				###Uncertain or conflicating analysis:
				PTCs_clinsig <- ( PTCs_gr[stringr::str_detect(PTCs_gr$clinsig, "(?i)conflicting_interpretations_of_pathogenicity") | stringr::str_detect(PTCs_gr$clinsig, "(?i)uncertain")] )

				#############################
				##All Rules NMDesc and NMDdeg
				#NMDesc
				PTCs_NMDesc_clinsig <- ( PTCs_NMDesc_gr[stringr::str_detect(PTCs_NMDesc_gr$clinsig, "(?i)conflicting_interpretations_of_pathogenicity") | stringr::str_detect(PTCs_NMDesc_gr$clinsig, "(?i)uncertain")] )
				#NMDdeg
				PTCs_NMDdeg_clinsig <- ( PTCs_NMDdeg_gr[stringr::str_detect(PTCs_NMDdeg_gr$clinsig, "(?i)conflicting_interpretations_of_pathogenicity") | stringr::str_detect(PTCs_NMDdeg_gr$clinsig, "(?i)uncertain")] )

				#############################
				##Canonical NMDesc and NMDdeg rules
				#NMDesc
				PTCs_NMDesc_CanRules_clinsig <- ( PTCs_NMDesc_CanRules_gr[stringr::str_detect(PTCs_NMDesc_CanRules_gr$clinsig, "(?i)conflicting_interpretations_of_pathogenicity") | stringr::str_detect(PTCs_NMDesc_CanRules_gr$clinsig, "(?i)uncertain")] ) 
				#NMDdeg
				PTCs_NMDdeg_CanRules_clinsig <- ( PTCs_NMDdeg_CanRules_gr[stringr::str_detect(PTCs_NMDdeg_CanRules_gr$clinsig, "(?i)conflicting_interpretations_of_pathogenicity") | stringr::str_detect(PTCs_NMDdeg_CanRules_gr$clinsig, "(?i)uncertain")] ) 

				#############################
				##Nonanonical NMDesc and NMDdeg rules
				#NMDesc
				PTCs_NMDesc_NonCanRules_clinsig <- ( PTCs_NMDesc_NonCanRules_gr[stringr::str_detect(PTCs_NMDesc_NonCanRules_gr$clinsig, "(?i)conflicting_interpretations_of_pathogenicity") | stringr::str_detect(PTCs_NMDesc_NonCanRules_gr$clinsig, "(?i)uncertain")] ) 
				#NMDdeg
				PTCs_NMDdeg_NonCanRules_clinsig <- ( PTCs_NMDdeg_NonCanRules_gr[stringr::str_detect(PTCs_NMDdeg_NonCanRules_gr$clinsig, "(?i)conflicting_interpretations_of_pathogenicity") | stringr::str_detect(PTCs_NMDdeg_NonCanRules_gr$clinsig, "(?i)uncertain")] )

			} else {stop("Invalid Clinsig Option given")}

			##Find the unique var numbers: 
			PTCs_clinsig_analysis <- PTC_analysis(PTCs_clinsig, NonPTC_gr)

			NMDesc_NMDdeg_clinsig_analysis_res <- NMDesc_NMDdeg_analysis(PTCs_clinsig, PTCs_NMDesc_clinsig, PTCs_NMDdeg_clinsig)

			NMDesc_NMDdeg_CanRules_clinsig_analysis_res <- NMDesc_NMDdeg_analysis(PTCs_clinsig, PTCs_NMDesc_CanRules_clinsig, PTCs_NMDdeg_CanRules_clinsig)

			NMDesc_NMDdeg_NonCanRules_clinsig_analysis_res <- NMDesc_NMDdeg_analysis(PTCs_clinsig, PTCs_NMDesc_NonCanRules_clinsig, PTCs_NMDdeg_NonCanRules_clinsig)

			######- >407bp:
			bp407_clinsig               <- PTCs_NMDesc_clinsig[PTCs_NMDesc_clinsig$res_aenmd$is_407plus == T]
			#Unique var - THIS HAS OVERLAP WITH ALL THE OTHER RULES
			num_uniq_bp407_clinsig      <- length(unique(bp407_clinsig$id))

			#########
			######- Start proximal:
			CSSprox_clinsig             <- PTCs_NMDesc_clinsig[PTCs_NMDesc_clinsig$res_aenmd$is_cssProximal == T]
			#Unique var - THIS HAS OVERLAP WITH ALL THE OTHER RULES
			num_uniq_CSSprox_clinsig    <- length(unique(CSSprox_clinsig$id))

			#########
			######- Single exon:
			singleExon_clinsig           <- PTCs_NMDesc_clinsig[PTCs_NMDesc_clinsig$res_aenmd$is_single == T]
			#Unique var - THIS HAS OVERLAP WITH ALL THE OTHER RULES
			num_uniq_singleExon_clinsig  <- length(unique(singleExon_clinsig$id))


			#Create the DF for this clinsig annotation
			Uni_Var_DF_byClinsig <- data.frame(functional_class = c("PTC", "Transcript_dependent_PTC", "NMD_escaping", "NMD_triggering", "Transcript_dependent", 
																                "NMD_escaping_canonical_rules", "NMD_triggering_canonical_rules", "Transcript_dependent_canonical_rules",
																                "NMD_escaping_noncanonical_rules", "NMD_triggering_noncanonical_rules", "Transcript_dependent_noncanonical_rules",
																                "NMD_escaping_>407bp_rule", "NMD_escaping_start_proximal_rule", "NMD_escaping_single_exon_rule"),
	                  		                   value          = format(c(PTCs_clinsig_analysis$num_uni_PTC, PTCs_clinsig_analysis$num_TxDepndnt_PTC, 
	                  		                   			 NMDesc_NMDdeg_clinsig_analysis_res$num_uni_NMDescPTC, NMDesc_NMDdeg_clinsig_analysis_res$num_uni_NMDdegPTC,
	                  		          					 		 NMDesc_NMDdeg_clinsig_analysis_res$num_TxDepndnt_NMDesc,
	                  		          							 NMDesc_NMDdeg_CanRules_clinsig_analysis_res$num_uni_NMDescPTC, NMDesc_NMDdeg_CanRules_clinsig_analysis_res$num_uni_NMDdegPTC, 
	                  		          							 NMDesc_NMDdeg_CanRules_clinsig_analysis_res$num_TxDepndnt_NMDesc,
	                  		          					 		 NMDesc_NMDdeg_NonCanRules_clinsig_analysis_res$num_uni_NMDescPTC, NMDesc_NMDdeg_NonCanRules_clinsig_analysis_res$num_uni_NMDdegPTC, 
	                  		          					 		 NMDesc_NMDdeg_NonCanRules_clinsig_analysis_res$num_TxDepndnt_NMDesc,
	                  		          					 		 num_uniq_bp407_clinsig, num_uniq_CSSprox_clinsig, num_uniq_singleExon_clinsig),
	                  		    					 scientific = F, digits = 2),
	                  		    					     percent       = format(c(PTCs_clinsig_analysis$num_uni_PTC / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100, PTCs_clinsig_analysis$num_TxDepndnt_PTC / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100, 
	                  		    					           NMDesc_NMDdeg_clinsig_analysis_res$num_uni_NMDescPTC / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100, NMDesc_NMDdeg_clinsig_analysis_res$num_uni_NMDdegPTC / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100,
	                  		    					           NMDesc_NMDdeg_clinsig_analysis_res$num_TxDepndnt_NMDesc / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100,
	                  		    					           NMDesc_NMDdeg_CanRules_clinsig_analysis_res$num_uni_NMDescPTC / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100, NMDesc_NMDdeg_CanRules_clinsig_analysis_res$num_uni_NMDdegPTC / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100, 
	                  		    					           NMDesc_NMDdeg_CanRules_clinsig_analysis_res$num_TxDepndnt_NMDesc / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100,
	                  		    					           NMDesc_NMDdeg_NonCanRules_clinsig_analysis_res$num_uni_NMDescPTC / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100, NMDesc_NMDdeg_NonCanRules_clinsig_analysis_res$num_uni_NMDdegPTC / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100, 
	                  		    					           NMDesc_NMDdeg_NonCanRules_clinsig_analysis_res$num_TxDepndnt_NMDes / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100,
	                  		    					           num_uniq_bp407_clinsig / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100, num_uniq_CSSprox_clinsig / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100, num_uniq_singleExon_clinsig / (PTCs_clinsig_analysis$num_uni_PTC + PTCs_clinsig_analysis$num_TxDepndnt_PTC) * 100),
	                  		    					scientific = F, digits = 2) )
			return(Uni_Var_DF_byClinsig)
		}
		
		#Function created, using it for analysis....
		##___________
		##___________
		##___________
		##___________
		###Benign:
		Clinsig_benign <- Clinvar_Clinsig_Analysis("benign")

		##___________
		##___________
		##___________
		##___________
		###Uncertain / Conflicting:
		Clinsig_confl <- Clinvar_Clinsig_Analysis("conflicting")

		##___________
		##___________
		##___________
		##___________
		###Pathogenic w/o Conflicting:
		Clinsig_patho <- Clinvar_Clinsig_Analysis("pathogenic")

		return(list(ext_clinvar_Uni_Var_DF = Uni_Var_DF, Clinsig_benign_DF = Clinsig_benign, Clinsig_confl_DF = Clinsig_confl, Clinsig_patho_DF =  Clinsig_patho))

	}

}