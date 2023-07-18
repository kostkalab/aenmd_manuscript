############____________________________________###########
############____________________________________###########
############____________________________________###########
## For reproducing and outputing all the random samples   ##
## used as examples for detailing discordant NMD escape   ##
## calls made by AENMD and VEP. Also outputs the results  ##
## from the automated checking of all the discordant      ##
## variant-transcript pairs for identifying the potential ##
## reason as to why the discordance is occurring.  				##
#
#' @param VEPxAENMD_Clinvar_ResultsMerged A data table created from an imported tab delimited text file with the results 
#'                  from merging AENMD and VEP results (see ./functions/Supp_Data_VEPxAENMD_OverlapAnalysis.R)
#' @param AENMD_base_transcriptset GRanges object. Granges for high-quality ENSTs analyzed by base AENMD with additional descriptive
#'                  information (see /functions/Retrieve_AENMD_BaseTxSetR).
#' @return List of 2 data frames. 
#'									[1] AnalyzedByHand_Var_AENMD_NMDesc_NotVEP , [2] AnalyzedByHand_VEP_NMDesc_notAENMD_SingleExonCheck;
#'                  Data frames containing a random sample of variant=-trainscript pairs with discordant NMDesc calls between
#'									VEP and AENMD that were analyzed by hand. 
#'										The - $Purpose - column tells user why the variant was included,
#'										The - $AENMD_NMDesc_VEPrules - this column joins the results from AENMD for the 4 NMDesc rules used by both AENMD
#'									and VEP (100 bp from coding start site, 50bp from the end of the penultimate coding exon boundary, final
#'									exon, single exon), 
#'										The - $VEP_NMDesc - column communicates the NMDesc result from VEP, and 
#'										The - $QC_Result - which communicates the result of the analysis done by a human.
#' @output 5 tab delimited data tables output to ./sup_data/raw_datatables. 
#'                  [3] AENMD_NMDesc_NotVEP_posStrand_LastExonCheck [4] AENMD_NMDesc_NotVEP_negStrand_LastExonCheck
#'									[5] AENMD_NMDesc_NotVEP_posStrand_PenultExonCheck [6] AENMD_NMDesc_NotVEP_negStrand_PenultExonCheck
#'									[7] VEP_NMDesc_notAENMD_cssProxCheck 
#'									[8] VEP_NMDesc_notAENMD_pos_LastExonCheck [9] VEP_NMDesc_notAENMD_neg_LastExonCheck 
#'									[10] VEP_NMDesc_notAENMD_pos_PenultExonCheck [11] VEP_NMDesc_notAENMD_neg_PenultExonCheck 
#'									[12] VEP_NMDesc_notAENMD_SingleExonCheck 
#'									are data frames containing analysis of ALL variant-transcript pairs with discordant NMDesc calls between VEP
#'									and AENMD using automated methods exhaustively checking for explainations for the discordant calls 
#'									(check_var_final_exon(), check_var_penult_exon(), and a couple lines checking for the single exon rule). Each
#'									automated method uses transcript builds from GRCh38p13 that are created by Retrieve_AENMD_DefaultTxSet.R
#'									of the aenmd.data.ensdb.v105 v0.2.1 package. For the CSS proxumial rule, the VEP 
#'									annotation for the variant's protein_position and cds_position is used. The title of the data.frame 
#'									communicates what is the discordance (VEP or AENMD calls the variant-transcript pair as NMDesc with the other
#'									not doing so), the strand of the variants -if necessary-, and the rule that was checked.
#'                            


Supp_Data_VEPxAENMD_DiffNMDescCalls <- function(VEPxAENMD_Clinvar_ResultsMerged, AENMD_base_transcriptset, exns_byTx) {
	##___________________________________________________________________________________
	##___________________________________________________________________________________
	##___________________________________________________________________________________
	##_______________________Starting the analysis of variants __________________________
	##________________that we discordantly called between AENMD and VEP__________________
	##___________________________________________________________________________________
	##___________________________________________________________________________________
	#Create the chr, pos, ref, alt columns from the uploaded variation column for the creation of 
	#the granges from the variants, which will be relevant later
	#Also seperate protein_position into start and end protein position check variants later
	VEPxAENMD_Clinvar_ResultsMerged <- VEPxAENMD_Clinvar_ResultsMerged %>% separate(uploaded_variation, into = c("chr", "pos", "ref", "alt"), sep = "_", remove = F) %>% 
																	 separate(protein_position, into = c("start_PP", "end_PP"), sep = "-" ,remove = F)

	VEPxAENMD_Clinvar_ResultsMerged$start_PP <- as.integer(VEPxAENMD_Clinvar_ResultsMerged$start_PP)
	VEPxAENMD_Clinvar_ResultsMerged$end_PP <- as.integer(VEPxAENMD_Clinvar_ResultsMerged$end_PP)

	#Creating a new column to communicate the purpose we are including the variant!!!!!:
	VEPxAENMD_Clinvar_ResultsMerged$Purpose <- NA
	VEPxAENMD_Clinvar_ResultsMerged$QC_Result <- NA

	##Isolate the variants that have different calls between AENMD and VEP:
	  #AENMD_NMDesc_notVEP
	VEPxAENMD_Clinvar_ResultsMerged_AENMD_NMDesc_notVEP <- VEPxAENMD_Clinvar_ResultsMerged %>% dplyr::filter(is.na(VEP_NMDesc) == T, AENMD_NMDesc_VEPrules == T) #(1,022)
	  #VEP_NMDesc_notAENMD
	VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD <- VEPxAENMD_Clinvar_ResultsMerged %>% dplyr::filter(is.na(VEP_NMDesc) == F, AENMD_NMDesc_VEPrules == F) #(31,206)

	#R

	##___________________
	##___________________
	##___________________
	##Start with the variants AENMD calls NMDesc but VEP does not.
	##___________________
	##___________________
	##___________________
	##Variants called cssProximal by AENMD but not NMDesc by VEP:
	Diff_cssProximal <- VEPxAENMD_Clinvar_ResultsMerged_AENMD_NMDesc_notVEP %>% dplyr::filter(res_aenmd.is_cssProximal == T) #(19)

	#random sample of 2 to process by hand to show examples. 
	#all of them are SNVs where VEP annotates them at cds position 102
	set.seed(321)
	index <- sample(1:nrow(Diff_cssProximal), 2)
	Diff_cssProximal_examples <- Diff_cssProximal[index, ]
	Diff_cssProximal_examples$Purpose <- "CSS proximal rule check."
	Diff_cssProximal_examples$QC_Result <- "Expected, within first 34 AA (AENMD) but not 100bp (VEP)"

	##___________________
	##___________________
	##___________________
	##Variants called in final coding exon by AENMD but not NMDesc by VEP:
	Diff_last <- VEPxAENMD_Clinvar_ResultsMerged_AENMD_NMDesc_notVEP %>% dplyr::filter(res_aenmd.is_last == T)

	#segmented by type:
	Diff_last_sbs <- Diff_last %>% dplyr::filter(type == "sbs") #(2)
	Diff_last_ins <- Diff_last %>% dplyr::filter(type == "ins") #(1)
	Diff_last_snv <- Diff_last %>% dplyr::filter(type == "snv") #(101)

	#Random sample for 1 SNV to be done by hand to show examples.
	#We then use code to analyze all SNVs altogether.
	#For sbs (2) and ins (1), since there is so few, we did it by hand.
	set.seed(321)
	index <- sample(1:nrow(Diff_last_snv), 1)
	Diff_last_snv_example <- Diff_last_snv[index, ]

	#join the data.frames from the 4 different checks:
	Diff_last_examples <- rbind(Diff_last_sbs,Diff_last_ins,Diff_last_snv_example)

	Diff_last_examples$Purpose <- "Final Exon rule check."
	Diff_last_examples$QC_Result <- "Expected, contained within the final exon."

	##-Analyzing the rest of the SNV calls using code:

	#function to check if the variant is in the final exon:
	check_var_final_exon <- function (variants_to_test, pos_or_neg, suppress_per_vartx_message = F) {

		results_final <- data.frame(uploaded_variation= c("0_0_0_X_X"),
	                                 feature = c("ENST"),
	                                 type = "", 
	                                 strand = "",
	                                 in_final = as.logical(NA),
	                                 distance_from_final_start = 0)
		#set up counter the number of var-tx pairs returned as not being in the final exon
		n = 0

		for(i in 1:nrow(variants_to_test)) {
		    test_gr <- GenomicRanges::GRanges( variants_to_test[i,]$chr,
		                                       IRanges::IRanges(as.double(variants_to_test[i,]$pos),
		                                                        as.double(variants_to_test[i,]$pos)))
		    #get the enst for the var-tx pair
		    test_enst <- variants_to_test[i,]$feature
		    #get the granges for the enst associated with the var-tx pair
		    reference <- exns_byTx[test_enst][[1]]
		    ind_reference_final <- (length(reference))
		    #check to see if the variant overlaps with the enst grange
		    ov <- GenomicRanges::findOverlaps(test_gr, reference)
		    #check to see if the variant overlaps with the final exon of the enst's grange
		    is_final <- S4Vectors::subjectHits(ov) == ind_reference_final
		    #output result, based on if the input variants were on the negative or postitive strand
		    if (pos_or_neg == "positive" | pos_or_neg == "pos") {
		    	if (is_final ==T ) {

		        	temp_df <- data.frame(uploaded_variation=variants_to_test[i,]$uploaded_variation,
	                                  feature = variants_to_test[i,]$feature,
	                                  type = variants_to_test[i,]$type,
	                                  strand = variants_to_test[i,]$strand,
	                                  in_final = TRUE,
	                                  distance_from_final_start = start(test_gr) - start(reference[ind_reference_final])
	                                  )
	       			 results_final <- rbind(results_final, temp_df)
		    	} else { 
		    		n = n + 1
		    		if (suppress_per_vartx_message == F) {message( paste0( "(", n, ")", variants_to_test[i,]$uploaded_variation, " is not in the final exon") )} 
		    	}

		    } else if (pos_or_neg == "negative" | pos_or_neg == "neg") {
		    	if (is_final ==T ) {

		    		temp_df <- data.frame(uploaded_variation=variants_to_test[i,]$uploaded_variation,
	                                  feature = variants_to_test[i,]$feature,
	                                  type = variants_to_test[i,]$type,
	                                  strand = variants_to_test[i,]$strand,
	                                  in_final = TRUE,
	                                  distance_from_final_start = end(reference[ind_reference_final]) - start(test_gr)
	                                  )
	       			 results_final <- rbind(results_final, temp_df)
	    		} else { 
	    			n = n + 1
	    			if (suppress_per_vartx_message == F) {message( paste0( "(", n, ")", variants_to_test[i,]$uploaded_variation, " is not in the final exon") )}
	    		}

		    } else { stop("ERROR: positive or negative are the options.") }
		}
		message( paste0( "(", n, ")", " out of ", nrow(variants_to_test),  " variants were found not in the final exon") )
		return(results_final)
	}

	#split by strand:

	#split the negative strand variants
	#then checking if the variant is in the final exon for negative strand:
	Diff_last_snv_negstrand <- Diff_last_snv %>% dplyr::filter(strand == "-1") #(51)
	neg_strand_last_res <- check_var_final_exon(Diff_last_snv_negstrand, "neg")

	neg_strand_last_res$Purpose <- "Final Exon rule check."
	neg_strand_last_res$QC_Result <- "Expected, contained within the final exon."

	#split the positive strand variants:
	#checking if the variant is in the final exon for positive strand:
	Diff_last_snv_posstrand <- Diff_last_snv %>% dplyr::filter(strand == "1") #(50)
	pos_strand_last_res <- check_var_final_exon(Diff_last_snv_posstrand, "pos")

	pos_strand_last_res$Purpose <- "Final Exon rule check."
	pos_strand_last_res$QC_Result <- "Expected, contained within the final exon."

	##___________________
	##___________________
	##___________________
	##Variants called in the final 50 bp of the penultimate exon by AENMD but not NMDesc by VEP:
	Diff_penultimate50 <- VEPxAENMD_Clinvar_ResultsMerged_AENMD_NMDesc_notVEP  %>% dplyr::filter(res_aenmd.is_penultimate == T)

	#segmented by type:
	Diff_penultimate50_snv <- Diff_penultimate50 %>% dplyr::filter(type == "snv") #(862)
	Diff_penultimate50_ins <- Diff_penultimate50 %>% dplyr::filter(type == "ins") #(28)
	Diff_penultimate50_del <- Diff_penultimate50 %>% dplyr::filter(type == "del") #(3)
	Diff_penultimate50_sbs <- Diff_penultimate50 %>% dplyr::filter(type == "sbs") #(6)

	#Random sample of 1-2 variants per variant type completed by hand to show examples. 
	#We then use code to the rest of the variants
	##########SNV
	set.seed(321)
	index <- sample(1:nrow(Diff_penultimate50_snv), 1)
	Diff_penultimate50_snv_example <- Diff_penultimate50_snv[index, ]

	##########INS
	set.seed(321)
	index <- sample(1:nrow(Diff_penultimate50_ins), 2)
	Diff_penultimate50_ins_example <- Diff_penultimate50_ins[index, ]

	##########DEL
	set.seed(321)
	index <- sample(1:nrow(Diff_penultimate50_del), 2)
	Diff_penultimate50_del_example <- Diff_penultimate50_del[index, ]

	##########SBS
	set.seed(321)
	index <- sample(1:nrow(Diff_penultimate50_sbs), 1)
	Diff_penultimate50_sbs_example <- Diff_penultimate50_sbs[index, ]

	Diff_penultimate50_examples <- rbind(Diff_penultimate50_snv_example, Diff_penultimate50_ins_example, Diff_penultimate50_del_example, Diff_penultimate50_sbs_example)
	Diff_penultimate50_examples$Purpose <- "50bp Penultimate Exon rule check."
	Diff_penultimate50_examples$QC_Result <- "Expected, contained final 50bp of the penultimate exon."

	##-Analyzing the rest of the variant calls using code:
	#function to check if the variant is in the final exon:
	check_var_penult_exon <- function (variants_to_test, pos_or_neg, suppress_per_vartx_message = F) {
		
		#setup the results table based on the strand of variants input
		if (pos_or_neg == "positive" | pos_or_neg == "pos") {

			results_penult <- data.frame(uploaded_variation = c("0_0_0_X_X"),
	                                 feature = c("ENST0"),
	                                 type = "", 
	                                 penultimate = as.logical(NA),
	                                 distance = 0,
	                                 adj_distance = 0,
	                                 amino_acids = "" )
		} else if (pos_or_neg == "negative" | pos_or_neg == "neg") {

			results_penult <- data.frame(uploaded_variation = c("0_0_0_X_X"),
	                                 			feature = c("ENST0"),
	                                 			type = "", 
	                                 			penultimate = as.logical(NA),
	                                 			distance = 0,
	                                 			amino_acids = "" )
		}

		#set up counter the number of var-tx pairs returned as not being in the final exon
		n = 0
		for(i in 1:nrow(variants_to_test)) {
			#create a gr for the variant-tx pair
		    test_gr <- GenomicRanges::GRanges( variants_to_test[i,]$chr,
		                                       IRanges::IRanges(as.double(variants_to_test[i,]$pos),
		                                                        as.double(variants_to_test[i,]$pos)))
		    #recover the ENST for the variant-tx pair
		    test_enst <- variants_to_test[i,]$feature
		    
		    #recover the gr for the ENST in question
		    reference <- exns_byTx[test_enst][[1]]
		    #find wthe index of the penultimate exon
		    ind_reference_penult <- (length(reference) - 1)
		    
		    #see if and where there is overlap between the var-tx pair and the tx
		    ov <- GenomicRanges::findOverlaps(test_gr, reference)
		    #check to see if the variant falls in the penultimate exon
		    is_penultimate <- S4Vectors::subjectHits(ov) == ind_reference_penult
		    
		    #output result, based on if the input variants were on the negative or postitive strand
		    if (pos_or_neg == "positive" | pos_or_neg == "pos") {
		    	if (is_penultimate ==T ) {

				    # If the variant is in the penultimate exon, return the variant, ENST, and how far the VARIANT (not PTC)
				    # is located from the end of the REFERENCE transcript (not alternative)
				    temp_df <- data.frame(uploaded_variation = variants_to_test[i,]$uploaded_variation,
				                            feature = variants_to_test[i,]$feature,
				                            type = variants_to_test[i,]$type,
				                            penultimate = TRUE,
				                            distance = end(reference[ind_reference_penult]) - start(test_gr),
				                            adj_distance = end(reference[ind_reference_penult]) - start(test_gr) - 
				                                  	stringr::str_length(variants_to_test[i,]$ref) + stringr::str_length(variants_to_test[i,]$alt),
				                            amino_acids =  variants_to_test[i,]$amino_acids
				                            ) #Here we include the adjusted distance metric to give insight into how AENMD sees deletion variants.
				        						#It does not give a good picture of how insertion and substitution variants, as the location of the 
				        						#PTC is not necessarily at the 5' most position. However, the real position of the insertion and 
				        						#substiution variants is likely between distance and adjusted distance.
				    results_penult <- rbind(results_penult, temp_df)

		    	} else { 
		    		n = n + 1
		    		if (suppress_per_vartx_message == F) { message( paste0( "(", n, ")", variants_to_test[i,]$uploaded_variation, " is not in the penultimate exon.") ) } 
		    	}

		    } else if (pos_or_neg == "negative" | pos_or_neg == "neg") {

		    	if (is_penultimate ==T ) {
		    	# If the variant is in the penultimate exon, return the variant, ENST, and how far the VARIANT (not PTC)
		    	# is located from the end of the REFERENCE transcript (not alternative)
		        temp_df <- data.frame(uploaded_variation = variants_to_test[i,]$uploaded_variation,
		                                  feature = variants_to_test[i,]$feature,
		                                  type = variants_to_test[i,]$type,
		                                  penultimate = TRUE,
		                                  distance = start(test_gr)  - start(reference[ind_reference_penult]),
		                                  amino_acids =  variants_to_test[i,]$amino_acids
		                                  ) #dont need an adjusted distance because the ref base pos of a variant on a negative strand
		       								#variant is the most 3' base of the cds, so all deletions delete bases more 5' of the pos
		       								#and the same goes for all insertions. Therefore, all the PTCs generated by deletions must
		       								#occour at that position while the alt pos of any insertions must be calculated; HOWEVER,
		       								#for insertions the distance, as calucated above, presents a "best case senerio" for VEP if
		       								#we are assuming VEP purely uses variant BP position as opposed to PTC variant position after
		       								#adjusting for bases in the insertion that might flank the PTC.
		        results_penult <- rbind(results_penult, temp_df)

	    		} else { 
	    			n = n + 1
	    			if (suppress_per_vartx_message == F) {message( paste0( "(", n, ")", variants_to_test[i,]$uploaded_variation, " is not in the penultimate exon.") )}
	    		}
		    } else { stop("ERROR: positive or negative are the options.") }
		}
		message( paste0( "(", n, ")", " out of ", nrow(variants_to_test),  " variants were found not in the penultimate exon.") )
		return(results_penult)
	}

	#split by strand -  NEGATIVE: 
	###!!! INSERTIONS AND SUBSTITUTIONS MUST BE ANALYZED SEPERATELY TO CHECK IF AENMD IS WORKING PROPERLY, BUT 
	###!!! GIVEN THAT VEP IS THEORETICALLY USING THE VARIANT POSITION AND NOT THE PTC POSITION, THE BELOW METHOD 
	###!!! SHOULD STILL TELL YOU WHETHER VEP SHOULD HAVE CALLED THE VARIANT AS NMD ESCAPING.
	Diff_penultimate50_neg <- Diff_penultimate50 %>% dplyr::filter(strand == "-1") #(856)
	neg_strand_penult50_res <- check_var_penult_exon(Diff_penultimate50_neg, "neg")

	#split by strand - POSITIVE:
	Diff_penultimate50_pos <- Diff_penultimate50 %>% dplyr::filter(strand == "1") #(43)
	pos_strand_penult50_res <- check_var_penult_exon(Diff_penultimate50_pos, "pos")

	#merge the data.tables of the variants for which was analyzed by hand
	#That AENMD has called NMDesc but VEP has not:
	AnalyzedByHand_Var_AENMD_NMDesc_NotVEP <- rbind(Diff_cssProximal_examples, Diff_last_examples, Diff_penultimate50_examples)

	##___________________
	##___________________
	##___________________
	##Next, the variants VEP calls NMDesc but AENMD does not.
	##___________________
	##___________________
	##___________________

	#Checking each:
	###1. Cds position within the first 100 or protein position within the 33 amino acids?
	###2. Final coding exon checking code.
	###3. Penultimate coding exon and distance from end of penultimate coding exon checking code.
	###4. Transcript a single exon transcript?

	#Seperate into 8 different groups, based on variant type and strand:
	VEPnotAENMD_SNV_pos <- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD %>% dplyr::filter(type == "snv", strand == "1") #(9)
	VEPnotAENMD_SNV_neg <- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD %>% dplyr::filter(type == "snv", strand == "-1") #(29732)

	VEPnotAENMD_ins_pos <- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD %>% dplyr::filter(type == "ins", strand == "1") #(7)
	VEPnotAENMD_ins_neg <- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD %>% dplyr::filter(type == "ins", strand == "-1") #(953)

	VEPnotAENMD_del_pos <- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD %>% dplyr::filter(type == "del", strand == "1") #(0)
	VEPnotAENMD_del_neg <- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD %>% dplyr::filter(type == "del", strand == "-1") #(237)

	VEPnotAENMD_sbs_pos <- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD %>% dplyr::filter(type == "sbs", strand == "1") #(0)
	VEPnotAENMD_sbs_neg <- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD %>% dplyr::filter(type == "sbs", strand == "-1") #(268)

	###1. Start proximal rule.
	VEP_NMDesc_notAENMD_cssProxCheck_AAbased <- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD  %>% dplyr::filter(start_PP <= 34) #(0 / 31,206 variants within the first 100bp or 33 amino acids)
	VEP_NMDesc_notAENMD_cssProxCheck_cdsBased<- VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD  %>% dplyr::filter(cds_position <= 100) #(0 / 31,206 variants within the first 100bp or 33 amino acids)

	VEP_NMDesc_notAENMD_cssProxCheck <- rbind(VEP_NMDesc_notAENMD_cssProxCheck_AAbased, VEP_NMDesc_notAENMD_cssProxCheck_cdsBased)

	###2. Final coding exon rule.

	#only SNV and ins on the pos strand
	pos_strand_SNV_last_res <- check_var_final_exon(VEPnotAENMD_SNV_pos, "pos") #(0 / 9 variants within the final exon)
	pos_strand_ins_last_res <- check_var_final_exon(VEPnotAENMD_ins_pos, "pos") #(0 / 7 variants within the final exon)
	#checking all on negative strand
	neg_strand_del_last_res <- check_var_final_exon(VEPnotAENMD_del_neg, "neg", suppress_per_vartx_message = T) #(0 / 237 variants within the final exon)
	neg_strand_sbs_last_res <- check_var_final_exon(VEPnotAENMD_sbs_neg, "neg", suppress_per_vartx_message = T) #(0 / 268 variants within the final exon)
	neg_strand_ins_last_res <- check_var_final_exon(VEPnotAENMD_ins_neg, "neg", suppress_per_vartx_message = T) #(0 / 953 variants within the final exon)
	neg_strand_SNV_last_res <- check_var_final_exon(VEPnotAENMD_SNV_neg, "neg", suppress_per_vartx_message = T) #(0 / 29,737 variants within the final exon)

	VEP_NMDesc_notAENMD_pos_LastExonCheck <- rbind(pos_strand_SNV_last_res, pos_strand_ins_last_res)
	VEP_NMDesc_notAENMD_neg_LastExonCheck <- rbind(neg_strand_del_last_res, neg_strand_sbs_last_res, neg_strand_ins_last_res, neg_strand_SNV_last_res)


	###3. Penultimate coding exon, final 50 bp

	#only SNV and ins on the pos strand
	pos_strand_SNV_penult_res <- check_var_penult_exon(VEPnotAENMD_SNV_pos, "pos") #(9 / 9 variants within the final exon)
	pos_strand_SNV_penult_res$Purpose <- "50bp Penultimate Exon rule check."
	pos_strand_SNV_penult_res$QC_Result <- "PTC does not overlap with the final 50bp."

	pos_strand_ins_penult_res <- check_var_penult_exon(VEPnotAENMD_ins_pos, "pos") #(7 / 7 variants within the final exon, all are expected given AENMD uses an alternative protein construction approach.
																														# This is because they are insertions that "push" the PTC out of range of the first 50bp of the penultimate
																														# coding exon when looking at teh alternative protein created by the insertion.)
	pos_strand_ins_penult_res$Purpose <- "50bp Penultimate Exon rule check."
	pos_strand_ins_penult_res$QC_Result <- "PTC is \"pushed\" out of penult rule by inserted bases."

	#checking all on negative strand
	neg_strand_del_penult_res <- check_var_penult_exon(VEPnotAENMD_del_neg, "neg", suppress_per_vartx_message = T) #(1 / 237 variants within the final exon)
	neg_strand_sbs_penult_res <- check_var_penult_exon(VEPnotAENMD_sbs_neg, "neg", suppress_per_vartx_message = T) #(9 / 268 variants within the final exon)
	neg_strand_ins_penult_res <- check_var_penult_exon(VEPnotAENMD_ins_neg, "neg", suppress_per_vartx_message = T) #(16 / 953 variants within the penultimate exon, 1 (10_90915591_C_CATCTTAT), expected given AENMD's alternative protein approach)
	neg_strand_SNV_penult_res <- check_var_penult_exon(VEPnotAENMD_SNV_neg, "neg", suppress_per_vartx_message = T) #(888 / 29,737 variants within the final exon)

	#binding all the neg strand results except for the one, 10_90915591_C_CATCTTAT, which we
	# expect VEP to call NMDesc but for which we do not expect AENMD to call NMDesc given
	# its alternative protein approach.
	VEP_NMDesc_notAENMD_PenultExonCheck_neg <- rbind(neg_strand_del_penult_res, neg_strand_sbs_penult_res, 
																								neg_strand_ins_penult_res %>% dplyr::filter( uploaded_variation != "10_90915591_C_CATCTTAT" ), 
																								neg_strand_SNV_penult_res)
	VEP_NMDesc_notAENMD_PenultExonCheck_neg$Purpose <- "50bp Penultimate Exon rule check."
	VEP_NMDesc_notAENMD_PenultExonCheck_neg$QC_Result <- "PTCs are not within 51bp."

	expected_neg_strand_ins_penult_res <- neg_strand_ins_penult_res %>% dplyr::filter( uploaded_variation == "10_90915591_C_CATCTTAT" )
	expected_neg_strand_ins_penult_res$Purpose <- "50bp Penultimate Exon rule check."
	expected_neg_strand_ins_penult_res$QC_Result <- "PTC is \"pushed\" out of penult rule by inserted bases."


	VEP_NMDesc_notAENMD_pos_PenultExonCheck <- rbind(pos_strand_SNV_penult_res, pos_strand_ins_penult_res)
	VEP_NMDesc_notAENMD_neg_PenultExonCheck <- rbind(expected_neg_strand_ins_penult_res, VEP_NMDesc_notAENMD_PenultExonCheck_neg)

	
	###4. Single exon
	#We do not need to variants that fall on possible strand transcripts because those have been solved.
	#First, we need to remove the variants we solved from the negative strand variant sets.
	VEPnotAENMD_ins_neg_penult_res_rmvd <- VEPnotAENMD_ins_neg %>% dplyr::filter( uploaded_variation != "10_90915591_C_CATCTTAT" )

	#Now, we will check 11 variants by hand. after that, we will check by code.
	#neg strand del
	set.seed(321)
	index <- sample(1:nrow(VEPnotAENMD_del_neg), 1)
	VEPnotAENMD_del_neg_example <- VEPnotAENMD_del_neg[index, ]

	#neg strand sbs
	set.seed(321)
	index <- sample(1:nrow(VEPnotAENMD_sbs_neg), 2)
	VEPnotAENMD_sbs_neg_example <- VEPnotAENMD_sbs_neg[index, ]

	#neg strand ins
	set.seed(321)
	index <- sample(1:nrow(VEPnotAENMD_ins_neg), 3)
	VEPnotAENMD_ins_neg_example <- VEPnotAENMD_ins_neg[index, ]

	#neg strand SNV
	set.seed(321)
	index <- sample(1:nrow(VEPnotAENMD_SNV_neg), 4)
	VEPnotAENMD_SNV_neg_example <- VEPnotAENMD_SNV_neg[index, ]

	#Merge Results:
	AnalyzedByHand_VEP_NMDesc_notAENMD_SingleExonCheck <- rbind(VEPnotAENMD_del_neg_example, VEPnotAENMD_sbs_neg_example, VEPnotAENMD_ins_neg_example, VEPnotAENMD_SNV_neg_example)
	AnalyzedByHand_VEP_NMDesc_notAENMD_SingleExonCheck$Purpose <- "Single Exon rule check."
	AnalyzedByHand_VEP_NMDesc_notAENMD_SingleExonCheck$QC_Result <- "Not in a single exon transcript."


	####Checking all to see if the variants are within single exon or single coding exon transcripts with code:
	#Get all the ENSTs:
	VEPnotAENMD_AllUniENSTs <- unique(VEPxAENMD_Clinvar_ResultsMerged_VEP_NMDesc_notAENMD$feature) #(2847)
	#Get the GRanges for all the ENSTs, which include information on the number of coding exons and number of exons for each transcript.
	#AENMD_base_transcriptset is generated in the R script: /functions/Retrieve_AENMD_BaseTxSet.R
	txs_VEPnotAENMD <- AENMD_base_transcriptset[VEPnotAENMD_AllUniENSTs]
	#check to see if any of the transcripts have a single coding exon or single exon.
	VEP_NMDesc_notAENMD_SingleExonCheck <- txs_VEPnotAENMD[(txs_VEPnotAENMD$num_exons_cds == 1 | txs_VEPnotAENMD$num_exons_txs == 1)] #(0)

	##___________________
	##___________________
	##___________________
	##Output:
	##Save the raw data used to generate data as tab separated tables:
  
	#AENMD NMDesc call, not called NMDesc by VEP
  write.table(neg_strand_last_res, file="./sup_data/raw_datatables/AENMD_NMDesc_NotVEP_negStrand_LastExonCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  write.table(pos_strand_last_res, file="./sup_data/raw_datatables/AENMD_NMDesc_NotVEP_posStrand_LastExonCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  
  write.table(neg_strand_penult50_res, file="./sup_data/raw_datatables/AENMD_NMDesc_NotVEP_negStrand_PenultExonCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  write.table(pos_strand_penult50_res, file="./sup_data/raw_datatables/AENMD_NMDesc_NotVEP_posStrand_PenultExonCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)

  #VEP NMDesc call, not called NMDesc by AENMD
  write.table(VEP_NMDesc_notAENMD_cssProxCheck, file="./sup_data/raw_datatables/VEP_NMDesc_notAENMD_cssProxCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  
  write.table(VEP_NMDesc_notAENMD_pos_LastExonCheck, file="./sup_data/raw_datatables/VEP_NMDesc_notAENMD_pos_LastExonCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  write.table(VEP_NMDesc_notAENMD_neg_LastExonCheck, file="./sup_data/raw_datatables/VEP_NMDesc_notAENMD_neg_LastExonCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  
  write.table(VEP_NMDesc_notAENMD_pos_PenultExonCheck, file="./sup_data/raw_datatables/VEP_NMDesc_notAENMD_pos_PenultExonCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  write.table(VEP_NMDesc_notAENMD_neg_PenultExonCheck, file="./sup_data/raw_datatables/VEP_NMDesc_notAENMD_neg_PenultExonCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)
  
  write.table(VEP_NMDesc_notAENMD_SingleExonCheck, file="./sup_data/raw_datatables/VEP_NMDesc_notAENMD_SingleExonCheck.txt", sep="\t", row.names=F, col.names = T,  quote = F)

  #Return the data.frames as well:
  return( list( AnalyzedByHand_Var_AENMD_NMDesc_NotVEP = AnalyzedByHand_Var_AENMD_NMDesc_NotVEP, AnalyzedByHand_VEP_NMDesc_notAENMD_SingleExonCheck = AnalyzedByHand_VEP_NMDesc_notAENMD_SingleExonCheck,
              AENMD_NMDesc_NotVEP_posStrand_LastExonCheck = pos_strand_last_res, AENMD_NMDesc_NotVEP_negStrand_LastExonCheck = neg_strand_last_res,
							AENMD_NMDesc_NotVEP_posStrand_PenultExonCheck = pos_strand_penult50_res, AENMD_NMDesc_NotVEP_negStrand_PenultExonCheck = neg_strand_penult50_res,
							VEP_NMDesc_notAENMD_cssProxCheck = VEP_NMDesc_notAENMD_cssProxCheck, 
							VEP_NMDesc_notAENMD_pos_LastExonCheck = VEP_NMDesc_notAENMD_pos_LastExonCheck, VEP_NMDesc_notAENMD_neg_LastExonCheck = VEP_NMDesc_notAENMD_neg_LastExonCheck,
							VEP_NMDesc_notAENMD_pos_PenultExonCheck = VEP_NMDesc_notAENMD_pos_PenultExonCheck, VEP_NMDesc_notAENMD_neg_PenultExonCheck = VEP_NMDesc_notAENMD_neg_PenultExonCheck,
							VEP_NMDesc_notAENMD_SingleExonCheck = VEP_NMDesc_notAENMD_SingleExonCheck
							) ) 
}