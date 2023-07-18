############____________________________________#########
############____________________________________#########
############____________________________________#########
#Code for annotating all clinvar variants and analyzing the data
#
#
##Library aenmd into R
library(aenmd.data.ensdb.v105)
library(aenmd)
library(GenomicRanges)
library(tidyr)

##Read in the downloaded and BCFtools processed Clinvar data
clinvarClinsig <- read.delim("./analysis/clinvar/resetID_clinvar_20221211_Clnsig.txt")

#_____________________________________________________________________________
###- 1. First we are going to retrieve how many total variants there are, then  
### how many total variants overlap with the coding regions of AENMD's default 
### transcript set.
#_____________________________________________________________________________
#Get total variants in Clinvar:
Variants_in_the_original_dataset <- nrow(clinvarClinsig)

#Make a GRange object of all the clinvar variants without filtering.
clinvar_gr_total  <- GenomicRanges::GRanges(clinvarClinsig$CHR,
                                            IRanges::IRanges(clinvarClinsig$POS,
                                                             clinvarClinsig$POS))
clinvar_gr_total$ref     <- clinvarClinsig$REF |> Biostrings::DNAStringSet()
clinvar_gr_total$alt     <- clinvarClinsig$ALT |> Biostrings::DNAStringSet()
clinvar_gr_total$id      <- clinvarClinsig$ID
clinvar_gr_total$clinsig      <- clinvarClinsig$CLNSIG

#Retrieve the base AENMD transcript set information: 
AENMD_base_exonset <- future::value(aenmd:::._EA_exn_grl)

#unlist and sort
AENMD_base_exonset_unlist_srtd <- sort(GenomeInfoDb::sortSeqlevels(unlist(AENMD_base_exonset)))

#find overlap between clinvar variants and the exons of the transcript set:
ov <- GenomicRanges::findOverlaps(clinvar_gr_total, AENMD_base_exonset_unlist_srtd)

#find the number of unique hits, such that a variant can only be counted once:
Variants_coding_our_tx_set <- length(unique(S4Vectors::queryHits(ov)))

#make a data.frame:
Total_var_stats <- data.frame(functional_class = c("Variants_in_the_original_dataset" , "Variants_coding_our_tx_set"),
                              value = c(Variants_in_the_original_dataset, Variants_coding_our_tx_set), percent = c(NA, Variants_coding_our_tx_set/Variants_in_the_original_dataset*100)
                              )
#_____________________________________________________________________________
#_____________________________________________________________________________

#_____________________________________________________________________________
###- 2. Next, we are going to run AENMD on the Clinvar data
#_____________________________________________________________________________
#get NA indexes and indexes where the alt is not a standard base (A, G, T, C)
NA_ind <- which(clinvarClinsig$ALT != "." & !stringr::str_detect(clinvarClinsig$ALT, "R|Y|S|W|K|M|B|D|H|V|N"))

#now put together the GRanges without NA to put into AENMD
#we do it this way to include the clinsig annotation
clinvar_gr  <- GenomicRanges::GRanges(clinvarClinsig$CHR[NA_ind],
                                       IRanges::IRanges(clinvarClinsig$POS[NA_ind],
                                                        clinvarClinsig$POS[NA_ind]))
clinvar_gr$ref     <- clinvarClinsig$REF[NA_ind] |> Biostrings::DNAStringSet()
clinvar_gr$alt     <- clinvarClinsig$ALT[NA_ind] |> Biostrings::DNAStringSet()
clinvar_gr$id      <- clinvarClinsig$ID[NA_ind]
clinvar_gr$clinsig      <- clinvarClinsig$CLNSIG[NA_ind]

#Process variants
clinvar_proc <- process_variants(clinvar_gr)

##- Annotate variants
#To compare to VEP annotation
res_clinvar_VEP  <- aenmd::annotate_nmd(clinvar_proc, css_prox_dist = 100)
#Our default setting annotation
res_clinvar_default <- aenmd::annotate_nmd(clinvar_proc)

#save results
saveRDS(res_clinvar_VEP, file= "./functions/accessory_files/aenmdVEPsettings_20221211clinvar105_results.rds")
saveRDS(res_clinvar_default, file= "./functions/accessory_files/aenmdDefaultSettings_20221211clinvar105_results.rds")

#_____________________________________________________________________________
#_____________________________________________________________________________

#_____________________________________________________________________________
###- 3. After that we get, create, and output the basic data frames for clinvar. 
#_____________________________________________________________________________

# Source and run /functions/accessory_files/Get_Canonical_Transcripts.R
source('./functions/Retrieve_CanonTxSet.R')
canonical_tx <- Get_Canonical_Transcripts()
# source and run /functions/Mutual_DataAnalysis.R
source('./functions/General_DataAnalysis.R')
#Run Basic_Data_Table()
Clinvar_basic_DFs <- Basic_Data_Table(res_clinvar_default, canonical_tx, Variants_coding_our_tx_set)

#Join the output Uni_Var_DF with the data frame created in the first step
DT_6_Clinvar_BasicStats_UniqueVar <- rbind(Total_var_stats, Clinvar_basic_DFs$Uni_Var_DF)

#Write the results:
readr::write_csv(Clinvar_basic_DFs$Var_Tx_DF, file="./sup_data/DT_5_Clinvar_BasicStats_VarTxPairs.csv")
readr::write_csv(DT_6_Clinvar_BasicStats_UniqueVar, file="./sup_data/DT_6_Clinvar_BasicStats_UniqueVar.csv")

#Get the additional data frames for clinvar - 
Clinvar_extended_DFs <- Basic_Data_Table(res_clinvar_default, canonical_tx, Variants_coding_our_tx_set, clinvar_extended_analysis = T)

#Write the results:
readr::write_csv(Clinvar_extended_DFs$ext_clinvar_Uni_Var_DF, file="./sup_data/DT_7_Clinvar_ClinsigFiltered_BasicStats_UniqueVar.csv")
readr::write_csv(Clinvar_extended_DFs$Clinsig_benign_DF, file="./sup_data/DT_8_Clinvar_Benign_BasicStats_VarTxPairs.csv")
readr::write_csv(Clinvar_extended_DFs$Clinsig_confl_DF, file="./sup_data/DT_9_Clinvar_Conflicting_VarTxPairs.csv")
readr::write_csv(Clinvar_extended_DFs$Clinsig_patho_DF, file="./sup_data/DT_10_Clinvar_Pathogenic_VarTxPairs.csv")

###- Creating plots from supplementary figure 1:
#Create GR for the entire, unfiltered clinvar data:
clinvar_full_gr  <- GenomicRanges::GRanges(clinvarClinsig$CHR,
                                      IRanges::IRanges(clinvarClinsig$POS,
                                                       clinvarClinsig$POS))
clinvar_full_gr$ref     <- clinvarClinsig$REF |> Biostrings::DNAStringSet()
clinvar_full_gr$alt     <- clinvarClinsig$ALT |> Biostrings::DNAStringSet()
clinvar_full_gr$id      <- clinvarClinsig$ID
clinvar_full_gr$clinsig      <- clinvarClinsig$CLNSIG

#execute plot function, note that the compute time is ~30 MINUTES,
#as such, the default is not to run the command (execute_analysis_and_plot = F)
#set execute_analysis_and_plot = T to have the command run.
#plots are auto saved into ./plot
Sup_fig_1_plots <- plot_clinvar_stack(clinvar_full_gr, res_clinvar_default, execute_analysis_and_plot = F)
