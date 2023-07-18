##Library aenmd and the data package into R
library(aenmd.data.ensdb.v105)
library(aenmd)

##Import VCF of variants that passed QC filter and with new unique identifiers
vcf_rng <-aenmd:::parse_vcf_vcfR("./analysis/gnomAD/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz")

#Make a function that pads the start position and remove duplicates, because this is gnomAD liftover from GRCh37/hg19 - 
#      such that some SNPs map to the same location- saving the file as you go:
step2_pad_and_removeDups <- function(vcf_rng, outputname) {
	starts <- GenomicRanges::start(vcf_rng) |> stringr::str_pad(9L, pad="0")
	kd <- vcf_rng$id[vcf_rng$id %in% vcf_rng$id[duplicated(vcf_rng$id)]] |> unique()
	ru <- vcf_rng[!(vcf_rng$id %in% kd)]
	rd <- vcf_rng[(vcf_rng$id %in% kd)]
	saveRDS(ru, file= paste0(outputname, '.rds')) 
	return(ru)
}
##Run the function:
ru_vcf_rng <- step2_pad_and_removeDups(vcf_rng$vcf_rng, "vcf_rng")

#_____________________________________________________________________________
###- 1. First we are going to retrieve how many total variants there are, then  
### how many total variants overlap with the coding regions of AENMD's default 
### transcript set.
#_____________________________________________________________________________
#Get total variants in gnomAD:
Variants_in_the_original_dataset <- length(ru_vcf_rng)

#Retrieve the base AENMD transcript set information: 
AENMD_base_exonset <- future::value(aenmd:::._EA_exn_grl)

#unlist and sort
AENMD_base_exonset_unlist_srtd <- sort(GenomeInfoDb::sortSeqlevels(unlist(AENMD_base_exonset)))

#find overlap between clinvar variants and the exons of the transcript set:
ov <- GenomicRanges::findOverlaps(ru_vcf_rng, AENMD_base_exonset_unlist_srtd)

#find the number of unique hits, such that a variant can only be counted once:
Variants_coding_our_tx_set <- length(unique(S4Vectors::queryHits(ov)))

#make a data.frame:
Total_var_stats <- data.frame(functional_class = c("Variants_in_the_original_dataset" , "Variants_coding_our_tx_set"),
                              value = c(Variants_in_the_original_dataset, Variants_coding_our_tx_set) , percent = c(NA, Variants_coding_our_tx_set/Variants_in_the_original_dataset*100)
                              )
#_____________________________________________________________________________
#_____________________________________________________________________________

#_____________________________________________________________________________
###- 2. Next, we are going to run AENMD on the gnomad data
#_____________________________________________________________________________

#Process the variants using AENMD:
vcf_rng_proc <- process_variants(ru_vcf_rng)

Run_Annotate_and_Save <- function(vcf_rng_proc, outputname) {
	vcf_rng_an  <- annotate_nmd(vcf_rng_proc)
	saveRDS(vcf_rng_an, file= paste0('./functions/accessory_files/', outputname, '.rds'))
	return(vcf_rng_an)
}

gnomad_an <- Run_Annotate_and_Save(vcf_rng_proc, "gnomad_an")

#_____________________________________________________________________________
#_____________________________________________________________________________

#_____________________________________________________________________________
###- 3. After that we get, create, and output the basic data frames for gnomAD 
#_____________________________________________________________________________
###- Get the basic data frames for clinvar - 
# Source and run /functions/accessory_files/Get_Canonical_Transcripts.R
source('./functions/Retrieve_CanonTxSet.R')
canonical_tx <- Get_Canonical_Transcripts()
# source and run /functions/Mutual_DataAnalysis.R
source('./functions/General_DataAnalysis.R')
# function Basic_Data_Table() can be found in /functions/Mutual_DataAnalysis.R
gnomad_basic_DFs <- Basic_Data_Table(gnomad_an, canonical_tx, Variants_coding_our_tx_set)

#Join the output Uni_Var_DF with the data frame created in the first step
DT_4_gnomad_BasicStats_UniqueVar <- rbind(Total_var_stats, gnomad_basic_DFs$Uni_Var_DF)

#Write the results:
readr::write_csv(gnomad_basic_DFs$Var_Tx_DF, file="./sup_data/DT_3_gnomAD_BasicStats_VarTxPairs.csv")
readr::write_csv(DT_4_gnomad_BasicStats_UniqueVar, file="./sup_data/DT_4_gnomad_BasicStats_UniqueVar.csv")