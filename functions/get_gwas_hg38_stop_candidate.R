library(GenomicRanges)
library(rsnps)

gwas_hg38 = readRDS("./functions/accessory_files/gwas_hg38.RDS")
exon_term = c("missense_variant", "synonymous_variant", "frameshift_variant", 
              "stop_gained", "stop_lost", "start_lost", 
              "inframe_deletion", "inframe_insertion",
              "coding_sequence_variant", "protein_altering_variant") # double check this


gwas_coding_hg38 = gwas_hg38[gwas_hg38$CONTEXT %in% exon_term]

# get how many unique snps
unique(gwas_coding_hg38) # 6676/6702
# how many traits
unique(gwas_coding_hg38$MAPPED_TRAIT) # 1769 traits

# context
table(gwas_coding_hg38$CONTEXT)

# take out frameshift and stopgain
stop_candidate = gwas_coding_hg38[gwas_coding_hg38$CONTEXT %in% c("stop_gained", "frameshift_variant")]
stop_candidate_id = unique(stop_candidate$SNPS)
length(stop_candidate_id) # 333 entries


# store out stop_candidate, it is important
saveRDS(stop_candidate, "./functions/accessory_files/gwas_coding_stopcandidate.RDS")



