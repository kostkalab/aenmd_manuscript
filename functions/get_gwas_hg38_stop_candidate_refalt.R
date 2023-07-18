library(GenomicRanges)
library(rsnps)

stop_candidate = readRDS("./functions/accessory_files/gwas_coding_stopcandidate.RDS")
stop_candidate_id = unique(stop_candidate$SNPS)

# get ref and alt alleles
ref_alt = ncbi_snp_query(stop_candidate_id) # 332 entries all here

saveRDS(ref_alt, "./functions/accessory_files/gwas_coding_stop_ref_alt.RDS")


