library(aenmd)
gwas_vcf = "./functions/accessory_files/gwas_stopcandidate.vcf"
vcf_rng  = aenmd:::parse_vcf_vcfR(gwas_vcf)[[1]]

vcf_rng  = process_variants(vcf_rng, check_ref = TRUE, verbose = TRUE)
vcf_rng |> print()

vcf_annotated <- annotate_nmd(vcf_rng) 
names(vcf_annotated) = NULL
vcf_annotated = as.list(vcf_annotated)

vcf_annotated2 = do.call("c", vcf_annotated)

# vcf_annotated$which_rule = 

saveRDS(vcf_annotated2, "./functions/accessory_files/gwas_vcf_nmd_annotated.RDS")