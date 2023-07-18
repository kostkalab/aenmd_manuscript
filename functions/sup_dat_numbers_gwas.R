library(stringr)
library(GenomicRanges)
library(aenmd)

# number of all variants
gwas_hg38 <- readRDS("./functions/accessory_files/gwas_hg38.RDS")
n_all = length(unique(gwas_hg38))

# number of coding variants
exon_term = c("missense_variant", "synonymous_variant", "frameshift_variant", 
              "stop_gained", "stop_lost", "start_lost", 
              "inframe_deletion", "inframe_insertion",
              "coding_sequence_variant", "protein_altering_variant") # double check this
gwas_coding_hg38 = gwas_hg38[gwas_hg38$CONTEXT %in% exon_term]
n_coding = length(unique(gwas_coding_hg38)) # 6676/6702

# number of coding variants in our transcript sets
AENMD_base_exonset <- future::value(aenmd:::._EA_exn_grl)
#unlist and sort
exn_rng_used <- sort(GenomeInfoDb::sortSeqlevels(unlist(AENMD_base_exonset)))
strand(exn_rng_used) = "*"

seqlevelsStyle(exn_rng_used) = "ucsc"
f = findOverlaps(gwas_coding_hg38, exn_rng_used)
g_used = gwas_coding_hg38[unique(queryHits(f))]
n_coding_ourtx = length(unique(g_used))

# number of ptc causing variants
# by variants
gwas_vcf_nmd_annotated = readRDS("./functions/accessory_files/gwas_vcf_nmd_annotated.RDS")
# add a column indicate if it is a nmd_all_rules
gwas_vcf_nmd_annotated$res_aenmd = cbind(gwas_vcf_nmd_annotated$res_aenmd, 
                                 nmd_escaping = apply(gwas_vcf_nmd_annotated$res_aenmd[,c(2,3,4,5,6)], 1, any))
# canonical rules
gwas_vcf_nmd_annotated$res_aenmd = cbind(gwas_vcf_nmd_annotated$res_aenmd, 
                                         nmd_escaping_canonical = apply(gwas_vcf_nmd_annotated$res_aenmd[,c(2,3)], 1, any))
# non-canonical rules
gwas_vcf_nmd_annotated$res_aenmd = cbind(gwas_vcf_nmd_annotated$res_aenmd, 
                                         nmd_escaping_noncanonical = apply(gwas_vcf_nmd_annotated$res_aenmd[,c(4,5,6)], 1, any))



unique_var = unique(str_extract(gwas_vcf_nmd_annotated$id, "rs[0-9]*"))
gwas_vcf_nmd_annotated$rsid = str_extract(gwas_vcf_nmd_annotated$id, "rs[0-9]*")
get_T_F = function(x){
  if(all(x)){
    return(T)
  }else if(!any(x)){
    return(F)
  }else{
    return("dependent")
  }
}

by_variant = lapply(unique_var, function(x){
  g = gwas_vcf_nmd_annotated[gwas_vcf_nmd_annotated$rsid == x]
  g_ptc = g[!g$res_aenmd$is_ptc == "FALSE",]
  if(length(g_ptc) == 0){
    df = data.frame(rsid = x, 
                    # is_ptc = any(g$res_aenmd$is_ptc),
                    is_ptc = get_T_F(g$res_aenmd$is_ptc),
                    is_nmd_escape = NA,
                    is_nmd_escape_canonical = NA,
                    is_nmd_escape_noncanonical = NA,
                    rule407 = NA,
                    rule_start = NA,
                    rule_single = NA)
  }else{
    df = data.frame(rsid = x, 
                    # is_ptc = any(g$res_aenmd$is_ptc),
                    is_ptc = get_T_F(g$res_aenmd$is_ptc),
                    is_nmd_escape = get_T_F(g_ptc$res_aenmd$nmd_escaping),
                    is_nmd_escape_canonical = get_T_F(g_ptc$res_aenmd$nmd_escaping_canonical),
                    is_nmd_escape_noncanonical = get_T_F(g_ptc$res_aenmd$nmd_escaping_noncanonical),
                    rule407 = get_T_F(g_ptc$res_aenmd$is_407plus),
                    rule_start = get_T_F(g_ptc$res_aenmd$is_cssProximal),
                    rule_single = get_T_F(g_ptc$res_aenmd$is_single))
  }
  
  return(df)
})

by_variant = do.call("rbind", by_variant)

n_ptc                = table(by_variant$is_ptc)["TRUE"]
n_ptc_dependent      = table(by_variant$is_ptc)["dependent"]
by_variant_ptc       = by_variant[!by_variant$is_ptc == "FALSE",]
n_escaping           = table(by_variant_ptc$is_nmd_escape)["TRUE"]
n_triggering         = table(by_variant_ptc$is_nmd_escape)["FALSE"]
n_dependent          = table(by_variant_ptc$is_nmd_escape)["dependent"]
n_escaping_cano      = table(by_variant_ptc$is_nmd_escape_canonical)["TRUE"]
n_triggering_cano    = table(by_variant_ptc$is_nmd_escape_canonical)["FALSE"]
n_dependent_cano     = table(by_variant_ptc$is_nmd_escape_canonical)["dependent"]
n_escaping_noncano   = table(by_variant_ptc$is_nmd_escape_noncanonical)["TRUE"]
n_triggering_noncano = table(by_variant_ptc$is_nmd_escape_noncanonical)["FALSE"] # add
n_dependent_noncano  = table(by_variant_ptc$is_nmd_escape_noncanonical)["dependent"] # add
n_407rule            = table(by_variant_ptc$rule407)["TRUE"]
n_start              = table(by_variant_ptc$rule_start)["TRUE"]
n_singleexon         = table(by_variant_ptc$rule_single)["TRUE"]

gwas_ptc = gwas_vcf_nmd_annotated[gwas_vcf_nmd_annotated$res_aenmd$nmd_escaping]$res_aenmd
n_annotation = apply(gwas_ptc, 1, function(x){
  sum(x[2:6])
}) # for each nmd escaping variant, how many rules are there? 
n_annotation = round(mean(n_annotation), 2)

data_suppl_gwas1 = data.frame(functional_class = c("all", 
                                                   # "coding", 
                                                   "coding_our_tx_set",
                                                   "ptc_causing", "transcript_dependent_ptc",
                                                   "nmd_escaping","nmd_triggering", "transcript_alt_allele_dependent",
                                                   "nmd_escaping_canonical", "nmd_triggering_canonical", "transcript_andor_alt_allele_dependent_canonical",
                                                   "nmd_escaping_noncanonical", "nmd_triggering_noncanonical", "transcript_andor_alt_allele_dependent_noncanonical",
                                                   "nmd_escaping_407rule",
                                                   "nmd_escaping_start_proximal_rule",
                                                   "nmd_escaping_single_exon_rule",
                                                   "nmd_escaping_annotation_per_nmd_escape_variant"),
                              number = c(n_all, 
                                         # n_coding, 
                                         n_coding_ourtx,
                                         n_ptc,n_ptc_dependent,
                                         n_escaping, n_triggering, n_dependent, 
                                         n_escaping_cano, n_triggering_cano, n_dependent_cano,
                                         n_escaping_noncano, n_triggering_noncano, n_dependent_noncano,
                                         n_407rule,
                                         n_start,
                                         n_singleexon,
                                         n_annotation))
data_suppl_gwas1$percent = NA
data_suppl_gwas1$percent[2] = data_suppl_gwas1$number[2]/data_suppl_gwas1$number[1]
data_suppl_gwas1$percent[3:4] = data_suppl_gwas1$number[3:4]/data_suppl_gwas1$number[2]
data_suppl_gwas1$percent[5:16] = data_suppl_gwas1$number[5:16]/(data_suppl_gwas1$number[3]+data_suppl_gwas1$number[4])
data_suppl_gwas1$percent = round(data_suppl_gwas1$percent*100,2)

vcf_ptc = gwas_vcf_nmd_annotated[gwas_vcf_nmd_annotated$res_aenmd$is_ptc]
n_pair_ptc      = length(vcf_ptc)
n_pair_escaping = table(vcf_ptc$res_aenmd$nmd_escaping)["TRUE"]
n_pair_trigger  = table(vcf_ptc$res_aenmd$nmd_escaping)["FALSE"]
n_tx_per_var = round(length(gwas_vcf_nmd_annotated)/length(unique(gwas_vcf_nmd_annotated$id)), digits = 2)

n_pair_escaping_cano        = table(vcf_ptc$res_aenmd$nmd_escaping_canonical)["TRUE"]
n_pair_triggering_cano      = table(vcf_ptc$res_aenmd$nmd_escaping_canonical)["FALSE"]
n_pair_escaping_noncano     = table(vcf_ptc$res_aenmd$nmd_escaping_noncanonical)["TRUE"]
n_pair_triggering_noncano   = table(vcf_ptc$res_aenmd$nmd_escaping_noncanonical)["FALSE"]
n_pair_escaping_407         = table(vcf_ptc$res_aenmd$is_407plus)["TRUE"]
n_pair_escaping_start       = table(vcf_ptc$res_aenmd$is_cssProximal)["TRUE"]
n_pair_escaping_single_exon = table(vcf_ptc$res_aenmd$is_single)["TRUE"]

# mean number of rules
get_mean_n_rules = function(x){ # x is nmd predicting table
  escaping = x[x$nmd_escaping,]
  rules = escaping[,c(2,3,4,6)]
  a = apply(rules, 1, sum)
  return(mean(a))
}
n_number_rules = round(get_mean_n_rules(gwas_vcf_nmd_annotated$res_aenmd),2)


# get canonical transcript
source("./functions/Retrieve_CanonTxSet.R")
canonical_tx <- Get_Canonical_Transcripts()
canonical_tx <- readRDS("~/project/aenmd_manuscript/data/canonical_tx.RDS")
canonical_tx = canonical_tx$tx_id

gwas_vcf_nmd_annotated_cntx = gwas_vcf_nmd_annotated[gwas_vcf_nmd_annotated$tx_id %in% canonical_tx$tx_id]
vcf_ptc_cntx = gwas_vcf_nmd_annotated_cntx[gwas_vcf_nmd_annotated_cntx$res_aenmd$is_ptc]

n_ptc_causing_cano_tx    = length(vcf_ptc_cntx)
n_nmd_escaping_cano_tx   = table (vcf_ptc_cntx$res_aenmd$nmd_escaping)["TRUE"]
n_nmd_triggering_cano_tx = table (vcf_ptc_cntx$res_aenmd$nmd_escaping)["FALSE"]
n_tx_per_var_cano_tx     = round (length(gwas_vcf_nmd_annotated_cntx)/length(unique(gwas_vcf_nmd_annotated_cntx$id)), digits = 2)



data_suppl_gwas2 = data.frame(functional_class = c("ptc_causing_pairs", "nmd_escaping_pairs",
                                                   "nmd_triggering_pairs", "transcript_per_variant",
                                                   "nmd_escaping_canon_rules_pairs",
                                                   "nmd_triggering_canon_rules_pairs",
                                                   "nmd_escaping_noncanon_rules_pairs",
                                                   "nmd_triggering_noncanonical_rules_pairs",
                                                   "nmd_escaping_407_rule_pairs",
                                                   "nmd_escaping_start_rule_pairs",
                                                   "nmd_escaping_single_exon_rule_pairs",
                                                   "nmd_escaping_annotation_per_nmd_escape_variant_tx_pair",
                                                   "ptc_causing_canon_tx_pairs",
                                                   "nmd_escaping_canon_tx_pairs",
                                                   "nmd_triggering_canon_tx_pairs",
                                                   "transcript_per_variant_canon_tx"),
                              number = c(n_pair_ptc, n_pair_escaping, n_pair_trigger, n_tx_per_var,
                                         n_pair_escaping_cano,
                                         n_pair_triggering_cano,
                                         n_pair_escaping_noncano,
                                         n_pair_triggering_noncano,
                                         n_pair_escaping_407,
                                         n_pair_escaping_start,
                                         n_pair_escaping_single_exon,
                                         n_number_rules,
                                         n_ptc_causing_cano_tx,
                                         n_nmd_escaping_cano_tx,
                                         n_nmd_triggering_cano_tx,
                                         n_tx_per_var_cano_tx))

data_suppl_gwas2$percent = NA
data_suppl_gwas2$percent[c(2:3,5:11)] = data_suppl_gwas2$number[c(2:3,5:11)]/data_suppl_gwas2$number[1]
data_suppl_gwas2$percent[14:15] = data_suppl_gwas2$number[14:15]/data_suppl_gwas2$number[13]
data_suppl_gwas2$percent = round(data_suppl_gwas2$percent*100,2)

readr::write_csv(data_suppl_gwas1, file="./sup_data/DT_2_sup_data_gwas_by_variant.csv") 
readr::write_csv(data_suppl_gwas2, file="./sup_data/DT_1_sup_data_gwas_by_variant_transcript_pair.csv") 
