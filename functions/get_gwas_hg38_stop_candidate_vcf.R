library(GenomicRanges)


ref_alt = readRDS("./functions/accessory_files/gwas_coding_stop_ref_alt.RDS") # there are 333 unique rsID, with 332 unique snps, since one is merged
stop_candidate = readRDS("./functions/accessory_files/gwas_coding_stopcandidate.RDS")

# delete duplicates
# which one is merged? 
a = ref_alt[c(which(duplicated(ref_alt$rsid, fromLast = T)), which(duplicated(ref_alt$rsid))),]$query
stop_candidate[stop_candidate$SNPS %in% a]

# get gwas risk allele from stop_candidate, for 333 variants
gwas_risk_allele = sapply(ref_alt$query, function(x){
  s1 = stop_candidate[stop_candidate$SNPS == x]
  risk_allele = unique(s1$STRONGEST.SNP.RISK.ALLELE) # different entries can have different risk allele.
  # should note this in future study, those different risk allele can probably lead to different traits
  # for now only concentrate on snps, not traits
  risk_allele = gsub(paste0(x,"-"),"", risk_allele)
  risk_allele = gsub("\\?", NA, risk_allele) # remove ? 
  risk_allele = risk_allele[!is.na(risk_allele)]
  risk_allele = paste(risk_allele, collapse = "|")
  
  return(risk_allele)
})

# get a dataset from ref_alt containing useful information, for 333 variants
df = data.frame(rsid = ref_alt$rsid,
                chr = ref_alt$chromosome,
                pos = ref_alt$bp,
                class = ref_alt$class,
                ancestral_allele = ref_alt$ancestral_allele, 
                # ref_seq = ref_alt$ref_seq, 
                variation_allele = ref_alt$variation_allele, 
                # alleles = ref_alt$alleles,
                gwas_risk_allele, 
                minor = ref_alt$minor, 
                maf = round(ref_alt$maf,3))

# remove ? in the gwas risk alleles

## 1 look at snv
snp   = df[df$class == "snv",]

# get vcf file
vcf_snp = lapply(snp$rsid, function(x){
  x1 = snp[snp$rsid == x,]
  ref = x1$ancestral_allele
  alt = stringr::str_split(x1$variation_allele,",")[[1]]
  class = x1$class
  gwas_risk = stringr::str_split(x1$gwas_risk_allele,"\\|")[[1]]
  if((gwas_risk[1] !="") & (!all(gwas_risk %in% c(ref, alt)))){
    print("not in it")
    alt2 = gwas_risk[!gwas_risk %in% ref]
    alt = c(alt,alt2)
  }
  
  
  s1_vcf = data.frame(chr = paste0("chr",x1$chr), 
                      pos = x1$pos, 
                      id = paste(x1$rsid, letters[1:length(alt)], class, sep = "-"),
                      ref = x1$ancestral_allele,
                      alt = alt)
  
  
  return(s1_vcf)
})
vcf_snp = do.call("rbind", vcf_snp)


## step 2 look at ins, del and delins
delins = df[df$class == "delins",]
ins = df[df$class == "ins",]
del = df[df$class == "del",]

library(BSgenome.Hsapiens.NCBI.GRCh38)
# seq1 = getSeq(Hsapiens, "17", start =76387026, end = 76387027)

# for ins
ins_vcf = lapply(ins$rsid, function(x){
  x1 = ins[ins$rsid == x,]
  class = x1$class
  pos = x1$pos-1
  seq1 = as.character(getSeq(Hsapiens, x1$chr, start =pos, end = pos))
  alt = stringr::str_split(x1$variation_allele,",")[[1]]
  alt = paste0(seq1, alt)
  vcf = data.frame(chr = paste0("chr",x1$chr), 
                   pos = pos,
                   id = paste(x1$rsid, letters[1:length(alt)],class, sep = "-"),
                   ref = seq1,
                   alt = alt)
})
ins_vcf = do.call("rbind", ins_vcf)

del_vcf = lapply(del$rsid, function(x){
  x1 = del[del$rsid == x,]
  class = x1$class
  pos = x1$pos-1
  len = nchar(x1$ancestral_allele)
  seq1 = as.character(getSeq(Hsapiens, x1$chr, start =pos, end = pos+len))
  alt = substr(seq1, 1, 1)
  vcf = data.frame(chr = paste0("chr",x1$chr), 
                   pos = pos,
                   id = paste(x1$rsid, letters[1:length(alt)],class, sep = "-"),
                   ref = seq1,
                   alt = alt)
  return(vcf)
})
del_vcf = do.call("rbind", del_vcf)



# remove indel duplicates, should be careful

delins = delins[!duplicated(delins$rsid),]
delins_vcf = lapply(delins$rsid, function(x){
  x1 = delins[delins$rsid == x,]
  class = x1$class
  ref = x1$ancestral_allele
  alt = stringr::str_split(x1$variation_allele,",")[[1]]
  pos = x1$pos
  s1_vcf = data.frame(chr = paste0("chr",x1$chr), 
                      # pos = x1$pos,
                      pos = pos,
                      id = paste(x1$rsid, letters[1:length(alt)],class, sep = "-"),
                      ref = ref,
                      alt = alt)
  
  return(s1_vcf)
})

delins_vcf = do.call("rbind", delins_vcf)

vcf = rbind(vcf_snp, ins_vcf, del_vcf, delins_vcf)

# add last three columns
vcf$QUAL = "."
vcf$FILTER = "PASS"
vcf$INFO = "."

cat("##fileformat=VCFv4.2\n", 
    "##Genome build GRch38v12/hg38\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    sep = "", 
    file="./functions/accessory_files/gwas_stopcandidate.vcf")
write.table(vcf, "./functions/accessory_files/gwas_stopcandidate.vcf", 
            append = T, sep = "\t", dec = ".",
            row.names = F, col.names = F, quote = F)


saveRDS(vcf, "./functions/accessory_files/gwas_coding_stopcandidate_vcf.RDS")


