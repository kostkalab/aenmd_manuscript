library(GenomicRanges)

FILE_CATALOG="./analysis/gwas/alternative"
# FILE_CHAIN="/data/projects/annotation/UCSC/hg38ToHg19.over.chain"

GW = read.csv(FILE_CATALOG,sep="\t",stringsAsFactors=FALSE, quote = "")

#- make a GRanges object from it / remove multi-location entries (cheap)
seqnames    = paste("chr",GW$CHR_ID,sep="")
starts      = as.numeric(GW$CHR_POS)
#- remove multi-loctation entries
id          = !is.na(starts)
seqnames    = seqnames[id]  # delete id which is NA
starts      = starts[id]
ends        = starts
SNPS        = GRanges(seqnames,IRanges(starts,ends))
mcols(SNPS) = GW[id,]

saveRDS(SNPS, "./functions/accessory_files/gwas_hg38.RDS")
