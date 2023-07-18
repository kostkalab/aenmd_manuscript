####gnomAD v2 liftover downloaded from https://gnomad.broadinstitute.org/downloads
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz

##Then filter for variants that passed gnomAD's QC filters and give each variant a unique identifier that 
##helps track variant information:  chr, pos, ref, alt

bcftools view gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz -f PASS | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -o gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
