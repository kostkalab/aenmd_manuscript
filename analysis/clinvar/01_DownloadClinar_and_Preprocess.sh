##########################################################################################################
#### Clinvar downloaded from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ dated to 2022 12 11 ####
##########################################################################################################
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/clinvar_20221211.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/clinvar_20221211.vcf.gz.tbi

##1. Give each variant a unique identifier that helps track variant information:  chr, pos, ref, alt
##2. Then Isolate the variants with the main VCF columns (chr, pos, id, ref, alt) along with the clinvar annotation (Clinsig)

#1. The file being created here is what will be fed for VEP annotation later
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' clinvar_20221211.vcf.gz > resetID_clinvar_20221211.vcf
bgzip resetID_clinvar_20221211.vcf
tabix resetID_clinvar_20221211.vcf.gz

#2. The file beling created here is what will be fed into R for AENMD annotation
echo -e CHROM\\tPOS\\tID\\tREF\\tALT\\tCLNSIG > resetID_clinvar_20221211_Clnsig.txt
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\n' resetID_clinvar_20221211.vcf.gz >> resetID_clinvar_20221211_Clnsig.txt
