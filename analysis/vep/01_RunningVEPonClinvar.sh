############################################################
#### Next we take the clinvar vcf with the udpated      ####
#### informative ID numbers and annotate them using VEP ####   																			
############################################################

##Running VEP using script:
./analysis/vep/run_vep.sh ./analysis/clinvar/resetID_clinvar_20221211.vcf.gz

#get all the stop gain (this is what VEP annnotates) and frameshifts (Just to compare totals with AENMD)
zegrep "^#" ./analysis/clinvar/resetID_clinvar_20221211.vcf.gz.vep_out.txt.gz > ./functions/accessory_files/SGandFS_resetID_clinvar_20221211_vep_out.txt
zegrep "stop_gained|frameshift_variant" ./analysis/clinvar/resetID_clinvar_20221211.vcf.gz.vep_out.txt.gz >> ./functions/accessory_files/SGandFS_resetID_clinvar_20221211_vep_out.txt