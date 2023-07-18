# README
------------------------------------------

AENMD is a computational tool 
(available for [R package](https://github.com/kostkalab/aenmd) and [command line interface](https://github.com/kostkalab/aenmd_cli)) 
that annotate predicted escape from nonsense-mediated decay (NMD).


This repository stores the code and result for our [preprint](https://www.biorxiv.org/content/10.1101/2023.03.17.533185v1).



## 0. Environment ##

Install AENMD R package, see [github](https://github.com/kostkalab/aenmd) for instructions. 

Install AENMD data R packages, see [github](https://github.com/kostkalab/aenmd_data) for instructions. 

Set up an environment providing `R` (at least version 4.1.2 (2021-11-01) ) and ensembl VEP version 108.2, but using ensmbl data version 105.
Below we show how this can be set up using the [mamba package manager](https://github.com/mamba-org/mamba). 

```
#- setup the environment, install R
$ mamba create -n aenmd_vep_v108.2 python=3.10
$ mamba activate aenmd_vep_v108.2
$ mamba install r-base

#- install bioperl (needed by VEP)
$ mamba install -c bioconda perl-bioperl
#- vep has specific requirements about Zlib that we have to deal with
$ wget https://cpan.metacpan.org/authors/id/P/PM/PMQS/Compress-Raw-Zlib-2.202.tar.gz
$ tar xvfz Compress-Raw-Zlib-2.202.tar.gz
$ cd Compress-Raw-Zlib-2.202
$ perl Makefile.pl
#- this had the wrong sysroot LDDLFLAGS - needed to edit the Makefkle, might be specific for the mac
#  thmis will give the correct sysroot: xcrun --sdk macosx --show-sdk-path
#  for me that was: /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
$ make
$ make install

#- install VEP version 108.2 and download version 105 annotations
$ mamba install -c bioconda ensembl-vep=108.2
#- install VEP
$ vep_install --NO_HTSLIB
$ vep_install  -a cf -s homo_sapiens -y GRCh38 --CACHE_version 105 -g NMD

To ensure that v105 and everything was installed, there should be a file:
~/.vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

Next, you can clone this repository, and the specific version of `aenmd` with the right data package and then install them.

```
#- clone repositories
$ git clone https://github.com/kostkalab/aenmd_manuscript.git
$ git clone https://github.com/kostkalab/aenmd.data.ensdb.v105.git 
$ git clone https://github.com/kostkalab/aenmd.git 
#- checkout the versions used
$ cd aenmd
$ git checkout a640839 #- version 0.2.15
$ cd ../aenmd.data.ensdb.v105
$ git checkout 32f7f22 #- version 0.2.2

#- now build and install the aenmd and aenmd.data.ensdb.v105 packages.
```



## 1. GWAS catalog analysis 

Here we generate results for running `aenmd` on the GWAS catalog. They generates the number for Supplemental Table 3 (ST3). 


```
$ R -e "source('./functions/get_gwas_hg38.R')"
$ R -e "source('./functions/get_gwas_hg38_stop_candidate.R')"
$ R -e "source('./functions/get_gwas_hg38_stop_candidate_refalt.R')"
$ R -e "source('./functions/get_gwas_hg38_stop_candidate_vcf.R')"
$ R -e "source('./functions/annotate_gwas_stop_candidate.R')"
$ R -e "source('./functions/sup_dat_numbers_gwas.R')"
```
Results will be in

* ./sup_data/sup_data_gwas_by_variant.csv
* ./sup_data/sup_data_gwas_by_variant_transcript_pair.csv




## 2. gnomAD analysis 

----BCFTOOLS IS REQUIRED FOR THIS SECTION----

Here we generate results for running `aenmd` on gnomAD. They generates the number for Supplemental Table 1 (ST1). 


1. Download gnomAD data to ./ :

This data can be downloaded a number of ways, see ./analysis/gnomAD/01_DownloadGnomAD_and_PreProcess.sh

```
$ ls ./analysis/gnomAD/01_DownloadGnomAD_and_PreProcess.sh

```

2. Process the gnomAD file : ./gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz - BCFTOOLS IS REQUIRED FOR THIS STEP:

```
$ ./analysis/gnomAD/01_DownloadGnomAD_and_PreProcess.sh

```

3. Run the annotation and analysis:

```
$ R -e "source('./analysis/gnomad/02_Annotating_gnomAD_withAENMD_andAnalysis.R')"
```
Results will be in

* ./sup_data/DT_3_gnomAD_BasicStats_VarTxPairs.cvs
* ./sup_data/DT_4_gnomad_BasicStats_UniqueVar.csv




## 3. Clinvar analysis 

Here we generate results for running `aenmd` on clinvar. They generates the number for Supplemental Table 2 (ST2)

1. Download clinvar data to ./ and process:

```
$ ./analysis/clinvar/01_DownloadClinar_and_Preprocess.sh

```

2. Run the annotation and analysis:

```
$ R -e "source('./analysis/clinvar/02_Annotating_clinvar_withAENMD_andAnalysis.R')"
```
Results will be in

* ./sup_data/DT_5_Clinvar_BasicStats_VarTxPairs.csv
* ./sup_data/DT_6_Clinvar_BasicStats_UniqueVar.csv
* ./sup_data/DT_7_Clinvar_ClinsigFiltered_BasicStats_UniqueVar.csv
* ./sup_data/DT_8_Clinvar_Benign_BasicStats_VarTxPairs.csv
* ./sup_data/DT_9_Clinvar_Conflicting_VarTxPairs.csv
* ./sup_data/DT_10_Clinvar_Pathogenic_VarTxPairs.csv




## 4. Run VEP - 

----A CONDA ENV THAT CONTAINS ENSEMBL VEP V105 IS REQUIRED AND MUST BE ACTIVATED FOR THIS SECTION----
To run VEP, a conda environment with ensembl VEP v105 must be activated. 
The below script will run vep and process the file to create the input for downstream analysis.

```
$ ./analysis/vep/01_RunningVEPonClinvar.sh

```



## 5. VEP comparison to AENMD for NMDesc annotation of Clinvar  

----This section REQUIRES that section 3 and 4 are BOTH compeleted before executing this section----
Here we generate results for the comparison between VEP and AENMD NMD escape annotation of Clinvar variants, they are supplemental tables 11, 12, 13 and 14.


1. After running AENMD on Clinvar in 3 and VEP on Clinvar in 4, run the comparison:

```
$ R -e "source('./analysis/VEPxAENMD_Comparison/VEPxAENMD_Comparison.R')"

```

Results will be in

* ./sup_data/DT_11_AENMDxVEP_AllVarTxPairs.cvs
* ./sup_data/DT_12_AENMDxVEPmerged_StopGainVarTxPairs.cvs
* ./sup_data/DT_13_AnalyzedByHand_Var_AENMD_NMDesc_NotVEP.cvs
* ./sup_data/DT_14_AnalyzedByHand_VEP_NMDesc_notAENMD_SingleExonCheck.cvs