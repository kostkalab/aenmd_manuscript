#########_______________________________________________#########
#########_______________________________________________#########
#########__Retrieve all ensembl canonical transcripts___#########
#########__to select variants from AENMD analysis_______#########
#########__that overlap canonical transcripts.__________#########
#########_______________________________________________#########
############____________________________________________#########
#' @return canonical_tx - List object containing ENSTs for ensembl canonical 
#'       transcripts. This object is used in Mutual_DataAnalysis.R for comparing 
#'       AENMD results when using high quality vs canonical high quality 
#'       transcript set.
#' @output Dataframe of canonical_tx to 
#'       ./functions/accessory_files/canonical_tx_set.txt
#'       
Get_Canonical_Transcripts <- function() {
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38 )
  library(ensembldb)
  library(AnnotationHub)
  
  Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens #- BSgenome.Hsapiens.UCSC.hg38_1.4.4
  genome(Hsapiens) <- 'GRCh38' 
  
  AH_VERSION <- 'AH98047'
  
  #- Get the data
  ah   <- AnnotationHub::AnnotationHub()
  edb  <- ah[[AH_VERSION]]
  
  #- COMPILE TRANSCRIPT SET
  #------------------------
  
  #- all protein-coding txs on standard chromosomes
  flt_chr <- AnnotationFilter::SeqNameFilter(as.character(c(1:22,"X","Y","MT")))
  flt_tbt <- AnnotationFilter::TxBiotypeFilter("protein_coding")
  flt_lst <- AnnotationFilter::AnnotationFilterList(flt_chr, flt_tbt, logicOp = "&")
  
  txs <- ensembldb::transcripts(edb, filter = flt_lst,
                                columns = c( "tx_id_version","gene_id", "tx_support_level", "tx_is_canonical"))
  
  #- Grab just the Canonical transcripts
  ind <- txs$tx_is_canonical == 1
  canonical_tx     <- txs[sort(which(ind))]
  
  return( data.frame(tx_id = canonical_tx$tx_id) )
  #- Output just a list of ENSTs
  write.table(data.frame(tx_id = canonical_tx$tx_id), "./functions/accessory_files/canonical_tx_set.txt", sep="\t", row.names=F, col.names = T,  quote = F)
}