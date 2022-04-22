args <- commandArgs(trailingOnly=TRUE)

gtf <- args[1]
threads <- args[2]

files <- list.files(path=".", pattern=".bam")
  
rs <-
  Rsubread::featureCounts(files=files, annot.ext=gtf, isGTFAnnotationFile=TRUE,
                          GTF.featureType="exon", GTF.attrType="gene_id", GTF.attrType.extra=NULL,
                          chrAliases=NULL, useMetaFeatures=TRUE, allowMultiOverlap=FALSE, 
                          isPairedEnd=TRUE, nthreads=threads)

invisible(lapply(c("counts", "annotation", "stat"), function(x){
  
  wr <- rs[[x]]
  
  #/ add genes which are the rownames directly to the matrix
  if(x=="counts") wr <- data.frame(Gene=rownames(wr), wr)
  
  write.table(wr, file=paste0(x, ".txt"),
              col.names=TRUE, row.names=FALSE,
              quote=FALSE, sep="\t")
 
  NULL
  
}))

