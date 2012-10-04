library(rtracklayer)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/find-il-position.R IL1-1_SNP\n")
	quit("no")
}
raw <- args[1]
GFFFILE <- paste("data",raw,"gff.RData",sep="/")
load(GFFFILE)
gene.range.mrna <- gene.range[gene.range$type=="mRNA",]

# To create a table with position
inputDir <- paste("data",raw,"eqtl",sep="/")
outputDir <- "output/note"
for (i in list.files(inputDir,pattern="csv")) 
{
  x10.1 <- read.csv(paste(inputDir,i,sep="/"))
  x10.1 <- x10.1[-1,]
  x10.1 <- x10.1[x10.1[,1]!="",]
  x10.1.pos <- gene.range.mrna[gene.range.mrna$ID %in% paste("mRNA:",x10.1[,1],sep=""),]

  x10.1.2 <- data.frame(Gene=sub("mRNA:","",x10.1.pos$ID),note=unlist(x10.1.pos$Note))
  colnames(x10.1) <- c("Gene","meanA","meanB","ratio","padj","b1","b2","z1")
  x10.1.3 <- merge(x10.1, x10.1.2, sort = FALSE)
  if (nrow(x10.1) == nrow(x10.1.2)) {
	  write.csv(x10.1.3,file=paste(outputDir,i,sep="/"))
    cat(paste(i,",",nrow(x10.1)),"\n")
  } else {
    cat(paste("Error not all genes were found",i),"\n")
  }
}
cat(paste("Total number of mRNAs,",nrow(gene.range.mrna)),"\n")

