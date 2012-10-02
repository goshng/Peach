library(rtracklayer)

x.f <- function(x1) {
  nrow(d[d$gene==x1,])
}

# args <- "raw"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/ven-eqtl.R raw\n")
	quit("no")
}

raw <- args[1]
GFFFILE <- paste("data",raw,"gff.RData",sep="/")
load(GFFFILE)
gene.range.mrna <- gene.range[gene.range$type=="mRNA",]

# To create a table with position
inputDir <- paste("data",raw,"eqtl",sep="/")
outputDir <- "output/note"
x.gene <- c()
x.il <- c()
for (i in list.files(inputDir,pattern="IL")) 
{
  
  x10.1 <- read.csv(paste(inputDir,i,sep="/"))
  x10.1 <- x10.1[-1,]
  x10.1 <- x10.1[x10.1[,1]!="",]
  ilname <- sub(".csv","",i)
  x.gene <- c(x.gene,as.character(x10.1[,1]))
  x.il <- c(x.il, rep(ilname,nrow(x10.1)))
}
rm(i,x10.1,ilname,inputDir,gene.range)
d <- data.frame(gene=x.gene,il=x.il)

v <- sapply(unique(as.character(d$gene)),x.f)
v.sorted <- v[order(names(v))]
rm(v)

# This is WRONG!
# x <- x[match(x$ID,paste("mRNA:",gname,sep="")),]
# IMPORTANT: RangedData does not work with subsetting 
# chr1gene, chr2gene, chr1gene. The names have to be ordered in chromosomal
# indicies: chr1gene1, chr1gene1, chr2gene1, chr2gene2. 
#x <- gene.range.mrna[gene.range.mrna$ID %in% paste("mRNA:",names(v.sorted),sep=""),]
#y <- x[match(paste("mRNA:",names(v.sorted),sep=""),x$ID),]
# This seems to be a much better way.
stopifnot(length(unlist(gene.range.mrna$Note))==nrow(gene.range.mrna))
x.note <- unlist(gene.range.mrna$Note)[match(paste("mRNA:",names(v.sorted),sep=""),gene.range.mrna$ID)]

write.csv(data.frame(gene=names(v.sorted),il.count=v.sorted,note=x.note),
          row.names=FALSE,
          file="output/summary/ven-eqtl.csv")

