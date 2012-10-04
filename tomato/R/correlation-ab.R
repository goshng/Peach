library(GenomicRanges)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5)
{                                                                               
  cat ("Rscript R/correlation-ab.R raw IL-expressionforcorrelation.csv for-apriori-subnetwork-a.csv for-apriori-subnetwork-b.csv ab\n")
  quit("no")
}
raw <- args[1]
expGene <- args[2]
geneSetA <- args[3]
geneSetB <- args[4]
ab.file <- args[5]

raw <- "raw"
expGene <- "IL-expressionforcorrelation.csv"
geneSetA <- "for-apriori-subnetwork-a.csv"
geneSetB <- "for-apriori-subnetwork-b.csv"
ab.file <- "ab"

################################################################################
# Read total expression values
# The first three lines are Gene, M82, and LA716.
# z2 are T/F for expression levels (> 1) of the 74 IL lines for genes.
# Select genes with IL values that are larger than 1, and M82 are positive.
x <- read.csv(paste("data",raw,expGene,sep="/"))
z2 <- x[,-1:-3] > 1
w <- apply(z2,1,sum) > 0 & x[,2] > 0
gname <- x[w,1]
m82 <- x[w,2]
z3 <- x[w,-1:-3]
z4 <- data.matrix(log2(z3/m82))
z5 <- apply(z4,1:2,function(x) if (is.infinite(x)) NA else x)
rownames(z5) <- gname
z5 <- t(z5)
rm(x,w,z2,z3,z4,m82,gname)
# z5 is the log2(IL expression/M82 expression) for those genes.
################################################################################

################################################################################
# Read 
geneNameSetA <- scan(paste("data",raw,geneSetA,sep="/"),what="character")
geneNameSetB <- scan(paste("data",raw,geneSetB,sep="/"),what="character")
geneNameSetA <- geneNameSetA[geneNameSetA %in% colnames(z5)]
geneNameSetB <- geneNameSetB[geneNameSetB %in% colnames(z5)]

stopifnot(ncol(z5) > length(geneNameSetA))
stopifnot(sum(geneNameSetA %in% colnames(z5)) == length(geneNameSetA))
z5A <- z5[,match(geneNameSetA,colnames(z5))]
stopifnot(sum(colnames(z5A) == geneNameSetA) == length(geneNameSetA))

stopifnot(ncol(z5) > length(geneNameSetB))
stopifnot(sum(geneNameSetB %in% colnames(z5)) == length(geneNameSetB))
z5B <- z5[,match(geneNameSetB,colnames(z5))]
stopifnot(sum(colnames(z5B) == geneNameSetB) == length(geneNameSetB))
# z5A and z5B are now two sets of gene expression.
################################################################################

################################################################################
# Compute p-values and correlations in pairwise fashion.
m.pval <- c()
m.cor <- c()
for (i in seq(ncol(z5A))) {
  x <- c()
  y <- c()
  for (j in seq(ncol(z5B))) {
    a <- na.omit(cbind(z5A[,i],z5B[,j]))

    if (nrow(a) >= 5)
    {
      b <- cor.test(a[,1],a[,2])$p.value
      x <- c(x,b)
      b <- cor.test(a[,1],a[,2])$estimate
      y <- c(y,b)
    }
    else
    {
      x <- c(x,NA)
      y <- c(y,NA)
    }
  }
  m.pval <- cbind(m.pval, x)
  m.cor <- cbind(m.cor, y)
}
colnames(m.pval) <- colnames(z5A)
rownames(m.pval) <- colnames(z5B)
colnames(m.cor) <- colnames(z5A)
rownames(m.cor) <- colnames(z5B)
m.pval <- t(m.pval)
m.cor <- t(m.cor)
rm(x,y,a,b,i,j)
# m.pval and m.cor: row is gene set A, col is gene set B. 
################################################################################

################################################################################
# Append annotations to correlation and correlation's p-value.
GFFFILE <- paste("data",raw,"gff.RData",sep="/")
load(GFFFILE)
gene.range.mrna <- gene.range[gene.range$type=="mRNA",]
stopifnot(sum(gene.range.mrna$ID %in% paste("mRNA:",rownames(m.pval),sep="")) == nrow(m.pval))
stopifnot(length(unlist(gene.range.mrna$Note))==nrow(gene.range.mrna))
x.note <- unlist(gene.range.mrna$Note)[match(paste("mRNA:",rownames(m.pval),sep=""),gene.range.mrna$ID)]

b1 <- apply(m.pval < 0.1,1,sum,na.rm=T)
b2 <- apply(m.pval < 0.01,1,sum,na.rm=T)
b3 <- apply(m.pval < 0.001,1,sum,na.rm=T)
b4 <- apply(m.pval < 0.0001,1,sum,na.rm=T)

f <- data.frame(m.pval,n1=b1,n2=b2,n3=b3,n4=b4,note=x.note)
write.csv(f,file=paste("output/cor-ab/",ab.file,"-pval.csv",sep=""))
f <- data.frame(m.cor,n1=b1,n2=b2,n3=b3,n4=b4,note=x.note)
write.csv(f,file=paste("output/cor-ab/",ab.file,"-cor.csv",sep=""))

rm(gene.range,gene.range.mrna)

