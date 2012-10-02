library(GenomicRanges)
library(rtracklayer)
library(Hmisc)
library(corrgram)
source("R/rcorr.R")

#args <- c("raw", "IL_mat_64.txt")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
{                                                                               
  cat ("Rscript R/correlation-metabolite.R raw IL_mat_30.txt\n")
  quit("no")
}
raw <- args[1]
metabolite.file <- args[2]

################################################################################
# Read total expression values
# The first three lines are Gene, M82, and LA716.
# z2 are T/F for expression levels (> 1) of the 74 IL lines for genes.
# Select genes with IL values that are larger than 1, and M82 are positive.
x <- read.csv(paste("data",raw,"IL-expressionforcorrelation.csv",sep="/"))
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
# Read metabolite quantity
m <- read.table(paste("data",raw,"metabolite",metabolite.file,sep="/"),head=TRUE,sep="\t")
m[m <= 0] <- NA
m <- t(log2(m))
stopifnot(sum(rownames(m) %in% rownames(z5)) == nrow(m))
# This is WRONG!
# z6 <- z5[rownames(z5) %in% rownames(m),]
# z5 <- z6[match(rownames(z6),rownames(m)),]
# > match(c(2,4,5,1,3),1:5)
# [1] 2 4 5 1 3
# > match(1:5,c(2,4,5,1,3))
# [1] 4 1 5 2 3
# > c(2,4,5,1,3)[match(1:5,c(2,4,5,1,3))]
# [1] 1 2 3 4 5
#z6 <- z5[rownames(z5) %in% rownames(m),]
#z5 <- z6[match(rownames(m),rownames(z6)),]
#rm(z6)
z5 <- z5[match(rownames(m),rownames(z5)),]
stopifnot(sum(match(rownames(m),rownames(z5)) == seq(nrow(m))) == nrow(m))
# z5 (gene expression) and m (metabolite quantities) are matched in their IL.
################################################################################

################################################################################
# Compute p-values and correlations in pairwise fashion.
m.pval <- c()
m.cor <- c()
for (i in seq(ncol(z5))) {
  x <- c()
  y <- c()
  for (j in seq(ncol(m))) {
    a <- na.omit(cbind(z5[,i],m[,j]))

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
colnames(m.pval) <- colnames(z5)
rownames(m.pval) <- colnames(m)
colnames(m.cor) <- colnames(z5)
rownames(m.cor) <- colnames(m)
m.pval <- t(m.pval)
m.cor <- t(m.cor)
rm(x,y,a,b,i,j)
# m.pval and m.cor: row is gene, col is metabolite.
################################################################################

################################################################################
# Adjust p-values.
#m.pval <- p.adjust(m.pval, method = "BH")
#dim(m.pval) <- dim(m.cor)
#dimnames(m.pval) <- dimnames(m.cor)

################################################################################
# Append annotations to correlation and correlation's p-value.
GFFFILE <- paste("data",raw,"gff.RData",sep="/")
load(GFFFILE)
gene.range.mrna <- gene.range[gene.range$type=="mRNA",]
# This is WRONG!
# x <- x[match(x$ID,paste("mRNA:",rownames(m.pval),sep="")),]
# This may not be wrong, but it is not safe to subset RangedData object.
#x <- gene.range.mrna[gene.range.mrna$ID %in% paste("mRNA:",rownames(m.pval),sep=""),]
#x <- x[match(paste("mRNA:",rownames(m.pval),sep=""),x$ID),]
stopifnot(sum(gene.range.mrna$ID %in% paste("mRNA:",rownames(m.pval),sep="")) == nrow(m.pval))
stopifnot(length(unlist(gene.range.mrna$Note))==nrow(gene.range.mrna))
x.note <- unlist(gene.range.mrna$Note)[match(paste("mRNA:",rownames(m.pval),sep=""),gene.range.mrna$ID)]

b1 <- apply(m.pval < 0.1,1,sum,na.rm=T)
b2 <- apply(m.pval < 0.01,1,sum,na.rm=T)
b3 <- apply(m.pval < 0.001,1,sum,na.rm=T)
b4 <- apply(m.pval < 0.0001,1,sum,na.rm=T)

f <- data.frame(m.pval,n1=b1,n2=b2,n3=b3,n4=b4,note=x.note)
write.csv(f,file=paste("output/metabolite/",metabolite.file,"-pval.csv",sep=""))
f <- data.frame(m.cor,n1=b1,n2=b2,n3=b3,n4=b4,note=x.note)
write.csv(f,file=paste("output/metabolite/",metabolite.file,"-cor.csv",sep=""))
rm(gene.range,gene.range.mrna)

################################################################################
# Network file by p-value
m.t <- which(m.pval < 0.01,arr.in=TRUE)
write.table( 
  data.frame(a=rownames(m.t), 
             c=rep("cor",nrow(m.t)),
             b=colnames(m.pval)[m.t[,2]]),
  file=paste("output/metabolite/",metabolite.file,"-cor.sif",sep=""),
  row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t"
)

x5 <-rcorr(m)
# Ignore when n is smaller than 5.
x5$r[x5$n < 5] <- 0.0
x5$P[x5$n < 5] <- 1.0
colnames(x5$r) <- colnames(m)
rownames(x5$r) <- colnames(m)
colnames(x5$n) <- colnames(m)
rownames(x5$n) <- colnames(m)
colnames(x5$P) <- colnames(m)
rownames(x5$P) <- colnames(m)

b1 <- apply(x5$P < 0.1,1,sum,na.rm=T)
b2 <- apply(x5$P < 0.01,1,sum,na.rm=T)
b3 <- apply(x5$P < 0.001,1,sum,na.rm=T)
b4 <- apply(x5$P < 0.0001,1,sum,na.rm=T)

f <- data.frame(x5$r,n1=b1,n2=b2,n3=b3,n4=b4)
write.csv(f,file=paste("output/metabolite-cor/",metabolite.file,".csv",sep=""))
f <- data.frame(x5$n)
write.csv(f,file=paste("output/metabolite-cor/",metabolite.file,"-n.csv",sep=""))
f <- data.frame(x5$P,n1=b1,n2=b2,n3=b3,n4=b4)
write.csv(f,file=paste("output/metabolite-cor/",metabolite.file,"-pval.csv",sep=""))

################################################################################
# Network file by p-value
m.t <- which(x5$P < 0.01,arr.in=TRUE)
write.table( 
  data.frame(a=rownames(m.t), 
             c=rep("cor",nrow(m.t)),
             b=colnames(x5$P)[m.t[,2]]),
  file=paste("output/metabolite-cor/",metabolite.file,"-cor.sif",sep=""),
  row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t"
)



