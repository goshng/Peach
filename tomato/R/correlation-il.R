library(rtracklayer)
library(Hmisc)
library(corrgram)
source("R/rcorr.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/correlation-il.R raw\n")
  quit("no")
}
raw <- args[1]
#raw <- "raw"

GFFFILE <- paste("data",raw,"gff.RData",sep="/")
load(GFFFILE)
gene.range.mrna <- gene.range[gene.range$type=="mRNA",]

list.il <- scan(paste("data",raw,"order.csv",sep="/"),what="character")
# list.il <- c("IL1-1-3", "IL1-1", "IL3-4", "IL6-2", "IL8-2-1", "IL11-1", "IL12-1")
# list.il <- c("IL5-3")
if (raw == "smallraw") {
  list.il <- c("IL1-1")
}

for (ilname in list.il) {
print(ilname)
x <- read.csv(paste("data",raw,"IL-expressionforcorrelation.csv",sep="/"))
# Select eQTL genes in the IL
y <- read.csv(paste("data/",raw,"/eqtl/",ilname,".csv",sep=""))
y <- y[-1,]                                                           
y <- y[y[,1]!="",]
z <- x[x$Gene %in% y[,1],]
rm(y)
# z5 is eQTL's expression levels of the gene (by column)
z2 <- z[,-1:-3] > 1
w <- apply(z2,1,sum) > 0 & z[,2] > 0
gname <- z[w,1]
m82 <- z[w,2]
z3 <- z[w,-1:-3]
z4 <- data.matrix(log2(z3/m82))
z5 <- apply(z4,1:2,function(x) if (is.infinite(x)) NA else x)
rownames(z5) <- gname
z5 <- t(z5)
rm(x,w,z2,z3,z4,m82)

x5 <-rcorr(z5)
# Ignore when n is smaller than 5.
x5$r[x5$n < 5] <- 0.0
x5$P[x5$n < 5] <- 1.0

colnames(x5$r) <- gname
rownames(x5$r) <- gname
colnames(x5$n) <- gname
rownames(x5$n) <- gname
colnames(x5$P) <- gname
rownames(x5$P) <- gname

b1 <- apply(x5$P < 0.1,1,sum,na.rm=T)
b2 <- apply(x5$P < 0.01,1,sum,na.rm=T)
b3 <- apply(x5$P < 0.001,1,sum,na.rm=T)
b4 <- apply(x5$P < 0.0001,1,sum,na.rm=T)

# This is WRONG!
# x <- x[match(x$ID,paste("mRNA:",gname,sep="")),]
# This may not be wrong, but it is not safe to subset RangedData object.
# x <- gene.range.mrna[gene.range.mrna$ID %in% paste("mRNA:",gname,sep=""),]
# x <- x[match(paste("mRNA:",gname,sep=""),x$ID),]
stopifnot(length(unlist(gene.range.mrna$Note))==nrow(gene.range.mrna))
x.note <- unlist(gene.range.mrna$Note)[match(paste("mRNA:",gname,sep=""),gene.range.mrna$ID)]

f <- data.frame(x5$r,n1=b1,n2=b2,n3=b3,n4=b4,note=x.note)
write.csv(f,file=paste("output/cor/",ilname,".csv",sep=""))
f <- data.frame(x5$n,note=x.note)
write.csv(f,file=paste("output/cor/",ilname,"-n.csv",sep=""))
f <- data.frame(x5$P,n1=b1,n2=b2,n3=b3,n4=b4,note=x.note)
write.csv(f,file=paste("output/cor/",ilname,"-pval.csv",sep=""))

################################################################################
# Network file by p-value
m.t <- which(x5$P < 0.01,arr.in=TRUE)
write.table( 
  data.frame(a=rownames(m.t), 
             c=rep("cor",nrow(m.t)),
             b=colnames(x5$P)[m.t[,2]]),
  file=paste("output/cor/",ilname,"-cor.sif",sep=""),
  row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t"
)

pdf(file=paste("output/cor/",ilname,".pdf",sep=""))
corrgram(x5$r,order=TRUE, lower.panel=panel.shade,
   upper.panel=panel.pie, text.panel=panel.txt,
   main=paste("Correlation Test",ilname))
dev.off()
rm(b1,b2,b3,b4)

}
