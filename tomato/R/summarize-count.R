# Intra SNPs
x <- read.table("output/count/snp.txt")
x$V2 <- sub("input/raw/IL_SNP_gene/","",x$V2)
x$V2 <- sub("_SNP","",x$V2)
y <- read.csv("input/chr-il.txt",head=F)
z <- match(y$V1,x$V2)
colnames(x) <- c("chr","IL","M","I","D","Total")
write.csv(t(x[z,]),file="output/summary/snp.csv",quote=F)

# Inter SNPs
x <- read.table("output/count/othersnp.txt")
x$V2 <- sub("input/raw/IL_SNP_gene/","",x$V2)
x$V2 <- sub("_SNP","",x$V2)
z <- match(y$V1,x$V2)
colnames(x) <- c("chr","IL","M","I","D","Total")
write.csv(t(x[z,]),file="output/summary/othersnp.csv",quote=F)

# Both Intra and Inter SNPs
outputDir <- "output/IL_SNP_gene"
il.line <- c()
num.m <- c()
num.i <- c()
num.d <- c()
num.t <- c()
for (i in list.files(outputDir,pattern="IL")) 
{
  il.line <- c(il.line,i)
  x <- read.table(paste(outputDir,i,sep="/"),head=T)
  num.m <- c(num.m, sum(x$type=="M"))
  num.i <- c(num.i, sum(x$type=="I"))
  num.d <- c(num.d, sum(x$type=="D"))
  num.t <- c(num.t, nrow(x))
}
il.line <- sub("_SNP","",il.line)
x <- data.frame(il=il.line,
                m=num.m,
                i=num.i,
                d=num.d,
                t=num.t)
y <- read.csv("input/chr-il.txt",head=F)
z <- match(y$V1,x$il)
write.csv(t(x[z,]),file="output/summary/bothsnp.csv",quote=F)

i <- "pen_SNP"
x <- read.table(paste(outputDir,i,sep="/"),head=T)
num.m <- c(sum(x$type=="M"))
num.i <- c(sum(x$type=="I"))
num.d <- c(sum(x$type=="D"))
num.t <- c(nrow(x))
x <- data.frame(il=i,
                m=num.m,
                i=num.i,
                d=num.d,
                t=num.t)
write.csv(t(x),file="output/summary/bothsnp-pen.csv",quote=F)

# SNP types
x <- read.table("output/count/pen-snp-type.txt")
y <- data.frame(chr=x[x$V3=="gene",2],
                gene=x[x$V3=="gene",4],
                mrna=x[x$V3=="mRNA",4],
                cds=x[x$V3=="CDS",4],
                exon=x[x$V3=="exon",4],
                utr5=x[x$V3=="5utr",4],
                utr3=x[x$V3=="3utr",4],
                intron=x[x$V3=="intron",4])
write.csv(t(y),file="output/summary/pen-snp-type.csv",quote=F)

# cis/trans
outputDir <- "output/position3"
il.line <- c()
cis.pen <- c()
cis.notpen <- c()
trans.pen <- c()
trans.notpen <- c()
total <- c()
for (i in list.files(outputDir,pattern="IL")) 
{
  il.line <- c(il.line,i)
  x <- read.table(paste(outputDir,i,sep="/"),head=T)
  cis.pen <- c(cis.pen,sum(x$cis==TRUE & x$pen==TRUE))
  cis.notpen <- c(cis.notpen,sum(x$cis==TRUE & x$pen==FALSE))
  trans.pen <- c(trans.pen,sum(x$cis==FALSE & x$pen==TRUE))
  trans.notpen <- c(trans.notpen,sum(x$cis==FALSE & x$pen==FALSE))
  total <- c(total,nrow(x)) 
}

il.line <- sub(".csv","",il.line)
x <- data.frame(il=il.line,
                cis.pen=cis.pen,
                cis.notpen=cis.notpen,
                trans.pen=trans.pen,
                trans.notpen=trans.notpen,
                total=total)
y <- read.csv("input/chr-il.txt",head=F)
z <- match(y$V1,x$il)
write.csv(t(x[z,]),file="output/summary/position3.csv",quote=F)

x <- read.table("output/count/snp-uniq.txt")
y <- data.frame(chr=x[x$V2=="M",1],
                m=x[x$V2=="M",3], i=x[x$V2=="I",3], d=x[x$V2=="D",3],
                total=x[x$V2=="M",3] + x[x$V2=="I",3] + x[x$V2=="D",3])
write.csv(t(y),file="output/summary/snp-uniq.csv",quote=F)

