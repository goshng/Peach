args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/find-il-position.R IL1-1_SNP\n")
  quit("no")
}
raw <- args[1]

# Intra SNPs
x <- read.table("output/count/snp.txt")
x$V2 <- sub(paste("data",raw,"IL_SNP_gene/",sep="/"),"",x$V2)
x$V2 <- sub("_SNP","",x$V2)
y <- read.csv(paste("data",raw,"chr-il.txt",sep="/"),head=F)
z <- match(y$V1,x$V2)
colnames(x) <- c("chr","IL","M","MD","I","D","Total")
write.csv(t(x[z,]),file="output/summary/snp.csv",quote=F)

# Inter SNPs
x <- read.table("output/count/othersnp.txt")
x$V2 <- sub(paste("data",raw,"IL_SNP_gene/",sep="/"),"",x$V2)
x$V2 <- sub("_SNP","",x$V2)
z <- match(y$V1,x$V2)
colnames(x) <- c("chr","IL","M","MD","I","D","Total")
write.csv(t(x[z,]),file="output/summary/othersnp.csv",quote=F)

# Both Intra and Inter SNPs
outputDir <- "output/IL_SNP_gene"
il.line <- c()
num.m <- c()
num.md <- c()
num.i <- c()
num.d <- c()
num.t <- c()
for (i in list.files(outputDir,pattern="IL")) 
{
  il.line <- c(il.line,i)
  x <- read.table(paste(outputDir,i,sep="/"),head=T)
  print(paste(levels(x$type)))
  num.m <- c(num.m, sum(x$type=="M"))
  num.md <- c(num.md, sum(x$type=="MD"))
  num.i <- c(num.i, sum(x$type=="I"))
  num.d <- c(num.d, sum(x$type=="D"))
  num.t <- c(num.t, nrow(x))
}
il.line <- sub("_SNP","",il.line)
x <- data.frame(il=il.line,
                m=num.m,
                md=num.md,
                i=num.i,
                d=num.d,
                t=num.t)
y <- read.csv(paste("data",raw,"chr-il.txt",sep="/"),head=F)
z <- match(y$V1,x$il)
write.csv(t(x[z,]),file="output/summary/bothsnp.csv",quote=F)

i <- "pen_SNP"
x <- read.table(paste(outputDir,i,sep="/"),head=T)
print(paste(levels(x$type)))
num.m <- c(sum(x$type=="M"))
num.md <- c(sum(x$type=="MD"))
num.i <- c(sum(x$type=="I"))
num.d <- c(sum(x$type=="D"))
num.t <- c(nrow(x))
x <- data.frame(il=i,
                m=num.m,
                md=num.md,
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
                intron=x[x$V3=="intron",4],
                annotation=x[x$V3=="annotation",4],
                noannotation=x[x$V3=="noannotation",4],
                all=x[x$V3=="all",4])
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

  cis.pen <- c(cis.pen,sum(x$withinil == TRUE & x$samechr == TRUE & x$pen==TRUE))
  cis.notpen <- c(cis.notpen,sum(x$withinil == TRUE & x$samechr == TRUE & x$pen==FALSE))

  trans.pen <- c(trans.pen,sum((x$withinil == FALSE | x$samechr == FALSE) & x$pen==TRUE))
  trans.notpen <- c(trans.notpen,sum((x$withinil == FALSE | x$samechr == FALSE) & x$pen==FALSE))
  total <- c(total,nrow(x)) 
}

il.line <- sub(".csv","",il.line)
x <- data.frame(il=il.line,
                cis.pen=cis.pen,
                cis.notpen=cis.notpen,
                trans.pen=trans.pen,
                trans.notpen=trans.notpen,
                total=total)
y <- read.csv(paste("data",raw,"chr-il.txt",sep="/"),head=F)
z <- match(y$V1,x$il)
write.csv(t(x[z,]),file="output/summary/cistranspenlike.csv",quote=F)

x <- read.table("output/count/snp-uniq.txt")
y <- data.frame(chr=x[x$V2=="M",1],
                m=x[x$V2=="M",3], i=x[x$V2=="I",3], d=x[x$V2=="D",3],
                total=x[x$V2=="M",3] + x[x$V2=="I",3] + x[x$V2=="D",3])
write.csv(t(y),file="output/summary/snp-uniq.csv",quote=F)

