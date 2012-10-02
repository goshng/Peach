args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/find-il-position.R IL1-1_SNP\n")
  quit("no")
}
raw <- args[1]

## Intra SNPs
#x <- read.table("output/count/snp.txt")
#x$V2 <- sub(paste("data",raw,"IL_SNP_gene/",sep="/"),"",x$V2)
#x$V2 <- sub("_SNP","",x$V2)
#y <- read.csv(paste("data",raw,"chr-il.txt",sep="/"),head=F)
#z <- match(y$V1,x$V2)
#colnames(x) <- c("chr","IL","M","MD","I","D","Total")
#write.csv(t(x[z,]),file="output/summary/snp.csv",quote=F)
#
## Inter SNPs
#x <- read.table("output/count/othersnp.txt")
#x$V2 <- sub(paste("data",raw,"IL_SNP_gene/",sep="/"),"",x$V2)
#x$V2 <- sub("_SNP","",x$V2)
#z <- match(y$V1,x$V2)
#colnames(x) <- c("chr","IL","M","MD","I","D","Total")
#write.csv(t(x[z,]),file="output/summary/othersnp.csv",quote=F)

################################################################################
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
  # print(paste(levels(x$type)))
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
# PEN
i.pen <- "pen_SNP"
x.pen <- read.table(paste(outputDir,i,sep="/"),head=T)
# print(paste(levels(x.pen$type)))
num.m <- c(sum(x.pen$type=="M"))
num.md <- c(sum(x.pen$type=="MD"))
num.i <- c(sum(x.pen$type=="I"))
num.d <- c(sum(x.pen$type=="D"))
num.t <- c(nrow(x.pen))
x.pen <- data.frame(il=i.pen,
                m=num.m,
                md=num.md,
                i=num.i,
                d=num.d,
                t=num.t)
x <- rbind(x,x.pen)
y <- read.csv(paste("data",raw,"chr-il.txt",sep="/"),head=F)
z <- match(y$V1,x$il)
z <- c(z,length(z)+1)
write.csv(x[z,],file="output/summary/bothsnp.csv",quote=F)
# write.csv(t(x),file="output/summary/bothsnp-pen.csv",quote=F)

################################################################################
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

################################################################################
# cis/trans eQTL
outputDir <- "output/position4"
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

#  cis.pen <- c(cis.pen,sum(x$withinil == TRUE & x$samechr == TRUE & x$pen==TRUE))
#  cis.notpen <- c(cis.notpen,sum(x$withinil == TRUE & x$samechr == TRUE & x$pen==FALSE))
#  trans.pen <- c(trans.pen,sum((x$withinil == FALSE | x$samechr == FALSE) & x$pen==TRUE))
#  trans.notpen <- c(trans.notpen,sum((x$withinil == FALSE | x$samechr == FALSE) & x$pen==FALSE))

  cis.pen <- c(cis.pen,sum(x$cis == TRUE & x$pen==TRUE))
  cis.notpen <- c(cis.notpen,sum(x$cis == TRUE & x$pen==FALSE))
  trans.pen <- c(trans.pen,sum(x$cis == FALSE & x$pen==TRUE))
  trans.notpen <- c(trans.notpen,sum(x$cis == FALSE & x$pen==FALSE))
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
write.csv(x[z,],file="output/summary/cistranspenlike.csv",quote=F)
################################################################################
# cis/trans pen-like/not by chromosome
outputDir <- "output/position4"
il.line <- c()
cis.pen <- c()
cis.notpen <- c()
trans.pen <- c()
trans.notpen <- c()
notunique.count <- c()
total <- c()
list.il <- scan(paste("data",raw,"order.csv",sep="/"),what="character")
for (chrI in 1:12) {
  list.il.chr <- list.il[grep(paste("IL",chrI,"-",sep=""),list.il)]

  x.chr <- c()
  for (ilname in list.il.chr) {
    x <- read.table(paste(outputDir,"/",ilname,".csv",sep=""),head=T)
    x <- cbind(x,il=ilname)
    x.chr <- rbind(x.chr,x)
  }
   
  notunique.each <- 0
  for (i in levels(x.chr[,2])) {
    #stopifnot(length(unique(apply(x.chr[x.chr$Gene==i,c("Gene","pen","cis")],1,paste,collapse=""))) == 1)
    a = unique(apply(x.chr[x.chr$Gene==i,c("Gene","pen","cis")],1,paste,collapse=""))
    if (length(a) != 1) {
      print(paste(chrI,i,length(a),a))
      print(x.chr[x.chr$Gene==i,c("il","Gene","pen","cis")])
      cat("\n")
      notunique.each <- notunique.each + length(a) - 1
    }
  }
  cis.pen <- c(cis.pen,length(unique(x.chr[x.chr$cis==TRUE & x.chr$pen==TRUE,"Gene"])))
  cis.notpen <- c(cis.notpen,length(unique(x.chr[x.chr$cis==TRUE & x.chr$pen==FALSE,"Gene"])))
  trans.pen <- c(trans.pen,length(unique(x.chr[x.chr$cis==FALSE & x.chr$pen==TRUE,"Gene"])))
  trans.notpen <- c(trans.notpen,length(unique(x.chr[x.chr$cis==FALSE & x.chr$pen==FALSE,"Gene"])))
  notunique.count <- c(notunique.count, notunique.each) 
  total <- c(total,length(unique(x.chr[,"Gene"])))
}
stopifnot(sum(total + notunique.count == (cis.pen + cis.notpen + trans.pen + trans.notpen)) == 12)
x <- data.frame(chr=1:12,
                cis.pen=cis.pen,
                cis.notpen=cis.notpen,
                trans.pen=trans.pen,
                trans.notpen=trans.notpen,
                notunique=notunique.count,
                total=total)
write.csv(x,file="output/summary/cistranspenlike-bychromosome.csv",quote=F)

################################################################################
# Unique SNP types
x <- read.table("output/count/snp-uniq.txt")
y <- data.frame(chr=x[x$V2=="M",1],
                m=x[x$V2=="M",3], i=x[x$V2=="I",3], d=x[x$V2=="D",3],
                total=x[x$V2=="M",3] + x[x$V2=="I",3] + x[x$V2=="D",3])
write.csv(t(y),file="output/summary/snp-uniq.csv",quote=F)

