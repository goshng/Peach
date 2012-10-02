library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
#ilname <- "IL12-2"
#args <- c(12, 
#  "data/karyotype/karyotype.sol.txt", 
#  paste("output/position3/",ilname,".csv",sep=""), 
#  paste("output/IL_SNP_gene/",ilname,"_SNP",sep=""), 
#  "output/ilmap2.dat",
#  "output/cis2", "output/trans2", "output/links2", "output/ilname", 
#  "output/cissnpinsideil", "output/cissnpoutsideil", "output/transsnpil")

if (length(args) != 12)
{                                                                               
  cat ("Rscript karyotype/karyotype.sol.txt R/eqtl-cis-trans.R position2/IL1-1 IL1-1 ilmap2.dat o/cis2 o/trans2 o/links2 o/ilname o/snp\n")
  quit("no")
}

chrLength <- read.table(args[2],head=F)[seq(as.integer(args[1])+1),c(3,6)]
# eQTL
x <- read.table(args[3],head=T)
# SNP
y <- read.table(args[4],head=T)
# IL map
z <- read.csv(args[5],head=T)

m <- regexec("(IL[0-9\\-]+)\\.csv",args[3])
ilname <- regmatches(args[3],m)[[1]][2]
z <- z[z$il==ilname & z$type=="main",]
m <- regexec("IL([0-9]+)",ilname)
mainILchr <- paste("sl",regmatches(ilname,m)[[1]][2],sep="")

print(paste("#######", ilname, "#######"))

# GRanges for SNPs
gr <-
  GRanges(seqnames =
          paste("sl",as.integer(sub("SL2.40ch","",y$chromosome)),sep=""),
          ranges =
          IRanges(start=y$position,width=1),
          strand = strand("*"))
seqlengths(gr) <- chrLength[match(names(seqlengths(gr)), chrLength[,1]),2]
# seqlengths(gr) <- chrLength[chrLength[,1] %in% runValue(seqnames(gr)),2]
gr <- flank(gr,width=50000,both=TRUE)
grIL <- 
  GRanges(seqnames = mainILchr,
          ranges = IRanges(start=z$start,end=z$end),
          strand = strand("*"))
seqlengths(grIL) <- chrLength[match(names(seqlengths(grIL)), chrLength[,1]),2]
# seqlengths(grIL) <- chrLength[chrLength[,1] %in% runValue(seqnames(grIL)),2]
gr <- c(gr,grIL)
# GRanges for eQTL
mainILstart <- as.integer(round((z$start + z$end)/2))
mainILend <- mainILstart + 10
grGene <- 
  GRanges(seqnames =
          sprintf("sl%d",as.integer(sub("Solyc","",sub("g.+","",x$Gene)))),
          ranges =
          IRanges(start=x$start,end=x$end),
          strand = strand("*"))
seqlengths(grGene) <- chrLength[match(names(seqlengths(grGene)), chrLength[,1]),2]
# seqlengths(grGene) <- chrLength[chrLength[,1] %in% runValue(seqnames(grGene)),2]
ol <- unique(as.matrix(findOverlaps(grGene,gr))[,1])

################################################################################
# cis2
i <- x[ol,]
write.table( 
  data.frame(chr=sprintf("sl%d",as.integer(sub("Solyc","",sub("g.+","",i$Gene)))),
             start=i$start, end=i$end, ratio=i$ratio),
  file=sprintf("%s/%s",args[6],ilname),
  row.names=FALSE, col.names=FALSE, quote=FALSE
)
cis.pen <- rep(FALSE,nrow(x))
cis.pen[ol] <- TRUE
write.table(
  data.frame(x,cis=cis.pen),
  file=sprintf("output/position4/%s.csv",ilname),
)

################################################################################
# trans2
i <- x[setdiff(seq(nrow(x)),ol),]
write.table( 
  data.frame(chr=sprintf("sl%d",as.integer(sub("Solyc","",sub("g.+","",i$Gene)))),
             start=i$start, end=i$end, ratio=i$ratio),
  file=sprintf("%s/%s",args[7],ilname),
  row.names=FALSE, col.names=FALSE, quote=FALSE
)

# links2
write.table( 
  data.frame(chrS=rep(mainILchr,nrow(i)),
             startS=rep(mainILstart,nrow(i)), endS=rep(mainILend,nrow(i)),
             chrD=sprintf("sl%d",as.integer(sub("Solyc","",sub("g.+","",i$Gene)))),
             startD=i$start, endD=i$end),
  file=sprintf("%s/%s",args[8],ilname),
  row.names=FALSE, col.names=FALSE, quote=FALSE
)

# ilname
write.table( 
  data.frame(chrS=mainILchr,
             startS=mainILstart, endS=mainILend,
             name=ilname),
  file=sprintf("%s/%s",args[9],ilname),
  row.names=FALSE, col.names=FALSE, quote=FALSE
)

################################################################################
# cis snp inside (within 50kbp) IL
gr <-
  GRanges(seqnames =
          paste("sl",as.integer(sub("SL2.40ch","",y$chromosome)),sep=""),
          ranges =
          IRanges(start=y$position,width=1),
          strand = strand("*"))
seqlengths(gr) <- chrLength[match(names(seqlengths(gr)), chrLength[,1]),2]
# seqlengths(gr) <- chrLength[chrLength[,1] %in% runValue(seqnames(gr)),2]
gr <- flank(gr,width=50000,both=TRUE)
ol <- unique(as.matrix(findOverlaps(gr,grIL))[,1])
i <- y[ol,]
write.table( 
  data.frame(chr=sprintf("sl%d",as.integer(sub("SL2.40ch","",i$chromosome))),
             start=i$position, end=i$position+10,
             color=rep("color=white",nrow(i))),
  file=sprintf("%s/%s",args[10],ilname),
  row.names=FALSE, col.names=FALSE, quote=FALSE
)

# cis snp outside IL
i.y <- y[setdiff(seq(nrow(y)),ol),]
gr <-
  GRanges(seqnames =
          paste("sl",as.integer(sub("SL2.40ch","",i.y$chromosome)),sep=""),
          ranges =
          IRanges(start=i.y$position,width=1),
          strand = strand("*"))
seqlengths(gr) <- chrLength[match(names(seqlengths(gr)), chrLength[,1]),2]
# seqlengths(gr) <- chrLength[chrLength[,1] %in% runValue(seqnames(gr)),2]
gr <- flank(gr,width=50000,both=TRUE)

ol <- unique(as.matrix(findOverlaps(gr,grGene))[,1])
i <- i.y[ol,]

write.table( 
  data.frame(chr=sprintf("sl%d",as.integer(sub("SL2.40ch","",i$chromosome))),
             start=i$position, end=i$position+10,
             color=rep("color=white",nrow(i))),
  file=sprintf("%s/%s",args[11],ilname),
  row.names=FALSE, col.names=FALSE, quote=FALSE
)

# trans snp
i <- i.y[setdiff(seq(nrow(i.y)),ol),]
write.table( 
  data.frame(chr=sprintf("sl%d",as.integer(sub("SL2.40ch","",i$chromosome))),
             start=i$position, end=i$position+10,
             color=rep("color=white",nrow(i))),
  file=sprintf("%s/%s",args[12],ilname),
  row.names=FALSE, col.names=FALSE, quote=FALSE
)

