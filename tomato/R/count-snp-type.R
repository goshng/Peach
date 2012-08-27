library(rtracklayer)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
{                                                                               
  cat ("Rscript R/find-il-position.R IL1-1_SNP\n")
	quit("no")
}
raw <- args[1]
numchr <- as.numeric(args[2])
GFFFILE <- paste("data",raw,"ITAG2.3_gene_models.gff3",sep="/")

gene.range <- import.gff3(GFFFILE)
gene.range.mrna <- gene.range[gene.range$type=="mRNA",]
gene.range.cds <- gene.range[gene.range$type=="CDS",]
gene.range.exon <- gene.range[gene.range$type=="exon",]
gene.range.5utr <- gene.range[gene.range$type=="five_prime_UTR",]
gene.range.3utr <- gene.range[gene.range$type=="three_prime_UTR",]
gene.range.gene <- gene.range[gene.range$type=="gene",]
gene.range.intron <- gene.range[gene.range$type=="intron",]

snp.pen <- read.table("output/IL_SNP_gene/pen_SNP",head=T)
for (i in 0:numchr) {
  j <- sprintf("SL2.40ch%02d",i)
  x <- snp.pen[snp.pen[,2]==j,]
  query <- IRanges(x$position,x$position)

  subject <- ranges(gene.range)[[i+1,]]
  tree <- IntervalTree(subject)
  y <- countOverlaps(query, tree)
  cat(paste("Chr",i,"annotation",sum(y > 0)),"\n")
  cat(paste("Chr",i,"noannotation",sum(y==0)),"\n")
  cat(paste("Chr",i,"all",length(y)),"\n")

  subject <- ranges(gene.range.mrna)[[i+1,]]
  tree <- IntervalTree(subject)
  y <- countOverlaps(query, tree)
  cat(paste("Chr",i,"mRNA",sum(y > 0)),"\n")

  subject <- ranges(gene.range.cds)[[i+1,]]
  tree <- IntervalTree(subject)
  y <- countOverlaps(query, tree)
  cat(paste("Chr",i,"CDS",sum(y > 0)),"\n")

  subject <- ranges(gene.range.exon)[[i+1,]]
  tree <- IntervalTree(subject)
  y <- countOverlaps(query, tree)
  cat(paste("Chr",i,"exon",sum(y > 0)),"\n")

  subject <- ranges(gene.range.5utr)[[i+1,]]
  tree <- IntervalTree(subject)
  y <- countOverlaps(query, tree)
  cat(paste("Chr",i,"5utr",sum(y > 0)),"\n")

  subject <- ranges(gene.range.3utr)[[i+1,]]
  tree <- IntervalTree(subject)
  y <- countOverlaps(query, tree)
  cat(paste("Chr",i,"3utr",sum(y > 0)),"\n")

  subject <- ranges(gene.range.gene)[[i+1,]]
  tree <- IntervalTree(subject)
  y <- countOverlaps(query, tree)
  cat(paste("Chr",i,"gene",sum(y > 0)),"\n")

  subject <- ranges(gene.range.intron)[[i+1,]]
  tree <- IntervalTree(subject)
  y <- countOverlaps(query, tree)
  cat(paste("Chr",i,"intron",sum(y > 0)),"\n")
}

# This may be misleading.
#  subject <- IRanges(x$position,x$position)
#  tree <- IntervalTree(subject)
#
#  subject <- IRanges(x$position,x$position)
#  tree <- IntervalTree(subject)
#  y <- countOverlaps(ranges(gene.range.mrna)[[i+1,]], tree)
#  cat(paste("Chr",i,"mRNA",sum(y)),"\n")
#  y <- countOverlaps(ranges(gene.range.cds)[[i+1,]], tree)
#  cat(paste("Chr",i,"CDS",sum(y)),"\n")
#  y <- countOverlaps(ranges(gene.range.exon)[[i+1,]], tree)
#  cat(paste("Chr",i,"exon",sum(y)),"\n")
#  y <- countOverlaps(ranges(gene.range.5utr)[[i+1,]], tree)
#  cat(paste("Chr",i,"5utr",sum(y)),"\n")
#  y <- countOverlaps(ranges(gene.range.3utr)[[i+1,]], tree)
#  cat(paste("Chr",i,"3utr",sum(y)),"\n")
#  y <- countOverlaps(ranges(gene.range.gene)[[i+1,]], tree)
#  cat(paste("Chr",i,"gene",sum(y)),"\n")
#  y <- countOverlaps(ranges(gene.range.intron)[[i+1,]], tree)
#  cat(paste("Chr",i,"intron",sum(y)),"\n")
