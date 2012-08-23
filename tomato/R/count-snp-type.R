library(rtracklayer)
GFFFILE <- "input/ITAG2.3_gene_models.gff3.txt.no-negative"

gene.range <- import.gff3(GFFFILE)
gene.range.mrna <- gene.range[gene.range$type=="mRNA",]
gene.range.cds <- gene.range[gene.range$type=="CDS",]
gene.range.exon <- gene.range[gene.range$type=="exon",]
gene.range.5utr <- gene.range[gene.range$type=="five_prime_UTR",]
gene.range.3utr <- gene.range[gene.range$type=="three_prime_UTR",]
gene.range.gene <- gene.range[gene.range$type=="gene",]
gene.range.intron <- gene.range[gene.range$type=="intron",]

snp.pen <- read.table("output/IL_SNP_gene/pen_SNP",head=T)
for (i in 0:12) {
  j <- sprintf("SL2.40ch%02d",i)
  x <- snp.pen[snp.pen[,2]==j,]

  # Ch12 SNP in mRNA
  subject <- IRanges(x$position,x$position)
  tree <- IntervalTree(subject)
  y <- countOverlaps(ranges(gene.range.mrna)[[i+1,]], tree)
  cat(paste("Chr",i,"mRNA",sum(y)),"\n")
  y <- countOverlaps(ranges(gene.range.cds)[[i+1,]], tree)
  cat(paste("Chr",i,"CDS",sum(y)),"\n")
  y <- countOverlaps(ranges(gene.range.exon)[[i+1,]], tree)
  cat(paste("Chr",i,"exon",sum(y)),"\n")
  y <- countOverlaps(ranges(gene.range.5utr)[[i+1,]], tree)
  cat(paste("Chr",i,"5utr",sum(y)),"\n")
  y <- countOverlaps(ranges(gene.range.3utr)[[i+1,]], tree)
  cat(paste("Chr",i,"3utr",sum(y)),"\n")
  y <- countOverlaps(ranges(gene.range.gene)[[i+1,]], tree)
  cat(paste("Chr",i,"gene",sum(y)),"\n")
  y <- countOverlaps(ranges(gene.range.intron)[[i+1,]], tree)
  cat(paste("Chr",i,"intron",sum(y)),"\n")
}


