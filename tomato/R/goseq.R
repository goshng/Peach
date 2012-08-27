library(rtracklayer)
library(goseq)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/find-il-position.R IL1-1_SNP\n")
  quit("no")
}
raw <- args[1]
GFFFILE <- paste("data",raw,"ITAG2.3_gene_models.gff3",sep="/")

gene.range <- import.gff3(GFFFILE)
gene.range.mrna <- gene.range[gene.range$type=="mRNA",]
gene.range.mrna$score <- 0
colnames(gene.range.mrna)[colnames(gene.range.mrna) == "Name"] <- "name"

f <- "data/go/mrna.bed"
export(gene.range.mrna,f,format="bed")
smutans.feature.genes = read.table(file=f,head=F)
stopifnot(ncol(smutans.feature.genes)==6)

f <- "data/go/smutans.gene2go"
smutans.go.genes <- read.table(file=f,head=F)
smutans.go.genes <- smutans.go.genes[,c(1,4,5)]
stopifnot(ncol(smutans.go.genes)==3)

f <- "data/go/smutans.go2ngene"
smutans.cat.desc <- read.table(file=f,head=F,sep="\t",quote="")                  
stopifnot(ncol(smutans.cat.desc)==3)  

feature.genes <- smutans.feature.genes
go.genes <- smutans.go.genes
cat.desc <- smutans.cat.desc 
qval <- 0.05
length.genes = feature.genes$V3 - feature.genes$V2                          
assayed.genes = feature.genes$V4                                            

inputDir <- paste("data",raw,"eqtl",sep="/")
outputDir <- "output/go"
for (i in list.files(inputDir,pattern="csv")) 
{
  print(i)
  x10.1 <- read.csv(paste(inputDir,i,sep="/"))
  x10.1 <- x10.1[-1,]
  x10.1 <- x10.1[x10.1[,1]!="",]
  de.genes <- x10.1$Gene
                                                                                  
  gene.vector = as.integer(assayed.genes %in% de.genes)                       
  names(gene.vector) = assayed.genes                                          
  pwf = nullp( gene.vector, bias.data=length.genes, plot.fit=FALSE )          
                                                                                  
  go = goseq(pwf,gene2cat=go.genes,method="Wallenius") # Length bias correction - Approximation
  go.fdr <- data.frame(go, qval=p.adjust(go$over_represented_pvalue,method="BH"))
  go.fdr <- go.fdr[go.fdr$qval < qval,]

  for (j in go.fdr$category)
  {
    cat(file=paste(outputDir,i,sep="/"), j, ",", as.character(cat.desc$V2[cat.desc$V1==j]), ",", go.fdr$qval[go.fdr$category==j], ",", as.character(cat.desc$V3[cat.desc$V1==j]), append=T)
    cat(file=paste(outputDir,i,sep="/"), "\n", append=T)
    # cat(file=paste(outputDir,i,sep="/"), "\\\\\n", append=T)
  }
}
