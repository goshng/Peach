args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/find-il-position.R IL1-1_SNP\n")
	quit("no")
}
x <- sub("_SNP","",args[1])
x <- sub("snp","ilposition",x)
outfile <- paste(x,"csv",sep=".")

x <- read.table(args[1])

y <- c(min(x[,2]), max(x[,2]), trunc((max(x[,2]) + min(x[,2]))/2))
cat(y,"\n",file=outfile)
