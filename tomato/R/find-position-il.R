args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/find-il-position.R IL1-1_SNP\n")
	quit("no")
}
# x <- read.csv(paste("input/snp/chr",i,"SNPINDELS.csv",sep=""))
# x <- read.csv("input/snp/chr10SNPINDELS.csv")

x <- sub("_SNP","",args[1])
x <- sub("snp","ilposition",x)
outfile <- paste(x,"csv",sep=".")

x <- read.table(args[1])

y <- c(min(x[,2]), max(x[,2]), trunc((max(x[,2]) + min(x[,2]))/2))
cat(y,"\n",file=outfile)

# x1 <- c()
# x2 <- c()
# for (i in grep("IL",colnames(x)))
# {
  # x1 <- c(x1,colnames(x)[i])
  # p <- (max(x[,i],na.rm=T) + min(x[,i],na.rm=T))/2
  # x2 <- c(x2,p)
# }

# y <- data.frame(il=x1,pos=x2)
# write.table(y,file="output/10.il",quote=F,row.names=F,col.names=F)
