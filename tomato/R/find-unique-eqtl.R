args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/count-snp-type.R raw 12\n")
	quit("no")
}
#args <- "smallraw"
raw <- args[1]

list.il <- scan(paste("data",raw,"order.csv",sep="/"),what="character")

v <- c()
for (chrI in 1:12) {
  list.il.chr <- list.il[grep(paste("IL",chrI,"-",sep=""),list.il)]
  z <- c()
  for (ilname in list.il.chr) {
    y <- read.csv(paste("data/",raw,"/eqtl/",ilname,".csv",sep=""))
    y <- y[-1,]                                                           
    y <- y[y[,1]!="",]
    z <- c(z,as.character(y[,1]))
  }
  v <- c(v,length(unique(z)))
}
names(v) <- 1:12
write.csv(t(v),file="output/summary/eqtl-uniq.csv",row.names=F)
