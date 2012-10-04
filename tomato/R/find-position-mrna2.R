outputDir <- "output/position"
outputDir2 <- "output/position2"
for (i in list.files(outputDir,pattern="csv")) 
{
  x <- read.csv(paste(outputDir,i,sep="/"))
  x <- x[,c(2,5,10,11)]
	x[,2] <- log2(x[,2])
  write.table(x,file=paste(outputDir2,i,sep="/"),quote=F,sep="\t",row.names=F)
  cat(paste(i,",",nrow(x)),"\n")
}
