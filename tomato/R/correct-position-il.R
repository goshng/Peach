inputDir <- "input/raw/ilposition-revised"
outputDir <- "output/ilposition"

for (i in list.files(inputDir))
{
  x <- scan(paste(inputDir,i,sep="/"))
  y <- c(x[1], x[2], trunc(x[1] + x[2])/2)
  cat(y,"\n",file=paste(outputDir,i,sep="/"))
}
