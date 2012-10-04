outputDir <- "output/position"
outputDir3 <- "output/position3"
pen <- "output/position/PENeQTL.csv"
y <- read.csv(pen)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/find-position-mrna3.R output/ilmap2.dat\n")
	quit("no")
}

for (i in list.files(outputDir,pattern="IL")) 
{
  # Find IL name and its location
  m <- regexec("(IL[0-9\\-]+)\\.csv",i)
  ilname <- regmatches(i,m)[[1]][2]
  z <- read.csv(args[1],head=T)
  z <- z[z$il==ilname & z$type=="main",]

  x <- read.csv(paste(outputDir,i,sep="/"))
  # Get the chromosome ID from the csv file name.
  chrN <- sub("IL","",i)
  chrN <- sub("-.*","",chrN)
  chrN <- as.numeric(chrN)
  genename.pattern <- sprintf("Solyc%02d",chrN)

  # TRUE for withinIL. 
  within.il <- rep(FALSE,nrow(x))
  within.il[z$start <= x$start & x$start <= z$end] <- TRUE
  within.il[-grep(genename.pattern,x$Gene)] <- FALSE

  # TRUE for "in the same chromosome"
  same.chr <- rep(FALSE,nrow(x))
  same.chr[grep(genename.pattern,x$Gene)] <- TRUE

  # TRUE for pen: x is for IL, and y is for PEN.
  # inpen.gre is the T/F for IL's existence in PEN's.
  inpen.gre <- x$Gene %in% y$Gene
  pen.gre <- inpen.gre
  x1 <- x[inpen.gre,] 
  y1 <- y[y$Gene %in% x$Gene,] 
  x2 <- merge(x1, y1, by.x="Gene", by.y="Gene", sort = FALSE)
  # pen.gre is T/F for IL's existence in PEN's, 
  # and additionally T/F for the same pattern of expression.
  pen.gre[inpen.gre] <- (log(x2$ratio.x) * log(x2$ratio.y) > 0)

  # Log Ratio of PEN
  ratioPen <- rep(NA,nrow(x))
  ratioPen[inpen.gre] <- log(x2$ratio.y)
  x$ratio <- log2(x$ratio)
  
  x <- data.frame(x,withinil=within.il,samechr=same.chr,inpen=inpen.gre,ratiopen=ratioPen,pen=pen.gre)
  write.table(x,file=paste(outputDir3,i,sep="/"),quote=F,sep="\t",row.names=F)
}

