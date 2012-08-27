outputDir <- "output/position"
outputDir3 <- "output/position3"
pen <- "output/position/PENeQTL.csv"
y <- read.csv(pen)

ilpositionDir <- "data/raw/ilposition"
for (i in list.files(outputDir,pattern="IL")) 
{
  x <- read.csv(paste(outputDir,i,sep="/"))
  # Get the chromosome ID from the csv file name.
  chrN <- sub("IL","",i)
  chrN <- sub("-.*","",chrN)
  chrN <- as.numeric(chrN)
  genename.pattern <- sprintf("Solyc%02d",chrN)

  # TRUE for withinIL. 
  within.il <- rep(FALSE,nrow(x))
  x.ilposition <- scan(paste(ilpositionDir,i,sep="/"))
  within.il[x.ilposition[1] <= x$start & x$start <= x.ilposition[2]] <- TRUE

  # TRUE for cis.
  same.chr <- rep(FALSE,nrow(x))
  same.chr[grep(genename.pattern,x$Gene)] <- TRUE

  # TRUE for pen
  ratioPen <- rep(NA,nrow(x))
  inpen.gre <- x$Gene %in% y$Gene
  pen.gre <- inpen.gre
  x1 <- x[inpen.gre,] 
  y1 <- y[y$Gene %in% x$Gene,] 

  x2 <- merge(x1, y1, by.x="Gene", by.y="Gene", sort = FALSE)
  pen.gre[inpen.gre] <- (log(x2$ratio.x) * log(x2$ratio.y) > 0)
  ratioPen[inpen.gre] <- x2$ratio.y
  
  x <- data.frame(x,withinil=within.il,samechr=same.chr,inpen=inpen.gre,ratiopen=ratioPen,pen=pen.gre)
  write.table(x,file=paste(outputDir3,i,sep="/"),quote=F,sep="\t",row.names=F)
}

