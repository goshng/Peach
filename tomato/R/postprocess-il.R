args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{                                                                               
  cat ("Rscript R/correlation-metabolite.R raw IL_mat_30.txt\n")
  quit("no")
}
# xil <- read.csv("output/cor/IL12-4-1-pval.csv")
xil <- read.csv(args[1])

met30 <- read.csv("output/metabolite/IL_mat_30.txt-pval.csv")
met41 <- read.csv("output/metabolite/IL_mat_41.txt-pval.csv")
met62 <- read.csv("output/metabolite/IL_mat_62.txt-pval.csv")
total <- read.csv("output/gene2gene-pval.csv")

# Append total.n2
# Not all eQTL may not be in the total set.
# We add a row to the total for using NA.
# stopifnot(sum(xil[,1] %in% total[,1]) == nrow(xil))
total <- rbind(total, total[nrow(total),])
total[nrow(total),c("n2","n3","n4")] <- NA
m <- match(xil[,1],total[,1])
m[is.na(m)] <- nrow(total)
yil <- xil[,c("X","n1","n2","n3","n4")]
yil <- cbind(yil, total.n2=total[m,"n2"])

# Append met30.n2
met30 <- rbind(met30, met30[nrow(met30),])
met30[nrow(met30),c("n1","n2","n3","n4")] <- NA
m <- match(xil[,1],met30[,1])
m[is.na(m)] <- nrow(met30)
yil <- cbind(yil, met30.n2=met30[m,"n2"])

# Append met41.n2
met41 <- rbind(met41, met41[nrow(met41),])
met41[nrow(met41),c("n1","n2","n3","n4")] <- NA
m <- match(xil[,1],met41[,1])
m[is.na(m)] <- nrow(met41)
yil <- cbind(yil, met41.n2=met41[m,"n2"])

# Append met62.n2
met62 <- rbind(met62, met62[nrow(met62),])
met62[nrow(met62),c("n1","n2","n3","n4")] <- NA
m <- match(xil[,1],met62[,1])
m[is.na(m)] <- nrow(met62)
yil <- cbind(yil, met62.n2=met62[m,"n2"])

# Append annotation
yil <- cbind(yil, note=xil[,"note"])

# Append other p-values
yil <- cbind(yil, xil[,grep("Sol",colnames(xil))])

# Print the result in a csv file.
write.csv(yil,file=paste(sub(".csv","",args[1]),"-processed.csv",sep=""))

