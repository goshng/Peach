x <- read.table("output/count/ven-il-snp.txt")
y <- read.table("output/count/ven-pen-snp.txt")
x.unique <- unique(apply(x,1,paste,collapse=""))
y.unique <- unique(apply(y,1,paste,collapse=""))
f <- "output/summary/ven-snp.txt"
cat(paste("IL unique SNPs", length(x.unique), "\n"),file=f)
cat(paste("PEN unique SNPs", length(y.unique), "\n"),file=f,append=TRUE)
cat(paste("Intersection unique SNPs", sum(x.unique %in% y.unique), "\n"),file=f,append=TRUE)
