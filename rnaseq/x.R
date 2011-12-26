x <- read.table("x.tab")
y <- read.table("x-with-deseq.tab")
y2 <- read.table("x-with-deseq-edger.tab")
z <- y2[! y2$Package %in% x$Package,]
z.package <- paste(z$Package, z$Version, sep="_")
z.package <- paste(z.package, "tar.gz", sep=".")
z.wget <- paste("wget -N -q http://www.bioconductor.org/packages/2.9/bioc/src/contrib", z.package, sep="/")
z2.wget <- paste("wget -N -q http://cran.r-project.org/src/contrib", z.package, sep="/")
z.install <- paste("leopard-x86_64/R-2.14.0/bin/R CMD INSTALL /Users/goshng/Documents/Projects/peach/rnaseq/downloads", z.package, sep="/")
write.table(z.wget, file="download-packages.sh", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(z2.wget, file="download-packages.sh", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(z.install, file="install-packages.sh", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)


