args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
  cat ("example: Rscript R/x.R snowleopard-i386\n")
  quit("no")
}
buildDir = args[1]

x <- installed.packages()
tabfile <- paste(buildDir, "x-with-edger-deseq.tab", sep="/")
notyetfile <- paste(buildDir, "not-yet-installd.txt", sep="/")
write.table(x, file=tabfile, sep="\t")

y <- read.table(tabfile)
y.full <- read.table("R/x-with-deseq-edger.tab")
z <- y.full[! y.full$Package %in% y$Package,]
z.package <- paste(z$Package, z$Version, sep="_")
if ( length(z.package) > 0 )
{
  z.package <- paste(z.package, "tar.gz", sep=".")
  write.table(z.package, file=notyetfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
} else {
  file.create(notyetfile, showWarnings = TRUE)
}
quit("no")

