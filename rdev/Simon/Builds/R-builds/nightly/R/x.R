args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
  cat ("example: Rscript R/x.R snowleopard-i386\n")
  quit("yes")
}
buildDir = args[1]

x <- installed.packages()
tabfile <- paste(buildDir, "x-with-edger-deseq.tab", sep="/")
write.table(x, file=tabfile, sep="\t")
