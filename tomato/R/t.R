metabolite.files <- c("IL_mat_30.txt", "IL_mat_41.txt", "IL_mat_62.txt")
metabolite.file <- "IL_mat_30.txt"
m <- read.table(paste("data",raw,metabolite.file,sep="/"),head=T,sep="\t")

#
m <- read.table("IL_mat_30.txt",head=T)

#
m <- read.table("IL_mat_41.txt",head=T,sep="\t")

#
m <- read.table("IL_mat_62.txt",head=T)
m <- t(m)
