#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) != 2)
#{                                                                               
# cat ("Rscript R/correlation-total.R 100 1\n")
#	quit("no")
#}
#total.segment <- as.integer(args[1])
total.segment <- 800
#index.segment <- as.integer(args[2])

################################################################################
#raw <- "raw"
#x <- read.csv(paste("data",raw,"IL-expressionforcorrelation.csv",sep="/"))
#
#z <- x
#z2 <- apply(z[,-1:-3],c(1,2),function(x) x > 1)
#gname <- z[apply(z2,1,sum) >= 5,1]
#m82 <- z[apply(z2,1,sum) >= 5,2]
#z3 <- z[apply(z2,1,sum) >= 5,-1:-3]
#
#gname <- gname[m82>0]
#z3 <- z3[m82>0,]
#m82 <- m82[m82>0]
#
#z4 <- data.matrix(log2(z3/m82))
#z5 <- apply(z4,1:2,function(x) if (is.infinite(x)) NA else x)
#z5 <- t(z5)
################################################################################
load("cor.RData")
m.total.pval <- c()
m.total.cor <- c()
for (index.segment in seq(total.segment)) {
  load(paste("output/cor-",index.segment,".RData",sep="")) 
  m.total.pval <- cbind(m.total.pval,m.pval)  
  m.total.cor <- cbind(m.total.cor,m.cor)  
}
rownames(m.total.pval) <- gname
colnames(m.total.pval) <- gname
rownames(m.total.cor) <- gname
colnames(m.total.cor) <- gname

stopifnot(nrow(m.total.pval)==ncol(z5), 
          ncol(m.total.pval)==ncol(z5),
          nrow(m.total.cor)==ncol(z5),
          ncol(m.total.cor)==ncol(z5))

# adjust pvalue 
P <- m.total.pval
p <- P[lower.tri(P)]
adj.p <- p.adjust(p, method = "BH")
P[lower.tri(P)] <- adj.p
P[upper.tri(P)] <- 0
P <- P + t(P)
P <- ifelse(P < 1e-06, 0, P)
P <- format(round(P, 6))
P[grep("NA", P)] <- ""
m.total.pval <- P
rm(P)

# Network file by p-value
m.total.t <- which(m.total.pval < 0.01,arr.in=TRUE)
write.table( 
  data.frame(a=rownames(m.total.t), 
             c=rep("cor",nrow(m.total.t)),
             b=colnames(m.total.pval)[m.total.t[,2]]),
  file="gene_to_gene-pval.sif",
  row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t"
)

write.csv(m.total.cor,file="gene2gene-cor.csv")
write.csv(m.total.pval,file="gene2gene-pval.csv")


