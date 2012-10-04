args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
{                                                                               
 cat ("Rscript R/correlation-total.R 100 1\n")
	quit("no")
}
total.segment <- as.integer(args[1])
index.segment <- as.integer(args[2])

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

# start.index <- (index.segment - 1) * as.integer(ncol(z5)/total.segment) + 1
# end.index <- index.segment * as.integer(ncol(z5)/total.segment)
start.index <- 17600 + (index.segment - 801) * 53 + 1
end.index <- 17600 + (index.segment - 801 + 1) * 53
if (end.index > ncol(z5)) {
  end.index <- ncol(z5)
}

m.pval <- c()
m.cor <- c()
m.gname <- gname[start.index:end.index]
m.index <- start.index:end.index
for (i in start.index:end.index) {
  x <- c()
  y <- c()
  for (j in seq(ncol(z5))) {
    a <- na.omit(cbind(z5[,i],z5[,j]))
    if (nrow(a) > 5)
    {
      b <- cor.test(a[,1],a[,2])$p.value
      x <- c(x,b)
      b <- cor.test(a[,1],a[,2])$estimate
      y <- c(y,b)
    }
    else
    {
      x <- c(x,NA)
      y <- c(y,NA)
    }
  }
  m.pval <- cbind(m.pval, x)
  m.cor <- cbind(m.cor, y)
}
save.image(paste("output/cor-",index.segment,".RData",sep="")) 

# x5 <-rcorr.adjust(z5)
# stem(j2)
