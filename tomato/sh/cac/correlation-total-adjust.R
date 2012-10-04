load(file="gene2gene-cor.RData")
                                                                                
# adjust pvalue                                                                 
P <- m.total.pval                                                               
rm(m.total.pval)

p <- P[lower.tri(P)]                                                            
adj.p <- p.adjust(p, method = "BH")                                             
P[lower.tri(P)] <- adj.p                                                        
P[upper.tri(P)] <- 0                                                            
P <- P + t(P)                                                                   
P <- ifelse(P < 1e-06, 0, P)                                                    
P <- format(round(P, 6))                                                        
P[grep("NA", P)] <- ""                                                          

save(P,file="gene2gene-cor-adjust.RData")
q("no")
