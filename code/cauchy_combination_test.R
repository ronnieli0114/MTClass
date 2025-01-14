library(ACAT)
library(tidyverse)
library(parallel)

cores <- 10
df <- read.csv(gzfile("/mnt/c/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/brain_tissue/brain_tensorqtl.txt.gz"), sep='\t', header=TRUE)
colnames(df) <- gsub("\\.","-",colnames(df))
print('loaded data')

pval.list <- mclapply(c(1:nrow(df)), mc.cores=cores, FUN=function(i) {
  if (i%%10000==0) {print(sprintf('row %s', i))} 
  pvals <- df[i,4:ncol(df)] %>% as.numeric() %>% na.omit()
  pval.cct <- ACAT(pvals) 
  return(pval.cct)
  }) %>% unlist() %>% as.numeric()

df$pval_cauchy <- pval.list
final <- df[,c('gene','variant','pair','pval_cauchy')]
write.table(final, "/mnt/c/Users/ronni/Desktop/brain_cauchy.txt", sep='\t', quote=F, row.names=F, col.names=T)
print('all done!')