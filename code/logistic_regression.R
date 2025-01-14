### runs multivariate logistic regression model
### on MTClass data, similar to SCOPA/META-SCOPA

library(arrow)
library(tidyverse)
library(parallel)

data.dir <- file.path("/mnt/c/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/multi-tissue")
results.dir <- file.path("/mnt/c/Users/ronni/OneDrive/Documents/")
cores <- 12

exp <- read_feather(file.path(data.dir, "expression/brain_tissue_TPM_PMMimpute.ftr")) %>% as.data.frame()
if ("__index_level_0__" %in% colnames(exp)) { exp <- subset(exp, select=-`__index_level_0__`) }

RunLogReg <- function(c) {
  
  chrom.res.list <- list()
  gt <- read.csv(gzfile(file.path(data.dir, sprintf('genotypes/chr%s_genotypes_10kb_binary.txt.gz', c))), sep='\t', header=TRUE) %>% as.data.frame()
  colnames(gt) <- gsub(colnames(gt), pattern="\\.", replacement="-")
  
  # find common genes and donors
  genes <- intersect(gt$gene_name, exp$gene)
  donors <- intersect(colnames(gt)[4:ncol(gt)], exp$donor)
  
  X <- exp[exp$donor %in% donors & exp$gene %in% genes,]
  Y <- gt[,c("gene_name","ID",donors)]
  print(sprintf("There are %s genes and %s donors in chr%s", length(genes), length(donors), c))
  
  for (g in genes) {
    
    # subset genotype and phenotype data
    geno_data_gene <- subset(Y, gene_name==g) %>% subset(., select=-gene_name)
    pheno_data_gene <- subset(X, gene==g) %>% subset(., select=-gene)
    rownames(pheno_data_gene) <- pheno_data_gene$donor
    pheno_data_gene <- pheno_data_gene[,-1]
    pheno_data_gene <- pheno_data_gene[,colSums(is.na(pheno_data_gene)) != nrow(pheno_data_gene)]
    snps <- sort(geno_data_gene$ID)
    
    gene.res.list <- mclapply(c(1:length(snps)), mc.cores=cores, FUN=function(i) {
      res <- list()
      v <- snps[i]
      print(sprintf("%s;%s",g,v))
      geno_data_var <- subset(geno_data_gene, ID==v)
      geno_data_var <- geno_data_var[!duplicated(geno_data_var$ID),]
      geno_data_var <- subset(geno_data_var, select=-ID) %>% t() %>% as.data.frame()
      colnames(geno_data_var) <- "genotype"
      Xy <- merge(pheno_data_gene, geno_data_var, by='row.names') %>% column_to_rownames("Row.names")
      Xy <- Xy[!is.na(Xy$genotype),]
      Xy <- Xy[Xy$genotype != 9,]
      
      skip_to_next <- FALSE
      glm1 <- tryCatch(glm(genotype ~ ., family='binomial', data=Xy), error=function(e) {skip_to_next <<- TRUE})
      glm2 <- tryCatch(glm(genotype ~ 1., family='binomial', data=Xy), error=function(e) {skip_to_next <<- TRUE})
      if (skip_to_next) {
        res[["gene"]] <- g
        res[["variant"]] <- v
        res[["pval"]] <- "NA"
      } else {
        anova.res <- anova(glm2, glm1)
        pval <- anova.res[2,5]
        res[["gene"]] <- g
        res[["variant"]] <- v
        res[["pval"]] <- pval
      }
      return(data.frame(res))
    })
    gene.res <- do.call(rbind, gene.res.list)
    gene.res$pair <- paste0(gene.res$gene,";",gene.res$variant)
    chrom.res.list[[g]] <- gene.res
    gc()
  }
  chrom.res <- do.call(rbind, chrom.res.list)
  return(chrom.res)
}

#### Run for each chromosome
datalist <- list()
start.time <- Sys.time()
for (c in c(1:22)) {
  chrom.res <- RunLogReg(c)
  datalist[[c]] <- chrom.res
}
end.time <- Sys.time()
print(end.time-start.time)

final = do.call(rbind, datalist)
print(head(final))

# write to file
write.csv(final, file.path(results.dir, "brain_LogReg.txt"), quote=F, sep='\t', row.names=F, col.names=T)
message("=== done ===")