library(MultiPhen)
library(arrow)
library(tidyverse)

experiment <- "multi-exon"
tissue <- "Whole_Blood"

data.dir <- file.path(sprintf("E:/lab_data/mtclass_eqtl/%s", experiment))
results.dir <- file.path("C:/Users/ronni/Desktop")

# data.dir <- file.path(sprintf("/projects/compbio/users/yli459/mtclass_eqtl/%s", experiment))
# results.dir <- file.path("/home/yli459/")

exp <- read_feather(file.path(data.dir, sprintf("tissue_level_expression/%s.ftr", tissue))) %>% as.data.frame()

Run.mPhen <- function(c) {
  
  chrom.res.list <- list()
  gt <- read.csv(gzfile(file.path(data.dir, 'genotypes', sprintf('chr%s_genotypes_10kb_binary.txt.gz', c))), sep='\t', header=TRUE) %>% as.data.frame()
  colnames(gt) <- gsub(colnames(gt), pattern="\\.", replacement="-")
  
  # find common genes and donors
  genes <- intersect(gt$gene_name, exp$gene_name)
  donors <- intersect(colnames(gt)[4:ncol(gt)], exp$donor)
  
  X <- exp[exp$donor %in% donors & exp$gene_name %in% genes,]
  Y <- gt[,c("gene_name","ID",donors)]
  message(sprintf("There are %s genes and %s donors in chr%s", length(genes), length(donors), c))
  
  for (gene in genes) {
    
    # subset genotype and phenotype data
    geno_data_gene <- subset(Y, gene_name==gene) %>% subset(., select=-gene_name)
    pheno_data_gene <- subset(X, gene_name==gene) %>% subset(., select=-gene_name)
    row.names(pheno_data_gene) <- pheno_data_gene$donor
    pheno_data_gene <- subset(pheno_data_gene, select=-donor)
    pheno_data_gene <- pheno_data_gene[,colSums(is.na(pheno_data_gene)) != nrow(pheno_data_gene)]
    snps <- sort(geno_data_gene$ID)
    
    gene.res.list <- lapply(c(1:length(snps)), FUN=function(i) {
      res <- list()
      snp <- snps[i]
      message(sprintf("%s-%s",gene,snp))
      geno_data_var <- subset(geno_data_gene, ID==snp)
      geno_data_var <- geno_data_var[!duplicated(geno_data_var$ID),]
      geno_data_var <- subset(geno_data_var, select=-ID) %>% t()
      colnames(geno_data_var) <- "genotype"
      Xy <- merge(pheno_data_gene, geno_data_var, by='row.names') %>% column_to_rownames("Row.names")
      Xy <- Xy[!is.na(Xy$genotype),]
      
      geno.data <- Xy[,ncol(Xy)] %>% as.matrix()
      pheno.data <- Xy[,c(1:ncol(Xy)-1)] %>% as.matrix()
      row.names(geno.data) <- row.names(pheno.data)
      num.exons <- ncol(pheno.data)
      
      skip_to_next <- FALSE
      results <- tryCatch(mPhen(geno.data, pheno.data), error=function(e) {skip_to_next <<- TRUE})
      if (skip_to_next) {
        res[["gene"]] <- gene
        res[["variant"]] <- snp
        res[["pval"]] <- "NA"
      } else {
        results <- results$Results
        pval <- results[,,,2][num.exons+1]
        res[["gene"]] <- gene
        res[["variant"]] <- snp
        res[["pval"]] <- pval
      }
      return(data.frame(res))
    })
    gene.res <- do.call(rbind, gene.res.list)
    gene.res$pair <- paste0(gene.res$gene,";",gene.res$variant)
    print(head(gene.res))
    chrom.res.list[[gene]] <- gene.res
    gc()
  }
  chrom.res <- do.call(rbind, chrom.res.list)
  return(chrom.res)
}

#### Run for each chromosome
datalist <- list()
start.time <- Sys.time()
for (c in c(1:22)) {
  chrom.res <- Run.mPhen(c)
  datalist[[c]] <- results
}
end.time <- Sys.time()
print(end.time-start.time)

final = do.call(rbind, datalist)
print(head(final))

# write to file
write.csv(final, file.path(results.dir, sprintf("%s_multiphen.txt", tissue)), quote=F, sep='\t', row.names=F, col.names=T)
message("=== done ===")
