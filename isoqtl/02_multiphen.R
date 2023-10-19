library(MultiPhen)
library(arrow)

computer <- 'pc'
out.filename <- 'isoqtl_multiphen_binary'
chr.range <- c(1:22)

if (computer=='pc') {
  data_dir <- file.path("E:/lab_data/isoqtl/for_mtclass")
  results_dir <- file.path("G:/My Drive/Lab/lab_projects/isoqtl/results")
} else if (computer=='cluster') {
  data_dir <- file.path("/projects/compbio/users/yli459/isoqtl/data")
  results_dir <- file.path("/projects/compbio/users/yli459/isoqtl/results")
}

start_time <- Sys.time()
datalist <- list()

for (c in chr.range) {
  
  gt <- read.table(gzfile(file.path(data_dir,'genotypes',sprintf('chr%s_genotypes_binary.txt.gz', c))), sep='\t', header=TRUE)
  gt = as.data.frame(gt)
  
  if ('__index_level_0__' %in% colnames(gt)) {
    gt <- subset(gt, select=-c(`__index_level_0__`))
  }
  
  colnames(gt) <- gsub('X','',colnames(gt), fixed=TRUE)
  
  exp <- arrow::read_feather(file.path(data_dir,'all_expression.ftr'))
  exp = as.data.frame(exp)

  if ('__index_level_0__' %in% colnames(exp)) {
    exp <- subset(exp, select=-c(`__index_level_0__`))
  }
    
  # only keep genes on the chromosome
  exp <- subset(exp, gene_name %in% gt$gene_name)
  exp <- subset(exp, genotypingID %in% colnames(gt))

  num_genes <- length(unique(exp$gene_name))
  num_donors <- length(unique(exp$genotypingID))
  
  print('=== Loaded data ===')
  print(sprintf('There are %s donors and %s genes on chr%s', num_donors, num_genes, c))

  genes <- as.character(unique(gt$gene_name))
  
  results_list <- list(gene=c(), variant=c(), pval=c())
  
  for (gene in genes) {
    
    print(sprintf('Now on gene %s', gene))
    
    # subset genotype and phenotype data
    geno_data_gene <- subset(gt, gene_name==gene)
    pheno_data_gene <- subset(exp, gene_name==gene)
    
    row.names(pheno_data_gene) <- pheno_data_gene$genotypingID
    pheno_data_gene <- subset(pheno_data_gene, select=-c(gene_name, gene_id, individualID, specimenID, genotypingID))
    pheno_data_gene <- pheno_data_gene[,colSums(is.na(pheno_data_gene))==0] # drop missing isoforms
    pheno_data_gene <- as.matrix(pheno_data_gene)
    
    snps <- as.character(geno_data_gene$ID)
    
    for (snp in snps) {
      
      geno_data_var <- subset(geno_data_gene, ID==snp)
      geno_data_var <- subset(geno_data_var, select=-c(gene_name, gene_id, ID))
      geno_data_var <- t(geno_data_var)
      colnames(geno_data_var) <- snp
      geno_data_var <- as.matrix(geno_data_var)
      
      Xy <- merge(pheno_data_gene, geno_data_var, by='row.names')
      Xy <- na.omit(Xy)
      
      geno_data <- Xy[,c(1,dim(Xy)[2])]
      row.names(geno_data) <- geno_data$Row.names
      geno_data <- subset(geno_data, select=-Row.names)
      
      pheno_data <- Xy[,c(1:dim(Xy)[2]-1)]
      row.names(pheno_data) <- pheno_data$Row.names
      pheno_data <- subset(pheno_data, select=-Row.names)
      
      geno_data <- as.matrix(geno_data)
      pheno_data <- as.matrix(pheno_data)
      
      skip_to_next <- FALSE
      
      results <- tryCatch(mPhen(geno_data, pheno_data), error=function(e) {skip_to_next <<- TRUE})
      if (skip_to_next) {
        results_list$gene <- append(results_list$gene, gene)
        results_list$variant <- append(results_list$variant, snp)
        results_list$pval <- append(results_list$pval, 'NA')
        next
      }
      
      results <- results$Results
      res <- results[,,,2]
      pval <- res[length(res)] # last entry is JointModel
      print(sprintf("SNP %s    p-value=%s", snp, pval))
      
      results_list$gene <- append(results_list$gene, gene)
      results_list$variant <- append(results_list$variant, variant)
      results_list$pval <- append(results_list$pval, pval)
    }
  }
  
  results_df <- data.frame(results_list)
  datalist[[c]] <- results_df
}

end_time <- Sys.time()
print(end_time - start_time)

final_df = do.call(rbind, datalist)
print(head(final_df))
arrow::write_feather(final_df, file.path(results_dir, sprintf('%s.ftr', out.filename)))
