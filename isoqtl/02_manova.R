# MANOVA script on isoQTL data

computer <- 'pc'
out.filename <- 'isoqtl_manova'
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
      
      row.names(Xy) <- Xy$Row.names
      Xy <- subset(Xy, select=-Row.names)
      
      independent_var <- Xy[,snp]
      dependent_vars <- as.matrix(Xy[,c(1:dim(Xy)[2]-1)])
      
      skip_to_next <- FALSE
      
      results <- tryCatch(manova(dependent_vars ~ independent_var, data=Xy), error=function(e) {skip_to_next <<- TRUE})
      pval <- tryCatch(summary(results, tol=0)$stats[1,6], error=function(e) {skip_to_next <<- TRUE})
      
      if (skip_to_next) {
        results_list$gene <- append(results_list$gene, gene)
        results_list$variant <- append(results_list$variant, snp)
        results_list$pval <- append(results_list$pval, 'NA')
        next
      }
  
      if (is.na(pval) || length(pval)==0) {
        results_list$gene <- append(results_list$gene, gene)
        results_list$variant <- append(results_list$variant, snp)
        results_list$pval <- append(results_list$pval, 'NA')
        next
      }
      
      print(sprintf("SNP %s    p-value=%s", variant, pval))
      results_list$gene <- append(results_list$gene, gene)
      results_list$variant <- append(results_list$variant, snp)
      results_list$pval <- append(results_list$pval, pval)
    }
  }
  
  results_df <- data.frame(results_list)
  return(results_df)
}

end_time <- Sys.time()
print(end_time - start_time)

final_df = do.call(rbind, datalist)
print(head(final_df))
arrow::write_feather(final_df, file.path(results_dir, sprintf('%s.ftr', out.filename)))