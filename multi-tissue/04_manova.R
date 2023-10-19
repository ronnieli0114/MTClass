library(arrow)

computer <- 'pc'
experiment <- '9_tissue'
out.filename <- '9tiss_manova'
chr.range <- c(1:22)

if (experiment == '9_tissue') {
  suffix <- '9_tissue'
} else if (experiment == 'brain') {
  suffix <- 'brain_tissue'
}

if (computer=='pc') {
  data_dir <- file.path("E:/lab_data/gtex_eqtl")
  results_dir <- file.path(sprintf("G:/My Drive/Lab/lab_projects/gtex_eqtl/results/%s", suffix))
} else if (computer=='cluster') {
  data_dir <- file.path("/projects/compbio/users/yli459/gtex_eqtl")
  results_dir <- file.path("/projects/compbio/users/yli459/gtex_eqtl", suffix)
}

datalist <- list()
start_time <- Sys.time()

for (c in chr.range) {
  
  gt <- read.table(gzfile(file.path(data_dir, 'genotypes', sprintf('chr%s_genotypes_10kb_binary.txt.gz', c))), sep='\t', header=TRUE)
  
  if (experiment == '9_tissue') {
    counts <- arrow::read_feather(file.path(data_dir, suffix, 'counts.ftr'))
  } else {
    counts <- arrow::read_feather(file.path(data_dir, suffix, 'PMM_imputed_counts.ftr'))
  }
  
  counts <- data.frame(counts)
  gt <- data.frame(gt)
  print('===== Loaded all data =====')
  
  if ('X__index_level_0__' %in% colnames(counts)) {
    counts <- subset(counts, select=-X__index_level_0__)
  }
  
  if ('X__index_level_0__' %in% colnames(gt)) {
    gt <- subset(gt, select=-X__index_level_0__)
  }
  
  # replace '.' with '-' for column names
  names(gt) <- gsub(names(gt), pattern="\\.", replacement="-")
  
  # only keep genes in top counts
  gt <- subset(gt, gene_name %in% counts$gene)
  
  # only keep donors in the genotype matrix
  counts <- subset(counts, donor %in% names(gt))
  num_donors <- length(unique(counts$donor))
  print(sprintf('There are %s donors in the %s case', num_donors, experiment))
  
  keep_cols <- c('gene_name', 'ID', sort(unique(counts$donor)))
  gt <- gt[, colnames(gt) %in% keep_cols]
  
  genes <- sort(unique(gt$gene_name))
  num_genes <- length(genes)
  print(sprintf('There are %s genes in chr%s in the %s case', num_genes, c, experiment))
  
  results_list <- list(gene=c(), variant=c(), pval=c())
  
  for (gene_id in genes) {
    
    # subset genotype and phenotype data
    geno_data_gene <- subset(gt, gene_name==gene_id)
    pheno_data_gene <- subset(counts, gene==gene_id)
    
    row.names(pheno_data_gene) <- pheno_data_gene$donor
    pheno_data_gene <- subset(pheno_data_gene, select=-c(gene, donor))
    pheno_data_gene <- as.matrix(pheno_data_gene)
    
    variants <- sort(geno_data_gene$ID)
    
    for (variant in variants) {
      
      geno_data_var <- subset(geno_data_gene, ID==variant)
      geno_data_var <- subset(geno_data_var, select=-c(gene_name, ID))
      geno_data_var <- t(geno_data_var)
      colnames(geno_data_var) <- variant
      geno_data_var <- as.matrix(geno_data_var)
      
      Xy <- merge(pheno_data_gene, geno_data_var, by='row.names')
      Xy <- na.omit(Xy)
      
      row.names(Xy) <- Xy$Row.names
      Xy <- subset(Xy, select=-Row.names)
      
      dep_vars <- as.matrix(Xy[,1:dim(Xy)[2]-1])
      indep_vars <- Xy[,dim(Xy)[2]]

      skip_to_next <- FALSE
      
      res <- tryCatch(manova(dep_vars ~ indep_vars, na.action=na.exclude), error=function(e) {skip_to_next <<-TRUE})
      pval <- tryCatch(summary(res)$stats[1,6], error=function(e) {skip_to_next <<- TRUE})
      if (skip_to_next) {
        results_list$gene <- append(results_list$gene, gene_id)
        results_list$variant <- append(results_list$variant, variant)
        results_list$pval <- append(results_list$pval, 'NA')
        next
      }
      
      results_list$gene <- append(results_list$gene, gene_id)
      results_list$variant <- append(results_list$variant, variant)
      
      if (is.na(pval) || length(pval)==0) {
        results_list$pval <- append(results_list$pval, 'NA')
      } else {
        results_list$pval <- append(results_list$pval, pval)
      
      print(sprintf("Done with %s, %s", gene_id,))
      }
    }
  }
  results_df <- data.frame(results_list)
  datalist[[c]] <- results_df
  # write.csv(results_df, file.path(results_dir, sprintf("chr%s_manova.csv", c)), row.names=FALSE)
  # print(sprintf("Wrote results for chr%s to CSV", c))
}

end_time <- Sys.time()
print(end_time - start_time)

final_df = do.call(rbind, datalist)
print(head(final_df))
arrow::write_feather(final_df, file.path(results_dir, sprintf("%s.ftr", out.filename)))
