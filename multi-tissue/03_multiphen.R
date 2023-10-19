library(MultiPhen)
library(arrow)

computer <- 'cluster'
experiment <- 'brain'
out.filename <- 'brain_multiphen'
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
  results_dir <- file.path(sprintf("/projects/compbio/users/yli459/gtex_eqtl/%s", suffix))
}

run_mphen <- function(c) {
  
  gt <- read.table(gzfile(file.path(data_dir, 'genotypes', sprintf('chr%s_genotypes_10kb_binary.txt.gz', c))), sep='\t', header=TRUE)
  
  if (experiment == '9_tissue') {
    counts <- arrow::read_feather(file.path(data_dir, suffix, 'counts.ftr'))
  } else if (experiment == 'brain') {
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
  gt <- gt[, names(gt) %in% keep_cols]
  
  genes <- sort(unique(gt$gene_name))
  num_genes <- length(genes)
  print(sprintf('There are %s genes in chr%s in the %s case', num_genes, c, experiment))
  
  results_list <- list(gene=c(), variant=c(), pval=c())
  
  for (gene_id in genes) {
    
    print(sprintf('Now on gene %s', gene_id))
    
    # subset genotype and phenotype data
    geno_data_gene <- subset(gt, gene_name==gene_id)
    pheno_data_gene <- subset(counts, gene==gene_id)
    
    row.names(pheno_data_gene) <- pheno_data_gene$donor
    pheno_data_gene <- subset(pheno_data_gene, select=-c(gene, donor))
    pheno_data_gene <- as.matrix(pheno_data_gene)
    
    variants <- sort(geno_data_gene$ID)
    
    for (variant in variants) {
    
      # print(sprintf('Now on %s', variant))
      
      geno_data_var <- subset(geno_data_gene, ID==variant)
      geno_data_var <- subset(geno_data_var, select=-c(gene_name, ID))
      geno_data_var <- t(geno_data_var)
      colnames(geno_data_var) <- variant
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
        
        results_list$gene <- append(results_list$gene, gene_id)
        results_list$variant <- append(results_list$variant, variant)
        results_list$pval <- append(results_list$pval, "NA")
        next
      }
      
      results <- results$Results
      
      if (experiment=='9_tissue') {
        pval <- results[,,,2][10] # because JointModel is last one after 9 tissues
      } else if (experiment == 'brain') {
        pval <- results[,,,2][14] # because there are 13 tissues
      }

      results_list$gene <- append(results_list$gene, gene_id)
      results_list$variant <- append(results_list$variant, variant)
      results_list$pval <- append(results_list$pval, pval)
    }
  }

  results_df <- data.frame(results_list)
  return(results_df)
}

datalist <- list()
start_time <- Sys.time()

for (c in chr.range) {
  results <- run_mphen(c)
  datalist[[c]] <- results
}

end_time <- Sys.time()
print(end_time - start_time)

final_df = do.call(rbind, datalist)
print(head(final_df))
#write.csv(final_df, file.path(results_dir, 'brain_multiphen_additive.csv'))
arrow::write_feather(final_df, file.path(results_dir, sprintf("%s.ftr", out.filename)))
