library(MultiPhen)
library(arrow)

computer <- 'pc'

if (computer=='pc') {
  counts_dir <- file.path('E:/lab_data/gtex_exoneqtl/tissue_level_expression')
  gt_dir <- file.path('E:/lab_data/gtex_exoneqtl/genotypes')
  results_dir <- file.path('G:/My Drive/Lab/lab_projects/gtex_exoneqtl/results')
} else if (computer=='cluster') {
  counts_dir <- file.path('/projects/compbio/users/yli459/gtex_exoneqtl/tissue_level_expression')
  gt_dir <- file.path('/projects/compbio/users/yli459/gtex_exoneqtl/genotypes')
  results_dir <- file.path('/projects/compbio/users/yli459/gtex_exoneqtl/manova_multiexon')
}

tissues <- list.files(counts_dir)
brain_tissues <- tissues[grepl("Brain_", tissues, fixed=TRUE) & !grepl("_gene_list.txt", tissues)]
brain_tissues <- unlist(sapply(brain_tissues, function(x) strsplit(x, "_counts.ftr")))

for (tissue in brain_tissues) {
  
  print(sprintf("Now on tissue %s", tissue))
  
  datalist <- list()
  start_time <- Sys.time()
  chr.range <- c(1:22)
  
  for (c in chr.range) {
    
    gt <- read.table(gzfile(file.path(gt_dir, sprintf('chr%s_genotypes_10kb_binary.txt.gz', c))), sep='\t', header=TRUE)
    counts <- arrow::read_feather(file.path(counts_dir, sprintf('%s_counts.ftr', tissue)))
    
    counts <- data.frame(counts)
    gt <- data.frame(gt)
    
    if ("X__index_level_0__" %in% colnames(counts)) {
      counts <- subset(counts, select=-X__index_level_0__)
    }
    if ("X__index_level_0__" %in% colnames(gt)) {
      gt <- subset(gt, select=-X__index_level_0__)
    }
    
    # only keep genes in low-proportion zeros gene list
    gene_list <- c(read.table(file.path(counts_dir, sprintf("%s_gene_list.txt", tissue)), header=TRUE)$gene_name)
    counts <- counts[counts$gene %in% gene_list, ]
    
    # make sure genotypes and counts have the same genes
    gt <- subset(gt, gene_name %in% counts$gene)
    
    # replace '.' with '-' for column names
    colnames(gt) <- gsub(colnames(gt), pattern="\\.", replacement="-")
    
    # only keep donors in the genotype matrix
    counts <- subset(counts, donor %in% names(gt))
    num_donors <- length(unique(counts$donor))
    print(sprintf('There are %s donors for tissue %s', num_donors, tissue))
    
    keep_cols <- c('gene_name', 'ID', sort(unique(counts$donor)))
    gt <- gt[, colnames(gt) %in% keep_cols]
    print("=== Loaded data ===")
    
    genes <- unlist(unique(gt$gene_name))
    num_genes <- length(genes)
    print(sprintf('There are %s genes in chr%s', num_genes, c))
    
    results_list <- list(gene=c(), variant=c(), pval=c())
    
    for (gene_id in genes) {
      
      print(sprintf('Now on gene %s', gene_id))
      
      # subset genotype and phenotype data
      geno_data_gene <- subset(gt, gene_name==gene_id)
      pheno_data_gene <- subset(counts, gene==gene_id)
      
      # drop NAs in columns of counts data (missing exons)
      pheno_data_gene <- pheno_data_gene[ , colSums(is.na(pheno_data_gene))==0]
      
      row.names(pheno_data_gene) <- pheno_data_gene$donor
      pheno_data_gene <- subset(pheno_data_gene, select=-c(gene, donor))
      pheno_data_gene <- as.matrix(pheno_data_gene)
      
      variants <- sort(geno_data_gene$ID)
      
      for (variant in variants) {
        
        geno_data_var <- subset(geno_data_gene, ID==variant)
        geno_data_var <- geno_data_var[!duplicated(subset(geno_data_var, select=c(gene_name, ID))),]
        geno_data_var <- subset(geno_data_var, select=-c(gene_name, ID))
        geno_data_var <- t(geno_data_var)
        colnames(geno_data_var) <- variant
        geno_data_var <- as.matrix(geno_data_var)
        
        Xy <- merge(pheno_data_gene, geno_data_var, by='row.names')
        
        geno_data <- Xy[,c(1,dim(Xy)[2])]
        row.names(geno_data) <- geno_data$Row.names
        geno_data <- subset(geno_data, select=-Row.names)
        
        pheno_data <- Xy[,c(1:dim(Xy)[2]-1)]
        row.names(pheno_data) <- pheno_data$Row.names
        pheno_data <- subset(pheno_data, select=-Row.names)
        
        geno_data <- as.matrix(geno_data)
        pheno_data <- as.matrix(pheno_data)
        
        skip_to_next <- FALSE
        
        res <- tryCatch(mPhen(geno_data, pheno_data), error=function(e) {skip_to_next <<- TRUE})
        result <- tryCatch(res$Results[,,,2], error=function(e) {skip_to_next <<- TRUE})
        pval <- tryCatch(as.numeric(result[length(result)]), error=function(e) {skip_to_next <<- TRUE})
        
        if (skip_to_next) {
          results_list$gene <- append(results_list$gene, gene_id)
          results_list$variant <- append(results_list$variant, variant)
          results_list$pval <- append(results_list$pval, 'NA')
          next
        }
        
        print(sprintf("SNP %s has p = %s", variant, pval))
        
        results_list$gene <- append(results_list$gene, gene_id)
        results_list$variant <- append(results_list$variant, variant)
        if (is.na(pval) || length(pval)==0) {
          results_list$pval <- append(results_list$pval, 'NA')
        } else {
          results_list$pval <- append(results_list$pval, pval)
        }
      }
    }
    
    results_df <- data.frame(results_list)
    datalist[[c]] <- results_df
  }
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  final_df <- do.call(rbind, datalist)
  final_df <- final_df[!duplicated(subset(final_df, select=c(gene,variant))),]
  print(head(final_df))
  
  final_path <- file.path(results_dir, sprintf('%s', tissue))
  if (!file.exists(final_path)) {
    dir.create(file.path(final_path))
  }
  
  arrow::write_feather(final_df, file.path(final_path, "multiphen.ftr"))
  
}