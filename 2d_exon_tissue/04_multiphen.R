library('MultiPhen')

for (c in c(7:9)) {
  
  data.dir <- 'E:/lab_data'
  # data.dir <- '/projects/compbio/users/yli459'
  
  results_list <- list(gene=c(), variant=c(), pval=c())
  
  # Load data
  exp <- read.table(file.path(data.dir, 'gtex_exoneqtl', '2d_expression', sprintf('chr%s_2d_expression.txt.gz', c)), sep='\t', header=TRUE)
  gt <- read.table(file.path(data.dir, 'gtex_exoneqtl', 'genotypes', sprintf('chr%s_genotypes_10kb_binary.txt.gz',c)), sep='\t', header=TRUE)
  print('===== Loaded data =====')
  
  for (gene in unlist(unique(exp$gene_id))) {
    
    # subset gene expression data
    gene.exp <- subset(exp, gene_id==gene)
    row.names(gene.exp) <- gene.exp$donor
    gene.exp <- subset(gene.exp, select=-c(donor,gene_id))
    gene.exp <- gene.exp[,colSums(is.na(gene.exp))==0]
    
    # subset genotypes
    gene.gt <- subset(gt, gene_id==gene)
    colnames(gene.gt) <- gsub('\\.', '-', colnames(gene.gt))
    
    for (snp in gene.gt$ID) {
      
      snp.gt <- subset(gene.gt, ID==snp)
      snp.gt <- snp.gt[!duplicated(snp.gt$ID),] # remove duplicate SNPs
      snp.gt <- subset(snp.gt, select=-c(gene_id, ID))
      snp.gt <- t(snp.gt)
      colnames(snp.gt) <- c('genotype')
      
      Xy <- merge(gene.exp, snp.gt, by='row.names')
      donors <- as.character(Xy$Row.names)
      
      geno.data <- as.matrix(Xy$genotype)
      colnames(geno.data) <- c(snp)
      row.names(geno.data) <- donors
      
      pheno.data <- subset(Xy, select=-c(Row.names,genotype))
      pheno.data <- as.matrix(pheno.data)
      row.names(pheno.data) <- donors
      
      error_occur <- FALSE
      results <- tryCatch(mPhen(geno.data, pheno.data), error=function(e) {error_occur <<- TRUE})
      if (error_occur) {
        
        results_list$gene <- append(results_list$gene, gene)
        results_list$variant <- append(results_list$variant, snp)
        results_list$pval <- append(results_list$pval, 'NA')
      
      } else {
        
        results <- results$Results
        num.features <- dim(pheno.data)[2]
        pval <- results[,,,2][num.features+1]
        
        results_list$gene <- append(results_list$gene, gene)
        results_list$variant <- append(results_list$variant, snp)
        results_list$pval <- append(results_list$pval, pval)
      }
    }
  }
  results.df <- data.frame(results_list)
  write.table(results.df, file=sprintf('G:/My Drive/Lab/lab_projects/gtex_exoneqtl/results/2d_exon_tissue/chr%s_2d_multiphen.txt',c), sep='\t', row.names=FALSE, quote=FALSE)
  print(sprintf('Wrote chr%s results to file',c))
}

