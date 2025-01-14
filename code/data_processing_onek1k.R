# load and pre-process downloaded data
library(Seurat)
library(pbapply)
library(parallel)
library(dplyr)

num.cores <- 6

# https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1
data.path <- "/mnt/c/Users/ronni/Desktop/OneK1K/data/"
sobj <- readRDS(file.path(data.path, "Yazar2022.rds"))
message("loaded data")

# find top variable features
sobj <- FindVariableFeatures(sobj, assay='RNA', nfeatures=10000)
var.features <- sobj@assays$RNA@var.features
sobj <- sobj[var.features,]

# normalize the data
sobj@assays$RNA@counts <- sobj@assays$RNA@data
sobj <- NormalizeData(sobj, normalization.method='LogNormalize', scale.factor=10000)

DefaultAssay(sobj) <- "RNA"
Idents(sobj) <- "predicted.celltype.l2"
all.celltypes <- sort(as.character(unique(sobj@meta.data$predicted.celltype.l2)))
counts <- sobj@assays$RNA@data
meta <- sobj[[]][,c('donor_id','predicted.celltype.l2')]
rm(sobj) # free memory
gc() # free memory

FeatureList <- mclapply(c(1:length(var.features)), function(i) {
  
  gene <- var.features[i]
  
  message(sprintf("Processing counts for #%s gene %s...", i, gene))
  counts.all.barcodes <- as.data.frame(counts[row.names(counts)==gene,])
  merged <- merge(counts.all.barcodes, meta, by='row.names')
  colnames(merged) <- c("barcode","counts","donor","celltype")
  merged <- merged[,-1]
  
  grouped <- merged %>% group_by(donor, celltype) %>% summarise(sum_counts=sum(counts)) %>% as.data.frame()
  grouped <- reshape(grouped, idvar='donor', timevar='celltype', direction='wide')
  rownames(grouped) <- grouped$donor
  grouped <- grouped[,-1]
  colnames(grouped) <- gsub('sum_counts.','', gsub('\\.','_', colnames(grouped)))
  grouped <- subset(grouped, select=-Doublet)
  grouped$gene <- gene
  grouped$donor <- rownames(grouped)
  
  return(grouped)
  
}, mc.cores=num.cores)

feature.mtx <- do.call(rbind, FeatureList)
new.cols <- c(c("gene","donor"), setdiff(colnames(feature.mtx), c("gene","donor")))
feature.mtx <- feature.mtx[, new.cols]

saveRDS(feature.mtx, file.path(data.path, "01_features.rds"))
# arrow::write_feather(feature.mtx, file.path(data.path, "01_features.ftr"))
message("=== done ===")
