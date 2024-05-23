dir= "~/Desktop/DF/DFCI_Paweletz/2024_Daiichi_DXD/"
obj.srt = readRDS(paste0(dir,"rds/P30348.24.05.08.rds"))

obj.srt@meta.data[1:3,]

obj.srt@meta.data$orig.ident %>% table()


obj.srt = subset(obj.srt, orig.ident %in% c("P30348_D", "P30348_A"))
obj.srt %>% saveRDS(paste0(dir,"rds/P30348_DA.24.05.22.rds"))


## perform default analysis
perform_default_analysis <- function(obj.srt, n_features = 2000, n_pcs = 20, 
                                     dims_for_neighbors = 1:20, 
                                     resolutions = c(0.2,0.4), 
                                     umap_dims = 1:20) {
  # Step 1: Find variable features
  obj.srt <- FindVariableFeatures(obj.srt, 
                                  selection.method = 'vst', 
                                  nfeatures = n_features)
  
  # Step 2: Scale and normalize data
  all_genes <- rownames(obj.srt)
  obj.srt <- NormalizeData(obj.srt)
  obj.srt <- ScaleData(obj.srt, features = all_genes)
  
  # Step 3: Run PCA
  obj.srt <- RunPCA(obj.srt, 
                    features = VariableFeatures(object = obj.srt), npcs = n_pcs)
  
  # Step 4: Find neighbors
  obj.srt <- FindNeighbors(obj.srt, dims = dims_for_neighbors)
  
  # Step 5: Find clusters
  obj.srt <- FindClusters(obj.srt, resolution = resolutions)
  
  # Step 6: Run UMAP
  obj.srt <- RunUMAP(obj.srt, dims = umap_dims)
  
  # Return the Seurat object with analysis results
  return(obj.srt)
}

# apply
obj.srt <- perform_default_analysis(obj.srt, umap_dims = 1:10)


DimPlot(obj.srt)
DimPlot(obj.srt, group.by = "orig.ident")
DimPlot(obj.srt, group.by = "RNA_snn_res.0.2")


obj.srt %>% saveRDS("~/Desktop/P30348.sample.rds")
