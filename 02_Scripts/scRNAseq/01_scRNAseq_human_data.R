library(Seurat)
library(tidyverse)
library(SCpubr)

here::here()
 
obj <- readRDS("pk_all.rds")
cells_to_remove <- rownames(obj@meta.data[obj@meta.data$Cell_type %in% "Endocrine cell",])
obj <- obj[, !colnames(obj) %in% cells_to_remove]

# # store mitochondrial percentage in object meta data
# obj  <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
# 
# obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# 
# # run sctransform
# obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# obj  <- NormalizeData(obj)
# obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(obj)
# obj <- ScaleData(obj, features = all.genes)
# 
# # These are now standard steps in the Seurat workflow for visualization and clustering
# obj <- RunPCA(obj, features = VariableFeatures(object = obj))
# obj <- RunUMAP(obj, dims = 1:30)
# 
# obj <- FindNeighbors(obj, dims = 1:30)
# obj <- FindClusters(obj, resolution = 0.5)

DimPlot(obj, label=TRUE)
FeaturePlot(obj, "SMPD1", raster=FALSE)

library(msigdbr)

kegg_list <- msigdbr(species = "Homo sapiens") %>%
  distinct(gs_name, gene_symbol) %>%
  nest(genes = c(gene_symbol)) %>%
  mutate(genes = purrr::map(genes, purrr::compose(as_vector, unname))) %>%
  deframe()

kegg_ok <- kegg_list[names(kegg_list) %in% grep(paste(c("_SPHINGOMYELIN", "_CERAMIDE"), collapse = "|"), names(kegg_list), value = TRUE)]

for (i in names(kegg_ok)) {
  print(i)
  
  #if (i %in% pathway_to_include){
  genes <- list(kegg_ok[[i]])
  
  obj <- AddModuleScore(
    object = obj,
    features = genes,
    name = i
  )
  
  p <- FeaturePlot(obj, features = paste0(i, "1"), keep.scale = "all")
  
  print(p)
}

p <- SCpubr::do_DotPlot(sample = obj, 
                        group.by = "Cell_type",
                        features = paste0(names(kegg_ok), "1"),
                        sequential.palette = "RdBu",
                        sequential.direction = -1,
                        scale= TRUE,
                        scale.by = "size",
                        flip=TRUE, 
                        cluster = TRUE, 
                        dot.scale = 8) 

print(p)

ggsave(paste0("module_score_dot_plot.pdf"),
       p,
       width = 13,
       height = 6,
       dpi=300)
