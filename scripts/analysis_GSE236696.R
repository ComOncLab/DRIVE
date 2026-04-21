library(dplyr)
library(Seurat)
all_samples <- list.files("~/data/GSE236696/",full.names = T)
# for (i in 1:length(all_samples)){
#   tt <- list.files(all_samples[i],full.names = T)
#   for (j in 1:length(tt)){
#     file.rename(tt[j],
#                 paste0(gsub("/[^/]*$","",tt[j]),"/",gsub("^[^.]+\\.", "", basename(tt[j]))))
#   }
# }
obj <- Read10X(all_samples)
obj <- CreateSeuratObject(obj)
obj$orig.ident <- as.character(obj$orig.ident)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
saveRDS(obj, "~/data/obj.rds")

###不进行整合
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:20, reduction = "pca")
obj <- FindClusters(obj, resolution = 0.3, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:20,
               reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))

###普通整合
rm(list = ls())
gc()
obj <- readRDS("~/data/obj.rds")
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
ElbowPlot(obj,ndims = 50)
obj <- IntegrateLayers(object = obj, method = CCAIntegration, 
                       orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)
# re-join layers after integration
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
saveRDS(obj,"~/data/obj_cca.rds")
##探索分辨率
library(clustree)
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:25)
sce <- FindClusters(
  object = obj,
  resolution = c(seq(0,1.5,.2)) #探究0~1.5，间隔为0.2
)
clustree(sce@meta.data, prefix = "RNA_snn_res.")

###用较优的分辨率进行聚类
obj <- FindClusters(obj, resolution = 0.4)
obj <- RunUMAP(obj, dims = 1:25, reduction = "integrated.cca")
DimPlot(obj, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))

raw_counts <- GetAssayData(obj, assay="RNA", slot='counts')
library(cellmarkeraccordion)
data(accordion_marker)##里面有细胞类型
library(data.table)
clusters <- data.table(cell = rownames(obj@meta.data), 
                       cluster = obj@meta.data$RNA_snn_res.0.4)
output <- accordion(raw_counts, assay ="RNA", 
                    species ="Human",  
                    cluster_info = clusters, annotation_resolution= "cluster", 
                    max_n_marker = 30, 
                    CL_celltypes = c("naive B cell","epithelial cell",
                                     "natural killer cell",
                                     "CD4-positive, alpha-beta T cell",
                                     "CD8-positive, alpha-beta T cell",
                                     "myeloid cell",
                                     "plasma cell",
                                     "endothelial cell", "fibroblast"),
                    include_detailed_annotation_info = TRUE, plot = TRUE)

ori_meta <- obj@meta.data
acc_meta <- output$cluster_annotation
ori_meta$cell <- rownames(ori_meta)
ori_meta <- ori_meta %>% 
  rename(cluster = RNA_snn_res.0.4) %>% 
  left_join(.,acc_meta %>% as.data.frame())
rownames(ori_meta) <- ori_meta$cell

obj@meta.data <- ori_meta
DimPlot(obj, reduction = "umap", group.by = c("accordion_per_cluster"),
        label = T)
saveRDS(obj,"~/data/obj_after_anno.rds")

######marker绘图



