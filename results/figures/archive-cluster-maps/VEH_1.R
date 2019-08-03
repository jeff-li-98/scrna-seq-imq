library(BiocManager)
library(Seurat)
library(cowplot)
library(dplyr)

library(reticulate)
conda_create("r-reticulate")
use_python("~/anaconda3/bin/python")
py_install("umap-learn")
# py_config()

### Load in VEH-1 data set ###
VEH1.data <-Read10X(data.dir = "m_VEH_7d_2_VEH2_mex")
VEH1 <- CreateSeuratObject(counts =VEH1.data, project = "1", min.cells = 3, min.features = 200)
VEH1@meta.data$stim <- "1 Epcam Sort"
VEH1 <- subset(x = VEH1, subset = nFeature_RNA > 500)
VEH1 <- NormalizeData(object = VEH1, verbose = FALSE)
VEH1 <- FindVariableFeatures(object = VEH1, selection.method = "vst", nfeatures = 2000)

### Load in VEH-2 data set ###
VEH2.data <-Read10X(data.dir = "m_VEH_7d_3_VEH3_mex")
VEH2 <- CreateSeuratObject(counts =VEH2.data, project = "2", min.cells = 3, min.features = 200)
VEH2@meta.data$stim <- "2 Epcam Sort"
VEH2 <- subset(x = VEH2, subset = nFeature_RNA > 500)
VEH2 <- NormalizeData(object = VEH2, verbose = FALSE)
VEH2 <- FindVariableFeatures(object = VEH2, selection.method = "vst", nfeatures = 2000)

### Load in VEH-3/ASUS/Desktop")
VEH3 <- CreateSeuratObject(counts =VEH2.data, project = "3", min.cells = 3, min.features = 200)
VEH3@meta.data$stim <- "3 Epcam Sort"
VEH3 <- subset(x = VEH3, subset = nFeature_RNA > 500)
VEH3 <- NormalizeData(object = VEH3, verbose = FALSE)
VEH3 <- FindVariableFeatures(object = VEH3, selection.method = "vst", nfeatures = 2000)

### Use Seurat v3 Integration Feature to merge data sets ###
VEH.combined <- FindIntegrationAnchors(object.list = list(VEH1, VEH2, VEH3), dims = 1:20)
VEH.combined <- IntegrateData(anchorset = VEH.combined, dims = 1:20)
VEH.combined <- ScaleData(object = VEH.combined, verbose = FALSE)
VEH.combined <- RunPCA(object = VEH.combined, npcs = 30, verbose = FALSE)
VEH.combined <- RunUMAP(object = VEH.combined, reduction = "pca", dims = 1:20)
VEH.combined <- RunTSNE(object = VEH.combined, reduction = "pca", dims = 1:20)
VEH.combined <- FindNeighbors(object = VEH.combined, reduction = "pca", dims = 1:20)
VEH.combined <- FindClusters(VEH.combined, resolution = 1.2)
save(VEH.combined, file ="Mouse VEH all Merged.Rdata")

### Dimension Reduction Visualization ###
p1 <- DimPlot(object = VEH.combined, reduction = "tsne", group.by = "stim") + NoLegend()
p2 <- DimPlot(object = VEH.combined, reduction = "tsne", label = TRUE) + NoLegend()
plot_grid(p1, p2)

p2 <- DimPlot(object = VEH.combined, reduction = "umap", group.by = "stim", label = FALSE) + NoLegend()
p1 <- DimPlot(object = VEH.combined, reduction = "umap", label = TRUE) + NoLegend()
plot_grid(p1, p2)
