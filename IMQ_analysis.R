#######################
#  Setup
# ---------------------
# IMPORTANT!!!! ONLY RUN THIS ONCE THE VERY FIRST TIME YOU OPEN THIS PROJECT
source("install_packages.R")
#######################

# Run these lines every time you open the project
source("config.R")
# Load packages
library(here)
library(BiocManager)
library(Seurat)
library(cowplot)
library(dplyr)
library(reticulate)

######################
# # Set up an Anaconda environment
# # IMPORTANT!!!! ONLY RUN THIS ONCE THE VERY FIRST TIME YOU OPEN THIS PROJECT
# # IMPORTANT!!!! ONLY RUN THIS ONCE THE VERY FIRST TIME YOU OPEN THIS PROJECT
# # IMPORTANT!!!! ONLY RUN THIS ONCE THE VERY FIRST TIME YOU OPEN THIS PROJECT
# conda_create("r-reticulate")
# conda_install("r-reticulate", "umap-learn")
# py_config()
#######################


####### IMQ+VEH Combined 

# Load in data using the Read10X function
IMQ1.data <-Read10X(data.dir = paste0(data_dir, "m_IMQ_7d_1_IMQ1_mex"))
IMQ1 <- CreateSeuratObject(counts =IMQ1.data, project = "1", min.cells = 3, min.features = 200)
IMQ1@meta.data$replicate <- "IMQ1" # 1 Epcam Sort
IMQ1 <- subset(x = IMQ1, subset = nFeature_RNA > 500)
IMQ1 <- NormalizeData(object = IMQ1, verbose = FALSE)
IMQ1 <- FindVariableFeatures(object = IMQ1, selection.method = "vst", nfeatures = 2000)

IMQ3.data <-Read10X(data.dir = paste0(data_dir, "m_IMQ_7d_3_IMQ3_mex"))
IMQ3 <- CreateSeuratObject(counts =IMQ3.data, project = "2", min.cells = 3, min.features = 200)
IMQ3@meta.data$replicate <- "IMQ3" # 2 Epcam Sort
IMQ3 <- subset(x = IMQ3, subset = nFeature_RNA > 500)
IMQ3 <- NormalizeData(object = IMQ3, verbose = FALSE)
IMQ3 <- FindVariableFeatures(object = IMQ3, selection.method = "vst", nfeatures = 2000)

IMQ6.data <-Read10X(data.dir = paste0(data_dir, "m_IMQ_7d_6_IMQ6_mex"))
IMQ6 <- CreateSeuratObject(counts =IMQ6.data, project = "3", min.cells = 3, min.features = 200)
IMQ6@meta.data$replicate <- "IMQ6" # 3 Epcam Sort
IMQ6 <- subset(x = IMQ6, subset = nFeature_RNA > 500)
IMQ6 <- NormalizeData(object = IMQ6, verbose = FALSE)
IMQ6 <- FindVariableFeatures(object = IMQ6, selection.method = "vst", nfeatures = 2000)

VEH2.data <-Read10X(data.dir = paste0(data_dir, "m_VEH_7d_2_VEH2_mex"))
VEH2 <- CreateSeuratObject(counts =VEH2.data, project = "4", min.cells = 3, min.features = 200)
VEH2@meta.data$replicate <- "VEH2" # 4 Epcam Sort
VEH2 <- subset(x = VEH2, subset = nFeature_RNA > 500)
VEH2 <- NormalizeData(object = VEH2, verbose = FALSE)
VEH2 <- FindVariableFeatures(object = VEH2, selection.method = "vst", nfeatures = 2000)

VEH3.data <-Read10X(data.dir = paste0(data_dir, "m_VEH_7d_3_VEH3_mex"))
VEH3 <- CreateSeuratObject(counts =VEH3.data, project = "5", min.cells = 3, min.features = 200)
VEH3@meta.data$replicate <- "VEH3" # 5 Epcam Sort
VEH3 <- subset(x = VEH3, subset = nFeature_RNA > 500)
VEH3 <- NormalizeData(object = VEH3, verbose = FALSE)
VEH3 <- FindVariableFeatures(object = VEH3, selection.method = "vst", nfeatures = 2000)

VEH4.data <-Read10X(data.dir = paste0(data_dir, "m_VEH_7d_4_VEH4_mex"))
VEH4 <- CreateSeuratObject(counts =VEH4.data, project = "6", min.cells = 3, min.features = 200)
VEH4@meta.data$replicate <- "VEH4" # 6 Epcam Sort
VEH4 <- subset(x = VEH4, subset = nFeature_RNA > 500)
VEH4 <- NormalizeData(object = VEH4, verbose = FALSE)
VEH4 <- FindVariableFeatures(object = VEH4, selection.method = "vst", nfeatures = 2000)

# Combine the runs and run PCA
all.combined <- FindIntegrationAnchors(object.list = list(IMQ1, IMQ3, IMQ6, VEH2, VEH3, VEH4), dims = 1:20)
all.combined <- IntegrateData(anchorset = all.combined, dims = 1:20)
all.combined <- ScaleData(object = all.combined, verbose = FALSE)
all.combined <- RunPCA(object = all.combined, npcs = 30, verbose = FALSE)

# Run UMAP & t-SNE 
all.combined <- RunUMAP(object = all.combined, reduction = "pca", dims = 1:20)
all.combined <- RunTSNE(object = all.combined, reduction = "pca", dims = 1:20)
all.combined <- FindNeighbors(object = all.combined, reduction = "pca", dims = 1:20)
all.combined <- FindClusters(all.combined, resolution = 0.47) #20 clusters

# # Group all IMQ together, all VEH together in stim var
all.combined@meta.data = all.combined@meta.data %>% mutate(stim = case_when(replicate == "IMQ1" ~ "IMQ",
                                                                            replicate == "IMQ3" ~ "IMQ",
                                                                            replicate == "IMQ6" ~ "IMQ",
                                                                            replicate == "VEH2" ~ "VEH",
                                                                            replicate == "VEH3" ~ "VEH",
                                                                            replicate == "VEH4" ~ "VEH"))

# Plot t-SNE by cluster 
all_p1 <- DimPlot(object = all.combined, reduction = "tsne", group.by = "replicate") + NoLegend() # what does group.by do????
all_p2 <- DimPlot(object = all.combined, reduction = "tsne", label = TRUE) + NoLegend()
pdf("tSNE_imq_veh_by_cluster.pdf")
plot_grid(all_p1, all_p2)
dev.off()

# Plot UMAP by cluster 
all_p3 <- DimPlot(object = all.combined, reduction = "umap", group.by = "replicate", label = FALSE, pt.size = 0.003) + NoLegend()
all_p4 <- DimPlot(object = all.combined, reduction = "umap", label = TRUE, pt.size = 0.003) + NoLegend()
pdf("UMAP_imq_veh_by_cluster.pdf")
plot_grid(all_p3, all_p4)
dev.off()

# Plot UMAP by 
all.combined$cellType <- Idents(all.combined)
Idents(all.combined) <- all.combined$replicate
DimPlot(all.combined, reduction = "umap")
Idents(all.combined) <- all.combined@meta.data$seurat_clusters

# Tabulate number of cells in each cluster by replicate
cluster_by_replicate_breakdown <- table(all.combined@meta.data$seurat_clusters, all.combined@meta.data$replicate) 
write.csv(cluster_by_replicate_breakdown, "results/data/num_cells_per_cluster_by_replicate.csv") # path is from cwd

# Identify top 10 differentially expressed genes between stim of the same cluster for all clusters (total = 10genes x num clusters)
all.combined@meta.data$cluster_stim <- paste(all.combined@meta.data$seurat_clusters, all.combined@meta.data$stim, sep = "_")
Idents(all.combined) <- all.combined@meta.data$cluster_stim

all.combined.deg <- FindMarkers(all.combined, ident.1 = "0_IMQ", ident.2 = "0_VEH", verbose = FALSE)
top200_deg <- head(all.combined.deg, n = 10)
top200_deg <- top200_deg %>% mutate(gene = rownames(top200_deg), cluster = 0)

for (x in 1:20) {
  top20_deg_curr_cluster <- head(FindMarkers(all.combined,
                                             ident.1 = paste0(x, "_IMQ"), 
                                             ident.2 = paste0(x, "_VEH"), 
                                             verbose = FALSE), 
                                 n=10)
  top20_deg_curr_cluster <- top20_deg_curr_cluster %>% 
    mutate(gene = rownames(top20_deg_curr_cluster), cluster = x)
  top200_deg <- rbind(top200_deg, top20_deg_curr_cluster)
}
rownames(top200_deg) <- top200_deg$gene
write.csv(top200_deg, "results/data/top10_deg_per_cluster.csv") # path is from cwd

# Identify top 200 DEG across all clusters to plot 
all.combined <- FindVariableFeatures(all.combined, selection.method = "vst", nfeatures = 2000)
all_top200 <- head(VariableFeatures(all.combined), 200)
# all_plot1 <- VariableFeaturePlot(all.combined)
# all_plot2 <- LabelPoints(plot = all_plot1, points = all_top200, repel = TRUE)
# pdf("combined_plot.pdf")
# CombinePlots(plots = list(all_plot1, all_plot2))
# dev.off()


### Downsample the clusters to a maximum of 200 cells each (makes the heatmap easier to see for small clusters)
all.small <- subset(all.combined, downsample = 200)
all.small <- ScaleData(all.small, features = rownames(all.small))

# # Find all variable gene markers across all clusters 
# rna.markers <- FindAllMarkers(all.small, assay = "RNA", only.pos = TRUE)
# write.csv(rna.markers, "results/data/all_rna_markers.csv") # path is from cwd

# # Draw heatmap for all variable gene markers; code below 
# DoHeatmap(all.small, features = unique(rna.markers$gene), angle = 90) + NoLegend()

# Draw heatmap for top 200 genes
all.heatmap <- DoHeatmap(all.small, features = all_top200, angle = 90)
ggsave2(all.heatmap, "combined_heatmap_21_clusters.pdf")
# pdf("combined_plot.pdf")
# CombinePlots(plots = list(all_plot1, all_plot2))
# dev.off()

# Draw heatmap for top 200 deg obtained from top 10 from each cluster  
all.heatmap.10per <- DoHeatmap(all.small, features = top200_deg$gene, angle = 90)
all.heatmap.10per


# IMQ Only 
# Combine IMQ Trials 
imq.combined <- FindIntegrationAnchors(object.list = list(IMQ1, IMQ3, IMQ6), dims = 1:20)
imq.combined <- IntegrateData(anchorset = imq.combined, dims = 1:20)
imq.combined <- ScaleData(object = imq.combined, verbose = FALSE)
imq.combined <- RunPCA(object = imq.combined, npcs = 30, verbose = FALSE)

# Run UMAP & t-SNE 
imq.combined <- RunUMAP(object = imq.combined, reduction = "pca", dims = 1:20)
imq.combined <- RunTSNE(object = imq.combined, reduction = "pca", dims = 1:20)
imq.combined <- FindNeighbors(object = imq.combined, reduction = "pca", dims = 1:20)
imq.combined <- FindClusters(imq.combined, resolution = 1.08) #20 clusters

# Plot t-SNE
imq_p1 <- DimPlot(object = imq.combined, reduction = "tsne", group.by = "stim") + NoLegend()
imq_p2 <- DimPlot(object = imq.combined, reduction = "tsne", label = TRUE) + NoLegend()
plot_grid(imq_p1, imq_p2)

# Plot UMAP
imq_p3 <- DimPlot(object = imq.combined, reduction = "umap", group.by = "stim", label = FALSE, pt.size = 0.003) + NoLegend()
imq_p4 <- DimPlot(object = imq.combined, reduction = "umap", label = TRUE, pt.size = 0.003) + NoLegend()
plot_grid(imq_p3, imq_p4)

# Identify top 200 genes to plot 
imq.combined <- FindVariableFeatures(imq.combined, selection.method = "vst", nfeatures = 2000)
imq_top200 <- head(VariableFeatures(imq.combined), 200)
imq_plot1 <- VariableFeaturePlot(imq.combined)
imq_plot2 <- LabelPoints(plot = imq_plot1, points = imq_top200, repel = TRUE)
CombinePlots(plots = list(imq_plot1, imq_plot2))

# Downsample the clusters to a maximum of 200 cells each (makes the heatmap easier to see for
# small clusters)
imq.small <- subset(imq.combined, downsample = 200)
imq.small <- ScaleData(imq.small, features = rownames(all.small))

# Draw heatmap with top 200 genes 
imq.heatmap <- DoHeatmap(imq.small, features = imq_top200, angle = 90)
imq.heatmap

#### VEH Only 
#Combine VEH Trials 
veh.combined <- FindIntegrationAnchors(object.list = list(VEH2, VEH3, VEH4), dims = 1:20)
veh.combined <- IntegrateData(anchorset = veh.combined, dims = 1:20)
veh.combined <- ScaleData(object = veh.combined, verbose = FALSE)
veh.combined <- RunPCA(object = veh.combined, npcs = 30, verbose = FALSE)

# Run UMAP & t-SNE 
veh.combined <- RunUMAP(object = veh.combined, reduction = "pca", dims = 1:20)
veh.combined <- RunTSNE(object = veh.combined, reduction = "pca", dims = 1:20)
veh.combined <- FindNeighbors(object = veh.combined, reduction = "pca", dims = 1:20)
veh.combined <- FindClusters(veh.combined, resolution = 0.9) #21 clusters

# Plot t-SNE
veh_p1 <- DimPlot(object = veh.combined, reduction = "tsne", group.by = "stim") + NoLegend()
veh_p2 <- DimPlot(object = veh.combined, reduction = "tsne", label = TRUE) + NoLegend()
plot_grid(veh_p1, veh_p2)

# Plot UMAP
veh_p3 <- DimPlot(object = veh.combined, reduction = "umap", group.by = "stim", label = FALSE, pt.size = 0.003) + NoLegend()
veh_p4 <- DimPlot(object = veh.combined, reduction = "umap", label = TRUE, pt.size = 0.003) + NoLegend()
plot_grid(veh_p3, veh_p4)

# Identify top 200 genes to plot 
veh.combined <- FindVariableFeatures(veh.combined, selection.method = "vst", nfeatures = 2000)
veh_top200 <- head(VariableFeatures(veh.combined), 200)
veh_plot1 <- VariableFeaturePlot(veh.combined)
veh_plot2 <- LabelPoints(plot = veh_plot1, points = veh_top200, repel = TRUE)
CombinePlots(plots = list(veh_plot1, veh_plot2))

# Downsample the clusters to a maximum of 200 cells each (makes the heatmap easier to see for
# small clusters)
veh.small <- subset(veh.combined, downsample = 200)
veh.small <- ScaleData(veh.small, features = rownames(all.small))

# Draw heatmap with top 200 genes 
veh.heatmap <- DoHeatmap(veh.small, features = veh_top200, angle = 90)
veh.heatmap

