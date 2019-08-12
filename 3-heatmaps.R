##########################################
# scrna seq imq visualizations
# deg analysis
#
# jeff li (jeffli@berkeley.edu)
# nolan pokpongkiat (nolanpokpongkiat@berkeley.edu)
#
# 1) generate heatmap of DEG
# 2) 
#
#
# inputs:
# "data_by_replicate.RDS"
# "data_by_condition.RDS"

# outputs: 

##########################################


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

# Load the raw data
all.combined.replicate <- readRDS(file = paste0(cache_dir, "data_by_replicate.RDS"))
all.combined.condition <- readRDS(file = paste0(cache_dir, "data_by_condition.RDS"))

# # Plot UMAP by stim
# all.combined$cellType <- Idents(all.combined)
# Idents(all.combined) <- all.combined$stim
# DimPlot(all.combined, reduction = "umap")
# Idents(all.combined) <- all.combined@meta.data$seurat_clusters

# # Tabulate number of cells in each cluster by stim
# cluster_by_stim_breakdown <- table(all.combined@meta.data$seurat_clusters, all.combined@meta.data$stim) 
# write.csv(cluster_by_stim_breakdown, "results/data/num_cells_per_cluster_by_stim.csv") # path is from cwd


####################### 
#    MAKE HEATMAPS 
#
#######################

### Downsample the clusters to a maximum of 200 cells each (makes the heatmap easier to see for small clusters)
all.small <- subset(all.combined.condition, downsample = 200)
all.small <- ScaleData(all.small, features = rownames(all.small))

# # Clean
# all.small@meta.data$seurat_clusters <- Idents(all.small)
# all.small@meta.data$orig.ident <- Idents(all.small)
# all.small@meta.data$stim <- Idents(all.small)
# all.small@meta.data$condition <- Idents(all.small)
# all.small@meta.data$cluster_condition <- Idents(all.small)
# all.small@meta.data$integrated_snn_res.0.47 <- Idents(all.small)

# ---------
# Heatmap of all genes (accuracy needs verification)
# --------

# Find all variable gene markers across all clusters
rna.markers <- FindAllMarkers(all.small, assay = "RNA", only.pos = TRUE)
write.csv(rna.markers, paste0(results_data_dir_cwd, "all_rna_markers.csv")) # path is from cwd
saveRDS(rna.markers, file = paste0(cache_dir, "DE_gene_markers.RDS"))

# Plot heatmap for all gene markers
all.heatmap.alldeg <- DoHeatmap(all.small, features = unique(rna.markers$gene), angle = 90) + NoLegend()
all.heatmap.alldeg # plot and save manually w/ size = 5000, 5000

# ---------
# Heatmap of top 200 intercluster DEG (accuracy needs verification)
# --------

# Find the top 200 DEG across all clusters to plot 
all.combined <- FindVariableFeatures(all.combined, selection.method = "vst", nfeatures = 2000)
all_top200 <- head(VariableFeatures(all.combined), 200)

# Plot heatmap for top 200 genes
all.heatmap.degtop200 <- DoHeatmap(all.small, features = all_top200, angle = 90)
all.heatmap.degtop200 # plot and save manually w/ size = 5000, 5000

# # Plot of unknown use???
# all_plot1 <- VariableFeaturePlot(all.combined)
# all_plot2 <- LabelPoints(plot = all_plot1, points = all_top200, repel = TRUE)
# CombinePlots(plots = list(all_plot1, all_plot2)) 

# ---------
# Heatmap of top 10 DEG from each cluster
# --------

# Identify top 10 differentially expressed genes between conditions of the same cluster for all clusters (total = 10genes x num clusters)
# New col name "$cluster_$condition" (1_IMQ, 3_VEH...)
all.combined.condition@meta.data$cluster_condition <- paste(all.combined.condition@meta.data$seurat_clusters, all.combined.condition@meta.data$stim, sep = "_")
Idents(all.combined.condition) <- all.combined.condition@meta.data$cluster_condition

all.combined.deg <- FindMarkers(all.combined.condition, ident.1 = "0_IMQ", ident.2 = "0_VEH", verbose = FALSE)
top10_deg <- head(all.combined.deg, n = 10)
top10_deg <- top10_deg %>% mutate(gene = rownames(top10_deg), cluster = 0)

for (x in 1:20) {
  top10_deg_per_cluster <- head(FindMarkers(all.combined.condition,
                                            ident.1 = paste0(x, "_IMQ"), 
                                            ident.2 = paste0(x, "_VEH"), 
                                            verbose = FALSE), 
                                n=10)
  top10_deg_per_cluster <- top10_deg_per_cluster %>% 
    mutate(gene = rownames(top10_deg_per_cluster), cluster = x)
  top10_deg <- rbind(top10_deg, top10_deg_per_cluster)
}
write.csv(top10_deg, paste0(results_data_dir_cwd, "top10_deg_per_cluster.csv")) # path is from cwd
saveRDS(top10_deg, file = paste0(cache_dir, "top10_deg_from_each_cluster.RDS"))

#Clean up to reset Idents back to clusters 
Idents(all.combined.condition) <- all.combined.condition@meta.data$seurat_clusters

# Draw heatmap for top 200 deg obtained from top 10 from each cluster  
all.heatmap.10per <- DoHeatmap(all.small, features = top10_deg$gene, angle = 90, group.by = "stim") + NoLegend()
all.heatmap.10per

# ---------
# Heatmap of DEG p<0.05 from each cluster
# --------
top10_deg <- readRDS(file = paste0(cache_dir, "top10_deg_from_each_cluster.RDS"))
degp_under05 <- top10_deg %>% filter(p_val_adj < 0.05) %>% group_by(cluster)
write.csv(degp_under05, paste0(results_data_dir_cwd, "degp_under05.csv")) # path is from cwd
hm_under05 <- DoHeatmap(all.small, features = degp_under05$gene, angle = 90, group.by = "seurat_clusters") + NoLegend()
hm_under05

# # IMQ Only 
# # Combine IMQ Trials 
# imq.combined <- FindIntegrationAnchors(object.list = list(IMQ1, IMQ3, IMQ6), dims = 1:20)
# imq.combined <- IntegrateData(anchorset = imq.combined, dims = 1:20)
# imq.combined <- ScaleData(object = imq.combined, verbose = FALSE)
# imq.combined <- RunPCA(object = imq.combined, npcs = 30, verbose = FALSE)
# 
# # Run UMAP & t-SNE 
# imq.combined <- RunUMAP(object = imq.combined, reduction = "pca", dims = 1:20)
# imq.combined <- RunTSNE(object = imq.combined, reduction = "pca", dims = 1:20)
# imq.combined <- FindNeighbors(object = imq.combined, reduction = "pca", dims = 1:20)
# imq.combined <- FindClusters(imq.combined, resolution = 1.08) #20 clusters
# 
# # Plot t-SNE
# imq_p1 <- DimPlot(object = imq.combined, reduction = "tsne", group.by = "stim") + NoLegend()
# imq_p2 <- DimPlot(object = imq.combined, reduction = "tsne", label = TRUE) + NoLegend()
# plot_grid(imq_p1, imq_p2)
# 
# # Plot UMAP
# imq_p3 <- DimPlot(object = imq.combined, reduction = "umap", group.by = "stim", label = FALSE, pt.size = 0.003) + NoLegend()
# imq_p4 <- DimPlot(object = imq.combined, reduction = "umap", label = TRUE, pt.size = 0.003) + NoLegend()
# plot_grid(imq_p3, imq_p4)
# 
# # Identify top 200 genes to plot 
# imq.combined <- FindVariableFeatures(imq.combined, selection.method = "vst", nfeatures = 2000)
# imq_top200 <- head(VariableFeatures(imq.combined), 200)
# imq_plot1 <- VariableFeaturePlot(imq.combined)
# imq_plot2 <- LabelPoints(plot = imq_plot1, points = imq_top200, repel = TRUE)
# CombinePlots(plots = list(imq_plot1, imq_plot2))
# 
# # Downsample the clusters to a maximum of 200 cells each (makes the heatmap easier to see for
# # small clusters)
# imq.small <- subset(imq.combined, downsample = 200)
# imq.small <- ScaleData(imq.small, features = rownames(all.small))
# 
# # Draw heatmap with top 200 genes 
# imq.heatmap <- DoHeatmap(imq.small, features = imq_top200, angle = 90)
# imq.heatmap
# 
# #### VEH Only 
# # Combine VEH Trials 
# veh.combined <- FindIntegrationAnchors(object.list = list(VEH2, VEH3, VEH4), dims = 1:20)
# veh.combined <- IntegrateData(anchorset = veh.combined, dims = 1:20)
# veh.combined <- ScaleData(object = veh.combined, verbose = FALSE)
# veh.combined <- RunPCA(object = veh.combined, npcs = 30, verbose = FALSE)
# 
# # Run UMAP & t-SNE 
# veh.combined <- RunUMAP(object = veh.combined, reduction = "pca", dims = 1:20)
# veh.combined <- RunTSNE(object = veh.combined, reduction = "pca", dims = 1:20)
# veh.combined <- FindNeighbors(object = veh.combined, reduction = "pca", dims = 1:20)
# veh.combined <- FindClusters(veh.combined, resolution = 0.9) #21 clusters
# 
# # Plot t-SNE
# veh_p1 <- DimPlot(object = veh.combined, reduction = "tsne", group.by = "stim") + NoLegend()
# veh_p2 <- DimPlot(object = veh.combined, reduction = "tsne", label = TRUE) + NoLegend()
# plot_grid(veh_p1, veh_p2)
# 
# # Plot UMAP
# veh_p3 <- DimPlot(object = veh.combined, reduction = "umap", group.by = "stim", label = FALSE, pt.size = 0.003) + NoLegend()
# veh_p4 <- DimPlot(object = veh.combined, reduction = "umap", label = TRUE, pt.size = 0.003) + NoLegend()
# plot_grid(veh_p3, veh_p4)
# 
# # Identify top 200 genes to plot 
# veh.combined <- FindVariableFeatures(veh.combined, selection.method = "vst", nfeatures = 2000)
# veh_top200 <- head(VariableFeatures(veh.combined), 200)
# veh_plot1 <- VariableFeaturePlot(veh.combined)
# veh_plot2 <- LabelPoints(plot = veh_plot1, points = veh_top200, repel = TRUE)
# CombinePlots(plots = list(veh_plot1, veh_plot2))
# 
# # Downsample the clusters to a maximum of 200 cells each (makes the heatmap easier to see for
# # small clusters)
# veh.small <- subset(veh.combined, downsample = 200)
# veh.small <- ScaleData(veh.small, features = rownames(all.small))
# 
# # Draw heatmap with top 200 genes 
# veh.heatmap <- DoHeatmap(veh.small, features = veh_top200, angle = 90)
# veh.heatmap

