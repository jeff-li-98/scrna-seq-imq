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

# install packages:
# source("install_packages.R")
##########################################

# Setup
source("0-config.R")

# Load the processed data
all.combined.replicate <- readRDS(file = paste0(data_dir, "data_by_replicate.RDS"))
all.combined.condition <- readRDS(file = paste0(data_dir, "data_by_condition.RDS"))

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
# Heatmap of all genes (check this -NP)
# --------

# Find all variable gene markers across all clusters
rna.markers <- FindAllMarkers(all.small, assay = "RNA", only.pos = TRUE)
write.csv(rna.markers, paste0(results_data_dir_cwd, "all_rna_markers.csv")) # path is from cwd
saveRDS(rna.markers, file = paste0(data_dir, "DE_gene_markers.RDS"))

# Plot heatmap for all gene markers
all.heatmap.alldeg <- DoHeatmap(all.small, features = unique(rna.markers$gene), angle = 90) + NoLegend()
all.heatmap.alldeg # plot and save manually w/ size = 5000, 5000

# ---------
# Heatmap of top 200 intercluster DEG (check this -NP)
# --------

# Find the top 200 DEG across all clusters to plot 
all.combined <- FindVariableFeatures(all.combined, selection.method = "vst", nfeatures = 2000)
all_top200 <- head(VariableFeatures(all.combined), 200)

# Plot heatmap for top 200 genes
all.heatmap.degtop200 <- DoHeatmap(all.small, features = all_top200, angle = 90)
all.heatmap.degtop200 # plot and save manually w/ size = 5000, 5000

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
saveRDS(top10_deg, file = paste0(data_dir, "top10_deg_from_each_cluster.RDS"))

#Clean up to reset Idents back to clusters 
Idents(all.combined.condition) <- all.combined.condition@meta.data$seurat_clusters

# Draw heatmap for top 200 deg obtained from top 10 from each cluster  
all.heatmap.10per <- DoHeatmap(all.small, features = top10_deg$gene, angle = 90, group.by = "stim") + NoLegend()
all.heatmap.10per

# ---------
# Heatmap of significant DEG between IMQ v. VEH from each cluster
# --------
top10_deg <- readRDS(file = paste0(data_dir, "top10_deg_from_each_cluster.RDS"))
degp_under05 <- top10_deg %>% filter(p_val_adj < 0.05) %>% group_by(cluster)
write.csv(degp_under05, paste0(results_data_dir_cwd, "degp_under05.csv")) # path is from cwd
hm_under05 <- DoHeatmap(all.small, features = degp_under05$gene, angle = 90, group.by = "seurat_clusters") + NoLegend()
hm_under05


