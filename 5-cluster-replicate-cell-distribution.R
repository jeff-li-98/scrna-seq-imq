##########################################
# scrna seq imq visualizations
# deg analysis
#
# jeff li (jeff.li@berkeley.edu)
# nolan pokpongkiat (nolanpokpongkiat@berkeley.edu)
#
# 1) generate distribution of cells from each replicate, per cluster 
#
#
# inputs:
# "data_by_replicate.RDS"
# "data_by_condition.RDS"

# outputs: 
# "results/data/num_cells_per_cluster_by_stim.csv"

# install packages:
# source("1-install-packages.R")
##########################################

# Setup
source("0-config.R")

# Load the processed data
all.combined.replicate <- readRDS(file = paste0(data_dir, "data_by_replicate.RDS"))
all.combined.condition <- readRDS(file = paste0(data_dir, "data_by_condition.RDS"))


####################### 
#    MAKE DISTRIBUTION 
#
#######################

# Plot UMAP by stim
all.combined$cellType <- Idents(all.combined)
Idents(all.combined) <- all.combined$stim
DimPlot(all.combined, reduction = "umap")
Idents(all.combined) <- all.combined@meta.data$seurat_clusters

# Tabulate number of cells in each cluster by stim
cluster_by_stim_breakdown <- table(all.combined@meta.data$seurat_clusters, all.combined@meta.data$stim) 
write.csv(cluster_by_stim_breakdown, "results/data/num_cells_per_cluster_by_stim.csv") # path is from cwd