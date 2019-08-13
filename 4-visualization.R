##########################################
# scrna seq imq
# deg analysis
#
# jeff li (jeffli@berkeley.edu)
# nolan pokpongkiat (nolanpokpongkiat@berkeley.edu)
#
# 1) generate heatmap of deg by cluster
#
#
# inputs: 
# "data_by_replicate.RDS"
# "data_by_condition.RDS"

# outputs: 
# "results/figures/tSNE_imq_veh.pdf"
# "results/figures/UMAP_imq_veh.pdf"
#
# install packages:
# source("install_packages.R")
##########################################

# Setup
source("0-config.R")

# Load the processed data
all.combined.replicate <- readRDS(file = paste0(data_dir, "data_by_replicate.RDS"))
all.combined.condition <- readRDS(file = paste0(data_dir, "data_by_condition.RDS"))

# Plot condition t-SNE vs cluster t-SNE
all_p1 <- DimPlot(object = all.combined.condition, reduction = "tsne", group.by = "stim", pt.size = 0.003) 
all_p2 <- DimPlot(object = all.combined.condition, reduction = "tsne", label = TRUE, pt.size = 0.003) + NoLegend()
pdf(paste0(fig_dir_cwd, "tSNE_imq_veh.pdf")) # takes relative path from current dir
plot_grid(all_p1, all_p2)
dev.off()

# Plot condition UMAP vs cluster UMAP
all_p3 <- DimPlot(object = all.combined.condition, reduction = "umap", group.by = "stim", label = FALSE, cols = c("blue", "orange"), order = "IMQ", pt.size = 0.003)
all_p4 <- DimPlot(object = all.combined.condition, reduction = "umap", group.by = "stim", label = FALSE, cols = c("orange", "blue"), order = "VEH", pt.size = 0.003)
all_p5 <- DimPlot(object = all.combined.condition, reduction = "umap", label = TRUE, pt.size = 0.003) + NoLegend()
pdf(paste0(fig_dir_cwd, "UMAP_imq_veh.pdf")) # takes relative path from current dir
plot_grid(all_p3, all_p4, all_p5)
dev.off()


