##########################################
# scrna seq imq
# deg analysis
#
# jeff li (jeffli@berkeley.edu)
# nolan pokpongkiat (nolanpokpongkiat@berkeley.edu)
#
# load raw rnaseq data into R format and save for further analysis
#
# inputs: data/.

# outputs: 
# "data_by_replicate.RDS"
# "data_by_condition.RDS"
##########################################


###### IMQ/VEH separated by replicate

# Load in data using the Read10X function
IMQ1.data <-Read10X(data.dir = paste0(data_dir, "m_IMQ_7d_1_IMQ1_mex"))
IMQ1 <- CreateSeuratObject(counts =IMQ1.data, project = "1", min.cells = 3, min.features = 200)
IMQ1@meta.data$stim <- "IMQ1" 
IMQ1 <- subset(x = IMQ1, subset = nFeature_RNA > 500)
IMQ1 <- NormalizeData(object = IMQ1, verbose = FALSE)
IMQ1 <- FindVariableFeatures(object = IMQ1, selection.method = "vst", nfeatures = 2000)

IMQ3.data <-Read10X(data.dir = paste0(data_dir, "m_IMQ_7d_3_IMQ3_mex"))
IMQ3 <- CreateSeuratObject(counts =IMQ3.data, project = "2", min.cells = 3, min.features = 200)
IMQ3@meta.data$stim <- "IMQ3" 
IMQ3 <- subset(x = IMQ3, subset = nFeature_RNA > 500)
IMQ3 <- NormalizeData(object = IMQ3, verbose = FALSE)
IMQ3 <- FindVariableFeatures(object = IMQ3, selection.method = "vst", nfeatures = 2000)

IMQ6.data <-Read10X(data.dir = paste0(data_dir, "m_IMQ_7d_6_IMQ6_mex"))
IMQ6 <- CreateSeuratObject(counts =IMQ6.data, project = "3", min.cells = 3, min.features = 200)
IMQ6@meta.data$stim <- "IMQ6" 
IMQ6 <- subset(x = IMQ6, subset = nFeature_RNA > 500)
IMQ6 <- NormalizeData(object = IMQ6, verbose = FALSE)
IMQ6 <- FindVariableFeatures(object = IMQ6, selection.method = "vst", nfeatures = 2000)

VEH2.data <-Read10X(data.dir = paste0(data_dir, "m_VEH_7d_2_VEH2_mex"))
VEH2 <- CreateSeuratObject(counts =VEH2.data, project = "4", min.cells = 3, min.features = 200)
VEH2@meta.data$stim <- "VEH2" 
VEH2 <- subset(x = VEH2, subset = nFeature_RNA > 500)
VEH2 <- NormalizeData(object = VEH2, verbose = FALSE)
VEH2 <- FindVariableFeatures(object = VEH2, selection.method = "vst", nfeatures = 2000)

VEH3.data <-Read10X(data.dir = paste0(data_dir, "m_VEH_7d_3_VEH3_mex"))
VEH3 <- CreateSeuratObject(counts =VEH3.data, project = "5", min.cells = 3, min.features = 200)
VEH3@meta.data$stim <- "VEH3" 
VEH3 <- subset(x = VEH3, subset = nFeature_RNA > 500)
VEH3 <- NormalizeData(object = VEH3, verbose = FALSE)
VEH3 <- FindVariableFeatures(object = VEH3, selection.method = "vst", nfeatures = 2000)

VEH4.data <-Read10X(data.dir = paste0(data_dir, "m_VEH_7d_4_VEH4_mex"))
VEH4 <- CreateSeuratObject(counts =VEH4.data, project = "6", min.cells = 3, min.features = 200)
VEH4@meta.data$stim <- "VEH4" 
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
saveRDS(all.combined, file = paste0(cache_dir, "data_by_replicate.RDS"))

###### IMQ/VEH separated by condition

# Group all IMQ together, all VEH together in condition var
all.combined@meta.data = all.combined@meta.data %>% mutate(stim = case_when(stim == "IMQ1" ~ "IMQ",
                                                                                 stim == "IMQ3" ~ "IMQ",
                                                                                 stim == "IMQ6" ~ "IMQ",
                                                                                 stim == "VEH2" ~ "VEH",
                                                                                 stim == "VEH3" ~ "VEH",
                                                                                 stim == "VEH4" ~ "VEH"))
saveRDS(all.combined, file = paste0(cache_dir, "data_by_condition.RDS"))
