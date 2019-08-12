#######################
#  Config
# ---------------------
# Load packages
library(BiocManager)
library(Seurat)
library(cowplot)
library(dplyr)
library(reticulate)

# Define directories
data_dir = here("data/")
data_dir_cwd = "data/"
results_dir = here("1-results")
results_dir_cwd = "1-results/"
fig_dir = here(results_dir, "figures/")
fig_dir_cwd = paste0(results_dir_cwd, "figures/")
results_data_dir = here(results_dir, "data/")
results_data_dir_cwd = paste0(results_dir_cwd, "data/")
cache_dir = here("0-cache/")

# Define project functions
rbind_rnames <- function(datalist) {
  require(plyr)
  temp <- rbind.fill(datalist)
  rownames(temp) <- unlist(lapply(datalist, row.names))
  temp
}
