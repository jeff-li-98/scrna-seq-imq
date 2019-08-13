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
data_dir = here("0-data/")
results_dir = here("6-results")
results_dir_cwd = "6-results/"
fig_dir = here(results_dir, "figures/")
fig_dir_cwd = paste0(results_dir_cwd, "figures/")
results_data_dir = here(results_dir, "data/")
results_data_dir_cwd = paste0(results_dir_cwd, "data/")


# Define project functions
rbind_rnames <- function(datalist) {
  require(plyr)
  temp <- rbind.fill(datalist)
  rownames(temp) <- unlist(lapply(datalist, row.names))
  temp
}
