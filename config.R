#######################
#  Config
# ---------------------
# Set up directories
library(here)

# Define directories
data_dir = here("data/")
results_dir = here("results")
fig_dir = here(results_dir, "figures/")
results_data_dir = here(results_dir, "data/")
cache_dir = here("cache/")

# Define project functions
rbind_rnames <- function(datalist) {
  require(plyr)
  temp <- rbind.fill(datalist)
  rownames(temp) <- unlist(lapply(datalist, row.names))
  temp
}
