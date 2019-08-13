# UCSF Cheng Lab Psoriasis Analysis 

## Overview

This project attempts to explore differences in gene expression across various immune cell types affected in psoriasis, using imiquimod as a model.

## Setup

Most of the analysis can be easily cloned and run on your local computer, with the exception of the `RunUMAP()` function in [2-preprocess_data.R](https://github.com/jeff-li-98/ucsf-cheng-psoriasis/blob/master/2-preprocess_data.R)., which requires setting up a Python environment within your RStudio. See below for details if you'd like to run this. Otherwise, the setup is relatively simple:

1. Create a folder on your computer where you want this project to live. (e.g. `/imq-analysis`)
2. Open a new terminal in that folder
3. Clone this repo into that folder by entering this into terminal: `git clone https://github.com/jeff-li-98/ucsf-cheng-psoriasis`

This should pull down all the analysis files onto your local computer. Then, fetch the data:

1. Request access to the data from jeff.li@berkeley.edu.
2. Download the `0-data.zip` file  [here](https://berkeley.box.com/s/pvu598x3zkq40rwsn8pszg7bhcrduz2q).
3. Move the `0-data.zip` file to your repo (e.g. `imq-analysis/ucsf-cheng-psoriasis/`) and unzip the file.
4. The zip can now be deleted.

Your repo is now set up!

## Project Structure

At a high level, the project is hierarchical and uses prefix numbering to indicate the workflow through the project. The files are divided as follows:

  - `0-config.R`: File paths and packages are loaded in here. All analyses source this file in the pre-amble. 

  - `0-data.R`: Processed data is stored here.

  - `1-install_packages.R`: All dependencies are installed here. All analyses include `source(1-install_packages.R)` in the pre-amble as a comment. When running this project for the first time, run this line to install packages. After the initial install, there is no need to run this again.

  - `2-preprocess_data.R`: This script loads in the rnaseq data from the 10XGenomics formats `.mtx` (fill here) using `Read10X()` and creates Seurat objects. Integration anchors are found for the data and we run PCA, UMAP, and t-SNE in this script. We then find clusters and save this object as an RDS for use in other scripts. Two RDS files are generated from this script: `data_by_replicate.RDS` (grouping cells by replicate e.g. IMQ1, IMQ3, IMQ5) and `data_by_condition.RDS` (grouping cells by condition e.g. IMQ/VEH). 

  - `3-heatmaps.R`: As stated, this script generates various heatmaps plotting genes on the y-axis and cells on the x-axis. We first downsample the data to 200 cells max per cluster, and plot three different heatmaps:
  
\space\space\space\space 1. All differentially expressed genes
\space\space\space\space 2. Top 10 differentially expressed genes between IMQ v. VEH, per cluster
\space\space\space\space 3. Significantly differentially expressed genes (p < 0.05) between IMQ v. VEH, per cluster

  - `4-visualization.R`: This script generates a UMAP and t-SNE plot for high-level visualization of clusters. 

  - `5-results.R`: All figures and data generated from the analysis scripts are saved here.


Created by [Jeffrey Li](https://github.com/jeff-li-98), [Nolan Pokpongkiat](https://github.com/nolanpokpongkiat).
