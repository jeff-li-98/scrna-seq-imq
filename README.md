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

`0-config.R`: File paths and packages are loaded in here. All analyses source this file in the pre-amble. 
`1-install_packages.R`: All dependencies are installed here. All analyses include `source(1-install_packages.R)` in the pre-amble as a comment. When running this project for the first time, run this line to install packages. After the initial install, there is no need to run this again.
`2-preprocess_data.R`: This script loads in the rnaseq data from the 10XGenomics formats `.mtx` (fill here) using `Read10X()` and creates Seurat objects. Integration anchors are found for the data and we run PCA, UMAP, and t-SNE in this script. We then find clusters and save this object as an RDS for use in other scripts. Two RDS files are generated from this script: `data_by_replicate.RDS` (grouping cells by replicate e.g. IMQ1, IMQ3, IMQ5) and `data_by_condition.RDS` (grouping cells by condition e.g. IMQ/VEH). 


As stated, necessary dependencies are Raw data import, management, and cleaning is handled by scripts in `0-data-prep/`, with intermediaries saved as R binaries. These intermediate binaries are then ingested by scripts in `1-power/` and `3-analysis/` to perform statistical computation, and the results are formatted and saved. The saved results are then ingested by scripts in `4-figures/` and `5-tables/` to complete the reproducible pipeline.

For those attempting to replicate results with posession of the raw data, a utility bash script (`0-run-project.sh`) has been included. This calls other bash scripts throughout the project (which sometimes call more bash scripts), so it is recommended that you navigate the tree of calls to understand how the scripts are being used. You will likely need to respecify directory paths in the bash scripts and `0-config.R` to fit your specific file structure and needs.

Finally, you'll notice a Python script, called `runFileSaveLogs`. This is an improved version of an `R CMD BATCH` call, and can be used to run any R Script and save the log files (`*.Rout`) produced to a directory of your choosing.
  - **Usage**: `./runFileSaveLogs [-h] [-p] [-l LOGDIRPREFIX] [-i IDENTIFIER] filepaths [filepaths ...]`
  - **Description**: Runs the argument R script(s) - in parallel if specified - and moves the
subsequent generated .Rout log files to a time-stamped, user-stamped, and
optionally, identifier-stamped directory.
  - **Options**
    - `-h` => help flag: show a help message and exits
    - `-p` => parallelization flag: runs the argument R scripts in parallel if specified (current version of `runFileSaveLogs` does not allow for log export if this flag is specified)
    - `-l` => logDirPrefix optional argument: The directory in which log files will be saved to (defaults to `/data/flu/flu-logs/`)
    - `-i` => identifier optional argument: Adds an identifier to the directory name where the logs are saved
    - `filepaths` => filepaths required argument(s): The filepath(s) of the R scripts to be run
  - **Example Usage**: `./runFileSaveLogs -l ~/flu_logs -i run_two_scripts 0-data-prep/1-prep-cdph-fluseas.R 3-analysis/1-primary/5a-absentee-p2-glm-unadj.R`
    - (will run specified files and place logs in a directory formatted at `~/flu_logs/[current_time]-[current_user]-run_two_scripts`)
    

Created by [Jeffrey Li](https://github.com/jeff-li-98), [Nolan Pokpongkiat](https://github.com/nolanpokpongkiat).
