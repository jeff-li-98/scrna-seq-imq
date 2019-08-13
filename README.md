# UCSF Cheng Lab Psoriasis Analysis 

## Overview

This project attempts to explore differences in gene expression across various immune cell types affected in psoriasis, using imiquimod as a model.

## Setup

Most of the analysis can be easily cloned and run on your local computer, with the exception of the `RunUMAP()` function in [2-preprocess_data.R](https://github.com/jeff-li-98/ucsf-cheng-psoriasis/blob/master/2-preprocess_data.R)., which requires setting up a Python environment within your RStudio. See below for details if you'd like to run this. Otherwise, the setup is relatively simple:

1. Create a folder on your computer where you want this project to live. (e.g. `/imq-analysis`)
2. Open a new terminal in that folder
3. Clone this repo into that folder by entering this into terminal: `git clone https://github.com/jeff-li-98/ucsf-cheng-psoriasis`
4. 

<<<<<<< HEAD
=======
1. hi
2. Measure the total effect and indirect effect of Shoo the Flu on laboratory-confirmed influenza
3. Measure the total effect and indirect effect of Shoo the Flu on influenza-related hospitalizationfdfd
>>>>>>> 202a8ff74f74ce53c88ef2487742ee5a74a37e40

## Project Structure

https://berkeley.box.com/s/pvu598x3zkq40rwsn8pszg7bhcrduz2q

At a high level, the project is hierarchical and uses prefix numbering to indicate how data "flows" through the project. Raw data import, management, and cleaning is handled by scripts in `0-data-prep/`, with intermediaries saved as R binaries. These intermediate binaries are then ingested by scripts in `1-power/` and `3-analysis/` to perform statistical computation, and the results are formatted and saved. The saved results are then ingested by scripts in `4-figures/` and `5-tables/` to complete the reproducible pipeline.

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
