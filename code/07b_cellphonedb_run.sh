#! /bin/bash

### INSTALLING CellPhoneDB ####################################################
# NOTE: this works with Python 3.5 or greater.
# 1. Create python virtual environment ----------------------------------------
cd /Users/Arun/research/covid-placenta/results/06_cellphonedb
python -m venv cpdb-venv

# 2. Activate virtualenv
source cpdb-venv/bin/activate

# 3. Install cellphoneDB
pip install cellphonedb

### RUNNING CellPhoneDB #######################################################
## covid samples --------------------------------------------------------------
cellphonedb method statistical_analysis ./01_input/covid_meta.txt ./01_input/covid_count.txt --project-name=covid --counts-data=gene_name --subsampling --subsampling-log true --subsampling-num-cells 5000

# Plotting statistical method results
cellphonedb plot dot_plot --means-path=./out/covid/means.txt --pvalues-path=./out/covid/pvalues.txt --output-path=./out/covid

cellphonedb plot heatmap_plot ./01_input/covid_meta.txt --pvalues-path=./out/covid/pvalues.txt --output-path=./out/covid 

## control samples ------------------------------------------------------------
cellphonedb method statistical_analysis ./01_input/cntrl_meta.txt ./01_input/cntrl_count.txt --project-name=cntrl --counts-data=gene_name --subsampling --subsampling-log true --subsampling-num-cells 5000

# Plotting statistical method results
cellphonedb plot dot_plot --means-path=./out/cntrl/means.txt --pvalues-path=./out/cntrl/pvalues.txt --output-path=./out/cntrl

cellphonedb plot heatmap_plot ./01_input/cntrl_meta.txt --pvalues-path=./out/cntrl/pvalues.txt --output-path=./out/cntrl
 