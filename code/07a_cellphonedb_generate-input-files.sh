#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=generate_cellphonedb_input
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load R

cd /home/arc78/scratch60/covid-placenta/code

Rscript  --no-save 07a_cellphonedb_generate-input-files.R
