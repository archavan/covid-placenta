#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=write_count_matrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load R

cd /home/arc78/scratch60/covid-placenta/code

Rscript  --no-save 99_paper_write-count-matrix.R

# gzip files

cd /home/arc78/scratch60/covid-placenta/data/count_matrix
gzip counts.tsv
gzip metadata.tsv
