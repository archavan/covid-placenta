#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=findmarkers
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load R

cd /home/arc78/scratch60/covid-placenta

Rscript  --no-save 03a_FindAllMarkers.R
