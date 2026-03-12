#!/bin/bash

#SBATCH --time=120:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --job-name="plot_data"

module load any/jdk/1.8.0_265
module load nextflow/22.04.3
module load any/singularity/3.7.3
module load squashfs/4.4
module load tabix

nextflow info

nextflow run main.nf -profile tartu_hpc -resume \
  --studyFile input/GTEx_V10_all/GTEx_V10_inputs_tx_all.tsv\
  --outdir /gpfs/helios/projects/eQTLCatalogue/coverage_plots/GTEx_V10_all_tx
