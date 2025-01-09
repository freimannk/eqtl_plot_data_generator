#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --job-name="plot_data"

module load any/jdk/1.8.0_265
module load nextflow
module load any/singularity/3.7.3
module load squashfs/4.4
module load tabix

nextflow run main.nf -profile tartu_hpc -resume \
  --studyFile /gpfs/space/projects/otar2069/eQTL_plot_data_testdata/testdata/GEUVADIS_testdata_studyfile.tsv\
  --outdir testdata_results