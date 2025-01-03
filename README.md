# eQTL-Catalogue/qtlmap output plot data generator workflow

The workflow builds upon [ Nurlan Kerimov's previous coverage plot workflow](https://github.com/kerimoff/coverage_plot) and incorporates the logic from Uku Raudvere's eqtl_plot_parsing workflow.


## Running the pipeline
```bash
nextflow run main.nf -profile tartu_hpc -resume\
  --studyFile ../GEUVADIS_studyfile.tsv\
  --outdir ../GEUVADIS_output
```
studyFile has to contain columns: 
* quant_method - ge/exon/tx/txrev/leafcutter
* dataset_id
* qtl_group	-  qtl_group in the study
* credible_sets_file	- File to the credible_sets.parquet file from the qtlmap workflow (./susie/*credible_sets.parquet)
* sample_meta	- Sample metadata file. Tab separated file
* vcf_file	
* bigwig_path	- Path to the bigwig files
* usage_matrix_norm	- Path to the normalised usage matrix
* exon_summ_stats_files	- Path to the file that contains full exon nominal summary statistics file paths (without header)
* all_summ_stats_files	- Path to the file that contains full gene nominal summary statistics file paths (without header)
* pheno_meta - Phenotype metadata file. Tab separated file
* scaling_factors - Path to the scaling_factors file

