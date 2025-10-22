# eQTL-Catalogue/qtlmap output plot data generator workflow

The workflow builds upon [ Nurlan Kerimov's previous coverage plot workflow](https://github.com/kerimoff/coverage_plot) and incorporates the logic from Uku Raudvere's eqtl_plot_parsing workflow.


## Running the pipeline
```bash
nextflow run main.nf -profile tartu_hpc -resume\
  --studyFile ../GEUVADIS_studyfile.tsv\
  --outdir ../GEUVADIS_output
```
studyFile has to contain columns: 
* dataset_id (unique)
* study_id
* quant_method - ge/exon/tx/txrev/leafcutter
* qtl_group	-  qtl_group in the study
* study_name
* credible_sets_file	- File to the credible_sets.parquet file from the qtlmap workflow (./susie/*credible_sets.parquet)
* sample_meta	- Sample metadata file. Tab separated file
* vcf_file	
* bigwig_path	- Path to the bigwig files
* usage_matrix_norm	- Path to the normalised usage matrix
* tpm_matrix - Path to the TPM matrix
* exon_summ_stats_files	- Path to the file that contains full exon nominal summary statistics file paths (without header). From qtlmap wf output: ../sumstats/{dataset_id}/all/ path.
* all_summ_stats_files	- Path to the file that contains full gene nominal summary statistics file paths (without header) From qtlmap wf output: ../sumstats/{dataset_id}/all/ path.
* pheno_meta - Phenotype metadata file. Tab separated file
* scaling_factors - Path to the scaling_factors file


In nextflow.config file: 
* vcf_sample_names_correction — default: false.
Set to true to modify sample names in VCF files.
* vcf_samples_old_string_part — default: "".
The substring to replace in VCF sample names.
* vcf_samples_new_string_part — default: "".
The replacement string for correcting VCF sample names.

