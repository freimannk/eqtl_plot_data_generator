# eQTL-Catalogue/qtlmap output plot data generator workflow

Workflow uses previous [ Nurlan Kerimov's coverage plot workflow](https://github.com/kerimoff/coverage_plot) and Uku Raudvere's eqtl_plot_parsing workflow logic.


## Running the pipeline
```bash
nextflow run main.nf -profile tartu_hpc -resume\
  --studyFile ../GEUVADIS_studyfile.tsv\
  --outdir ../GEUVADIS_output
```
studyFile has to contain columns: quant_method	name_of_study	qtl_group	susie_purity_filtered	sample_meta	vcf_file	bigwig_path	usage_matrix_norm	exon_summ_stats_files	all_summ_stats_files	pheno_meta	scaling_factors
