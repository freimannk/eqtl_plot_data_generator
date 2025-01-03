
process generate_plot_ge_data { 
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    label "process_medium"
    container "quay.io/kfkf33/coverage_plot:v3"


    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path("output_dir_*"),  emit: plot_data_directories

    path "${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log", emit: log_file

    script:
    """
    Rscript $projectDir/bin/generate_plot_ge_data.R \
        --name_of_study $dataset_id \
        --qtl_group $qtl_group \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --vcf_file $vcf_file \
        --bigwig_path $bigwig_path \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --usage_matrix_norm $usage_matrix_norm

    cp .command.log ${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}