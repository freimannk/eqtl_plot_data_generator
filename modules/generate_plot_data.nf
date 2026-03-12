
process generate_plot_ge_data { 
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    label "generate_plot_data"
    container "quay.io/kerimoff/coverage_plot:v6"


    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), file(sample_meta), file(coverage_parquet), file(usage_matrix_norm), file(tpm_matrix), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path("output_dir_*"),  emit: plot_data_directories

    path "${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log", emit: log_file

    script:
    vcf_sample_names_correction_bad_symbol = params.vcf_sample_names_correction ? "--vcf_sample_bad_symbol {$params.vcf_samples_old_string_part}" : ""
    vcf_sample_names_correction_replaced_string = params.vcf_sample_names_correction ? "--vcf_sample_replacement_symbol {$params.vcf_samples_new_string_part}" : ""
    """
    Rscript $projectDir/bin/generate_plot_ge_data.R \
        --qtl_group $qtl_group \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --vcf_file $vcf_file \
        --coverage_parquet $coverage_parquet \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --usage_matrix_norm $usage_matrix_norm \
        --tpm_matrix $tpm_matrix $vcf_sample_names_correction_bad_symbol $vcf_sample_names_correction_replaced_string

    cp .command.log ${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}


process generate_plot_exon_data { 
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    label "generate_plot_data"
    container "quay.io/kerimoff/coverage_plot:v6"


    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), file(sample_meta), file(coverage_parquet), file(usage_matrix_norm), file(tpm_matrix), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path("output_dir_*"),  emit: plot_data_directories

    path "${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log", emit: log_file

    script:
    vcf_sample_names_correction_bad_symbol = params.vcf_sample_names_correction ? "--vcf_sample_bad_symbol {$params.vcf_samples_old_string_part}" : ""
    vcf_sample_names_correction_replaced_string = params.vcf_sample_names_correction ? "--vcf_sample_replacement_symbol {$params.vcf_samples_new_string_part}" : ""
    """
    Rscript $projectDir/bin/generate_plot_exon_data.R \
        --qtl_group $qtl_group \
        --phenotype_meta $phenotype_meta \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --vcf_file $vcf_file \
        --coverage_parquet $coverage_parquet \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --usage_matrix_norm $usage_matrix_norm \
        --tpm_matrix $tpm_matrix $vcf_sample_names_correction_bad_symbol $vcf_sample_names_correction_replaced_string

    cp .command.log ${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}

process generate_plot_leafcutter_data { 
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    label "generate_plot_data"
    container "quay.io/kerimoff/coverage_plot:v6"


    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), file(sample_meta), file(coverage_parquet), file(usage_matrix_norm), file(tpm_matrix), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path("output_dir_*"),  emit: plot_data_directories

    path "${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log", emit: log_file

    script:
    vcf_sample_names_correction_bad_symbol = params.vcf_sample_names_correction ? "--vcf_sample_bad_symbol {$params.vcf_samples_old_string_part}" : ""
    vcf_sample_names_correction_replaced_string = params.vcf_sample_names_correction ? "--vcf_sample_replacement_symbol {$params.vcf_samples_new_string_part}" : ""
    """
    Rscript $projectDir/bin/generate_leafcutter_plot_data.R \
        --qtl_group $qtl_group \
        --phenotype_meta $phenotype_meta \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --vcf_file $vcf_file \
        --coverage_parquet $coverage_parquet \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --usage_matrix_norm $usage_matrix_norm \
        --tpm_matrix $tpm_matrix $vcf_sample_names_correction_bad_symbol $vcf_sample_names_correction_replaced_string

    cp .command.log ${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}

process generate_plot_tx_data {
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    label "generate_plot_data"
    container "quay.io/kerimoff/coverage_plot:v6"


    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), file(sample_meta), file(coverage_parquet), file(usage_matrix_norm), file(tpm_matrix), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file
    path tx_gtf_file

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path("output_dir_*"),  emit: plot_data_directories
    path "${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log", emit: log_file

    script:
    vcf_sample_names_correction_bad_symbol = params.vcf_sample_names_correction ? "--vcf_sample_bad_symbol {$params.vcf_samples_old_string_part}" : ""
    vcf_sample_names_correction_replaced_string = params.vcf_sample_names_correction ? "--vcf_sample_replacement_symbol {$params.vcf_samples_new_string_part}" : ""
    """
    Rscript $projectDir/bin/generate_tx_plot_data.R \
        --qtl_group $qtl_group \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --phenotype_meta $phenotype_meta \
        --vcf_file $vcf_file \
        --coverage_parquet $coverage_parquet \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --tx_gtf_file $tx_gtf_file \
        --usage_matrix_norm $usage_matrix_norm \
        --tpm_matrix $tpm_matrix $vcf_sample_names_correction_bad_symbol $vcf_sample_names_correction_replaced_string

    cp .command.log ${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}

process generate_plot_txrev_data {
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    label "generate_plot_data"
    container "quay.io/kerimoff/coverage_plot:v6"


    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), file(sample_meta), file(coverage_parquet), file(usage_matrix_norm), file(tpm_matrix), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), file(susie_purity_filtered)
    path mane_transcript_gene_map
    path mane_gtf_file
    path txrev_gtf_file

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path("output_dir_*"),  emit: plot_data_directories
    path "${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log", emit: log_file

    script:
    vcf_sample_names_correction_bad_symbol = params.vcf_sample_names_correction ? "--vcf_sample_bad_symbol {$params.vcf_samples_old_string_part}" : ""
    vcf_sample_names_correction_replaced_string = params.vcf_sample_names_correction ? "--vcf_sample_replacement_symbol {$params.vcf_samples_new_string_part}" : ""


    """
    Rscript $projectDir/bin/generate_txrev_plot_data.R \
        --qtl_group $qtl_group \
        --finemap_susie $susie_purity_filtered \
        --sample_meta $sample_meta \
        --phenotype_meta $phenotype_meta \
        --vcf_file $vcf_file \
        --coverage_parquet $coverage_parquet \
        --mane_transcript_gene_map $mane_transcript_gene_map \
        --gtf_file $mane_gtf_file \
        --div_scaling_factors $scaling_factors \
        --txrev_gtf_file $txrev_gtf_file \
        --tpm_matrix $tpm_matrix \
        --usage_matrix_norm $usage_matrix_norm $vcf_sample_names_correction_bad_symbol $vcf_sample_names_correction_replaced_string


    cp .command.log ${dataset_id}_${qtl_group}_${quant_method}_${susie_purity_filtered.simpleName}.log
    """
}