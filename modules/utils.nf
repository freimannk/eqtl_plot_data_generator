process tabix_index {
    tag "${dataset_id}"
    storeDir "${projectDir}/vcfs"
    container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"

    input:
    tuple val(dataset_id), file(vcf_file)

    output:
    tuple val(dataset_id), file(vcf_file), file("${vcf_file}.tbi")

    script:
    """
    tabix $vcf_file
    """
}

process prepare_batches {
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    container "quay.io/kfkf33/polars"
    label "process_low"

    input:
    tuple val(dataset_id), val(study_id), val(quant_method), val(qtl_group), val(study_name), file(susie_purity_filtered), file(sample_meta), file(coverage_parquet), file(usage_matrix_norm),file(tpm_matrix), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index)

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), file(sample_meta), file(coverage_parquet), file(usage_matrix_norm) ,file(tpm_matrix), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), emit: study_tsv_inputs_ch
    tuple val(dataset_id), val(quant_method), val(qtl_group), file("${dataset_id}_${qtl_group}_${quant_method}*.parquet"), emit: susie_batches
    tuple val(study_id), val(dataset_id), val(quant_method), val(qtl_group), file("${study_id}_${dataset_id}_${qtl_group}_${quant_method}.parquet"), emit: dataset_id_credible_set_tables

    script:
    """
    $projectDir/bin/prepare_batches.py \
        -a $study_name \
        -n $all_summ_stats_files \
        -e $exon_summ_stats_files \
        -s $susie_purity_filtered \
        -p $phenotype_meta \
        -c ${params.chunk_size} \
        -d $dataset_id \
        -i $study_id \
        -q $quant_method \
        -g $qtl_group
        
    """
}

process writeFileFromChannel {
    tag "${dataset_id}_${qtl_group}_${quant_method}"

    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), val(files)

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path("${dataset_id}_${qtl_group}_${quant_method}_dir_paths.txt")

    script:
    """
    echo "${files.join(" ")}" | grep -o '[^ ]*' > ${dataset_id}_${qtl_group}_${quant_method}_dir_paths.txt
    """
}

process generate_dataset_ids_sqlites {
    tag "${dataset_id}"
    container "quay.io/kfkf33/eqtl_ploting:v2"
    publishDir "${params.outdir}/${dataset_id}_${qtl_group}_${quant_method}", mode: 'copy', overwrite: true, pattern: "*.sqlite"


    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path(directory_paths)

    output:
    path("*.sqlite")

    script:
    """
    $projectDir/bin/generate_datasets_sqlites.py \
        -d $dataset_id \
        -s $directory_paths 
    """
}

process generate_credible_sets_db {
    container = 'quay.io/kfkf33/duckdb_sqlite_with_path:v1'
    publishDir "${params.outdir}/${study_id}", mode: 'copy', overwrite: true, pattern: "*.sqlite"


    input:
    tuple val(study_id), val(datasets_cs)
    

    output:
    path("*.sqlite")

    script:
    """
    generate_credible_set_db.py -f ${datasets_cs.join(' ')} -o ${study_id}.sqlite

    """
}

process convert_parquet_format {
    container = 'quay.io/kfkf33/duckdb_sqlite_with_path:v1'

    input:
    tuple val(dataset_id), val(study_id), val(quant_method), val(qtl_group), val(study_name), file(susie_purity_filtered), file(sample_meta), file(coverage_parquet), file(usage_matrix_norm),file(tpm_matrix), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index)
    

    output:
    tuple val(dataset_id), val(study_id), val(quant_method), val(qtl_group), val(study_name), file(susie_purity_filtered), file(sample_meta), file(coverage_parquet), file("${usage_matrix_norm.simpleName}.parquet"),file("${usage_matrix_norm.simpleName}_TPM.parquet"), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), emit: converted_study_input

    script:
    """
    convert_to_parquet.py -i $usage_matrix_norm -t $tpm_matrix -o ${usage_matrix_norm.simpleName}

    """
}