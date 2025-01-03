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

process split_into_batches {
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    container "quay.io/kfkf33/polars"
    label "process_low"

    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), file(susie_purity_filtered), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index)

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), file(sample_meta), file(bigwig_path), file(usage_matrix_norm), file(exon_summ_stats_files), file(all_summ_stats_files), file(phenotype_meta), file(scaling_factors), file(vcf_file), file(vcf_file_index), emit: study_tsv_inputs_ch
    tuple val(dataset_id), val(quant_method), val(qtl_group), file("${dataset_id}_${qtl_group}_${quant_method}*.parquet"), emit: susie_batches

    script:
    """

    $projectDir/bin/prepare_batches.py \
        -n $all_summ_stats_files \
        -e $exon_summ_stats_files \
        -s $susie_purity_filtered \
        -p $phenotype_meta \
        -c ${params.chunk_size} \
        -o ${dataset_id}_${qtl_group}_${quant_method}
    """
}

process writeFileFromChannel {
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    label "process_low"

    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), val(files)

    output:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path("${dataset_id}_${qtl_group}_${quant_method}_dir_paths.txt")

    script:
    """
    for file in ${files.join(" ")};do
        echo \$file >> ${dataset_id}_${qtl_group}_${quant_method}_dir_paths.txt
    done
    """
}

process generate_sqlite {
    tag "${dataset_id}_${qtl_group}_${quant_method}"
    container "quay.io/kfkf33/eqtl_ploting:v2"
    publishDir "${params.outdir}/${dataset_id}_${qtl_group}_${quant_method}", mode: 'copy', overwrite: true, pattern: "*.sqlite"


    input:
    tuple val(dataset_id), val(quant_method), val(qtl_group), path(directory_paths)

    output:
    path("*.sqlite")

    script:
    """
    $projectDir/bin/generate_sqlites.py \
        -d $dataset_id \
        -s $directory_paths 
    """
}