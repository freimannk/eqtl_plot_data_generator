#!/usr/bin/env nextflow
nextflow.enable.dsl=2

Channel
    .fromPath(params.mane_transcript_gene_map, checkIfExists: true)
    .set { mane_transcript_gene_map_ch }

Channel
    .fromPath(params.mane_gtf_file, checkIfExists: true)
    .set { mane_gtf_file_ch }

include { generate_plot_exon_data } from  '../modules/generate_plot_data'
include { prepare_batches } from  '../modules/utils'
include { writeFileFromChannel } from '../modules/utils'
include { generate_dataset_ids_sqlites } from '../modules/utils'

workflow recap_plot_exon {
    take: 
    study_tsv_inputs_ch
    
    main:
    prepare_batches(study_tsv_inputs_ch)

    generate_plot_exon_data(
        prepare_batches.out.study_tsv_inputs_ch.combine(prepare_batches.out.susie_batches.transpose(), by:[0,1,2]),
        mane_transcript_gene_map_ch.collect(),
        mane_gtf_file_ch.collect()
    )
    groupped_ge_dir_paths = generate_plot_exon_data.out.plot_data_directories.groupTuple(by: [0,1,2])
    writeFileFromChannel_ge_input_ch = groupped_ge_dir_paths.map{tuple -> 
       def name = tuple[0]
        def method = tuple[1]
        def qtl_group = tuple[2]
        def directories = tuple[3].flatten()  // flatten the directories
        return [name, method, qtl_group, directories]}
    writeFileFromChannel(writeFileFromChannel_ge_input_ch)
    generate_dataset_ids_sqlites(writeFileFromChannel.out)
 
    emit:
    dataset_id_credible_set_tables = prepare_batches.out.dataset_id_credible_set_tables

}
