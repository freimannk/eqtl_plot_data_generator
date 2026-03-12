#!/usr/bin/env nextflow

nextflow.enable.dsl=2


custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input files
 */ 

Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.dataset_id, row.study_id, row.quant_method, row.qtl_group, row.study_name, file(row.credible_sets_file), file(row.sample_meta), file(row.coverage_path), file(row.usage_matrix_norm), file(row.tpm_matrix), file(row.exon_summ_stats_files), file(row.all_summ_stats_files), file(row.pheno_meta), file(row.scaling_factors) ]}
    .branch {
        ge: it[2] == "ge"
        exon: it[2] == "exon"
        tx: it[2] == "tx"
        txrev: it[2] == "txrev"
        leafcutter: it[2] == "leafcutter"
    }
    .set { study_file_ch }

Channel.fromPath(params.studyFile)
  .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
  .splitCsv(header: true, sep: '\t', strip: true)
  .map{row -> [ row.dataset_id, file(row.vcf_file) ]}
  .distinct()
  .set { vcf_file_ch }

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

eQTL-Catalogue/qtlmap v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']        = 'eQTL-Catalogue/recap_plot'
summary['Pipeline Version']     = workflow.manifest.version
summary['Run Name']             = custom_runName ?: workflow.runName
summary['Study file']           = params.studyFile

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="

include { recap_plot_ge } from './workflows/recap_plot_ge_wf'
include { recap_plot_exon } from './workflows/recap_plot_exon_wf'
include { recap_plot_leafcutter } from './workflows/recap_plot_leafcutter_wf'
include { recap_plot_tx } from './workflows/recap_plot_tx_data_wf'
include { recap_txrev_plot_data } from './workflows/recap_plot_txrev_data_wf'
include { tabix_index } from './modules/utils'
include { generate_credible_sets_db } from './modules/utils'





workflow {
    tabix_index(vcf_file_ch)
    recap_plot_ge(study_file_ch.ge.combine(tabix_index.out, by: 0).distinct())
    recap_plot_exon(study_file_ch.exon.combine(tabix_index.out, by: 0).distinct())
    recap_plot_leafcutter(study_file_ch.leafcutter.combine(tabix_index.out, by: 0).distinct())
    recap_plot_tx(study_file_ch.tx.combine(tabix_index.out, by: 0).distinct())
    recap_txrev_plot_data(study_file_ch.txrev.combine(tabix_index.out, by: 0).distinct())

    concated_study_id_datasets_ch_oututs = recap_plot_ge.out.dataset_id_credible_set_tables
    .concat(recap_plot_exon.out.dataset_id_credible_set_tables)
    .concat(recap_plot_leafcutter.out.dataset_id_credible_set_tables)
    .concat(recap_plot_tx.out.dataset_id_credible_set_tables)
    .concat(recap_txrev_plot_data.out.dataset_id_credible_set_tables)

    concated_study_id_datasets_ch_oututs_collated = concated_study_id_datasets_ch_oututs.collate(5).collect().flatMap()
    study_id_datasets = concated_study_id_datasets_ch_oututs_collated.groupTuple(by: 0).map { study_id, datasets, quant_methods, qtl_groups, parquets -> [study_id, parquets] }
    generate_credible_sets_db(study_id_datasets)

}

workflow.onComplete {
    log.info "Pipeline Complete"
}
