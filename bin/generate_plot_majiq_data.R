message(" ## Loading libraries: optparse")
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  make_option(c("-f", "--finemap_susie"), type="character", default=NULL,
              help="Purity filtered susie output. Tab separated file", metavar = "type"),
  make_option(c("-s", "--sample_meta"), type="character", default=NULL,
              help="Sample metadata file. Tab separated file", metavar = "type"),
  make_option(c("-q", "--qtl_group"), type="character", default=NULL,
              help="The selected qtl_group in the study", metavar = "type"),
  make_option(c("-p", "--phenotype_meta"), type="character", default=NULL,
              help="Phenotype metadata file. Tab separated file", metavar = "type"),
  make_option(c("-n", "--name_of_study"), type="character", default=NULL,
              help="Name of the study. Optional", metavar = "type"),
  make_option(c("-v", "--vcf_file"), type="character", default=NULL,
              help="TPM quantile TSV file with phenotype_id column", metavar = "type"),
  make_option(c("-c", "--coverage_parquet"), type="character", default=NULL,
              help="Path to the coverage parquet file", metavar = "type"),
  make_option(c("-m", "--mane_transcript_gene_map"), type="character", default=NULL,
              help="Path to the MANE transcripts of genes map file", metavar = "type"),
  make_option(c("-g", "--gtf_file"), type="character", default=NULL,
              help="Path to the GTF file to get exons of transcripts", metavar = "type"),
  make_option(c("--div_scaling_factors"), type="character", default=NULL,
              help="Path to the scaling_factors file", metavar = "type"),
  make_option(c("-u", "--usage_matrix_norm"), type="character", default=NULL,
              help="Path to the MAJIQ normalised PSI values matrix", metavar = "type"),
  make_option(c("--tpm_matrix"), type="character", default=NULL,
              help="Path to the MAJIQ raw PSI values matrix", metavar = "type"),
  make_option(c("--vcf_sample_bad_symbol"), type = "character", default = NULL,
              help = "Bad symbol to replace in VCF sample names (default: none)", metavar = "type"),
  make_option(c("--vcf_sample_replacement_symbol"), type = "character", default = NULL,
              help = "Replacement symbol for bad VCF sample names (default: none)", metavar = "type")
)

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))
library(seqminer)
library(stringr)
library(rlist)
library("wiggleplotr")
library(cowplot)
library(purrr)
library(arrow)
library(GenomicRanges)

susie_file_path = opt$f
sample_meta_path = opt$s
qtl_group_in = opt$q
phenotype_meta_path = opt$p
study_name = opt$n
vcf_file_path = opt$v
coverage_parquet_path = opt$c
mane_transcript_gene_map_path = opt$m
gtf_file_path = opt$g
scaling_factors_path = opt$div_scaling_factors
tpm_matrix_path = opt$tpm_matrix
usage_matrix_path = opt$usage_matrix_norm
vcf_sample_bad_symbol = opt$vcf_sample_bad_symbol
vcf_sample_replacement_symbol = opt$vcf_sample_replacement_symbol

conf.level <- 0.95
ci.value <- -qnorm( ( 1 - conf.level ) / 2 )

make_transcript_exon_granges <- function(gff, transcript_ids) {
  exon_list <- list()
  for (transcript_id in transcript_ids) {
    transcript_exons_temp <- gff[(base::gsub("\\..*","",SummarizedExperiment::elementMetadata(gff)[,"transcript_id"]) == transcript_id)]
    gene_id <- transcript_exons_temp$gene_name[1]
    exon_list[[paste0("GENE:", gene_id, ":", transcript_id)]] <- transcript_exons_temp
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

make_transcript_exon_granges_ccds <- function(gff, transcript_ids) {
  exon_list <- list()
  for (transcript_id in transcript_ids) {
    transcript_exons_temp <- gff[(base::gsub("\\..*","",SummarizedExperiment::elementMetadata(gff)[,"transcript_id"]) == transcript_id & !is.na(SummarizedExperiment::elementMetadata(gff)[,"ccds_id"]))]
    gene_id <- transcript_exons_temp$gene_name[1]
    exon_list[[paste0("GENE:", gene_id, ":", transcript_id)]] <- transcript_exons_temp
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

make_pseudo_exons <- function(df_introns, chr, ps_exon_len = 10){
  split_coords <- do.call(rbind, strsplit(df_introns$junction_coord, '-'))

  df_introns$intron_start <- as.integer(split_coords[, 1])
  df_introns$intron_end <- as.integer(split_coords[, 2])
  uniq_intron_ids <- df_introns %>% dplyr::pull(intron_id)

  introns_oi_exon1 <- df_introns %>% dplyr::mutate(ps_exon_start = intron_start - ps_exon_len, ps_exon_end = intron_start)
  introns_oi_exon2 <- df_introns %>% dplyr::mutate(ps_exon_start = intron_end, ps_exon_end = intron_end + ps_exon_len)
  df_introns <- BiocGenerics::rbind(introns_oi_exon1, introns_oi_exon2)
  exon_list <- list()
  for (int_id in uniq_intron_ids) {
    df_introns_sub <- df_introns %>% dplyr::filter(intron_id == int_id)
    if (all(grepl('intron', df_introns_sub$intron_id))) {
      df_introns_sub <- df_introns_sub[1,]
      exon_list[[int_id]] <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(start = df_introns_sub$intron_start, end = df_introns_sub$intron_end),
        strand = df_introns_sub$strand,
        mcols = data.frame(intron_id = df_introns_sub$intron_id,
                           group_id = df_introns_sub$quant_id,
                           event_id = df_introns_sub$event_id,
                           gene_id = df_introns_sub$gene_id,
        )
    )
    } else {
      exon_list[[int_id]] <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(start = df_introns_sub$ps_exon_start, end = df_introns_sub$ps_exon_end),
        strand = df_introns_sub$strand,
        mcols = data.frame(intron_id = df_introns_sub$intron_id,
                           group_id = df_introns_sub$quant_id,
                           event_id = df_introns_sub$event_id,
                           gene_id = df_introns_sub$gene_id)
    )
    }
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

prepareTranscriptStructureForPlotting <- function(exon_ranges, cds_ranges, transcript_annotations){
  #Combine exon_ranges and cds_ranges into a single data.frame that also contains transcript rank

  #Convert exon ranges into data.frame and add transcript rank
  exons_df = purrr::map_df(exon_ranges, data.frame, .id = "transcript_id")
  exons_df = dplyr::mutate(exons_df, transcript_rank = as.numeric(factor(exons_df$transcript_id)), type = "")
  transcript_rank = unique(exons_df[,c("transcript_id", "transcript_rank", "type")])

  #Convert CDS ranges into a data.frame
  cds_df = purrr::map_df(cds_ranges, data.frame, .id = "transcript_id")
  cds_df = dplyr::left_join(cds_df, transcript_rank, by = "transcript_id") #Add matching transcript rank

  #Join exons and cdss together
  exons_df = dplyr::mutate(exons_df, feature_type = "exon")
  cds_df = dplyr::mutate(cds_df, feature_type = "cds")
  transcript_struct = rbind(exons_df, cds_df)

  #Add transcript label to transcript structure
  transcript_struct = dplyr::left_join(transcript_struct, transcript_annotations, by = "transcript_id")
  return(transcript_struct)
}

format_trait_matrix <- function(trait_matrix_oi, column_name,value_type_id) {
  trait_matrix_oi <- tibble::column_to_rownames(.data = trait_matrix_oi,var = "phenotype_id")
  trait_matrix_oi <- trait_matrix_oi %>% base::t() %>%
    GenomicRanges::as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_id")
  trait_matrix_oi <- trait_matrix_oi %>%
    tidyr::pivot_longer(cols = -sample_id, names_to=value_type_id, values_to = column_name)
  return(trait_matrix_oi)
}

read_and_filter_parquet <- function(file_list, variant_to_match, phenotype_id) {
  if (!is.vector(file_list) || length(file_list) == 0) {
    stop("file_list must be a non-empty vector of file names.")
  }
  total_files <- length(file_list)
  i <- 1
  for (file_name in file_list) {
    if (!file.exists(file_name)) {
      warning(paste("File not found:", file_name))
      next
    }

    dataset <- open_dataset(file_name)
    filtered_data <- dataset %>%
      filter(variant == variant_to_match & gene_id == phenotype_id) %>%
      collect()

    if (nrow(filtered_data) > 0) {
      return(filtered_data)
    }
    if (i == total_files) {
      return(filtered_data)
    }
    i <- i +1
  }
}

ir_df <- read.csv(phenotype_meta_path, sep='\t', comment.char = "#")
highest_pip_vars_per_cs <- read_parquet(susie_file_path)
highest_pip_vars_per_cs$nominal_exon_cc_path <- lapply(highest_pip_vars_per_cs$nominal_exon_cc_path, function(x) {
  strsplit(gsub("[\\[\\]'\" ]", "", x), ",")
})
highest_pip_vars_per_cs$nominal_cc_path <- lapply(highest_pip_vars_per_cs$nominal_cc_path, function(x) {
  strsplit(gsub("[\\[\\]'\" ]", "", x), ",")
})
sample_metadata <- readr::read_tsv(sample_meta_path)
scaling_factor_data <- readr::read_tsv(scaling_factors_path)
gtf_ref <- rtracklayer::import(gtf_file_path,
                               colnames = c("type", "gene_id", "gene_name", "gene_biotype",
                                            "transcript_id", "transcript_name","transcript_biotype", "exon_number", "exon_id", "ccds_id"),
                               feature.type = c("exon"))

mane_transcript_gene_map <- readr::read_tsv(mane_transcript_gene_map_path)
message(" ## Reading TMP matrix")

tpm_exp_df <- read_parquet(tpm_matrix_path)
colnames(tpm_exp_df)[-1] <- gsub('\\.', '-', colnames(tpm_exp_df)[-1])

tpm_exp_df <- tpm_exp_df %>%
  mutate(
    phenotype_id = str_replace(phenotype_id, '^gene:', ''),
    quant_id = str_replace(str_replace(phenotype_id, ':\\d+-\\d+$', ''), '^gene:', ''),
    junction_coord = str_extract(phenotype_id, '\\d+-\\d+$'),
    gene_id = str_extract(quant_id, '(?<=gene:)[^:]+'),
         # ToDo: check if strand info is important here
    strand = rep('+', nrow(tpm_exp_df))

  ) %>%
  filter(quant_id != "") %>%
  distinct(quant_id, junction_coord, .keep_all = TRUE)

message(" ## Reading normalised usage matrix")

norm_exp_df <- read_parquet(usage_matrix_path, sep='\t')
colnames(norm_exp_df)[-1] <- gsub('\\.', '-', colnames(norm_exp_df)[-1])

norm_exp_df <- norm_exp_df %>%
  mutate(
    phenotype_id = str_replace(phenotype_id, '^gene:', ''),
    quant_id = str_replace(str_replace(phenotype_id, ':\\d+-\\d+$', ''), '^gene:', ''),
    junction_coord = str_extract(phenotype_id, '\\d+-\\d+$'),
    gene_id = str_extract(quant_id, '(?<=gene:)[^:]+'),
    # ToDo: check if strand info is important here
    strand = rep('+', nrow(norm_exp_df))
  ) %>%
  filter(quant_id != "") %>%
  distinct(quant_id, junction_coord, .keep_all = TRUE)

highest_pip_vars_per_cs <- highest_pip_vars_per_cs %>%
    mutate(quant_id = paste0(sapply(strsplit(molecular_trait_id, ':'), function(x) paste(head(x, -1), collapse = ':'))))

variant_regions_vcf <- highest_pip_vars_per_cs %>%
  dplyr::mutate(region = str_remove(region, '^chr')) %>%
  dplyr::mutate(region = str_replace_all(region, ':-\\d+', ':1'))

message(" ## Reading all variants from VCF_file")
snps_all <- seqminer::tabix.read.table(vcf_file_path, variant_regions_vcf$region)

# Apply replacement if vcf sample_bad_symbol is not empty
if (!is.null(vcf_sample_bad_symbol)) {
  # Replace bad symbol with the replacement symbol
  names(snps_all) <- gsub(pattern = vcf_sample_bad_symbol,
                          replacement = vcf_sample_replacement_symbol,
                          x = names(snps_all),
                          fixed = TRUE)
}
message(paste0("## Number of credible sets to plot: ", as.character(length(unique(highest_pip_vars_per_cs$cs_id)))))

for (index in 1:nrow(highest_pip_vars_per_cs)) {
  ss_oi <- highest_pip_vars_per_cs[index,]

  variant_regions_vcf <- ss_oi %>%
    dplyr::mutate(region = str_remove(region, '^chr')) %>%
    dplyr::mutate(region = str_replace_all(region, ':-\\d+', ':1'))

  snps_filt <- snps_all %>%
    dplyr::filter(ID %in% variant_regions_vcf$variant) %>%
    dplyr::arrange(CHROM, POS)
    if (nrow(snps_filt) == 0) {
      next
    }

  var_genotype <- snps_filt %>%
    dplyr::select(-c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")) %>%
    base::t() %>%
    BiocGenerics::as.data.frame() %>%
    dplyr::rename("GT_DS" = "V1") %>%
    dplyr::mutate(GT = gsub(pattern = "\\:.*", replacement = "", x = GT_DS)) %>%
    dplyr::mutate(REF = gsub(pattern = "\\|.*", replacement = "", x = GT)) %>%
    dplyr::mutate(ALT = gsub(pattern = ".*\\|", replacement = "", x = GT)) %>%
    dplyr::mutate(DS = as.numeric(REF) + as.numeric(ALT)) %>%
    dplyr::mutate(genotype_id = BiocGenerics::rownames(.)) %>%
    dplyr::mutate(genotype_id = gsub(pattern = "\\.", replacement = "-", x = genotype_id))

  sample_meta_clean = sample_metadata %>%
    dplyr::mutate(
        rna_qc_passed = as.logical(rna_qc_passed),
        genotype_qc_passed = as.logical(genotype_qc_passed)
    ) %>%
    dplyr::filter(rna_qc_passed, genotype_qc_passed) %>%
    dplyr::select(sample_id, genotype_id, qtl_group)

  var_genotype <- sample_meta_clean %>%
    dplyr::left_join(var_genotype, by = "genotype_id")

  track_data_study <- var_genotype %>%
    dplyr::select(sample_id, qtl_group, DS) %>%
    dplyr::left_join(scaling_factor_data) %>%
    dplyr::rename(scaling_factor = scaling_factors) %>%
    dplyr::mutate(track_id = qtl_group) %>%
    dplyr::mutate(colour_group = as.character(DS)) %>%
    dplyr::mutate(bigWig = "") %>%
    dplyr::select(sample_id, scaling_factor, bigWig, track_id, colour_group, qtl_group) %>%
    dplyr::filter(qtl_group==qtl_group_in)

  message(" ## Extracting junctions")
  tryCatch ({
    cluster_introns <- ir_df %>% dplyr::filter(quant_id %in% ss_oi$quant_id)

    message(cluster_introns$quant_id[1])
    message(ss_oi$cs_id[1])
    message(ss_oi$variant[1])
    chr <- str_remove(unlist(strsplit(ss_oi[1,]$region, ':'))[1], '^chr')
    ps_exons <- make_pseudo_exons(df_introns = cluster_introns, chr = chr)
    ps_exons_cdss <- make_pseudo_exons(df_introns = cluster_introns, chr = chr)

    error = function(e) {
        message(" ## Problem with generating pseudo exons")
        message(e)
        message(cluster_introns)
        return(NULL)
        }
    })

  exons_to_plot = ps_exons
  exons_cdss_to_plot = ps_exons_cdss

  message(" ## Extracting MANE transcripts")

  MANE_transcript_oi <- mane_transcript_gene_map %>%
      dplyr::filter(gene_id %in% ss_oi$gene_id) %>%
      dplyr::pull(transcript_id)

  if (length(MANE_transcript_oi) > 0) {
    mane_transcript_exons <-  make_transcript_exon_granges(gff = gtf_ref, transcript_ids = MANE_transcript_oi)
    exons_to_plot <- append(exons_to_plot, mane_transcript_exons)
    exons_cdss_to_plot <- append(exons_cdss_to_plot, mane_transcript_exons)
  }
  message(" ## Extracting exon-level summstats")

  nom_exon_cc_sumstats_variant_phenotype_id <- read_and_filter_parquet(
    file_list = ss_oi$nominal_exon_cc_path[[1]],
    variant_to_match = ss_oi$variant,
    phenotype_id=ss_oi$gene_id
  )
  if(!is.null(nom_exon_cc_sumstats_variant_phenotype_id)) {
    # Extract the QTLs of exons according to gene and variant of interest
    nom_exon_cc_sumstats_filt <- nom_exon_cc_sumstats_variant_phenotype_id %>%
      dplyr::group_by(molecular_trait_id, variant) %>%
      dplyr::filter(!is.na(rsid) | all(is.na(rsid))) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(exon_end = as.numeric(gsub(pattern = ".*\\_", replacement = "", x = molecular_trait_id))) %>%
      dplyr::mutate(exon_start = gsub(pattern = "_[^_]+$", replacement = "", x = molecular_trait_id)) %>%
      dplyr::mutate(exon_start = as.numeric(gsub(pattern = ".*\\_", replacement = "", x = exon_start))) %>%
      dplyr::mutate(exon_center = round((exon_start + exon_end) / 2)) %>%
      dplyr::mutate(exon_length = abs(exon_start - exon_end) + 1) %>%
      dplyr::mutate(exon_row_num = dplyr::row_number()) %>%
      dplyr::mutate(se_top = beta + se) %>%
      dplyr::mutate(se_bottom = beta - se) %>%
      dplyr::mutate(interval = ci.value * se) %>%
      dplyr::mutate(p_fdr = p.adjust(pvalue, method = "fdr"))

    if (nrow(nom_exon_cc_sumstats_filt) > 0) {
      nom_exon_granges <- list()
      nom_exon_granges[[paste0("GENE:",ss_oi$gene_name)]] = GenomicRanges::GRanges(
        seqnames = nom_exon_cc_sumstats_filt$chromosome,
        ranges = IRanges::IRanges(start = nom_exon_cc_sumstats_filt$exon_start, end = nom_exon_cc_sumstats_filt$exon_end),
        strand = ifelse(test = ss_oi$strand == 1, yes = "+", no = "-"),
        mcols = data.frame(exon_id = nom_exon_cc_sumstats_filt$molecular_trait_id,
                           gene_id = nom_exon_cc_sumstats_filt$gene_id))

      exons_to_plot <- append(exons_to_plot, nom_exon_granges)
    }
  }

  message(" ## Extracting coverage data")
  coverage_data_list = tryCatch(wiggleplotr::extractCoverageData(exons = exons_to_plot,
                                                                 cdss = exons_cdss_to_plot, # Does var will be default NULL if not stated?
                                                                 plot_fraction = 0.2,
                                                                 track_data = track_data_study,
                                                                 coverage_ranges_pq = coverage_parquet_path),

  error = function(e) {
      message(" ## Problem with generating coverage_data wiggleplotr")
      message(e)
    })
  if (!exists("coverage_data_list")) {
      message(" ERROR: !exists")
      next
  }
  if (all(is.na(coverage_data_list)) | length(coverage_data_list) == 0) {
      message(" ERROR: is.na(coverage_data_list) | length(coverage_data_list)")
      next
  }

  tx_structure_df = prepareTranscriptStructureForPlotting(exon_ranges = coverage_data_list$tx_annotations$exon_ranges,
                                                    cds_ranges = coverage_data_list$tx_annotations$cds_ranges,
                                                    transcript_annotations = coverage_data_list$plotting_annotations)
  # BOX PLOTS
  tpm_exp_df_oi <- tpm_exp_df %>% dplyr::filter(quant_id %in% cluster_introns$quant_id)
  if (nrow(tpm_exp_df_oi) == 0) {
    message('Skipping box plot - tpm info missing')
    next
  }
  tpm_exp_df_oi <- format_trait_matrix(tpm_exp_df_oi, "tpm_exp","intron_id")

  norm_exp_df_oi <- norm_exp_df %>% dplyr::filter(quant_id %in% cluster_introns$quant_id)
  if (nrow(norm_exp_df_oi) == 0) {
    message('Skipping box plot - norm info missing')
    next
  }
  norm_exp_df_oi <- format_trait_matrix(norm_exp_df_oi, "norm_exp","intron_id")

  track_data_study_box <- track_data_study %>%
    dplyr::mutate(genotype_text = as.factor(colour_group)) %>%
    dplyr::filter(qtl_group %in% qtl_group_in) %>%
    dplyr::mutate(condition_name = qtl_group) %>%
    dplyr::mutate(gene_name = ss_oi$gene_name) %>%
    dplyr::mutate(snp_id = ss_oi$variant)

  track_data_study_box <- norm_exp_df_oi %>%
    dplyr::left_join(track_data_study_box, by = "sample_id")
  track_data_study_box <- track_data_study_box %>%
    dplyr::left_join(tpm_exp_df_oi, by = c("sample_id", "intron_id"))

  nom_cc_sumstats_variant_phenotype_id <- read_and_filter_parquet(
    file_list = ss_oi$nominal_cc_path[[1]],
    variant_to_match = ss_oi$variant,
    phenotype_id=ss_oi$gene_id
  )

  # Keep only 1 rsid per variant per molecular_trait_id
  nom_cc_sumstats <- nom_cc_sumstats_variant_phenotype_id %>%
    dplyr::group_by(molecular_trait_id, variant) %>%
    dplyr::filter(!is.na(rsid) | all(is.na(rsid))) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  nom_cc_sumstats_filt <- nom_cc_sumstats %>%
    dplyr::filter(variant %in% ss_oi$variant) %>%
    dplyr::select(molecular_trait_id, pvalue, beta, se, maf) %>%
    dplyr::rename(intron_id = molecular_trait_id)

  track_data_study_box_wrap <- track_data_study_box %>%
    dplyr::left_join(nom_cc_sumstats_filt, by = "intron_id")

  track_data_study_box_wrap_for_RDS <- track_data_study_box_wrap %>%
    dplyr::select(genotype_text, tpm_exp, norm_exp, intron_id, pvalue, beta, se, snp_id, maf)

  print(track_data_study_box_wrap_for_RDS)
  track_data_study_box_wrap_for_RDS <- track_data_study_box_wrap_for_RDS[sample(nrow(track_data_study_box_wrap_for_RDS)),]
  tx_str_df <- tx_structure_df %>% dplyr::mutate(limit_max = max(coverage_data_list$limits))

  message('## Writing plot data to .pq files')
  signal_name <- gsub(pattern = "&", replacement = "\\&", x = paste0(gsub(pattern = ":", replacement = "_", x = ss_oi$molecular_trait_id), "__", ss_oi$variant, "__", ss_oi$gene_id))

  output_path <- paste0("output_dir_", signal_name)
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  arrow::write_parquet(coverage_data_list$coverage_df, file.path(output_path, paste0("coverage_df_", signal_name, ".parquet")))
  arrow::write_parquet(tx_str_df, file.path(output_path, paste0("tx_str_", signal_name, ".parquet")))
  arrow::write_parquet(track_data_study_box_wrap_for_RDS, file.path(output_path, paste0("box_plot_df_", signal_name, ".parquet")))
  arrow::write_parquet(ss_oi, file.path(output_path, paste0("ss_oi_df_", signal_name, ".parquet")))
  nom_exon_cc_sumstats_filt <- nom_exon_cc_sumstats_filt %>%
    mutate(rsid = stringr::str_trim(rsid))
  arrow::write_parquet(nom_exon_cc_sumstats_filt, file.path(output_path, paste0("nom_exon_cc_", signal_name, ".parquet")))
}


