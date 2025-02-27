message(" ## Loading libraries: optparse")
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  make_option(c("-f", "--finemap_susie"), type="character", default=NULL,
              help="Purity filtered susie output. Tab separated file", metavar = "type"),
  make_option(c("-s", "--sample_meta"), type="character", default=NULL,
              help="Sample metadata file. Tab separated file", metavar = "type"),
  make_option(c("-p", "--phenotype_meta"), type="character", default=NULL,
              help="Phenotype metadata file. Tab separated file", metavar = "type"),
  make_option(c("-q", "--qtl_group"), type="character", default=NULL,
              help="The selected qtl_group in the study", metavar = "type"),
  make_option(c("-n", "--name_of_study"), type="character", default=NULL,
              help="Name of the study. Optional", metavar = "type"),
  make_option(c("-v", "--vcf_file"), type="character", default=NULL,
              help="TPM quantile TSV file with phenotype_id column", metavar = "type"),
  make_option(c("-b", "--bigwig_path"), type="character", default=NULL,
              help="Path to the bigwig files", metavar = "type"),
  make_option(c("-m", "--mane_transcript_gene_map"), type="character", default=NULL,
              help="Path to the MANE transcripts of genes map file", metavar = "type"),
  make_option(c("-g", "--gtf_file"), type="character", default=NULL,
              help="Path to the GTF file to get exons of transcripts", metavar = "type"),
  make_option(c("-u", "--usage_matrix_norm"), type="character", default=NULL,
              help="Path to the normalised usage matrix", metavar = "type"),
  make_option(c("--div_scaling_factors"), type="character", default=NULL,
              help="Path to the scaling_factors file", metavar = "type")
)

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("seqminer"))
suppressPackageStartupMessages(library("arrow"))


make_transcript_exon_granges <- function(gff, transcript_ids) {
  exon_list <- list()
  for (transcript_id in transcript_ids) {
    transcript_exons_temp <- gff[(base::gsub("\\..*","",SummarizedExperiment::elementMetadata(gff)[,"transcript_id"]) == transcript_id)]
    gene_id = transcript_exons_temp$gene_name[1]
    exon_list[[paste0("GENE:", gene_id, ":", transcript_id)]] <- transcript_exons_temp
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

make_transcript_exon_granges_ccds <- function(gff, transcript_ids) {
  exon_list <- list()
  for (transcript_id in transcript_ids) {
    transcript_exons_temp <- gff[(base::gsub("\\..*","",SummarizedExperiment::elementMetadata(gff)[,"transcript_id"]) == transcript_id & !is.na(SummarizedExperiment::elementMetadata(gff)[,"ccds_id"]))]
    gene_id = transcript_exons_temp$gene_name[1]
    exon_list[[paste0("GENE:", gene_id, ":", transcript_id)]] <- transcript_exons_temp
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

make_pseudo_exons <- function(df_introns, ps_exon_len = 50){
  uniq_intron_ids <- df_introns %>% dplyr::pull(phenotype_id)
  introns_oi_exon1 <- df_introns %>% dplyr::mutate(ps_exon_start = intron_start - ps_exon_len, ps_exon_end = intron_start)
  introns_oi_exon2 <- df_introns %>% dplyr::mutate(ps_exon_start = intron_end, ps_exon_end = intron_end + ps_exon_len)
  df_introns <- BiocGenerics::rbind(introns_oi_exon1, introns_oi_exon2) %>% dplyr::arrange(phenotype_id)
  exon_list <- list()
  for (intron_id in uniq_intron_ids) {
    df_introns_sub <- df_introns %>% dplyr::filter(phenotype_id == intron_id)
    exon_list[[intron_id]] <- GenomicRanges::GRanges(
      seqnames = df_introns_sub$chromosome,
      ranges = IRanges::IRanges(start = df_introns_sub$ps_exon_start, end = df_introns_sub$ps_exon_end),
      strand = ifelse(test = df_introns_sub$strand == 1, yes = "+", no = "-"),
      mcols = data.frame(intron_id = df_introns_sub$phenotype_id, 
                         group_id = df_introns_sub$group_id, 
                         gene_id = df_introns_sub$gene_id)
    )
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

read_and_filter_parquet <- function(file_list, variant_to_match, phenotype_id) {
  if (!is.vector(file_list) || length(file_list) == 0) {
    stop("file_list must be a non-empty vector of file names.")
  }
  total_files <- length(file_list)
  i=1
  
  for (file_name in file_list) {
    if (!file.exists(file_name)) {
      warning(paste("File not found:", file_name))
      next
    }
    
    dataset <- open_dataset(file_name)
    
    filtered_data <- dataset %>%
      filter(variant == variant_to_match & molecular_trait_object_id == phenotype_id) %>%
      collect()
    
    if (nrow(filtered_data) > 0) {
      return(filtered_data)
    }
    if (i == total_files) {
      return(filtered_data)
    }
    i = i +1
  }
}

#Debugging
if (FALSE) {
  opt = list()
  opt$n = "Lepik_2017"
  opt$f = "/Users/kerimov/Work/temp_files/debug/c96dadbf3eb6600783be4ab63b8905/single_batch_debug_mode.tsv"
  opt$s = "/Users/kerimov/Work/temp_files/debug/c96dadbf3eb6600783be4ab63b8905/Lepik_2017.tsv"
  opt$p = "/Users/kerimov/Work/temp_files/debug/c96dadbf3eb6600783be4ab63b8905/leafcutter_metadata.txt.gz"
  opt$q = "blood"
  opt$v = "/Users/kerimov/Work/temp_files/debug/c96dadbf3eb6600783be4ab63b8905/Lepik_2017_GRCh38.filtered.vcf.gz"
  opt$b = "/Users/kerimov/Work/temp_files/debug/c96dadbf3eb6600783be4ab63b8905/bigwig/"
  opt$m = "/Users/kerimov/Work/temp_files/debug/c96dadbf3eb6600783be4ab63b8905/MANE_transcript_gene_map.txt"
  opt$g = "/Users/kerimov/Work/temp_files/debug/c96dadbf3eb6600783be4ab63b8905/Homo_sapiens.GRCh38.105.gtf"
  opt$div_scaling_factors = "/Users/kerimov/Work/temp_files/debug/c96dadbf3eb6600783be4ab63b8905/Lepik_2017.blood.scaling_factors.tsv.gz"
  opt$u = "/Users/kerimov/Work/temp_files/debug/c96dadbf3eb6600783be4ab63b8905/Lepik_2017.blood.leafcutter_CPM_norm.tsv.gz"
  index = 1
}

susie_file_path = opt$f
sample_meta_path = opt$s
phenotype_meta_path = opt$p
qtl_group_in = opt$q
study_name = opt$n
vcf_file_path = opt$v
bigwig_files_path = opt$b
mane_transcript_gene_map_file = opt$m
gtf_file_path = opt$g
norm_usage_matrix_path = opt$u
scaling_factors_path = opt$div_scaling_factors

message("######### Options: ######### ")
message("######### Working Directory  : ", getwd())
message("######### qtl_group          : ", qtl_group_in)
message("######### susie_file_path    : ", susie_file_path)
message("######### sample_meta_path   : ", sample_meta_path)
message("######### phenotype_meta_path: ", phenotype_meta_path)
message("######### opt_study_name     : ", study_name)
message("######### vcf_file_path      : ", vcf_file_path)
message("######### bigwig_files_path  : ", bigwig_files_path)
message("######### mane_map_file_path : ", mane_transcript_gene_map_file)
message("######### gtf_file_path      : ", gtf_file_path)
message("######### scaling_fct_path   : ", scaling_factors_path)
message("######### norm_usage_matrix  : ", norm_usage_matrix_path)


################## Global variable definitions ################
conf.level = 0.95
ci.value <- -qnorm( ( 1 - conf.level ) / 2 )

time_here <- function(prev_time, message_text = "Time in this point: "){
  message(message_text, Sys.time() - prev_time)
  return(Sys.time())
}
###############################################################
start_time <- Sys.time()
message(" ## Reading GTF file")
gtf_ref <- rtracklayer::import(gtf_file_path, 
                               colnames = c("type", "gene_id", "gene_name", "gene_biotype", 
                                            "transcript_id", "transcript_name","transcript_biotype", "exon_number", "exon_id", "ccds_id"),
                               feature.type = c("exon"))
message(" ## Reading sample_metadata file")
sample_metadata <- readr::read_tsv(sample_meta_path)

message(" ## Reading mane_transcript_gene_map file")
mane_transcript_gene_map <- readr::read_tsv(mane_transcript_gene_map_file)

message(" ## Reading susie_purity_filtered file")
highest_pip_vars_per_cs <- read_parquet(susie_file_path) 
highest_pip_vars_per_cs$nominal_exon_cc_path <- lapply(highest_pip_vars_per_cs$nominal_exon_cc_path, function(x) {
  strsplit(gsub("[\\[\\]'\" ]", "", x), ",")
})
highest_pip_vars_per_cs$nominal_cc_path <- lapply(highest_pip_vars_per_cs$nominal_cc_path, function(x) {
  strsplit(gsub("[\\[\\]'\" ]", "", x), ",")
})
message(" ## Reading normalised usage matrix")
norm_exp_df <- readr::read_tsv(norm_usage_matrix_path)

message(" ## Reading leafcutter metadata file")
leafcutter_metadata <- readr::read_tsv(phenotype_meta_path, col_types = "cccccddiccidddddddd") 

message(" ## Reading scaling_factors file")
scaling_factor_data <- readr::read_tsv(scaling_factors_path, col_types = "cd") 

start_time <- time_here(prev_time = start_time, message_text = " >> Read input TSVs in: ")
if (is.null(study_name)) { 
  assertthat::has_name(sample_metadata, "study" )
  study_name <- sample_metadata$study[1] 
}

if(assertthat::assert_that(all(!is.na(leafcutter_metadata$gene_id) & all(!is.na(leafcutter_metadata$gene_name))), 
                           msg = "All gene_id's and gene_name's in leafcutter_metadata should be non-NA")) {
  message("All the phenotypes in leafcutter_metadata has properly assigned to a certain gene!")
}

variant_regions_vcf <- highest_pip_vars_per_cs %>% 
  dplyr::select(variant, chromosome, position) %>% 
  dplyr::mutate(region = paste0(chromosome,":", position, "-", position))

message(" ## Reading all variants from VCF_file")
snps_all <- seqminer::tabix.read.table(vcf_file_path, variant_regions_vcf$region)
if (study_name == "Steinberg_2020") {
  names(snps_all) <- gsub(pattern = ".", replacement = ":", x = names(snps_all), fixed = T)
}
if (study_name == "Quach_2016") {
  names(snps_all) <- gsub(pattern = ".", replacement = "@", x = names(snps_all), fixed = T)
}
if (study_name %in% c("Schmiedel_2018", "Bossini-Castillo_2019", "ROSMAP", "iPSCORE")) {
  names(snps_all) <- gsub(pattern = "X", replacement = "", x = names(snps_all), fixed = T)  
}
if (study_name %in% c("iPSCORE")) {
  names(snps_all) <- gsub(pattern = ".", replacement = "-", x = names(snps_all), fixed = T)  
}
message(" ## Reading all variants from VCF_file complete")

start_time <- time_here(prev_time = start_time, message_text = " >> seqminer tabix took: ")
message(" ## Will plot batch of ", nrow(highest_pip_vars_per_cs)," highest pip per credible set signals.")
for (index in 1:nrow(highest_pip_vars_per_cs)) {
  start_time = Sys.time()
  ss_oi = highest_pip_vars_per_cs[index,]
  message("index: ", index, ", intron_id: ", ss_oi$molecular_trait_id, ", variant: ", ss_oi$variant)
  
  # get all the intons in leafcutter cluster
  cluster_introns <- leafcutter_metadata %>% dplyr::filter(group_id %in% ss_oi$group_id)
  ps_exons <- make_pseudo_exons(df_introns = cluster_introns)
  ps_exons_cdss <- make_pseudo_exons(leafcutter_metadata %>% dplyr::filter(phenotype_id %in% ss_oi$molecular_trait_id))
  
  exons_to_plot = ps_exons
  exon_cdss_to_plot = ps_exons_cdss
  
  variant_regions_vcf <- ss_oi %>% 
    dplyr::select(variant, chromosome, position) %>% 
    dplyr::mutate(region = paste0(chromosome,":", position, "-", position))
  
  snps_filt <- snps_all %>% 
    dplyr::filter(ID %in% variant_regions_vcf$variant) %>% 
    dplyr::arrange(CHROM, POS)
  
  if (study_name == "Lepik_2017") {
    var_genotype <- snps_filt %>% 
      dplyr::select(-c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")) %>% 
      base::t() %>% 
      BiocGenerics::as.data.frame() %>% 
      dplyr::rename("GT_DS" = "V1") %>% 
      dplyr::mutate(GT = gsub(pattern = "\\:.*", replacement = "", x = GT_DS)) %>% 
      dplyr::mutate(REF = gsub(pattern = "\\/.*", replacement = "", x = GT)) %>% 
      dplyr::mutate(ALT = gsub(pattern = ".*\\/", replacement = "", x = GT)) %>% 
      dplyr::mutate(DS = as.numeric(REF) + as.numeric(ALT)) %>% 
      dplyr::mutate(genotype_id = BiocGenerics::rownames(.)) %>% 
      dplyr::mutate(genotype_id = gsub(pattern = "\\.", replacement = "-", x = genotype_id))
  } else  {
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
  }
  
  sample_meta_clean = sample_metadata %>% 
    dplyr::filter(rna_qc_passed, genotype_qc_passed) %>%  
    dplyr::select(sample_id, genotype_id, qtl_group)
  
  var_genotype <- sample_meta_clean %>% 
    left_join(var_genotype, by = "genotype_id")
  
  bigwig_files <- list.files(bigwig_files_path, full.names = T)
  
  track_data_study <- data.frame(file_name = (gsub(pattern = paste0(bigwig_files_path, "/"), replacement = "", x = bigwig_files, fixed = T)),
                                 bigWig = bigwig_files)
  
  track_data_study <- track_data_study %>% 
    dplyr::mutate(sample_id = gsub(pattern = ".bigwig", replacement = "", x = file_name)) %>% 
    dplyr::filter(sample_id %in% sample_meta_clean$sample_id)
  
  track_data_study <-track_data_study %>% 
    dplyr::left_join(var_genotype %>% dplyr::select(sample_id, qtl_group, DS), by = c("sample_id")) %>% 
    dplyr::left_join(scaling_factor_data) %>% 
    dplyr::rename(scaling_factor = scaling_factors) %>% 
    dplyr::mutate(track_id = qtl_group) %>% 
    dplyr::mutate(colour_group = as.character(DS)) %>% 
    dplyr::select(sample_id, scaling_factor, bigWig, track_id, colour_group, qtl_group) %>% 
    dplyr::filter(qtl_group==qtl_group_in)
  
  # Generate the output path 
  signal_name <- paste0(gsub(pattern = ":", replacement = "_", x = ss_oi$molecular_trait_id), "__", ss_oi$variant, "__", ss_oi$gene_id)
  
  nom_exon_cc_sumstats_variant_phenotype_id <- read_and_filter_parquet( 
    file_list = ss_oi$nominal_exon_cc_path[[1]], 
    variant_to_match = ss_oi$variant,
    phenotype_id=ss_oi$gene_id
  )
  
  start_time <- time_here(prev_time = start_time, message_text = " >> prepared track_data_study in: ")
  if(!is.null(nom_exon_cc_sumstats_variant_phenotype_id)) {
    # Extract the QTLs of exons according to gene and variant of interest
    nom_exon_cc_sumstats_filt <- nom_exon_cc_sumstats_variant_phenotype_id %>% 
      dplyr::filter(rsid == rsid[1]) %>% # if variant has more than one rsid keep only the first unique rsid 
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
    start_time <- time_here(prev_time = start_time, message_text = " >> prepared nom_exon_cc_sumstats_filt in: ")
  }
  
  if (!is.null(mane_transcript_gene_map_file)) {
    MANE_transcript_oi <- mane_transcript_gene_map %>% 
      dplyr::filter(gene_id %in% ss_oi$gene_id) %>% 
      dplyr::pull(transcript_id)
    
    if (length(MANE_transcript_oi) > 0) {
      mane_transcript_exons <-  make_transcript_exon_granges(gff = gtf_ref, transcript_ids = MANE_transcript_oi)
      exons_to_plot <- append(exons_to_plot, mane_transcript_exons)
      exon_cdss_to_plot <- append(exon_cdss_to_plot, mane_transcript_exons)
    }
    start_time <- time_here(prev_time = start_time, message_text = " >> prepared MANE_transcript_oi in: ")
  }
    
  message(" ## Extracting coverage data")
  coverage_data_list = tryCatch(wiggleplotr::extractCoverageData(exons = exons_to_plot, 
                                                                 cdss = exon_cdss_to_plot, # Does var will be default NULL if not stated?
                                                                 plot_fraction = 0.2,
                                                                 track_data = track_data_study), 
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
  
  message(" ## Prepare transcript structure for plotting")
  tx_structure_df = prepareTranscriptStructureForPlotting(exon_ranges = coverage_data_list$tx_annotations$exon_ranges, 
                                                          cds_ranges = coverage_data_list$tx_annotations$cds_ranges, 
                                                          transcript_annotations = coverage_data_list$plotting_annotations)
  
  message(" ## Generate coverage data plots (wiggleplotr)")
  wiggle_plots <- wiggleplotr::plotCoverageData(coverage_data_list, alpha = 1, 
                                                fill_palette = wiggleplotr::getGenotypePalette(), 
                                                coverage_type = "line", return_subplots_list = TRUE, 
                                                show_legend = TRUE)
  
  start_time <- time_here(prev_time = start_time, message_text = " >> generated wiggleplots in: ")
  
  coverage_data_list$coverage_df <- coverage_data_list$coverage_df %>% dplyr::filter(!is.na(coverage))
  
  # permute the rows so that it becomes anonymous
  coverage_data_list$coverage_df <- coverage_data_list$coverage_df[sample(nrow(coverage_data_list$coverage_df)),]
  
  # BOXPLOTS START HERE
  message(" ## Prepare boxplot data")
  norm_exp_df_oi <- norm_exp_df %>% dplyr::filter(phenotype_id %in% cluster_introns$phenotype_id)
  norm_exp_df_oi <- tibble::column_to_rownames(.data = norm_exp_df_oi,var = "phenotype_id")
  norm_exp_df_oi <- norm_exp_df_oi %>% base::t() %>% 
    GenomicRanges::as.data.frame() %>% 
    tibble::rownames_to_column(var = "sample_id")
  norm_exp_df_oi <- norm_exp_df_oi %>% 
    tidyr::pivot_longer(cols = -sample_id, names_to="intron_id", values_to = "norm_exp")
  
  track_data_study_box <- track_data_study %>% 
    dplyr::mutate(genotype_text = as.factor(colour_group)) %>% 
    dplyr::filter(qtl_group %in% qtl_group_in) %>%
    dplyr::mutate(condition_name = qtl_group) %>% 
    dplyr::mutate(gene_name = ss_oi$gene_name) %>% 
    dplyr::mutate(snp_id = ss_oi$variant) 
  
  track_data_study_box <- norm_exp_df_oi %>%  
    dplyr::left_join(track_data_study_box, by = "sample_id") %>% 
    dplyr::mutate(is_significant = intron_id == ss_oi$molecular_trait_id)
  
  nom_cc_sumstats_variant_phenotype_id <- read_and_filter_parquet( 
    file_list = ss_oi$nominal_cc_path[[1]],
    variant_to_match = ss_oi$variant,
    phenotype_id=ss_oi$molecular_trait_id
  )
  
  
  #nom_cc_sumstats <- nom_cc_sumstats_variant_phenotype_id %>% # OK? Filtered by variant and molecular_trait_id
  #  dplyr::filter(variant %in% variant_regions_vcf$variant)
  
  # Keep only 1 rsid per variant per molecular_trait_id
  nom_cc_sumstats <- nom_cc_sumstats_variant_phenotype_id %>% 
    dplyr::group_by(molecular_trait_id, variant) %>% 
    dplyr::filter(rsid == rsid[1]) %>% 
    dplyr::ungroup()
  
  nom_cc_sumstats_filt <- nom_cc_sumstats %>% 
    dplyr::filter(variant %in% ss_oi$variant) %>% # do I need to filter even?
    dplyr::select(molecular_trait_id, pvalue, beta, se, maf) %>% 
    dplyr::rename(intron_id = molecular_trait_id)
  
  track_data_study_box_wrap <- track_data_study_box %>% 
    dplyr::left_join(nom_cc_sumstats_filt, by = "intron_id") %>% 
    dplyr::mutate(stats_text = paste0("Pval: ", pvalue, "			BETA: ", beta, 
                                      "\nSE: ", se)) %>% 
    dplyr::mutate(intron_id_with_stats = paste0(intron_id, "\n", stats_text))
  
  start_time <- time_here(prev_time = start_time, message_text = " >> until boxplot ready: ")
  track_data_study_box_wrap_for_RDS <- track_data_study_box_wrap %>%
    dplyr::select(genotype_text, norm_exp, is_significant, intron_id, pvalue, beta, se, snp_id, maf)
  
  track_data_study_box_wrap_for_RDS <- track_data_study_box_wrap_for_RDS[sample(nrow(track_data_study_box_wrap_for_RDS)),]
  
  tx_str_df <- tx_structure_df %>% dplyr::mutate(limit_max = max(coverage_data_list$limits))
  
  signal_name <- gsub(pattern = "&", replacement = "\\&", x = signal_name)

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
  start_time <- time_here(prev_time = start_time, message_text = " >> saving the plot data took: ")
}