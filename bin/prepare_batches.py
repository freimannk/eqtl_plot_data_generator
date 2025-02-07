#!/usr/bin/env python

import math
import polars as pl
import argparse



def parse_parquet_file_info_from_file(file_list_path):
    file_paths_df = pl.read_csv(file_list_path, has_header=False, new_columns=["file_path"])
    # Strip leading/trailing whitespace from the file_path column
    file_paths_df = file_paths_df.with_columns(
        file_path=file_paths_df["file_path"].str.strip_chars()
    )
    # Extract metadata using regex
    file_paths_df = file_paths_df.with_columns(
        chromosome=file_paths_df["file_path"].str.extract(r"_chr_(\d+|X)_", 1).alias("chromosome"),
        start=file_paths_df["file_path"].str.extract(r"_chr_(?:\d+|X)_(\d+)_", 1).cast(pl.Int64).alias("start"),
        end=file_paths_df["file_path"].str.extract(r"_chr_(?:\d+|X)_(?:\d+)_([\d]+)\.parquet", 1).cast(pl.Int64).alias("end")
    )
    return file_paths_df


def find_matching_files(chromosome, position, file_info_df):
    """
    Find all matching Parquet files for a given chromosome and position,
    and return them as a single list.
    """
    matches = file_info_df.filter(
        (pl.col("chromosome") == str(chromosome)) &
        (pl.col("start") <= position) &
        (pl.col("end") >= position)
    )
    return matches["file_path"].to_list()

def prepare_study_id_credible_set(credible_sets, highest_pip_vars_per_cs, study_id, dataset_id, qtl_group, study_name, quant_method):
    """
    Prepare study_id credible set file.
    """
    variant_plot_data = highest_pip_vars_per_cs.select("cs_id","variant")
    variant_plot_data= variant_plot_data.rename({"variant": "plot_variant"})
    credible_sets_w_plot_variant = credible_sets.join(variant_plot_data, on="cs_id", how="left")
    credible_sets_w_plot_variant = credible_sets_w_plot_variant.with_columns(
    pl.lit(study_id).alias("study_id"),
    pl.lit(study_name).alias("study_label"),
    pl.lit(dataset_id).alias("dataset_id"),
    pl.lit(quant_method).alias("quantification_method"),
    pl.lit(qtl_group).alias("dataset_label")
)
    credible_sets_w_plot_variant= credible_sets_w_plot_variant.rename({"cs_id": "credible_set","cs_size":"credible_set_size"})
    credible_sets_w_plot_variant = credible_sets_w_plot_variant.select(["study_id","study_label","dataset_id","molecular_trait_id", "gene_id","gene_name", "variant","rsid","quantification_method","credible_set", "credible_set_size","pip","pvalue","beta","se","dataset_label","plot_variant"])
    filename = f"{study_id}_{dataset_id}_{qtl_group}_{quant_method}.parquet"
    credible_sets_w_plot_variant.write_parquet(filename,compression="snappy")
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split susie output files into batches.")
    parser.add_argument('-a', '--study_name', required=True, help="Study name")
    parser.add_argument('-n', '--nominal_sumstats_files', required=True, help="File with nominal sumstats files (tsv format).")
    parser.add_argument('-d', '--dataset_id', required=True, help="Unique dataset_id.")
    parser.add_argument('-i', '--study_id', required=True, help="Study id.")
    parser.add_argument('-q', '--quant_method', required=True, help="Quantification method.")
    parser.add_argument('-g', '--qtl_group', required=True, help="Qtl group.")
    parser.add_argument('-e', '--nominal_sumstats_exon_files', required=True, help="File with nominal sumstats exon files (tsv format).")
    parser.add_argument('-s', '--susie_output_file', required=True, help="Purity filtered susie output (parquet format).")
    parser.add_argument('-p', '--phenotype_metadata', required=True, help="Phenotype metadata file. Tab separated file")
    parser.add_argument('-c', '--chunk_size', required=True,type=int,help="Perform analysis in chunks. Eg value 5 10 would indicate that phenotypes are split into 10 chunks and the 5th one of those will be processed.")


    args = parser.parse_args()


    # Read the file paths from the text file created by Nextflow
    nominal_sumstats_info = parse_parquet_file_info_from_file(args.nominal_sumstats_files)
    nominal_exon_sumstats_info = parse_parquet_file_info_from_file(args.nominal_sumstats_exon_files)
    df_credible_sets = pl.read_parquet(args.susie_output_file)
    phenotype_metadata_df = pl.read_csv(args.phenotype_metadata, separator='\t', schema_overrides={"chromosome": pl.Utf8})
    df_credible_sets = df_credible_sets.join(
    phenotype_metadata_df.select(["phenotype_id", "strand","gene_name","group_id","quant_id"]),  
    left_on="molecular_trait_id",  
    right_on="phenotype_id",       
    how="left"                    
)
    df_credible_sets = df_credible_sets.with_columns(
    (pl.lit(args.dataset_id+"_") + pl.col("cs_id")).alias("cs_id")
)
    original_columns = df_credible_sets.columns
    highest_pip_vars_per_cs = (
        df_credible_sets
        .filter(pl.col("position").is_not_null())
        .sort(["cs_id", "pip"], descending=[False, True])
        .group_by("cs_id")
        .agg(pl.all().first())
    )
    # keep order
    highest_pip_vars_per_cs = highest_pip_vars_per_cs.select(original_columns)
    prepare_study_id_credible_set(df_credible_sets, highest_pip_vars_per_cs, args.study_id, args.dataset_id, args.qtl_group, args.study_name, args.quant_method)
    # Add nominal_cc_path and nominal_exon_cc_path columns
    highest_pip_vars_per_cs = highest_pip_vars_per_cs.with_columns([
        # Find and combine all matching files for nominal_cc_path
        pl.struct(["chromosome", "position"]).map_elements(
            lambda x: find_matching_files(x['chromosome'], x["position"], nominal_sumstats_info),
            return_dtype=pl.List(pl.Utf8)  # Ensure return type is a list of strings
        ).alias("nominal_cc_path"),

        # Find and combine all matching files for nominal_exon_cc_path
        pl.struct(["chromosome", "position"]).map_elements(
            lambda x: find_matching_files(x['chromosome'], x["position"], nominal_exon_sumstats_info),
            return_dtype=pl.List(pl.Utf8)  # Ensure return type is a list of strings
        ).alias("nominal_exon_cc_path")
    ])
    chunk_output_prefix = args.dataset_id + "_" + args.qtl_group + "_" + args.quant_method
    chunk_size = args.chunk_size
    n_chunks = math.ceil(len(highest_pip_vars_per_cs) / chunk_size)
    for i in range(1, n_chunks + 1):
        start = (i - 1) * chunk_size
        end = min(i * chunk_size, len(highest_pip_vars_per_cs))
        chunk = highest_pip_vars_per_cs[start:end]
        filename = f"{chunk_output_prefix}_{i}_{n_chunks}.parquet"
        chunk.write_parquet(filename,compression="snappy")


