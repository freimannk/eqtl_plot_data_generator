#!/usr/bin/env python

import argparse
import sqlite3
import duckdb
import os

def create_db_table(sqlite_db:str, parquet_files:list):
    # Set DuckDB home directory and extension path
    os.environ["DUCKDB_HOME"] = "/opt/duckdb_extensions"
    os.environ["DUCKDB_EXTENSION_PATH"] = "/opt/duckdb_extensions/extensions"

    # Connect to DuckDB
    con = duckdb.connect()
    con.execute("SET home_directory = '/opt/duckdb_extensions';")
    con.execute("SET extension_directory = '/opt/duckdb_extensions/extensions';")
    con.execute("LOAD sqlite;")
    con.execute(f"ATTACH '{sqlite_db}' AS sqlite_db (TYPE SQLITE);")
    con.execute(f"""
        CREATE TABLE sqlite_db.credible_set_table (
            id INTEGER PRIMARY KEY NOT NULL,
            study_id TEXT,
            study_label TEXT,
            dataset_id TEXT,
            molecular_trait_id TEXT,
            gene_id TEXT,
            gene_name TEXT,
            variant TEXT,
            rsid TEXT,
            quantification_method TEXT,
            credible_set TEXT,
            credible_set_size INTEGER,
            pip FLOAT,
            pvalue FLOAT,
            beta FLOAT,
            se FLOAT,
            dataset_label TEXT,
            plot_variant TEXT
        );
    """)

    query = f"""
        INSERT INTO sqlite_db.credible_set_table (
              id, study_id, study_label, dataset_id, molecular_trait_id, gene_id, gene_name, 
            variant, rsid, quantification_method, credible_set, credible_set_size, 
            pip, pvalue, beta, se, dataset_label, plot_variant
        ) 
        SELECT ROW_NUMBER() OVER() AS id ,* FROM read_parquet({parquet_files});
    """
    con.execute(query)
    con.close()

def create_indexes(sqlite_db:str):
    conn = sqlite3.connect(sqlite_db)
    cursor = conn.cursor()

    indexes = [
        "CREATE INDEX idx_molecular_trait_id ON credible_set_table(molecular_trait_id);",
        "CREATE INDEX idx_gene_name ON credible_set_table(gene_name);",
        "CREATE INDEX idx_credible_set ON credible_set_table(credible_set);",
        "CREATE INDEX idx_variant ON credible_set_table(variant);",
        "CREATE INDEX idx_rsid ON credible_set_table(rsid);",
        "CREATE INDEX idx_dataset ON credible_set_table(study_label);",
        "CREATE INDEX idx_for_plotting ON credible_set_table(dataset_id, gene_id, molecular_trait_id, variant);"
    ]
    for index in indexes:
        cursor.execute(index)
    conn.commit()
    conn.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate multiple parquet files using DuckDB.")
    parser.add_argument('-f', '--files', nargs='+', required=True, help="List of parquet files to concatenate.")
    parser.add_argument('-o', '--output', required=True, help="Output sqlite file.")
    args = parser.parse_args()
    study_id_datasets_pqs= args.files
    sqlite_file_name = args.output
    create_db_table(sqlite_file_name, study_id_datasets_pqs)
    create_indexes(sqlite_file_name)