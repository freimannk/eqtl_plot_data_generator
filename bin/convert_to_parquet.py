#!/usr/bin/env python

import duckdb
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('-i', help='inverse-normal transformed matrix file', required=True)
argparser.add_argument('-t', help='TPM matrix', required=True)
argparser.add_argument('-o', help='The output pq file beginning', required=True)

args = argparser.parse_args()

norm_exp_file = args.i
tpm_exp_file = args.t
output = args.o

norm_file_output = output + ".parquet"
tpm_file_output = output + "_TPM.parquet"

def convert(txt_file,file_name):
    con = duckdb.connect()
    con.execute(f"""
        COPY (
            SELECT *
            FROM read_csv_auto('{txt_file}')) 
        TO '{file_name}' (FORMAT 'parquet');
    """)
    con.close()

convert(norm_exp_file,norm_file_output)
convert(tpm_exp_file,tpm_file_output)