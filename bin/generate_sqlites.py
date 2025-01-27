#!/usr/bin/env python


import argparse
import datetime
from typing import List, Tuple, Any,Optional
import io
import sqlite3
from functools import lru_cache
from math import isnan
from pathlib import Path
import pandas as pd
import pyarrow as pa
from pyarrow import parquet as pq
from pydantic import BaseModel

SerializedParquet = bytes
Pydict = dict[str, list[Any]]


def _pydict_to_bytes(data: Pydict) -> SerializedParquet:
    table = pa.Table.from_pydict(data)
    buffer = io.BytesIO()
    pq.write_table(table, buffer, compression='brotli')
    buffer.seek(0)
    serialized_bytes = buffer.read()
    return serialized_bytes


def _bytes_to_pydict(serialized_bytes: SerializedParquet) -> Pydict:
    buffer = io.BytesIO(serialized_bytes)
    buffer.seek(0)
    table = pq.read_table(buffer)
    data = table.to_pydict()
    return data


class SerializableDF(BaseModel):
    def _serialize_df(self) -> bytes:
        #data = self.dict()
        data = self.model_dump()
        serialized_bytes = _pydict_to_bytes(data)
        return serialized_bytes

    @classmethod
    def _from_bytes(cls, serialized_bytes: bytes):
        data = _bytes_to_pydict(serialized_bytes)
        return cls(**data)


class Metadata(BaseModel):
    gene_id: str
    credible_set_id: str
    molecular_trait_id: str
    variant: str
    gene_name: str
    x_limit: int
    y_limit: float

    dataset_id: Optional[str] = None
    study_label: Optional[str] = None
    tissue_label: Optional[str] = None
    condition_label: Optional[str] = None

    @classmethod
    def _from_directory(cls, path: Path) -> 'Metadata':  # , *, molecular_trait_id, variant, gene_id
        molecular_trait_id = get_df_from_directory(path, 'ss_oi_df').loc[0, 'molecular_trait_id']
        gene_name = get_df_from_directory(path, 'ss_oi_df').loc[0, 'gene_name']
        credible_set_id = get_df_from_directory(path, 'ss_oi_df').loc[0, 'cs_id']
        variant = get_df_from_directory(path, 'ss_oi_df').loc[0, 'variant']
        gene_id = get_df_from_directory(path, 'ss_oi_df').loc[0, 'gene_id']
        x_limit = get_df_from_directory(path, 'tx_str').loc[0, 'limit_max']

        # y_limit
        # extract the coverage plot and take the max value of the coverage field
        coverage_df = get_df_from_directory(path, 'coverage_df')
        y_limit = coverage_df['coverage'].max()
        return Metadata(gene_id=gene_id,
                        credible_set_id=credible_set_id,
                        molecular_trait_id=molecular_trait_id,
                        variant=variant,
                        gene_name=gene_name,
                        x_limit=x_limit,
                        y_limit=y_limit
                        )


class TranscriptAnno(SerializableDF):
    molecular_trait_ids: list[str]
    start: list[int]
    end: list[int]
    molecular_trait_ranks: list[int]
    feature_type: list[str]
    exon_index: list[Optional[int]]

    @staticmethod
    def _from_directory(path, gene_name: str):
        tx = get_df_from_directory(path, 'tx_str')
        filtered_tx = tx[tx['transcript_id'] == f'GENE:{gene_name}'].copy()
        # Add the 'exon_index' column to the filtered_tx dataframe
        filtered_tx['exon_index'] = range(1, len(filtered_tx) + 1)
        tx.loc[filtered_tx.index, 'exon_index'] = filtered_tx.exon_index
        # if exon_index is still missing, add it as NaN
        if 'exon_index' not in tx.columns:
            tx['exon_index'] = None

        tx = tx.astype({'feature_type': 'category'})
        return TranscriptAnno(
            molecular_trait_ids=tx.transcript_id.tolist(),
            start=tx.start.tolist(),
            end=tx.end.tolist(),
            molecular_trait_ranks=tx.transcript_rank.tolist(),
            feature_type=tx.feature_type.tolist(),
            exon_index=[int(i) if not isnan(i) else None for i in tx.exon_index.tolist()],
        )


class Coverage(SerializableDF):
    bins: list[int]
    cov_gt0: Optional[list[float]]
    cov_gt1: Optional[list[float]]
    cov_gt2: Optional[list[float]]

    @staticmethod
    def _from_directory(file_path: Path):
        coverage = get_df_from_directory(file_path, 'coverage_df')
        coverage = coverage[coverage.coverage.notnull()]  # drop rows where coverage is null
        new_coverage = pd.pivot_table(coverage, index='bins', columns='colour_group', values='coverage').reset_index()

        # rename columns
        # 0 -> cov_gt0
        # 1 -> cov_gt1
        # 2 -> cov_gt2

        new_coverage.rename(columns={"0": 'cov_gt0', "1": 'cov_gt1', "2": 'cov_gt2'}, inplace=True, errors='ignore')
        coverage = new_coverage.astype({'bins': int}, errors='raise')

        try:
            coverage = coverage.astype({'cov_gt0': float}, errors='ignore')
        except KeyError:
            pass

        try:
            coverage = coverage.astype({'cov_gt1': float}, errors='ignore')
        except KeyError:
            pass

        try:
            coverage = coverage.astype({'cov_gt2': float}, errors='ignore')
        except KeyError:
            pass

        return Coverage(
            bins=coverage.bins.tolist(),
            cov_gt0=None if 'cov_gt0' not in coverage.columns else coverage.cov_gt0.tolist(),
            cov_gt1=None if 'cov_gt1' not in coverage.columns else coverage.cov_gt1.tolist(),
            cov_gt2=None if 'cov_gt2' not in coverage.columns else coverage.cov_gt2.tolist(),
        )



    def _to_pydict(self) -> Pydict:
        data = {'bins': self.bins}

        if self.cov_gt0 is not None:
            data['cov_gt0'] = self.cov_gt0

        if self.cov_gt1 is not None:
            data['cov_gt1'] = self.cov_gt1

        if self.cov_gt2 is not None:
            data['cov_gt2'] = self.cov_gt2

        return data

    @classmethod
    def _from_pydict(cls, data: Pydict) -> 'Coverage':
        return cls(
            bins=data['bins'],
            cov_gt0=data.get('cov_gt0'),
            cov_gt1=data.get('cov_gt1'),
            cov_gt2=data.get('cov_gt2')
        )

    def _serialize(self) -> bytes:
        coverage_data = self._to_pydict()
        return _pydict_to_bytes(coverage_data)

    @classmethod
    def _from_bytes(cls, data: bytes) -> 'Coverage':
        coverage_data = _bytes_to_pydict(data)
        return cls._from_pydict(coverage_data)


class Exon_sumstats(SerializableDF):
    exon_index: list[int]
    molecular_trait_id: list[str]
    beta: list[float]
    ci: list[float]
    p_fdr: list[float]

    @staticmethod
    def _from_directory(file_path) -> Optional['Exon_sumstats']:
        nom_exon = get_df_from_directory(file_path, 'nom_exon_cc')
        if len(nom_exon) == 0:
            return None

        return Exon_sumstats(
            exon_index=nom_exon.exon_row_num.tolist(),
            molecular_trait_id=nom_exon.molecular_trait_id.tolist(),
            beta=nom_exon.beta.tolist(),
            ci=nom_exon.interval.tolist(),
            p_fdr=nom_exon.p_fdr.tolist(),
        )


class SingleBoxPlot(BaseModel):
    molecular_trait_id: str
    pvalue: float
    beta: float
    se: float
    maf: float
    variant: str

    abs_gt0: Optional[list[float]] = None
    abs_gt1: Optional[list[float]] = None
    abs_gt2: Optional[list[float]] = None

    norm_gt0: Optional[list[float]] = None
    norm_gt1: Optional[list[float]] = None
    norm_gt2: Optional[list[float]] = None

    def _to_pydict(self) -> Pydict:
        data = {
            'molecular_trait_id': [],
            'pvalue': [],
            'beta': [],
            'se': [],
            'maf': [],
            'variant': [],
            'genotype_text': [],
            'norm_exp': [],
        }

        genotypes = {'abs_gt0': 'genotype_0', 'abs_gt1': 'genotype_1', 'abs_gt2': 'genotype_2'}

        for genotype_key, genotype_text in genotypes.items():
            genotype_data = getattr(self, genotype_key)
            if genotype_data is not None:
                for i in genotype_data:
                    data['molecular_trait_id'].append(self.molecular_trait_id)
                    data['pvalue'].append(self.pvalue)
                    data['beta'].append(self.beta)
                    data['se'].append(self.se)
                    data['maf'].append(self.maf)
                    data['variant'].append(self.variant)
                    data['genotype_text'].append(genotype_text)
                    data['norm_exp'].append(i)

        return data


class Boxplot(BaseModel):
    boxplots: list[SingleBoxPlot]

    @classmethod
    def _from_directory(cls, file_path: Path) -> 'Boxplot':
        box_plot_df = get_df_from_directory(file_path, 'box_plot_df')
        # remove rows where all values are null
        box_plot_df = box_plot_df[box_plot_df.norm_exp.notnull()]
        boxplots = cls._df_to_boxplot_list(box_plot_df)
        return cls(boxplots=boxplots)

    @classmethod
    def _df_to_boxplot_list(cls, box_plot_df: pd.DataFrame) -> list[SingleBoxPlot]:
        # rename intron_id or tx_id to molecular_trait_id
        if 'intron_id' in box_plot_df.columns:
            box_plot_df = box_plot_df.rename(columns={'intron_id': 'molecular_trait_id'})
        elif 'tx_id' in box_plot_df.columns:
            box_plot_df = box_plot_df.rename(columns={'tx_id': 'molecular_trait_id'})
        intron_statistics = {}
        for molecular_trait_id, df in box_plot_df.groupby('molecular_trait_id'):
            intron_statistics[molecular_trait_id] = df.iloc[0][
                ['pvalue', 'beta', 'se', 'maf', 'snp_id', 'genotype_text', 'molecular_trait_id']].to_dict()
        intron_statistics = {k: v for k, v in intron_statistics.items()}
        # rename snp_id to variant
        for key in intron_statistics:
            intron_statistics[key]['variant'] = intron_statistics[key]['snp_id']
            del intron_statistics[key]['snp_id']
        box_plot_df = box_plot_df[['norm_exp', 'molecular_trait_id', 'genotype_text']]
        # we sort by norm_exp as well to make sure we don't accidentally leak the original order
        box_plot_df = box_plot_df.sort_values(by=['molecular_trait_id', 'genotype_text', 'norm_exp'])
        for group, df in box_plot_df.groupby(['molecular_trait_id', 'genotype_text'],observed=False):
            molecular_trait_id, genotype = group
            norm_exp_values = df['norm_exp'].tolist()
            if molecular_trait_id not in intron_statistics:
                continue
            # let's take the previously collected metadata and get a reference to it
            boxplot = intron_statistics[molecular_trait_id]
            # add the exp_values in the format we expect
            boxplot[f'abs_gt{genotype}'] = norm_exp_values
        boxplots = []
        for molecular_trait_id, boxplot in intron_statistics.items():
            boxplots.append(SingleBoxPlot(**boxplot))
        return boxplots

    def _to_pydict(self) -> Pydict:
        # calls the _to_pydict method of each SingleBoxPlot
        # and then combines the results into a single dict by extending the lists
        data = {
            'molecular_trait_id': [],
            'pvalue': [],
            'beta': [],
            'se': [],
            'maf': [],
            'variant': [],
            'genotype_text': [],
            'norm_exp': []
        }
        for boxplot in self.boxplots:
            boxplot_data = boxplot._to_pydict()
            for key, value in boxplot_data.items():
                data[key].extend(value)
        return data

    @classmethod
    def _from_pydict(cls, data: Pydict) -> 'Boxplot':
        '''
        reverses the _to_pydict method
        '''
        boxplots = []
        # Create a dictionary to store temporary boxplot data
        temp_boxplots = {}
        # Iterate through the data Pydict
        for idx in range(len(data['molecular_trait_id'])):
            molecular_trait_id = data['molecular_trait_id'][idx]

            # If the molecular_trait_id is not in the temp_boxplots dictionary, create a new entry
            if molecular_trait_id not in temp_boxplots:
                temp_boxplots[molecular_trait_id] = {
                    'molecular_trait_id': molecular_trait_id,
                    'pvalue': data['pvalue'][idx],
                    'beta': data['beta'][idx],
                    'se': data['se'][idx],
                    'maf': data['maf'][idx],
                    'variant': data['variant'][idx],
                    'genotype_0': [],
                    'genotype_1': [],
                    'genotype_2': []
                }

            # Add the norm_exp value to the corresponding genotype list
            genotype_key = data["genotype_text"][idx]
            temp_boxplots[molecular_trait_id][genotype_key].append(data['norm_exp'][idx])

        # Convert the temp_boxplots dictionary to a list of Boxplot objects
        for molecular_trait_id, boxplot_data in temp_boxplots.items():
            boxplot = SingleBoxPlot(
                molecular_trait_id=boxplot_data['molecular_trait_id'],
                pvalue=boxplot_data['pvalue'],
                beta=boxplot_data['beta'],
                se=boxplot_data['se'],
                maf=boxplot_data['maf'],
                variant=boxplot_data['variant'],
                abs_gt0=boxplot_data['genotype_0'] or None,
                abs_gt1=boxplot_data['genotype_1'] or None,
                abs_gt2=boxplot_data['genotype_2'] or None
            )
            boxplots.append(boxplot)

        # Instantiate a new Box_plot_df2 object with the list of Boxplot objects
        return cls(boxplots=boxplots)

    def _serialize(self):
        box_plot_df = self._to_pydict()
        return _pydict_to_bytes(box_plot_df)

    @classmethod
    def _from_bytes(cls, data: bytes):
        box_plot_df = _bytes_to_pydict(data)
        return cls._from_pydict(box_plot_df)


class SerializedApiData(BaseModel): 
    gene_id: str
    credible_set_id:str
    molecular_trait_id: str
    variant: str
    gene_name: str
    x_limit: int
    y_limit: float

    tx_structure_serialized: bytes
    coverage_serialized: bytes
    nom_exon_cc_serialized: Optional[bytes]
    box_plot_df_serialized: bytes

    def _deserialize(self) -> 'ApiData':
        return ApiData(
            meta=Metadata(
                gene_id=self.gene_id,
                molecular_trait_id=self.molecular_trait_id,
                variant=self.variant,
                gene_name=self.gene_name,
                x_limit=self.x_limit,
                y_limit=self.y_limit,
            ),
            transcript_anno=TranscriptAnno._from_bytes(self.tx_structure_serialized),
            coverage=Coverage._from_bytes(self.coverage_serialized),
            exon_sumstats=Exon_sumstats._from_bytes(
                self.nom_exon_cc_serialized) if self.nom_exon_cc_serialized else None,
            boxplot=Boxplot._from_bytes(self.box_plot_df_serialized),
        )


class ApiData(BaseModel):
    meta: Metadata
    transcript_anno: TranscriptAnno
    coverage: Coverage
    exon_sumstats: Optional[Exon_sumstats]
    boxplot: Boxplot

    @staticmethod
    def _from_directory(directory_path: Path) -> 'ApiData':

        meta = Metadata._from_directory(directory_path)
        transcript_anno = TranscriptAnno._from_directory(directory_path, meta.gene_name)

        coverage = Coverage._from_directory(directory_path)
        exon_sumstats = Exon_sumstats._from_directory(directory_path)
        box_plot = Boxplot._from_directory(directory_path)

        return ApiData(meta=meta, transcript_anno=transcript_anno, coverage=coverage, exon_sumstats=exon_sumstats, boxplot=box_plot)

    def serialize(self) -> SerializedApiData:

        tx_structure_serialized = self.transcript_anno._serialize_df()
        coverage_serialized = self.coverage._serialize()
        nom_exon_cc_serialized = self.exon_sumstats._serialize_df() if self.exon_sumstats else None
        box_plot_df_serialized = self.boxplot._serialize()

        return SerializedApiData(
            gene_id=self.meta.gene_id,
            credible_set_id= self.meta.credible_set_id,
            molecular_trait_id=self.meta.molecular_trait_id,
            variant=self.meta.variant,
            gene_name=self.meta.gene_name,
            x_limit=self.meta.x_limit,
            y_limit=self.meta.y_limit,
            tx_structure_serialized=tx_structure_serialized,
            coverage_serialized=coverage_serialized,
            nom_exon_cc_serialized=nom_exon_cc_serialized,
            box_plot_df_serialized=box_plot_df_serialized,
        )

    @classmethod
    @lru_cache(maxsize=10)
    def from_database(cls, database: str, gene_id: str, molecular_trait_id: str, variant: str) -> Optional['ApiData']:
        conn = sqlite3.connect(database)

        # Get the metadata from the metadata table
        cursor_meta = conn.execute("SELECT key, value FROM metadata")
        metadata = {row[0]: row[1] for row in cursor_meta}

        # Fetch the data from the main table based on gene_id, molecular_trait_id, and variant
        cursor = conn.execute('''SELECT gene_id, molecular_trait_id, variant, gene_name, x_limit, y_limit,
                                        tx_structure_serialized, coverage_serialized, nom_exon_cc_serialized, box_plot_df_serialized
                                  FROM main
                                  WHERE gene_id=? AND molecular_trait_id=? AND variant=?''',
                              (gene_id, molecular_trait_id, variant))

        # Get the first (and only) row of the result
        row = cursor.fetchone()
        if row is None:

            # Check if gene_id is present in the database
            cursor_gene = conn.execute("SELECT * FROM main WHERE gene_id=?", (gene_id,))
            if cursor_gene.fetchone() is None:
                raise AssertionError(f"gene_id={gene_id} not found")

            # Check if molecular_trait_id is present in the database
            cursor_phenotype = conn.execute("SELECT * FROM main WHERE molecular_trait_id=?", (molecular_trait_id,))
            if cursor_phenotype.fetchone() is None:
                raise AssertionError(f"molecular_trait_id={molecular_trait_id} not found")

            # Check if variant is present in the database
            cursor_variant = conn.execute("SELECT * FROM main WHERE variant=?", (variant,))
            if cursor_variant.fetchone() is None:
                raise AssertionError(f"variant={variant} not found")

            # If gene_id, molecular_trait_id, and variant are individually present, but the specific combination is not found
            raise AssertionError("combination of gene_id, molecular_trait_id, and variant not found")
        else:
            # You can now return the data or process it as needed
            pass

        # If a row is found, create a dictionary mapping column names to their respective values in the row
        if row:
            column_names = [desc[0] for desc in cursor.description]
            row_dict = dict(zip(column_names, row))

            apidata = SerializedApiData(
                gene_id=row_dict['gene_id'],
                molecular_trait_id=row_dict['molecular_trait_id'],
                variant=row_dict['variant'],
                gene_name=row_dict['gene_name'],
                x_limit=row_dict['x_limit'],
                y_limit=row_dict['y_limit'],
                tx_structure_serialized=row_dict['tx_structure_serialized'],
                coverage_serialized=row_dict['coverage_serialized'],
                nom_exon_cc_serialized=row_dict['nom_exon_cc_serialized'],
                box_plot_df_serialized=row_dict['box_plot_df_serialized']
            )._deserialize()

            # Update the meta attributes in ApiData instance
            apidata.meta.dataset_id = metadata.get('dataset_id')
            apidata.meta.study_label = metadata.get('study_label')
            apidata.meta.tissue_label = metadata.get('tissue_label')
            apidata.meta.condition_label = metadata.get('condition_label')

            return apidata

        else:
            return None

def get_database_stem(file_name: str):
    return Path(file_name).name.removeprefix('plot_data_').removesuffix('.tar.gz')


def get_df_from_directory(dir_path: Path, df_name: str):
    """
    Reads a Parquet file from a directory based on a given prefix (df_name).

    Args:
        dir_path (Path): Path to the directory containing the files.
        df_name (str): Prefix of the file to select.

    Returns:
        pd.DataFrame: The DataFrame loaded from the selected Parquet file.

    Raises:
        AssertionError: If the df_name is not valid or the file is not found.
    """
    # Validate the df_name
    valid_df_names = ['nom_exon_cc', 'box_plot_df', 'tx_str', 'ss_oi_df', 'coverage_df']
    assert df_name in valid_df_names, f"'{df_name}' is not a valid dataframe name. Must be one of {valid_df_names}"
    # Search for the file starting with df_name in the directory
    dir_path = Path(dir_path)
    matching_files = list(dir_path.glob(f"{df_name}_*.parquet"))
    if not matching_files:
        raise FileNotFoundError(f"No file found in directory '{dir_path}' starting with '{df_name}'")
    if len(matching_files) > 1:
        raise ValueError(f"Multiple files found in directory '{dir_path}' starting with '{df_name}': {matching_files}")
    parquet_file_path = matching_files[0]
    df = pd.read_parquet(parquet_file_path)
    return df

def initialize_database(conn_path: str | Path) -> None:
    """Initialize the database with main and metadata tables."""
    # Create a connection to the database
    conn = sqlite3.connect(conn_path)
    # Create the main table
    conn.execute('''CREATE TABLE main
                 (id INTEGER PRIMARY KEY,
                 
                 -- this quadruplet is the pseudokey for this table
                 -- dataset, molecular_trait, variant defines the row
                 -- credible_set is defined by the dataset and molecular_trait, so credible (credible_set, variant) also define the row
                 dataset_id TEXT NOT NULL, -- this duplicates the dataset filename, but I find it's cleaner to have it in the table
                 credible_set_id TEXT NOT NULL,
                 molecular_trait_id TEXT NOT NULL,
                 variant TEXT NOT NULL,
                 
                 -- these are metadata for the API endpoint
                 gene_id TEXT NOT NULL,
                 gene_name TEXT NOT NULL,
                 x_limit INTEGER NOT NULL,
                 y_limit INTEGER NOT NULL,
                 
                 
                 tx_structure_serialized BLOB NOT NULL,
                 coverage_serialized BLOB NOT NULL,
                 nom_exon_cc_serialized BLOB,
                 box_plot_df_serialized BLOB NOT NULL,
                 UNIQUE(gene_id, molecular_trait_id, variant),
                 UNIQUE(credible_set_id, variant)
                 )
                 ''')

    # Create the metadata table with database_generated timestamp
    conn.execute('''CREATE TABLE metadata
                    (key TEXT PRIMARY KEY,
                     value TEXT NOT NULL)''')

    # Add the database_generated timestamp to the metadata table
    now = datetime.datetime.now(datetime.UTC).isoformat()
    conn.execute('''INSERT INTO metadata (key, value)
                     VALUES (?, ?)''', ('database_generated_utc', now))
    # Commit the changes and close the connection
    conn.commit()
    conn.close()


def write_batch_data_to_main_table(conn: sqlite3.Connection, data_list: List[Tuple]) -> None:
    """Insert a batch of serialized API data into the main table."""
    # Insert the data into the main table using executemany
    conn.executemany('''INSERT INTO main 
    (dataset_id, credible_set_id, gene_id, molecular_trait_id, variant, gene_name, x_limit, y_limit, tx_structure_serialized, coverage_serialized, nom_exon_cc_serialized, box_plot_df_serialized) 
    VALUES
    (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) 
    ''', data_list)
    conn.commit()

def process_directory_files(directory_path: Path) -> SerializedApiData:
    """Process gzipped files and return serialized data."""
    try:
        data = ApiData._from_directory(directory_path)
        serialized = data.serialize()
        return serialized
    except Exception as e:
        print(f"Error occurred while processing '{directory_path}': {e}")
        raise

def process_and_write(directory_containing_pqs, conn, dataset_id):
    data = process_directory_files(directory_path=directory_containing_pqs)
    output_row = (
        dataset_id,
        data.credible_set_id,
        data.gene_id,
        data.molecular_trait_id,
        data.variant,
        data.gene_name,
        data.x_limit,
        data.y_limit,
        data.tx_structure_serialized,
        data.coverage_serialized,
        data.nom_exon_cc_serialized,
        data.box_plot_df_serialized
    )
    write_batch_data_to_main_table(conn, [output_row])


def main(dataset_id,dir_paths) -> None:
    dir_paths_df = pd.read_csv(dir_paths, header=None, names=["path"])
    TARGET_DB = f'{dataset_id}.sqlite'
    initialize_database(TARGET_DB)
    conn = sqlite3.connect(TARGET_DB)
    dir_paths_df['path'].apply(lambda x: process_and_write(x, conn, dataset_id))
    conn.execute('''CREATE INDEX idx_main_gene_pheno_var ON main(gene_id, molecular_trait_id, variant)''')
    conn.execute('analyze')
    conn.commit()
    conn.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--dataset_id', required=True, type=str, help="Dataset ID") 
    parser.add_argument('-s', '--source_root_file', required=True, type=str, help="File with  paths to the directories containing parquet files.")

    args = parser.parse_args() 

    dataset_id = args.dataset_id
    source_roots = args.source_root_file
    main(dataset_id, source_roots)

