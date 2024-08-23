#!/usr/bin/env python3

import pathlib
import pandas as pd
from logzero import logger
from k4neo.index.index import KmerIndex
from k4neo.database.database import DataBase
from k4neo.database.queries import Queries
from k4neo.annotator import EXPECTED_CTS_COLUMNS
import xxhash
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class FastaHandler:
    """
    Generic class to handle FASTA operations
    """
    @staticmethod
    def write_fasta(entries: pd.DataFrame, fasta_file: str):
        """
        Write query context sequences into fasta format
        """
        assert all([x in entries.columns for x in ['query_cts_id', 'query_sequence']]), "Columns are missing"
        sequences = entries[['query_cts_id', 'query_sequence']].drop_duplicates()
        fasta_entries = []
        with open(fasta_file, 'w') as file_handle:
            for row in sequences.itertuples(index=False):
                record = SeqRecord(
                    Seq(row.query_sequence),
                    id = row.query_cts_id,
                    name = row.query_cts_id,
                    description='')
                fasta_entries.append(record)
            SeqIO.write(fasta_entries, file_handle, "fasta")

class Annotator:
    """
    Annotator class
    """
    def __init__(self, working_dir: str, sequence_table: str, db_file: str, index_kmer_size: int = 21) -> None:
        self.working_dir = pathlib.Path(working_dir)
        self.sequence_table, self.non_queryable = self.read_context_seq(sequence_table, index_kmer_size)
        self.query_fasta = self.working_dir / "query.fa"
        self.db = DataBase(db_file)
        self.queries = Queries(self.db)

    @staticmethod
    def read_context_seq(sequence_table: pathlib.Path, index_kmer_size: int) -> tuple:
        """Read context sequenc table.

        Read context sequence table and filter sequences that do not fullfill the minimal
        query length requirement of Raptor. (query_len + 4)

        Args:
            sequence_table (pathlib.Path): Path to context sequence table
            index_kmer_size (int): K-mer size of index. Used to determine minimal query length.

        Returns:
            tuple: Returns two dataframes with sequences to process and sequences that can be queried.
        """
        with open(sequence_table, "r") as file_handle:
            seq = pd.read_csv(file_handle, sep="\t")
        assert all([x in seq.columns for x in EXPECTED_CTS_COLUMNS]), "Missing columns in context sequence table"
        # Sequences that are shorter than minimiser_size are excluded from search
        seq_to_short = seq[seq["cts_seq"].str.len() < index_kmer_size + 4]
        seq = seq[~seq['cts_id'].isin(seq_to_short['cts_id'])]
        return seq, seq_to_short

    def _generate_soi(self):
        """
        Generate target sequence of interest by extracting 
        subsequence according to position and query_length column
        """
        self.sequence_table['query_sequence'] = \
            self.sequence_table.apply(lambda x: self._generate_short_sequence(x.cts_seq, x.pos, x.query_length), axis=1)

    def _generate_soi_cts_id(self):
        """
        Generate cts_id of target sequence of interest. ID is hash value of the query sequence
        """
        self.sequence_table['query_cts_id'] = \
            self.sequence_table.apply(lambda x: xxhash.xxh64(x.query_sequence).hexdigest(), axis=1)

    @staticmethod
    def _generate_short_sequence(cts_seq: str, pos: int = None, length: int = None) -> str:
        """Generate query sequence

        Generate sequence of interest to query in kmer index e.g. the fusion breakpoint, splice junction. If length
        is null the original sequence is returned for search in index e.g. a full length transcript or SNP sequence.

        Args:
            cts_seq (str): The full-length context sequence.
            pos (int, optional): Position in context sequence to extract query seqeunce from. Defaults to None.
            length (int, optional): Length of query sequence. +/- length/2 around pos. Defaults to None.

        Returns:
            str: Context sequence string
        """
        if pd.isnull(pos) or pd.isnull(length):
            return cts_seq

        start = pos - round(length / 2)
        if start < 0:
            start = 0

        stop = start + length
        if stop > len(cts_seq) - 1:
            stop = len(cts_seq) - 1

        sequence_of_interest = cts_seq[start:stop]
        return sequence_of_interest

    def _write_to_fasta(self):
        """
        Wrapper to call FastaHandler class
        """
        logger.info(f"-> Writing context sequences to fasta file: {self.query_fasta}")
        FastaHandler.write_fasta(self.sequence_table, self.query_fasta)

    def prepare_cts(self):
        """
        Prepare context sequences provided by user for targeted search and
        dump into fasta format for query
        """
        logger.info("-> Generating breakpoint sequences.")
        self._generate_soi()
        logger.info("-> Generating cts_id for search sequences.")
        self._generate_soi_cts_id()
        logger.info("-> Dumping query sequences to disk.")
        self._write_to_fasta()

    def search_cts(self,
                   pipeline: pathlib.Path,
                   index_manifest: pathlib.Path,
                   kmer_ratio: float = 0.7,
                   slurm: bool = False,
                   cores: int = 16) -> pd.DataFrame:
        """Search context sequence in k-mer indices

        Search query fasta in k-mer indices of manifest with QueryPipeline.

        Args:
            pipeline (pathlib.Path): Path to SnakeMake pipeline (Snakefile)
            index_manifest (pathlib.Path): Path to k4neo index manifest (yaml)
            kmer_ratio (float, optional): Required fraction of shared k-mers between query and sample. Defaults to 0.7.
            slurm (bool, optional):  If True, QueryPipeline will submit jobs to slurm scheduler. Defaults to False.
            cores (int, optional): Number of cores for pipeline. Defaults to 16.

        Returns:
            pd.DataFrame: A pandas DataFrame with parsed results for each method from manifest file.
        """
        index = KmerIndex(pipeline=pipeline,
                          index_manifest=index_manifest,
                          kmer_ratio=kmer_ratio
                          )
        query_pipeline_results = index.search_index(self.query_fasta, self.working_dir, slurm=slurm, cores=cores)
        parsed_results = index.result_parser(query_pipeline_results=query_pipeline_results, cores=cores)
        method_results = {}
        for this_method, this_parsed_results in parsed_results.items():
            dict_to_pandas = []
            # This could be done more efficiently by zipping samples to cts
            for cts, samples in this_parsed_results.items():
                dict_to_pandas.extend([(cts, this_sample) for this_sample in samples])
            method_results[this_method] = pd.DataFrame(dict_to_pandas, columns=['cts_id', 'sample_name'])
        return method_results

    def _annotate_studies(self, parsed_results: pd.DataFrame) -> pd.DataFrame:
        """Annotate sample hits with study id

        Retrive sample to study mapping from database and merge with search results.

        Args:
            parsed_results (pd.DataFrame): Parsed results of k-mer pipeline.

        Returns:
            pd.DataFrame: Updated DataFrame with study column.
        """
        study_annotation = self.queries.get_sample_study()
        parsed_results = parsed_results.merge(study_annotation,
                                              how="left",
                                              on="sample_name")
        return parsed_results

    def _annotate_sample_metadata(self, parsed_results: pd.DataFrame) -> pd.DataFrame:
        """Annotate sample hits with metadata

        Each occurence of the query sequence in the k-mer index is annotated with
        the sample-level metadata from the database.
        Sample hits are aggregated for each study by tissue, developmental state and disease
        to determine the occurence of the sequence. 

        cts_1   5   healthy adult   liver   study1
        cts_1   4   healthy fetal   liver   study2

        Args:
            parsed_results (pd.DataFrame): Parsed results of k-mer pipeline.

        Returns:
            pd.DataFrame: Aggregated DataFrame with counts for each combination of tissue,
            developmental_stage and tissue per study.
        """
        parsed_results = \
            parsed_results.groupby('study_id', dropna=False).apply(
                lambda sub_df: self.queries.annotate_samples_of_project(sub_df)
            )
        parsed_results.reset_index(drop=True, inplace=True)
        # Subset to required columns
        parsed_results = parsed_results.loc[:,
                         ['cts_id', 'study_id', 'sample_name', 'tissue', 'developmental_stage', 'disease']]
        # Count combinations of tissue, disease and developmental stage
        parsed_results = parsed_results.groupby(
            ['cts_id', 'study_id', 'disease', 'developmental_stage', 'tissue'])[
            ['disease', 'developmental_stage', 'tissue']].size().to_frame('count').reset_index()

        df = parsed_results[['cts_id']].drop_duplicates()
        df = pd.merge(df, parsed_results, how="left")
        df['count'] = df['count'].fillna(0).astype('int')
        return df

    def _annotate_counts(self, parsed_results: pd.DataFrame) -> pd.DataFrame:
        """Add pre-computed tissue counts.

        Retrieve pre-computed total tissue counts from database and add to
        DataFrame.

        Args:
            parsed_results (pd.DataFrame): Parsed results of k-mer pipeline.

        Returns:
            pd.DataFrame: A DataFrame with total counts for each combination of tissue,
            developmental_stage and tissue per study.
        """
        parsed_results =\
            parsed_results.groupby('study_id', group_keys=False).apply(
                lambda sub_df: self.queries.annotate_tissue_counts(sub_df)
            )

        parsed_results.reset_index(drop=True, inplace=True)
        return parsed_results

    @staticmethod
    def _split_found(parsed_results: pd.DataFrame) -> pd.DataFrame:
        """Remove sequences not detected in index

        Sequences that were not detected in any sample of the index
        are removed from the result DataFrame to handle and annotate
        them separately.

        Args:
            parsed_results (pd.DataFrame): Parsed results of k-mer pipeline.

        Returns:
            pd.DataFrame: DataFrame of non-detected sequences in final output format.
        """
        
        not_expressed = \
            parsed_results.loc[
                parsed_results['sample_name'].isnull(), ["cts_id"]]
        not_expressed['count'] = 0
        not_expressed['total'] = 0
        not_expressed['disease'] = np.nan
        not_expressed['developmental_stage'] = np.nan
        not_expressed['tissue'] = np.nan
        return not_expressed

    def _calculate_sample_rate(self, parsed_results: pd.DataFrame):
        """
        High level aggregate. Sum all tissue hits of a developmental stage
        per cts_id an calculate sample rate. The number of samples per tissue
        containing the sequence of interest.
        """
        tissue_counts = self.queries.get_tissue_counts()
        tissue_counts = tissue_counts.groupby(['developmental_stage', 'tissue'])['total'].\
            sum().reset_index()
        parsed_results = parsed_results.groupby(
            ['cts_id', 'developmental_stage', 'tissue'])['count'].\
                sum().reset_index()
        df = parsed_results[['cts_id']].drop_duplicates()
        df = pd.merge(df, parsed_results, how="left")
        df['count'] = df['count'].fillna(0).astype('int')
        df = pd.merge(df, tissue_counts, how="left")
        df['sample_rate'] = df.apply(lambda row: round(row['count'] / row['total'], 2), axis=1)
        
        return df

    def annotate_cts(self, parsed_results: pd.DataFrame, annot_style: str = "normal"):
        """
        Given all hits in an index collect tissue and number of tissue samples in whole index.
        Combine annotated CTS with not expressed targets and return aggegrated table
        :param parsed_results:
        :param annot_style:
        :return:
        """
        logger.info("-> Annotating sample hits with corresponding study annotation.")
        # Select cts not found in index and append columns required to merge later with annotated results

        parsed_results = self._annotate_studies(parsed_results)

        # Group all indexing results by project_id. This allows us to query the database for each table once, regardless
        # of the query sequence. Annotation results are then merged back to the dataframe
        not_expressed = self._split_found(parsed_results)
        parsed_results = parsed_results.dropna()
        if len(parsed_results.index) == 0:
            logger.info('-> None of the queried sequences was found in index.')
            return not_expressed
        logger.info("-> Annotating sample hits with sample level metadata.")
        parsed_results = self._annotate_sample_metadata(parsed_results)
        logger.info("-> Annotating with pre-computed counts.")
        parsed_results = self._annotate_counts(parsed_results)
        parsed_results = pd.concat([parsed_results, not_expressed])
        parsed_results =\
            parsed_results[['cts_id', 'count', 'total', 'disease', 'developmental_stage', 'tissue', 'study_id']]
        return parsed_results

    def annotate_sequences(self, annotated_cts):
        """
        Merge annotated cts results with original table supplied by the user
        :param annotated_cts:
        :return: Annotated search table
        """
        df = pd.merge(self.sequence_table, annotated_cts, left_on="query_cts_id", right_on="cts_id")
        df.drop('cts_id_y', inplace=True, axis=1)
        df.rename(columns={'cts_id_x': 'cts_id'}, inplace=True)
        df["total"] = pd.to_numeric(df["total"])

        return df
    
    def annotate_sample_rate(self, annotated_cts, min_total=15):
        """
        Add sample rate to sequences
        """
        sample_rate_df = self._calculate_sample_rate(annotated_cts)
        df = pd.merge(self.sequence_table, sample_rate_df, left_on="query_cts_id", right_on="cts_id")
        df.drop('cts_id_y', inplace=True, axis=1)
        df.rename(columns={'cts_id_x': 'cts_id'}, inplace=True)
        df["sample_rate"] = pd.to_numeric(df["sample_rate"])
        df = df.loc[df['total'] >= min_total]

        return df



