#!/usr/bin/env python3

import pathlib
import pandas as pd
from logzero import logger
from neokant.index.index import KmerIndex
from neokant.database.database import DataBase
from neokant.database.queries import Queries
from neokant.annotator import EXPECTED_CTS_COLUMNS
import xxhash


class FastaHandler:
    @staticmethod
    def write_fasta(entries: pd.DataFrame, fasta_file: str):
        assert all([x in entries.columns for x in ['query_cts_id', 'query_sequence']]), "Columns are missing"
        with open(fasta_file, 'w') as file_handle:
            for row in entries.itertuples(index=False):
                file_handle.write(f">{row.query_cts_id}\n")
                file_handle.write(f"{row.query_sequence}\n")


class Annotator:
    def __init__(self, working_dir: str, sequence_table: str, db_file: str) -> None:
        self.working_dir = pathlib.Path(working_dir)
        self.sequence_table = self.read_context_seq(sequence_table)
        self.query_fasta = self.working_dir / "query.fa"
        self.db = DataBase(db_file)
        self.queries = Queries(self.db)

    @staticmethod
    def read_context_seq(sequence_table):
        """
        Read context sequence table into dataframe.

        :param sequence_table: Path to sequence table
        :return: pandas dataframe
        """
        with open(sequence_table, "r") as file_handle:
            seq = pd.read_csv(file_handle, sep="\t")
        assert all([x in seq.columns for x in EXPECTED_CTS_COLUMNS]), "Missing columns in context sequence table"
        return seq

    def _generate_soi(self):
        """
        Generate target sequence of interest by extractin subsequence sccording to position and query_length column
        :return:
        """
        self.sequence_table['query_sequence'] = \
            self.sequence_table.apply(lambda x: self._generate_short_sequence(x.cts_seq, x.pos, x.query_length), axis=1)

    def _generate_soi_cts_id(self):
        """
        Generate cts_id of target sequence of interest
        :return:
        """
        self.sequence_table['query_cts_id'] = \
            self.sequence_table.apply(lambda x: xxhash.xxh64(x.query_sequence).hexdigest(), axis=1)

    @staticmethod
    def _generate_short_sequence(cts_seq, pos, length):
        """
        Generate sequence of interest to query in kmer index e.g. the fusion breakpoint, splice junction. If length
        is null the original sequence is returned for search in index
        :param cts_seq: A context sequence
        :param pos:  Position of interest in cts to query
        :param length: The length of the query sequence
        :return:
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
        :return:
        """
        logger.info(f"Writing context sequences to fasta file {self.query_fasta}")
        FastaHandler.write_fasta(self.sequence_table, self.query_fasta)

    def prepare_cts(self):
        """
        Prepare context sequences provided by user for targeted search and
        dump into fasta format for query
        :return:
        """
        logger.info("Generating breakpoint sequences...")
        self._generate_soi()
        logger.info("Generating cts_id for search sequences")
        self._generate_soi_cts_id()
        logger.info("Dumping query sequences to disk")
        self._write_to_fasta()

    def search_cts(self,
                   pipeline: pathlib.Path,
                   index: pathlib.Path,
                   method: str = "raptor",
                   kmer_ratio: float = 0.5,
                   reindeer_sample_mapping: str = None,
                   raptor_sample_mapping: str = None,
                   kmindex_cutoff: float = 0.7,
                   slurm: bool = False):
        """
        Call query pipeline and parse results into pandas DataFrame
        """
        index = KmerIndex(pipeline=pipeline,
                          index=index,
                          method=method,
                          reindeer_sample_mapping=reindeer_sample_mapping,
                          raptor_sample_mapping=raptor_sample_mapping,
                          kmer_ratio=kmer_ratio,
                          kmindex_cutoff=kmindex_cutoff
                          )
        hits = index.search_index(self.query_fasta, self.working_dir, slurm=slurm)
        parsed_results = index.result_parser(hits)

        dict_to_pandas = []
        for cts, samples in parsed_results.items():
            for this_sample in samples:
                dict_to_pandas.append((cts, this_sample))
        parsed_results = pd.DataFrame(dict_to_pandas, columns=['cts_id', 'sample_name'])
        return parsed_results

    def annotate_cts(self, parsed_results: pd.DataFrame, annot_style: str = "normal"):
        """
        Given all hits in an index collect tissue and number of tissue samples in whole index
        :param parsed_results:
        :param annot_style:
        :return:
        """
        logger.info("Annotating sample hits with corresponding project id")
        parsed_results['study_id'] = parsed_results.apply(lambda row: self.queries.get_project_id(row.sample_name),
                                                            axis=1)
        # Group all indexing results by project_id. This allows us to query the database for each table once, regardless
        # of the query sequence. Annotation results are then merged back to the dataframe
        logger.info("Annotating sample hits with sample level metadata")
        parsed_results = parsed_results.groupby('study_id').apply(lambda df: self.queries.annotate_samples_of_project(df))
        # Subset to required columns
        parsed_results = parsed_results.loc[:, ['cts_id', 'sample_name', 'tissue', 'developmental_stage', 'disease']]
        # Count combinations of tissue, disease and developmental stage
        parsed_results = parsed_results.groupby(['cts_id', 'disease', 'developmental_stage', 'tissue'])[
            ['disease', 'developmental_stage', 'tissue']].size().to_frame('count').reset_index()

        df = parsed_results[['cts_id']].drop_duplicates()
        df = pd.merge(df, parsed_results, how="left")
        df['count'] = df['count'].fillna(0).astype('int')
        df = df[['cts_id', 'count', 'disease', 'developmental_stage', 'tissue']]
        return df

    def annotate_sequences(self, annotated_cts):
        """
        Merge annotated cts results with original table supplied by the user
        :param annotated_cts:
        :return:
        """
        df = pd.merge(self.sequence_table, annotated_cts, left_on="query_cts_id", right_on="cts_id")
        df.drop('cts_id_y', inplace=True, axis=1)
        return df


