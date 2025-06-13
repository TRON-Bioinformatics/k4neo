#!/usr/bin/env python3

import pathlib
import pandas as pd
from k4neo.index.index import KmerIndex
from k4neo.database.database import DataBase
from k4neo.database.queries import Queries
from k4neo.annotator import (
    EXPECTED_CTS_COLUMNS,
    NON_TUMOR_TISSUE,
    TUMOR_TISSUE,
    IMMUNO_PRIVILIGED_TISSUE,
)
from k4neo.helper.helper import FastaHandler, SequenceOperation, InputValidation
import numpy as np
from loguru import logger


class Annotator:
    """
    Annotator class
    """

    def __init__(
        self,
        working_dir: str,
        sequence_table: str,
        db_file: str,
        index_kmer_size: int = 21,
    ) -> None:
        self.working_dir = pathlib.Path(working_dir)
        self.sequence_table = self.read_context_seq(sequence_table)
        self.non_queryable = pd.DataFrame(
            columns=["cts_id", "cts_seq", "query_length", "pos", "query_sequence", "query_cts_id"]
        )
        # , self.non_queryable, index_kmer_size
        self.query_fasta = self.working_dir / "query.fa"
        self.db = DataBase(db_file)
        self.queries = Queries(self.db)
        self.index_kmer_size = index_kmer_size

    @staticmethod
    def read_context_seq(sequence_table: pathlib.Path) -> tuple:
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
        ret, missing_cols = InputValidation.columns_missing(seq, EXPECTED_CTS_COLUMNS)
        assert not ret, f"-> Missing columns: {missing_cols} in input table"

        return seq

    def _generate_target_sequence(self):
        """
        Generate target sequence of interest by extracting
        subsequence according to position and query_length column
        """
        self.sequence_table["query_sequence"] = self.sequence_table.apply(
            lambda x: SequenceOperation.subset_cts(x.cts_seq, x.pos, x.query_length),
            axis=1,
        )

        self.sequence_table["query_cts_id"] = self.sequence_table.apply(
            lambda x: SequenceOperation.cts_id(x.query_sequence), axis=1
        )

    def _filter_seq_to_short(self):
        """
        Filter sequences that are shorter than minimiser_size are excluded from search
        """
        self.non_queryable = pd.concat(
            [
                self.non_queryable,
                self.sequence_table[
                    self.sequence_table["query_sequence"].str.len() < self.index_kmer_size + 4
                ],
            ]
        )
        self.sequence_table = self.sequence_table[
            ~self.sequence_table["cts_id"].isin(self.non_queryable["cts_id"])
        ]

    def _write_to_fasta(self):
        """
        Wrapper to call FastaHandler class
        """
        if not self.query_fasta.exists():
            logger.info(f"-> Writing context sequences to fasta file: {self.query_fasta}")
            FastaHandler.write_fasta(self.sequence_table, self.query_fasta)
            return
        logger.info(
            f"-> Fasta file: {self.query_fasta} already exists in working directory. Re-using for annotation"
        )

    def prepare_cts(self):
        """
        Prepare context sequences provided by user for targeted search and
        dump into fasta format for query
        """
        logger.info("-> Generating breakpoint sequences.")
        self._generate_target_sequence()

        logger.info("-> Filtering sequences to short for query in k-mer index.")
        self._filter_seq_to_short()

        logger.info("-> Writing target sequences to disk.")
        self._write_to_fasta()

    def search_cts(
        self,
        pipeline: pathlib.Path,
        workflow_profile: pathlib.Path,
        index_manifest: pathlib.Path,
        kmer_ratio: float = 0.7,
        slurm: bool = False,
        cores: int = 16,
    ) -> pd.DataFrame:
        """Search context sequence in k-mer indices

        Search query fasta in k-mer indices of manifest with instance of QueryPipeline.

        Args:
            pipeline (pathlib.Path): Path to SnakeMake pipeline (Snakefile)
            index_manifest (pathlib.Path): Path to k4neo index manifest (yaml)
            kmer_ratio (float, optional): Required fraction of shared k-mers between query and sample. Defaults to 0.7.
            slurm (bool, optional):  If True, QueryPipeline will submit jobs to slurm scheduler. Defaults to False.
            cores (int, optional): Number of cores for pipeline. Defaults to 16.

        Returns:
            pd.DataFrame: A pandas DataFrame with parsed results for each method from manifest file.
        """
        index = KmerIndex(pipeline=pipeline, workflow_profile=workflow_profile, index_manifest=index_manifest, kmer_ratio=kmer_ratio)
        query_pipeline_results = index.search_index(
            self.query_fasta, self.working_dir, slurm=slurm, cores=cores
        )
        parsed_results = index.result_parser(
            query_pipeline_results=query_pipeline_results, cores=cores
        )
        method_results = {}
        for this_method, this_parsed_results in parsed_results.items():
            dict_to_pandas = []
            # This could be done more efficiently by zipping samples to cts
            for cts, samples in this_parsed_results.items():
                dict_to_pandas.extend([(cts, this_sample) for this_sample in samples])
            method_results[this_method] = pd.DataFrame(
                dict_to_pandas, columns=["cts_id", "sample_name"]
            )
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
        parsed_results = parsed_results.merge(study_annotation, how="left", on="sample_name")
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
        parsed_results = parsed_results.groupby("study_id", dropna=False).apply(
            lambda sub_df: self.queries.annotate_samples_of_project(sub_df)
        )
        parsed_results.reset_index(drop=True, inplace=True)
        # Subset to required columns
        parsed_results = parsed_results.loc[
            :,
            [
                "cts_id",
                "study_id",
                "sample_name",
                "tissue",
                "developmental_stage",
                "disease",
            ],
        ]
        # Count combinations of tissue, disease and developmental stage
        parsed_results = (
            parsed_results.groupby(
                ["cts_id", "study_id", "disease", "developmental_stage", "tissue"]
            )[["disease", "developmental_stage", "tissue"]]
            .size()
            .to_frame("count")
            .reset_index()
        )

        df = parsed_results[["cts_id"]].drop_duplicates()
        df = pd.merge(df, parsed_results, how="left")
        df["count"] = df["count"].fillna(0).astype("int")
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
        parsed_results = parsed_results.groupby("study_id", group_keys=False).apply(
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

        not_expressed = parsed_results.loc[parsed_results["sample_name"].isnull(), ["cts_id"]]
        not_expressed["count"] = 0
        not_expressed["total"] = 0
        not_expressed["disease"] = np.nan
        not_expressed["developmental_stage"] = np.nan
        not_expressed["tissue"] = np.nan
        return not_expressed

    def _calculate_sample_rate(self, parsed_results: pd.DataFrame):
        """
        High level aggregate. Sum all tissue hits of a developmental stage
        per cts_id an calculate sample rate. The number of samples per tissue
        containing the sequence of interest.
        """
        tissue_counts = self.queries.get_tissue_counts()
        tissue_counts = (
            tissue_counts.groupby(["developmental_stage", "tissue"])["total"].sum().reset_index()
        )
        parsed_results = (
            parsed_results.groupby(["cts_id", "developmental_stage", "tissue"])["count"]
            .sum()
            .reset_index()
        )
        df = parsed_results[["cts_id"]].drop_duplicates()
        df = pd.merge(df, parsed_results, how="left")
        df["count"] = df["count"].fillna(0).astype("int")
        df = pd.merge(df, tissue_counts, how="left")
        df["sample_rate"] = df.apply(lambda row: round(row["count"] / row["total"], 2), axis=1)

        return df

    @staticmethod
    def _calculate_healthy_sample_rate(parsed_results: pd.DataFrame, tissue_counts: pd.DataFrame):
        """
        High level aggregate. Sum all tissue hits of a developmental stage
        per cts_id an calculate sample rate. The number of samples per tissue
        containing the sequence of interest.
        """
        # Remove TCGA and other tumor tissues
        tissue_counts = tissue_counts.loc[tissue_counts["disease"].isin(NON_TUMOR_TISSUE)]
        # Get number of samples per tissue and developmental state
        tissue_counts = (
            tissue_counts.groupby(["developmental_stage", "tissue"])["total"].sum().reset_index()
        )
        tissue_counts = tissue_counts.rename(columns={"total": "samples_per_tissue"})
        # tissue_counts['samples_per_index'] = tissue_counts['samples_per_tissue'].sum()

        # Generate all tissue / cts combinations
        cts_tissue_comb = (
            parsed_results[["cts_id"]].drop_duplicates().merge(tissue_counts, how="cross")
        )
        # Count occurence of each cts per tissue
        parsed_counts = parsed_results.groupby(
            ["cts_id", "developmental_stage", "tissue"], as_index=False, dropna=False
        ).agg(count=("count", "sum"))

        # Merge and calculate sample rate
        count_table = cts_tissue_comb.merge(
            parsed_counts, how="left", on=["cts_id", "developmental_stage", "tissue"]
        ).fillna({"count": 0})
        count_table["sample_rate"] = count_table["count"] / count_table["samples_per_tissue"]
        # parsed_results['index_sample_rate'] = parsed_results['total_index_count'] / parsed_results['samples_per_index']

        return count_table

    @staticmethod
    def _calculate_tumor_sample_rate(parsed_results: pd.DataFrame, tissue_counts: pd.DataFrame):
        """
        High level aggregate. Sum all tissue hits of a developmental stage
        per cts_id an calculate sample rate. The number of samples per tissue
        containing the sequence of interest.
        """
        # Subset for valid TCGA cancers
        tumor_counts = tissue_counts.loc[tissue_counts["disease"].isin(TUMOR_TISSUE)]
        # Get number of samples per cancer entity
        tumor_counts = tumor_counts.groupby(["disease", "tissue"])["total"].sum().reset_index()
        tumor_counts = tumor_counts.rename(columns={"total": "index_count"})

        # Generate all tumor / cts combinations
        cts_tumor_comb = (
            parsed_results[["cts_id"]].drop_duplicates().merge(tumor_counts, how="cross")
        )
        # Count occurence of each cts per tumor
        parsed_counts = parsed_results.groupby(
            ["cts_id", "disease", "tissue"], as_index=False, dropna=False
        ).agg(count=("count", "sum"))

        # Merge and calculate tumor rate
        count_table = cts_tumor_comb.merge(
            parsed_counts, how="left", on=["cts_id", "disease", "tissue"]
        ).fillna({"count": 0})
        count_table["sample_rate"] = count_table["count"] / count_table["index_count"]

        return count_table

    def _annotate_tumor_specificity(healthy_ts_rate: pd.DataFrame):
        """
        Annotate tumor specificity criteria based on expression in healthy tissues
        """
        pass

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
            logger.info("-> None of the queried sequences was found in index.")
            return not_expressed

        logger.info("-> Annotating sample hits with sample level metadata.")
        parsed_results = self._annotate_sample_metadata(parsed_results)

        logger.info("-> Annotating with pre-computed counts.")
        parsed_results = self._annotate_counts(parsed_results)
        parsed_results = pd.concat([parsed_results, not_expressed])
        parsed_results = parsed_results[
            [
                "cts_id",
                "count",
                "total",
                "disease",
                "developmental_stage",
                "tissue",
                "study_id",
            ]
        ]
        return parsed_results

    def annotate_sequences(self, annotated_cts):
        """
        Merge annotated cts results with original table supplied by the user
        :param annotated_cts:
        :return: Annotated search table
        """
        df = pd.merge(
            self.sequence_table,
            annotated_cts,
            left_on="query_cts_id",
            right_on="cts_id",
        )
        df.drop("cts_id_y", inplace=True, axis=1)
        df.rename(columns={"cts_id_x": "cts_id"}, inplace=True)
        df["total"] = pd.to_numeric(df["total"])

        return df

    def annotate_sample_rate(self, annotated_cts, min_total=15):
        """
        Add sample rate to sequences
        """
        sample_rate_df = self._calculate_sample_rate(annotated_cts)
        df = pd.merge(
            self.sequence_table,
            sample_rate_df,
            left_on="query_cts_id",
            right_on="cts_id",
        )
        df.drop("cts_id_y", inplace=True, axis=1)
        df.rename(columns={"cts_id_x": "cts_id"}, inplace=True)
        df["sample_rate"] = pd.to_numeric(df["sample_rate"])
        df = df.loc[df["total"] >= min_total]

        return df

    def annotate_sample_rate2(self, annotated_cts, min_total=1):
        """
        Add sample rate to sequences
        """
        tissue_counts = self.queries.get_tissue_counts()
        healthy_sample_rate = self._calculate_healthy_sample_rate(annotated_cts, tissue_counts)
        healthy_sample_rate = pd.merge(
            self.sequence_table,
            healthy_sample_rate,
            left_on="query_cts_id",
            right_on="cts_id",
        )
        healthy_sample_rate.drop("cts_id_y", inplace=True, axis=1)
        healthy_sample_rate.rename(columns={"cts_id_x": "cts_id"}, inplace=True)
        healthy_sample_rate["sample_rate"] = pd.to_numeric(healthy_sample_rate["sample_rate"])
        healthy_sample_rate["count"] = healthy_sample_rate["count"].astype(int)
        healthy_sample_rate["samples_per_tissue"] = healthy_sample_rate["samples_per_tissue"].astype(int)
        healthy_sample_rate = healthy_sample_rate.loc[
            healthy_sample_rate["samples_per_tissue"] >= min_total
        ]

        tumor_sample_rate = self._calculate_tumor_sample_rate(annotated_cts, tissue_counts)
        tumor_sample_rate = pd.merge(
            self.sequence_table,
            tumor_sample_rate,
            left_on="query_cts_id",
            right_on="cts_id",
        )
        tumor_sample_rate.drop("cts_id_y", inplace=True, axis=1)
        tumor_sample_rate.rename(columns={"cts_id_x": "cts_id"}, inplace=True)
        tumor_sample_rate["sample_rate"] = pd.to_numeric(tumor_sample_rate["sample_rate"])
        tumor_sample_rate["count"] = tumor_sample_rate["count"].astype(int)
        tumor_sample_rate["index_count"] = tumor_sample_rate["index_count"].astype(int)

        tumor_sample_rate = tumor_sample_rate.loc[tumor_sample_rate["index_count"] >= min_total]

        return healthy_sample_rate, tumor_sample_rate
