#!/usr/bin/env python3

import pathlib
import pandas as pd
from k4neo.index.index_loader import load_metaindex_from_manifest
from k4neo.index.index_processor import KmerIndexProcessor
from k4neo.database_sqlite.queries import Queries
from k4neo.annotator import (
    EXPECTED_CTS_COLUMNS,
    NON_TUMOR_TISSUE,
    TUMOR_TISSUE,
    load_annotator_config,
)
from k4neo.helper.helper import FastaHandler, SequenceOperation, InputValidation, DiskIO
import numpy as np
from loguru import logger


class Annotator:
    """
    Annotator class
    """

    def __init__(
        self,
        config_yaml: pathlib.Path,
        index_kmer_size: int = 21,
        working_dir: pathlib.Path | None = None,
    ) -> None:
        
        self.config = load_annotator_config(config_yaml)
        self.sequence_table = self.read_context_seq(self.config.sequence_table_output)
        
        if working_dir is None:
            self.working_dir = self.config.working_dir
        else:
            self.working_dir = working_dir
        
        
        self.non_queryable = pd.DataFrame(
            columns=["cts_id", "cts_seq", "query_length", "pos", "query_sequence", "query_cts_id"]
        )
        # , self.non_queryable, index_kmer_size
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
        meta_index = load_metaindex_from_manifest(index_manifest)

        index_processor = KmerIndexProcessor(
            meta_index=meta_index, pipeline=pipeline, workflow_profile=workflow_profile
        )
        query_pipeline_results = index_processor.search_index(
            self.config.query_fasta, self.working_dir, slurm=slurm, cores=cores, kmer_ratio=kmer_ratio
        )
        parsed_results = index_processor.result_parser2(
            query_pipeline_results=query_pipeline_results, cores=cores, kmer_ratio=kmer_ratio
        )
        return parsed_results

    def _count_aggregation(self, parsed_results: pd.DataFrame) -> pd.DataFrame:
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
        not_expressed["study_id"] = np.nan

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
            ["cts_id", "developmental_stage", "tissue"], as_index=False, dropna=True
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
            ["cts_id", "disease", "tissue"], as_index=False, dropna=True
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

    def annotate_cts(
        self, parsed_results: pd.DataFrame, queries: Queries, annot_style: str = "normal"
    ):
        """
        Given all hits in an index collect tissue and number of tissue samples in whole index.
        Combine annotated CTS with not expressed targets and return aggegrated table
        :param parsed_results:
        :param annot_style:
        :return:
        """
        # Select cts not found in index and append columns required to merge later with annotated results

        # Group all indexing results by project_id. This allows us to query the database for each table once, regardless
        # of the query sequence. Annotation results are then merged back to the dataframe
        logger.debug("Removing sequences not detetcted in any sample of index")

        not_expressed = self._split_found(parsed_results)
        parsed_results.dropna(inplace=True, ignore_index=True)

        logger.debug("Annotating sample hits with corresponding study annotation.")
        study_annotation = queries.get_sample_study()
        parsed_results = parsed_results.merge(study_annotation, how="left", on="sample_name")

        if len(parsed_results.index) == 0:
            logger.warning("None of the queried sequences was found in index.")
            return not_expressed
        logger.debug("Annotating sample hits with sample level metadata.")
        parsed_results = parsed_results.groupby("study_id", dropna=False).apply(
            lambda sub_df: queries.annotate_samples_of_project(sub_df)
        )
        parsed_results.reset_index(drop=True, inplace=True)
        parsed_results = self._count_aggregation(parsed_results)

        logger.debug("Annotating with pre-computed counts from database.")
        parsed_results = parsed_results.groupby("study_id", group_keys=False).apply(
            lambda sub_df: queries.annotate_tissue_counts(sub_df)
        )
        parsed_results.reset_index(drop=True, inplace=True)

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

    def annotate_sample_rate2(self, annotated_cts, queries: Queries, min_total=1):
        """
        Add sample rate to sequences
        """
        tissue_counts = queries.get_tissue_counts()

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
        healthy_sample_rate["samples_per_tissue"] = healthy_sample_rate[
            "samples_per_tissue"
        ].astype(int)
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
