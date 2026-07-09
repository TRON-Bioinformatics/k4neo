#!/usr/bin/env python3

import pathlib
from typing import Dict
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
from k4neo.helper.helper import FastaHandler, SequenceOperation, InputValidation, DiskIO, Worker
from k4neo.helper.async_writer import AsyncDFWriter
from k4neo.parser.index_parser import IndexResultParser
import numpy as np
from joblib import Parallel, delayed
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
        
        if self.working_dir is None:
            raise ValueError(
                "Working directory is not set. Please provide it either via the CLI "
                "or in the annotator configuration (working_dir)."
            )

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
        sample_integer_encoding: Dict[str, int] | None = None,
        kmer_ratio: float = 0.7,
        slurm: bool = False,
        cores: int = 16,
    ) -> pd.DataFrame:
        """Search context sequence in k-mer indices

        Search query fasta in k-mer indices of manifest with instance of QueryPipeline.

        Args:
            pipeline (pathlib.Path): Path to SnakeMake pipeline (Snakefile)
            index_manifest (pathlib.Path): Path to k4neo index manifest (yaml)
            sample_integer_encoding (Dict[str, int] | None): Mapping of sample names to unique integers for memory optimization.
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
            self.config.query_fasta,
            self.working_dir,
            slurm=slurm,
            cores=cores,
            kmer_ratio=kmer_ratio,
        )

        parsed_results = index_processor.result_parser(
            query_pipeline_results=query_pipeline_results, cores=cores, kmer_ratio=kmer_ratio, sample_integer_encoding=sample_integer_encoding
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
                "sample_id",
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

        not_expressed = parsed_results.loc[parsed_results["sample_id"].isnull(), ["cts_id"]]
        not_expressed["count"] = 0
        not_expressed["total"] = 0
        not_expressed["disease"] = np.nan
        not_expressed["developmental_stage"] = np.nan
        not_expressed["tissue"] = np.nan
        not_expressed["study_id"] = np.nan

        return not_expressed
    
    def _calculate_index_sample_rate(self, parsed_results: pd.DataFrame, tissue_counts: pd.DataFrame):
        """Calculate sample rate of sequence of interest in whole index

        Args:
            parsed_results (pd.DataFrame): Parsed results of k-mer pipeline.
            tissue_counts (pd.DataFrame): DataFrame with tissue counts.

        Returns:
            pd.DataFrame: DataFrame with sample rates for each sequence in the index.
        """
        tissue_counts = (
            tissue_counts.loc[tissue_counts["disease"].isin(NON_TUMOR_TISSUE)]["total"].sum()
        )
        parsed_results = parsed_results.loc[parsed_results["disease"].isin(NON_TUMOR_TISSUE)]
        if tissue_counts == 0:
            logger.warning("Denominator (total sample count) is zero. Cannot calculate index sample rate.")
            return parsed_results.assign(index_sample_rate=np.nan)

        index_sample_rate = parsed_results.groupby("cts_id")["count"].sum().reset_index()
        index_sample_rate["index_sample_rate"] = index_sample_rate.apply(
            lambda row: row["count"] / tissue_counts, axis=1
        )
        return index_sample_rate

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
            ["cts_id", "developmental_stage", "tissue"]
        )["count"].sum().reset_index()

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
            ["cts_id", "disease", "tissue"]
        )["count"].sum().reset_index()

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
        parsed_results = parsed_results.merge(study_annotation, how="left", on="sample_id")

        if len(parsed_results.index) == 0:
            logger.warning("None of the queried sequences was found in index.")
            return not_expressed[
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


        logger.debug("Annotating sample hits with sample level metadata.")
        parsed_results = parsed_results.groupby("study_id", dropna=False)[['cts_id', 'sample_id', 'study_id']].apply(
            lambda sub_df: queries.annotate_samples_of_project(sub_df)
        )
        parsed_results.reset_index(drop=True, inplace=True)
        parsed_results = self._count_aggregation(parsed_results)

        logger.debug("Annotating with pre-computed counts from database.")
        parsed_results = parsed_results.groupby("study_id", group_keys=False)[['cts_id', 'count', 'study_id', 'tissue', 'developmental_stage', 'disease']].apply(
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
        # Use categorical data types to reduce memory usage and speed up processing
        parsed_results["disease"] = parsed_results["disease"].astype("category")
        parsed_results["developmental_stage"] = parsed_results["developmental_stage"].astype("category")
        parsed_results["tissue"] = parsed_results["tissue"].astype("category")
        parsed_results["study_id"] = parsed_results["study_id"].astype("category")

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

    def annotate_sample_rate(self, annotated_cts, queries: Queries, min_total=1):
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

        index_sample_rate = self._calculate_index_sample_rate(annotated_cts, tissue_counts)
        index_sample_rate = pd.merge(
            self.sequence_table,
            index_sample_rate,
            left_on="query_cts_id",
            right_on="cts_id",
        )
        index_sample_rate.drop("cts_id_y", inplace=True, axis=1)
        index_sample_rate.rename(columns={"cts_id_x": "cts_id"}, inplace=True)
        index_sample_rate["index_sample_rate"] = pd.to_numeric(index_sample_rate["index_sample_rate"])
        index_sample_rate["count"] = index_sample_rate["count"].astype(int)
        index_sample_rate["samples_per_index"] = tissue_counts.loc[tissue_counts["disease"].isin(NON_TUMOR_TISSUE)]["total"].sum()

        return healthy_sample_rate, tumor_sample_rate, index_sample_rate

    def write_non_queryable(self, output_prefix: str, compression: bool) -> None:
        """Write non-queryable sequences to disk if any exist.

        Sequences too short to be queried in the k-mer index are stored in
        ``self.non_queryable``.  This method writes them to a TSV file next to
        the other output files.  The file is skipped when the DataFrame is empty.

        Args:
            output_prefix (str): Output file path prefix (without extension).
            compression (bool): Write gzip-compressed output when True.
        """
        if len(self.non_queryable.index) > 0:
            logger.info("-> Writing non-queryable sequences to disk")
            output_non_queryable = (
                pathlib.Path(output_prefix + "_non_querable.tsv.gz")
                if compression
                else pathlib.Path(output_prefix + "_non_querable.tsv")
            )
            DiskIO.write_df(self.non_queryable, output_non_queryable, compression)

    @staticmethod
    def _write_annotation_batch(
        annot_writer,
        healthy_writer,
        tumor_writer,
        index_writer,
        sample_hits,
        healthy_sample_rate,
        tumor_sample_rate,
        index_sample_rate,
        first_chunk: bool,
    ) -> None:
        """Send one batch of annotation results to the async writer threads.

        Writes the three result DataFrames produced by a single annotation
        worker to their respective :class:`AsyncDFWriter` instances.  On the
        first call (``first_chunk=True``) the header row is written; on
        subsequent calls it is suppressed and the file is appended.

        Args:
            annot_writer (AsyncDFWriter): Writer for the per-sample annotation table.
            healthy_writer (AsyncDFWriter): Writer for healthy-tissue sample-rate table.
            tumor_writer (AsyncDFWriter): Writer for tumor sample-rate table.
            index_writer (AsyncDFWriter): Writer for index sample-rate table.
            sample_hits (pd.DataFrame): Annotated sample hits from one chunk.
            healthy_sample_rate (pd.DataFrame): Healthy sample rates for the chunk.
            tumor_sample_rate (pd.DataFrame): Tumor sample rates for the chunk.
            index_sample_rate (pd.DataFrame): Index sample rates for the chunk.
            first_chunk (bool): If True, write the header and open the file fresh.
        """
        annot_writer.write(
            sample_hits,
            ["cts_id", "count", "total", "disease", "developmental_stage", "tissue", "study_id"],
            append=not first_chunk,
            header=first_chunk,
        )
        healthy_writer.write(
            healthy_sample_rate,
            ["cts_id", "developmental_stage", "tissue", "sample_rate"],
            append=not first_chunk,
            header=first_chunk,
        )
        tumor_writer.write(
            tumor_sample_rate,
            ["cts_id", "disease", "tissue", "sample_rate"],
            append=not first_chunk,
            header=first_chunk,
        )
        index_writer.write(
            index_sample_rate,
            ["cts_id", "index_sample_rate", "samples_per_index"],
            append=not first_chunk,
            header=first_chunk,
        )

    def annotate_result_dict(
        self,
        result_dict: dict,
        database: pathlib.Path,
        output_prefix: str,
        cpu: int,
        chunk_size: int,
        compression: bool,
    ) -> None:
        """Annotate all methods in result_dict and write results to disk.

        Iterates over each method in result_dict, runs parallel annotation workers,
        and streams results to async writer threads.

        Args:
            result_dict (dict): Parsed k-mer index results keyed by method name.
            database (pathlib.Path): Path to SQLite annotation database.
            output_prefix (str): Output file prefix.
            cpu (int): Number of parallel workers.
            chunk_size (int): CTS IDs per DataFrame chunk.
            compression (bool): Compress output files with gzip.
        """
        for method_name in result_dict:
            logger.info(f"-> Annotating query results of method: {method_name}")

            ext = ".tsv.gz" if compression else ".tsv"
            output_annotated = pathlib.Path(output_prefix + f"_annotated_{method_name}{ext}")
            output_healthy_rate = pathlib.Path(output_prefix + f"_healthy_sample_rate_{method_name}{ext}")
            output_tumor_rate = pathlib.Path(output_prefix + f"_tumor_sample_rate_{method_name}{ext}")
            output_index_rate = pathlib.Path(output_prefix + f"_index_sample_rate_{method_name}{ext}")
            
            first_chunk = True

            # Start writer threads
            healthy_writer = AsyncDFWriter(output_healthy_rate, compression=compression)
            healthy_writer.start()

            tumor_writer = AsyncDFWriter(output_tumor_rate, compression=compression)
            tumor_writer.start()

            annot_writer = AsyncDFWriter(output_annotated, compression=compression)
            annot_writer.start()

            index_writer = AsyncDFWriter(output_index_rate, compression=compression)
            index_writer.start()


            results = Parallel(n_jobs=cpu, return_as="generator_unordered", pre_dispatch="n_jobs")(
                delayed(Worker.annotator_worker)(this_chunk, self, database)
                for _, _, this_chunk in IndexResultParser.generate_dataframe_in_batches(
                    {method_name: result_dict[method_name]}, batch_size=chunk_size
                )
            )

            # Get batch length and results from chunk and result tuple
            for sample_hits, healthy_sample_rate, tumor_sample_rate, index_sample_rate in results:
                self._write_annotation_batch(
                    annot_writer, 
                    healthy_writer, 
                    tumor_writer,
                    index_writer,
                    sample_hits, 
                    healthy_sample_rate, 
                    tumor_sample_rate,
                    index_sample_rate,
                    first_chunk,
                )
                first_chunk = False  # turn off headers after first write

            logger.info("Waiting for writer threads to finish")
            # Wait for writer threads to finish and close
            annot_writer.wait_until_done()
            healthy_writer.wait_until_done()
            tumor_writer.wait_until_done()
            index_writer.wait_until_done()

            annot_writer.stop()
            healthy_writer.stop()
            tumor_writer.stop()
            index_writer.stop()
