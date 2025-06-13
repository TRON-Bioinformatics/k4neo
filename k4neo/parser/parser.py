import csv
import pandas as pd
import pathlib
from multiprocessing import Pool
from collections import ChainMap, defaultdict
from functools import partial
from loguru import logger
import pyarrow.parquet as pq
from typing import Any
from k4neo.parser import (
    EXPECTED_SAMPLE_COLUMNS,
    EXPECTED_TISSUE_COLUMNS,
    SUPPORTED_TOOLS,
)


class Parser:
    """
    Class to provide generic parser functions
    """

    @staticmethod
    def _read_tissue_map(tissue_map: pathlib.Path) -> dict:
        """Read k4neo tissue map

        The k4neo index metadata library provides a mapping of public tissue
        identifiers to internal tissue definitions.

        Args:
            tissue_map (pathlib.Path): The path to the tissue mapping file

        Returns:
            dict: A mapping of public tissue identifiers to k4neo tissue definitions
        """
        tissue_mapping = []
        with open(tissue_map, "r") as file_handle:
            reader = csv.DictReader(file_handle, delimiter="\t")
            for line in reader:
                assert all(
                    [x in line.keys() for x in EXPECTED_TISSUE_COLUMNS]
                ), "Missing columns in input"
                tissue_mapping.append(
                    {
                        "tissue_public": line[
                            "tissue_description_found_in_public_data"
                        ],
                        "tissue": line["tissue"],
                        "subtissue": line["subtissue"],
                    }
                )
        return tissue_mapping

    @staticmethod
    def parse_tissuemap_into_document(tissue_map):
        """
        High level function to import tissue map
        """
        return Parser._read_tissue_map(tissue_map)

    @staticmethod
    def _read_sample_file(data_table: pathlib.Path) -> list[dict]:
        """
        Reader for k4neo annotation files. Make sure required fieldnames are present to ensure
        compatibility with document db.
        """
        file_content = []
        with open(data_table, "r") as file_handle:
            reader = csv.DictReader(file_handle, delimiter="\t")
            for line in reader:
                assert all(
                    [x in line.keys() for x in EXPECTED_SAMPLE_COLUMNS]
                ), "Missing columns in input file"
                file_content.append(line)

        return file_content

    @staticmethod
    def parse_sample_into_document(data_table: pathlib.Path):
        """
        High level function to import study tissue sheets.
        """
        return Parser._read_sample_file(data_table)


class IndexResultParser:
    """
    Class provides functions to parse table formats
    returned by kmer indexing pipeline and map to sample names
    if required.
    """

    def __init__(self, query_pipeline_results, cores: int = 8) -> None:

        self.query_tables = query_pipeline_results.query_path
        self.cores = cores

    def parse_results(self, kmer_ratio=0.7) -> dict:
        """Parse results of QueryPipeline

        Args:
            kmer_ratio (float, optional):
                Required fraction of shared k-mers between query and sample. Only required for kmindex results.
                Defaults to 0.7.

        Returns:
            dict: A dictionary containing parsed results of k-mer methods.
            For example:
                {'raptor': {cts: set(P1,P2,P3), cts_2: set(P1)},
                 'kimindex': {cts: set(P2,P3), cts_2: set(P1,P2)}
                }
        """
        result = {}
        logger.info("-> Parsing parquet file of k-mer query pipeline")
        for this_method, result_parquet in self.query_tables.items():
            logger.info(f"-> Parsing query results of method: {this_method}")
            parsed_result = self.parse_parquet(
                result_parquet, self.cores, method=this_method, kmer_ratio=kmer_ratio
            )
            result[this_method] = parsed_result
        return result

    @staticmethod
    def write_result(results: dict, out_file: str):
        """Write parsed results into tabular text format

        Args:
            results (dict): Parsed k-mer query results
            out_file (str): Path to output file
        """
        logger.info(f"-> Writing results to file: {out_file}")
        with open(out_file, "w") as file_handle:
            for query, samples in results.items():
                for sample in samples:
                    line = f"{query}\t{sample}\n"
                    file_handle.write(line)

    @staticmethod
    def parse_parquet(
        path: pathlib.Path,
        batch_size: int = 1000,
        cores: int = 8,
        method: str = "raptor",
        kmer_ratio: float = 0.7,
    ) -> dict:
        """Parse k-mer pipeline output

        The k-mer pipeline saves index hits in Apache Parquet format.
        Depending on the index size, it is not feasible to hold the dataframe
        in memory. Therefore, this function iterates in chunks over the parquet
        files and stores results in a dictionary mapping context sequences to samples.

        Args:
            path:
                File path of parquet file.
            batch_size:
                Iterate over parquet file yielding this many records at once.

        Returns:
            A dict mapping each search sequence (cts) to detected samples.

        """
        query_results = defaultdict(set)
        with Pool(processes=cores) as pool:
            for this_batch in IndexResultParser._iterate_parquet_in_batches(
                path, batch_size
            ):
                detected_samples_list = pool.map(
                    partial(
                        IndexResultParser._parse_table_row,
                        method=method,
                        kmer_ratio=kmer_ratio,
                    ),
                    this_batch,
                )
                for parse_dict in detected_samples_list:
                    for this_cts, this_sample_set in parse_dict.items():
                        query_results[this_cts].update(this_sample_set)
                #query_results.extend(detected_samples_list)
        
        #query_results = ChainMap(*query_results)
        return query_results

    @staticmethod
    def _parse_table_row(table_row: dict, method: str, kmer_ratio: float) -> dict:
        """Parse a dataframe row

        Parse a row of a dataframe stored in Apache Parquet file. Extract
        context sequence name and samples that matched the query in the k-mer index.

        Args:
            table_row (dict): A dict represenation of pyarrow RecordBatch.
            method (str): The k-mer method used to generate the table.
            kmer_ratio (float): Required fraction of shared k-mers between query and sample.

        Returns:
            dict: A dict containing all samples (column names) matching the query (based on k_mer ratio)
        """
        query_results = {}
        cts_id = ""
        detected_samples = set()
        for this_key, this_value in table_row.items():
            if this_key == "__index_level_0__":
                cts_id = this_value
            else:
                # If method kmindex and less than ratio k-mers are shared -> skip
                # If method raptor and value is None -> Skip
                if method == "kmindex":
                    if this_value < kmer_ratio:
                        continue
                elif method == "raptor":
                    if this_value is None:
                        continue
                detected_samples.add(this_key)
        if not detected_samples:
            query_results[cts_id] = set(
                [
                    None,
                ]
            )
        else:
            query_results[cts_id] = detected_samples

        return query_results

    @staticmethod
    def _iterate_parquet_in_batches(path: pathlib.Path, batch_size=1000) -> list:
        """ "Iterate over a parquet file in batches.

        Each row in the "dataframe" is represented by a dict.

        Args:
            path: File path of parquet file..
            batch_size:
                Iterate over parquet file yielding this many records at once.

        Yields:
            Rows of Table/RecordBatch as a list of dictionaries.
        """

        parquet_file = pq.ParquetFile(path)
        for batch in parquet_file.iter_batches(batch_size=batch_size):
            yield batch.to_pylist()
