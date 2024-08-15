import os
import sys
import csv
import pandas as pd
from logzero import logger
import pyarrow
import pyarrow.parquet as pq
from typing import Any
from collections.abc import Iterator
from k4neo.parser import EXPECTED_SAMPLE_COLUMNS, EXPECTED_TISSUE_COLUMNS, SUPPORTED_TOOLS


class Parser:
    @staticmethod
    def _read_tissue_map(tissue_map):
        """
        Read k4eno tissue map to homogenize tissue identifiers and make them compatible with GTEx.
        :return:
        """
        tissue_mapping = []
        with open(tissue_map) as file_handle:
            reader = csv.DictReader(file_handle, delimiter="\t")
            for line in reader:
                assert all([x in line.keys() for x in EXPECTED_TISSUE_COLUMNS]), "Missing columns"
                tissue_mapping.append({'tissue_public': line['tissue_description_found_in_public_data'],
                                       'tissue': line['tissue'],
                                       'subtissue': line['subtissue']})
        return tissue_mapping

    @staticmethod
    def parse_tissuemap_into_document(tissue_map):
        return Parser._read_tissue_map(tissue_map)

    @staticmethod
    def _read_sample_file(data_table):
        """
        Read neokant sample sheets
        :return: file_content
        """
        file_content = []
        with open(data_table, 'r') as file_handle:
            reader = csv.DictReader(file_handle, delimiter='\t')
            for line in reader:
                assert all([x in line.keys() for x in EXPECTED_SAMPLE_COLUMNS]), "Missing columns"
                file_content.append(line)

        return file_content


    @staticmethod
    def parse_sample_into_document(data_table):
        """
        Parse sample sheets into document format. Make sure required fieldnames are present to ensure
        compatibility with document db.
        :return: Set of dictionaries
        """
        return Parser._read_sample_file(data_table)

class IndexResultParser:
    """
    Class provides functions to parse table formats
    returned by kmer indexing tools and map to sample names
    if required.
    """
    def __init__(self, 
                 indexing_table: str, 
                 tool: str,
                 raptor_sample_mapping:str = None,
                 kmindex_cutoff: float = 0.7) -> None:

        self.indexing_table = indexing_table
        assert tool in SUPPORTED_TOOLS, f"Selected method '{tool}' not supported by neokant..."
        self.tool = tool
        # Parameters specific for prediction tools
        self.raptor_sample_mapping = raptor_sample_mapping
        if self.tool == "raptor":
            assert self.raptor_sample_mapping is not None and self.raptor_sample_mapping != "",\
                "Parsing Raptor results requires a sample/index mapping file"
        self.kmindex_cutoff = kmindex_cutoff

    def parse_results(self) -> dict:
        result = {}
        match self.tool:
            case "kmindex":
                logger.info("Parsing KMINDEX index query results...")
                result = self._parse_kmindex()
            case "raptor":
                logger.info("Parsing RAPTOR index query results...")
                result = self._parse_raptor()
            case _:
                logger.error("Tool is unknown. Cannot parse results")
        return result

    @staticmethod
    def write_result(results: dict, out_file: str):
        """
        Write parsed results into tabular format to be processed by user
        :param results: A dictionary qith cts as key and sample hits as value
        :param out_file: Output file
        :return:
        """
        with open(out_file, "w") as file_handle:
            for query, samples in results.items():
                for sample in samples:
                    line = f"{query}\t{sample}\n"
                    file_handle.write(line)
    
    def _parse_kmindex(self) -> dict:
        results = {}
        logger.info(f"Using {self.kmindex_cutoff} as cutoff to determine presence/absence of target sequences...")
        with open(self.indexing_table) as file_handle:
            reader = csv.DictReader(file_handle, delimiter='\t')
            for line in reader:
                cts_id = line["samples"].split(":")[1]
                results[cts_id] = []
                sample_count = 0
                for sample, prediction in line.items():
                    if sample == "samples":
                        continue
                    prediction = float(prediction)
                    if prediction >= self.kmindex_cutoff:
                        results[cts_id].append(sample)
                        sample_count += 1
                if sample_count == 0:
                    results[cts_id] = [None]

        logger.info(f"Parsed {len(results)} target sequences")
        return results

    def _parse_raptor(self) -> dict:
        """
        Parse raptor search results into dictionary with cts_ids as keys and samples
        as values.

        :return: Result mapping
        """
        dataset_mapping = {}
        sample_name_mapping = {}
        results = {}

        with open(self.raptor_sample_mapping, 'r') as file_handle:
            logger.info("Reading sample/minimiser mapping file to match raptor bin ids to sample_names")
            reader = csv.DictReader(file_handle, delimiter="\t")
            for row in reader:
                sample_name_mapping[row['minimiser_id']] = row['sample_name']

        with open(self.indexing_table) as file_handle:
            for line in file_handle:
                elements = line.rstrip().split("\t")
                # Skip config section returned in raptor output file
                if line.startswith("##"):
                    continue
                # Get key/sample mapping from header
                elif line.startswith("#"):
                    # Stop collecting samples when header of results section start
                    if elements[0] != "#QUERY_NAME":
                        dataset_mapping[int(elements[0][1:])] = \
                            os.path.basename(elements[1]).rstrip().rstrip('.minimiser')
                else:
                    elements = line.rstrip().split('\t')
                    cts_id = elements[0]
                    results[cts_id] = []
                    if len(elements) > 1:
                        for this_sample in elements[1].split(","):
                            results[cts_id].append(sample_name_mapping[dataset_mapping[int(this_sample)]])
                    else:
                        results[cts_id] = [None]

        logger.info(f"Parsed {len(results)} target sequences")
        return results

class IndexResultParser:
    """
    Class provides functions to parse table formats
    returned by kmer indexing tools and map to sample names
    if required.
    """
    def __init__(self, 
                 query_tables: dict[str]) -> None:

        self.query_table = query_tables

    def parse_results(self) -> dict:
        result = self.parse_parquet(self.query_table)
        return result

    @staticmethod
    def write_result(results: dict, out_file: str):
        """
        Write parsed results into tabular format to be processed by user
        :param results: A dictionary qith cts as key and sample hits as value
        :param out_file: Output file
        :return:
        """
        with open(out_file, "w") as file_handle:
            for query, samples in results.items():
                for sample in samples:
                    line = f"{query}\t{sample}\n"
                    file_handle.write(line)
    @staticmethod
    def parse_parquet():
        query_results = {}
        for this_batch in _iterate_parquet_in_batches(path, batch_size):
            cts_id = ""
            detected_samples = []
            for key, value in this_batch.items():
                if key == "":
                    cts_id = ""
                else:
                    if value is None:
                        continue
                    detected_samples.append(key)
            query_results[cts_id] = detected_samples
        return query_results


    @staticmethod
    def _iterate_parquet_in_batches(path: pathlib.Path, batch_size = 100):
        
        parquet_file = pq.ParquetFile(path)
        for batch in parquet_file.iter_batches(batch_size=batch_size):
            yield from batch.to_pylist()
    