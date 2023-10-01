import os
import sys
import csv
import pandas as pd
from logzero import logger
from neokant.parser import EXPECTED_SAMPLE_COLUMNS, EXPECTED_TISSUE_COLUMNS, SUPPORTED_TOOLS


class Parser:
    @staticmethod
    def _read_tissue_map(tissue_map):
        """
        Read neokant tissue map to homogenize tissue identifiers and make them compatible with ExBio.
        :return:
        """
        tissue_map = []
        with open(tissue_map) as file_handle:
            reader = csv.DictReader(file_handle, delimiter="\t")
            for line in reader:
                assert all([x in line.keys() for x in EXPECTED_TISSUE_COLUMNS]), "Missing columns"
                tissue_map.append({'tissue_public': line['tissue_public'],
                                  'tissue': line['tissue'],
                                  'subtissue': line['subtissue']})
        return tissue_map

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
                 reindeer_sample_mapping:str = None,
                 kmindex_cutoff: float = 0.7) -> None:

        self.indexing_table = indexing_table
        assert tool in SUPPORTED_TOOLS, f"Selected method '{tool}' not supported by neokant..."
        self.tool = tool
        # Parameters specific for prediction tools
        self.reindeer_sample_mapping = reindeer_sample_mapping
        if self.tool == "reindeer":
            assert self.reindeer_sample_mapping is not None and self.reindeer_sample_mapping != "",\
                "Parsing REINDEER results requires a sample/index mapping file..."
        self.kmindex_cutoff = kmindex_cutoff

    def parse_results(self) -> dict:
        result = pd.DataFrame()
        match self.tool:
            case "reindeer":
                logger.info("Parsing REINDEER index query...")
                result = self._parse_reindeer()
            case "cobs":
                logger.info("Parsing COBS index query...")
                result = self._parse_cobs()
            case "kmindex":
                logger.info("Parsing KMINDEX index query...")
                result = self._parse_kmindex()
            case "raptor":
                logger.info("Parsing RAPTOR index query...")
                result = self._parse_raptor()
            case _:
                print("Tool is unknown. Cannot parse results")
        return result

    def _parse_reindeer(self) -> dict:
        """
        Parse reindeer results into dictionary format showing mapping of
        
        """
        dataset_mapping = {}
        results = {}
        with open(self.reindeer_sample_mapping, 'r') as file_handle:
            logger.info("Reading sample index mapping table")
            for line in file_handle:
                elements = line.rstrip().split("\t")
                dataset_mapping[int(elements[0])] = elements[1]
        with open(self.indexing_table) as file_handle:
            for line in file_handle:
                elements = line.rstrip().split()
                query_name = elements[0].lstrip(">")
                found_list = []
                results[query_name] = []
                for i, ele in enumerate(elements[1:]):
                    for item in ele.split(","):
                        if item == "*":
                            break

                        key, val = item.split(":")
                        if val != "*":
                            found_list.append(dataset_mapping[i])
                            break
                results[query_name].extend(found_list)       
        logger.info(f"Parsed {len(results)} target sequences")
        return results
    
    def _parse_kmindex(self) -> dict:
        results = {}
        logger.info(f"Using {self.kmindex_cutoff} as cutoff to determine presence/absence of target sequences...")
        with open(self.indexing_table) as file_handle:
            reader = csv.DictReader(file_handle, delimiter='\t')
            for line in reader:
                cts_id = line["samples"].split(":")[1]
                results[cts_id] = []
                for sample, prediction in line.items():
                    if sample == "samples":
                        continue
                    prediction = float(prediction)
                    if prediction >= self.kmindex_cutoff:
                        found = True
                    else:
                        found = False
                    results[cts_id].append(sample)
        logger.info(f"Parsed {len(results)} target sequences")
        return results

    def _parse_raptor(self) -> dict:
        """
        Parse raptor search results into dictionary with cts_ids as keys and samples
        as values.

        :return: Result mapping
        """
        dataset_mapping = {}
        results = {}
        
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
                            results[cts_id].append(dataset_mapping[int(this_sample)])

        logger.info(f"Parsed {len(results)} target sequences")
        return results

    def _parse_cobs(self) -> dict:
        results = {}
        cts_id = ""
        with open(self.indexing_table) as file_handle:
            for line in file_handle:
                elements = line.rstrip().split('\t')
                if line.startswith("*"):
                    cts_id = elements[0].lstrip('*')
                    results[cts_id] = set()
                # If the line does not start with an asteriks it is a sample hit
                else:
                    results[cts_id].add([elements[0]])

        logger.info(f"Parsed {len(results)} target sequences")
        return results


                
