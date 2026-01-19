import os
import sys
from csv import DictReader
from collections import defaultdict
import pandas as pd
from loguru import logger
from multiprocessing import Pool
from typing import Dict, Set, Generator, Tuple
import itertools
from joblib import Parallel, delayed


class IndexResultParser2:
    """
    Class provides functions to parse table formats
    returned by kmer indexing pipeline and map to sample names
    if required.
    """

    def __init__(self, query_pipeline_results, cores: int = 8) -> None:

        self.query_tables = query_pipeline_results
        self.cores = cores

    @staticmethod
    def parse_results_of_kmer_search(
        this_method, this_index, this_sample_mapping, this_result_path, kmer_ratio
    ):

        kmer_parser = BinaryKmerIndexResultParser(
            search_results=this_result_path,
            method=this_method,
            raptor_sample_mapping=this_sample_mapping,
            kmer_ratio=kmer_ratio,
        )
        return this_method, kmer_parser.parse_results()

    def parse_result2(self, kmer_ratio=0.7) -> dict:
        """Parse results of QueryPipeline

        Args:
            kmer_ratio (float, optional):
                Required fraction of shared k-mers between query and sample.
                Only required for kmindex results. Defaults to 0.7.

        Returns:
            dict: A dictionary containing parsed results of k-mer methods.
            For example:
                {'raptor': {cts: set(P1,P2,P3), cts_2: set(P1)},
                 'kimindex': {cts: set(P2,P3), cts_2: set(P1,P2)}
                }
        """
        query_results = defaultdict(lambda: defaultdict(lambda: {None}))
        logger.debug("Parsing k-mer result files in parallel")

        results = Parallel(n_jobs=self.cores, backend="multiprocessing")(
            delayed(IndexResultParser2.parse_results_of_kmer_search)(m, i, s, r, kmer_ratio)
            for m, i, s, r in self.query_tables
        )

        for this_method, detected_samples in results:
            for this_cts, this_sample_set in detected_samples.items():
                IndexResultParser2.update_sample_set(
                    query_results[this_method][this_cts], this_sample_set
                )

        return query_results

    def parse_results(self, kmer_ratio=0.7) -> dict:
        """Parse results of QueryPipeline

        Args:
            kmer_ratio (float, optional):
                Required fraction of shared k-mers between query and sample.
                Only required for kmindex results. Defaults to 0.7.

        Returns:
            dict: A dictionary containing parsed results of k-mer methods.
            For example:
                {'raptor': {cts: set(P1,P2,P3), cts_2: set(P1)},
                 'kimindex': {cts: set(P2,P3), cts_2: set(P1,P2)}
                }
        """
        query_results = defaultdict(lambda: defaultdict(lambda: {None}))
        logger.debug("Parsing k-mer result files")
        # [("method", "subindex_name", "subindex_mapping", "result_path")]
        for this_method, this_index, this_sample_mapping, this_result_path in self.query_tables:
            logger.debug(f"Parsing {this_index} results obtained with {this_method}")
            kmer_parser = BinaryKmerIndexResultParser(
                search_results=this_result_path,
                method=this_method,
                raptor_sample_mapping=this_sample_mapping,
                kmer_ratio=kmer_ratio,
            )
            detected_samples = kmer_parser.parse_results()

            for this_cts, this_sample_set in detected_samples.items():
                query_results[this_method][this_cts].update(this_sample_set)

        return query_results

    @staticmethod
    def generate_dataframe_in_batches(
        parsed_results: Dict[str, Dict[str, Set[str]]], batch_size: int = 10000
    ) -> Generator[Tuple[str, int, pd.DataFrame], None, None]:
        for method_name, cts_dict in parsed_results.items():
            it = iter(cts_dict.items())
            while True:
                batch = list(itertools.islice(it, batch_size))
                if not batch:
                    break
                # Generate cts/sample generators for dataframe creation
                rows = ((cts, sample) for cts, samples in batch for sample in samples)
                df = pd.DataFrame.from_records(rows, columns=["cts_id", "sample_name"])
                yield method_name, len(batch), df

    @staticmethod
    def update_sample_set(target_set: set, new_set: set):
        """Remove placeholder None

        This method can be used to remove the classic “None” as a placeholder problem.
        Each k-mer index result is parsed idenpendently using None as placeholder if the cts was
        not detected in any sample of the subindex. If the sequence is missing in just one of the indices, the
        placeholder will also be present in the set of detected samples.

        Clean up the set of samples for each cts so that None is removed if there is at least one “real” sample name
        in the parsed results.

        Args:
            target_set (set): The set to update with sample names.
            new_set (set): A set with sample names from parsing.
        Returns:
            set: The modified target_set without placeholder if detected in at least one sample or a placeholder set.

        """
        pre_existing = len(target_set - {None})
        new_entries = new_set - {None}

        target_set.update(new_entries)

        if pre_existing == 0 and new_entries:
            target_set.discard(None)


class BinaryKmerIndexResultParser:
    """
    Class provides functions to parse table formats
    returned by presence/absence kmer indexing tools and map to sample names
    if required.
    """

    def __init__(
        self,
        search_results: str,
        method: str,
        raptor_sample_mapping: str = None,
        kmer_ratio: float = 0.7,
    ) -> None:

        self.search_results = search_results
        self.method = method
        # Parameters specific for prediction tools
        self.raptor_sample_mapping = raptor_sample_mapping
        if self.method == "raptor":
            assert (
                self.raptor_sample_mapping is not None and self.raptor_sample_mapping != ""
            ), "Parsing Raptor results requires a sample/index mapping file"
        self.kmer_ratio = kmer_ratio

    def parse_results(self) -> pd.DataFrame:
        """
        Choose parsing method for specified method
        """
        result = {}
        match self.method:
            case "kmindex":
                logger.debug("-> Parsing KMINDEX index query results...")
                result = self._parse_kmindex()
            case "raptor":
                logger.debug("-> Parsing RAPTOR index query results...")
                result = self._parse_raptor()
            case _:
                logger.error(f"Tool {self.method} is unknown. Cannot parse results")
        return result

    def _parse_kmindex(self) -> pd.DataFrame:
        """
        Parse tabular output format of kmindex
        """
        results = {}
        with open(self.search_results) as file_handle:
            reader = DictReader(file_handle, delimiter="\t")
            for line in reader:
                cts_id = line["samples"].split(":")[1]
                detected_samples = set()
                for sample, prediction in line.items():
                    if sample == "samples":
                        continue
                    prediction = round(float(prediction), 2)
                    # k-mer fraction smaller than minimum.
                    if prediction < self.kmer_ratio:
                        continue
                    detected_samples.add(sample)
                if not detected_samples:
                    results[cts_id] = set(
                        [
                            None,
                        ]
                    )
                else:
                    results[cts_id] = detected_samples

        logger.debug(f"Parsed {len(results.keys())} query sequences from kmindex output.")
        return results

    def _parse_raptor(self) -> pd.DataFrame:
        """
        Parse raptor search results

        :return: Result mapping
        """
        dataset_mapping = {}
        sample_name_mapping = {}
        results = {}

        with open(self.raptor_sample_mapping, "r") as file_handle:
            logger.debug(
                "Reading minimiser2sample mapping file to match raptor bin ids to sample names"
            )
            reader = DictReader(file_handle, delimiter="\t")
            for row in reader:
                sample_name_mapping[row["minimiser_id"]] = row["sample_name"]
        with open(self.search_results) as file_handle:
            for line in file_handle:
                elements = line.rstrip().split("\t")
                # Skip config section returned in raptor output file -> Starts with '##'
                if line.startswith("##"):
                    continue
                # Get key/sample mapping from header
                elif line.startswith("#"):
                    # Stop collecting samples when header of results section start
                    if elements[0] != "#QUERY_NAME":
                        dataset_mapping[int(elements[0][1:])] = elements[1].rstrip()
                # Parse index hits -> Each query has a number of bins assigned Q1   3,4,5,6
                else:
                    elements = line.rstrip().split("\t")
                    cts_id = elements[0]
                    results[cts_id] = {}
                    detected_samples = set()
                    # If at least one bin is reported by raptor
                    if len(elements) > 1:
                        # Iterate over samples bins
                        for this_sample in elements[1].split(","):
                            # Save detected bins and directly translate into sample identifier
                            detected_samples.add(
                                sample_name_mapping[dataset_mapping[int(this_sample)]]
                            )
                    if not detected_samples:
                        results[cts_id] = set(
                            [
                                None,
                            ]
                        )
                    else:
                        results[cts_id] = detected_samples

        logger.info(f"Parsed {len(results.keys())} query sequences from raptor output")
        return results


class QuantitativeKmerIndexParser:
    """
    Class provides functions to parse results from quantitative indices
    Currently limited to JellyFish
    """

    def __init__(self, search_results: str, method: str) -> None:

        self.search_results = search_results
        self.method = method

    def parse_jellyfish(self) -> dict:
        cts_kmer_count = defaultdict(list)
        result = pd.DataFrame(
            columns=[
                "cts_id",
                "median_kmer_count",
                "mean_kmer_count",
                "max_kmer_count",
                "min_kmer_count",
                "rate_non_zero_kmers",
                "rate_zero_kmers",
                "variance",
                "cv",
            ]
        )

        with open(self.search_results, "r") as file_handle:
            for this_line in file_handle:
                cts_id, _, count = this_line.rstrip().split("\t")
                cts_kmer_count[cts_id].append(int(count))

        return cts_kmer_count
