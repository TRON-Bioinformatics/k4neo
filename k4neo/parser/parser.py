import csv
import pathlib
from k4neo.parser import (
    EXPECTED_SAMPLE_COLUMNS,
    EXPECTED_TISSUE_COLUMNS,
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
                        "tissue_public": line["tissue_description_found_in_public_data"],
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
