import pathlib
from typing import List, Tuple
from k4neo.parser.parser import Parser
from tinydb import TinyDB, Query
from tinydb.storages import MemoryStorage
import pandas as pd
from loguru import logger


class DataBase:
    """
    Generic class representing the k4neo database
    """

    def __init__(self, db_file: pathlib.Path, test: bool = False):
        """Parameter initialization

        Args:
            db_file (pathlib.Path):  A path to a TinyDB database file.
            test (bool, optional): Establish TinyDB database in-memory for testing. Defaults to False.
        """
        self.db_file = db_file
        if test:
            self.database = TinyDB(storage=MemoryStorage)
        else:
            self.database = TinyDB(self.db_file, sort_keys=True, indent=4)


class CreateDataBase(DataBase):
    def __init__(self, db_file, data_set_file, tissue_map, test: bool = False):
        """
        Database initialization
        """
        super().__init__(db_file, test=test)
        self.data_set_file = data_set_file
        self.tissue_map = tissue_map

    def _parse_study_table(self):
        """
        Parse study table and yield study specific arguments
        """
        with open(self.data_set_file, "r") as file_handle:
            for line in file_handle:
                elements = line.rstrip().split("\t")
                yield elements[0], elements[1], elements[2]

    def _sample_study_table(self, study_id: str, study_annot: dict):
        """
        Create an table mapping sample_names to study_ids. This table
        is required in the annotation step to query the correct table for
        each found sample. TinyDB has table support however we cannot connect
        the tables with keys.
        """
        sample_study_query = Query()
        sample_study_table = self.database.table("sample_study_table")
        for this_sample in study_annot:
            exists = sample_study_table.contains(
                (sample_study_query.sample_name == this_sample["sample_name"])
                & (sample_study_query.study_id == study_id)
            )
            if not exists:
                sample_study_table.insert(
                    {"sample_name": this_sample["sample_name"], "study_id": study_id}
                )

    def _add_samples(self, study_id: str, study_annot: dict, sample_count: str):
        """
        Check consistency of metadata database. Compare in-memory database to sample tables shipped
        with k4neo. Add or update documents in database.

        Args:
            study_id (str): A unique study identifier for samples to be added.
            study_annot (dict): A document representation of the sample metadata.
            sample_count (str): Number of samples in study.
        Returns:
            None
        ToDo:
            - Simplify arguments
            - Simplify function
        """

        study_query = Query()
        study_table = self.database.table(study_id)
        # Check if sampe table is already initialiazed
        if len(study_table) == sample_count:
            logger.info(f"-> Sample table {study_id} already in database")
            return
        if len(study_table) == 0:
            logger.info(f"-> Initializing table for study {study_id}")
            study_table.insert_multiple(study_annot)
            logger.info(f"-> Loaded {len(study_annot)} documents into study table")
        # Query DB for missing entries
        elif len(study_table) != sample_count:
            counter = 0
            for element in study_annot:
                exists = study_table.contains(study_query.sample_name == element["sample_name"])
                if not exists:
                    study_table.insert(element)
                    counter += 1
            logger.info(f"-> Loaded {counter} missing documents into {study_id} table")
        else:
            logger.error("-> Don't know what to do here")

    def _add_tissues(self):
        """
        Add a general tissue map to database to match public tissue identifiers to
        k4neo tissue terminology. This mapping is currently defined by the k4neo
        index data table. However, this should be shipped as part of the package.
        """
        counter = 0
        tissue_table = self.database.table("tissue_map")
        tissue_query = Query()
        available_tissues = Parser.parse_tissuemap_into_document(self.tissue_map)
        for element in available_tissues:
            exists = tissue_table.contains(
                (tissue_query.tissue_public == element["tissue_public"])
                & (tissue_query.tissue == element["tissue"])
                & (tissue_query.subtissue == element["subtissue"])
            )
            if not exists:
                tissue_table.insert(element)
                counter += 1
        logger.info(f"-> Loaded {counter} tissue map documents into tissue_map table")

    def _update_sample_document_with_tissue(self, sample: dict) -> Tuple[dict, bool]:
        """

        Update the document representation of a sample to use the k4neo tissue identifiers
        instead of public tissue descriptions.

        Args:
            sample (dict): A document representation of a sample.

        Returns:
            Tuple[dict, bool]: A tuple containing the updated dict with new keys for tissue and subtissue
                and a boolean indicting if a tissue match was detected/the dict was updated.
        """
        ## Query database for tissue
        ## Match if not leave out
        tissue_query = Query()
        tissue_map = self.database.table("tissue_map")
        if len(tissue_map) == 0:
            logger.warning("-> Tissue map is not initialized. Can not update tissue identifier")
            return sample, False

        tissue_public = sample["tissue"]
        # Find tissue match of public tissue identifier
        tissue_match = tissue_map.get(tissue_query.tissue_public == tissue_public)
        if not tissue_match:
            logger.warning(
                f"-> Could not find for sample {sample['sample_name']} a tissue match. Ignoring for annotation"
            )
            return sample, False

        sample["tissue"] = tissue_match.get("tissue")
        sample["subtissue"] = tissue_match.get("subtissue")
        return sample, True

    def setup_db(self):
        """
        Initialize when establishing database handle
        """
        logger.info("-> Adding tissue mapping into database")
        self._add_tissues()
        logger.info("-> Adding samples into database")
        for study_id, study_annot, sample_count in self._parse_study_table():
            # If study is not in database parse table into document format
            study_elements = Parser.parse_sample_into_document(study_annot)
            # Update samples with tissue mapping and add subtissue section
            study_elements = [self._update_sample_document_with_tissue(x) for x in study_elements]
            # Drop samples without a tissue match
            study_elements = [x[0] for x in study_elements if x[1]]
            self._sample_study_table(study_id, study_elements)
            self._add_samples(study_id, study_elements, sample_count)

    def precomputations(self):
        """
        Contains precomputations that would be an unnecessary overhead when compuzted always on the fly. Should be run
        after inserting all samples
        :return:
        """
        tissue_count_table = self.database.table("tissue_counts")
        tissue_count_query = Query()
        for study_id, _, _ in self._parse_study_table():
            exists = tissue_count_table.contains(tissue_count_query.study_id == study_id)
            if exists:
                logger.info(
                    f"-> Precomputed counts for {study_id} already in database. Skipping calculation"
                )
                continue

            table = pd.DataFrame(self.database.table(study_id))
            table["study_id"] = study_id
            table = table[["tissue", "developmental_stage", "disease", "study_id"]].value_counts()
            # Returns record in document format
            tissue_counts = table.to_frame().reset_index().to_dict(orient="records")
            tissue_count_table.insert_multiple(tissue_counts)
            logger.info(
                f"-> Added {len(tissue_counts)} precomputed tissue counts for {study_id} into database"
            )
