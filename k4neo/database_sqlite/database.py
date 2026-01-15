"""k4neo metadata database

This module aims to replace the TinyDB based metadata database implementation
used in k4neo annotation steps. The reason for this module is to enable parallelization
of the annotation class. As TinyDB does not support multiple processes accessing the
database, we decided to port this to SQLITE3. The new database module will not use ORM
and will only be based on standard sqlite3 library. We plan to provide a thin wrapper on
top that somehow acts as the old implementation and provides API compatibility for existing methods.

Todo:
    Better error and value handling
"""

import pathlib
from typing import List, Tuple
import sqlite3
from k4neo.parser.parser import Parser
from k4neo.database_sqlite.validate import validate_tissue_record, validate_sample_record
import pandas as pd
from loguru import logger


class DataBase:
    """
    Generic class representing the k4neo metadata database
    """

    def __init__(self, db_file: pathlib.Path, test: bool = False, timeout: float = 30.0):
        """Parameter initialization

        Args:
            db_file (pathlib.Path):  A path to a database file.
            test (bool, optional): Establish database in-memory for testing. Defaults to False.
            timeout(float, optional): How many seconds should we wait before raising error that table is locked. Defaults to 30s.
        """
        self.db_file = db_file
        self.connection = None
        self.timeout = timeout
        self.test = test

    def connect(self):
        """Establish connection and apply PRAGMA statements"""
        if self.test:
            self.connection = sqlite3.connect(
                ":memory:", timeout=self.timeout, check_same_thread=False
            )
            logger.info("Established in-memory database")
        else:
            self.connection = sqlite3.connect(
                self.db_file, timeout=self.timeout, check_same_thread=False
            )
            logger.info(f"Established connection to: {self.db_file}")
        self.connection.execute("PRAGMA journal_mode=WAL;")
        self.connection.execute("PRAGMA foreign_keys=ON;")

    def close(self):
        if not self.connection is None:
            self.connection.close()
            self.connection = None

    # Required for usage with context manager
    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class CreateDataBase(DataBase):
    """
    Class to construct k4neo metadata database from k4neo index data.
    """

    def __init__(
        self,
        db_file: pathlib.Path,
        data_set_file: pathlib.Path,
        tissue_map: pathlib.Path,
        test: bool = False,
        timeout: float = 30.0,
    ):
        """Generate metadata database

        Args:
            db_file (pathlib.Path): A path to a database file. Can be empty string, when test equals to true.
            data_set_file (pathlib.Path): A file listing documents to insert into database.
            tissue_map (pathlib.Path): A file with accepted tissue identifiers.
            test (bool, optional): Establish database in-memory for testing. Defaults to False.
            timeout (float, optional): Seconds before table is locked error is raised. Defaults to 30.0.
        """
        super().__init__(db_file, test=test, timeout=timeout)
        self.data_set_file = data_set_file
        self.tissue_map = tissue_map

    def create_static_tables(self):
        """
        Create metadata database tables. Currently the database consists of 4 static tables

        * sample_study_mapping: A table mapping samples to study ids
        * tissue_map: A table mapping public tissue identifiers to k4neo tissue identifiers.
        * samples: A table describing indexed samples.
        * tissue_counts: A table with tissue counts per study, diasease and developmental_stage
        * index_information:

        """
        cursor = self.connection.cursor()
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS sample_study_mapping (
                sample_name TEXT,
                study_id TEXT,
                PRIMARY KEY (sample_name, study_id)
            );
            """
        )

        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS tissue_map (
                tissue_public TEXT PRIMARY KEY,
                tissue TEXT NOT NULL,
                subtissue TEXT
            );
            """
        )

        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS samples (
                sample_name TEXT,
                study_id TEXT, 
                runs TEXT NOT NULL,
                tissue TEXT NOT NULL,
                developmental_stage TEXT NOT NULL,
                disease TEXT NOT NULL,
                sex TEXT,
                PRIMARY KEY (sample_name, study_id),
                FOREIGN KEY (tissue)
                    REFERENCES tissue_map(tissue_public)
                    ON UPDATE CASCADE
                    ON DELETE RESTRICT
            );
            """
        )

        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS tissue_counts (
                tissue TEXT NOT NULL,
                study_id TEXT NOT NULL,
                disease TEXT NOT NULL,
                developmental_stage TEXT NOT NULL,
                count INTEGER NOT NULL,
                PRIMARY KEY (tissue, study_id, disease, developmental_stage)
            );
            """
        )

        self.connection.commit()
        cursor.close()

    def insert_sample_table(self, study_records: dict):
        """Insert records into sample table

        Method to insert sample records batchwise into database

        Returns:
            None

        """
        cursor = self.connection.cursor()
        cursor.executemany(
            """
            INSERT INTO samples (sample_name, study_id, runs, tissue, developmental_stage, disease, sex) 
            VALUES (:sample_name, :study_id, :runs, :tissue, :developmental_stage, :disease, :sex)
            """,
            study_records,
        )
        self.connection.commit()
        cursor.close()

    def insert_tissue_table(self, tissue_records: dict):
        """Insert records into tissue table

        Method to insert tissue records batchwise into database

        Returns:
            None

        """
        cursor = self.connection.cursor()
        cursor.executemany(
            """
            INSERT INTO tissue_map (tissue_public, tissue, subtissue) 
            VALUES (:tissue_public, :tissue, :subtissue)
            """,
            tissue_records,
        )
        self.connection.commit()
        cursor.close()

    def _parse_study_table(self) -> Tuple[str, str, int]:
        """
        Parse study table and yield study specific arguments

        Yields:
            Tuple (str): Study name, Path to study annotation table, Number of samples
        """
        with open(self.data_set_file, "r") as file_handle:
            for line in file_handle:
                elements = line.rstrip().split("\t")
                yield elements[0], elements[1], int(elements[2])

    def _add_samples(
        self, study_id: str, study_annot: dict, sample_count: str, update: bool = False
    ):
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
        cursor = self.connection.cursor()

        # Returns a list of
        sample_table = cursor.execute(
            "SELECT * from samples WHERE study_id == ?", (study_id,)
        ).fetchall()
        # Check if sampe table is already initialiazed
        if len(sample_table) == sample_count:
            logger.info(f"-> Samples of {study_id} are already in database")
            return

        if len(sample_table) == 0:
            logger.info(f"-> Inserting samples of study {study_id}")

            self.insert_sample_table(study_annot)

            logger.info(f"-> Loaded {len(study_annot)} documents into study table")

        # Query DB for missing entries
        elif len(sample_table) != sample_count:
            counter = 0
            for element in study_annot:
                exists = cursor.execute(
                    "SELECT * from ? WHERE sample_name == ?", (study_id, element["sample_name"])
                )
                if not exists:
                    cursor.execute(
                        """
                        INSERT INTO samples (sample_name, study_id, runs, tissue, developmental_stage, diesease, age, sex) 
                        VALUES (:sample_name, :study_id, :runs, :tissue, :developmental_stage, :disease, :age, :sex)
                        """,
                        element,
                    )
                    counter += 1
            logger.info(f"-> Loaded {counter} missing documents from {study_id} into sample table")
        else:
            logger.error("-> Don't know what to do here")
        cursor.close()

    def setup_db(self):
        """
        Initialize when establishing database handle
        """
        self.create_static_tables()
        logger.info("-> Adding tissue mapping into database")
        available_tissues = Parser.parse_tissuemap_into_document(self.tissue_map)
        available_tissues = validate_tissue_record(available_tissues)
        self.insert_tissue_table(available_tissues)

        logger.info("-> Adding samples into database")
        for study_id, study_annot, sample_count in self._parse_study_table():
            # If study is not in database parse table into document format
            study_elements = Parser.parse_sample_into_document(study_annot)
            for this_element in study_elements:
                this_element["study_id"] = study_id
            # Update samples with tissue mapping and add subtissue section
            study_elements = validate_sample_record(study_elements)
            if len(study_elements) != sample_count:
                logger.debug(
                    f"Dropped {sample_count - len(study_elements)} samples because of validation"
                )
                sample_count = sample_count - len(study_elements)

            self._add_samples(study_id, study_elements, sample_count)

    def precomputations(self):
        """
        Contains precomputations that would be an unnecessary overhead when computed always on the fly. Should be run
        after inserting all samples
        """
        table = pd.read_sql(
            """SELECT s.study_id, s.disease, s.developmental_stage,
            t.tissue AS tissue
            FROM samples s
            JOIN tissue_map t
            ON s.tissue = t.tissue_public
            """,
            self.connection,
        )
        table = (
            table[["tissue", "developmental_stage", "disease", "study_id"]]
            .value_counts()
            .to_frame()
            .reset_index()
        )
        # Returns record in document format
        table.to_sql(name="tissue_counts", con=self.connection, if_exists="replace", index=False)
        logger.info(f"-> Added {len(table)} precomputed tissue counts for into database")
