import pandas as pd
from tinydb import Query
from k4neo.database2.database import DataBase
import numpy as np


class Queries:
    """
    Generic query class that provides methods to retrieve specific data
    SQL metadata database. It provides methods to retrieve data as pd.DataFrames
    and defines which methods the k4neo.Annotator requires for annotation.
    """

    def __init__(self, db: DataBase):
        """Parameter initialization

        Args:
            db (DataBase): An instance of the k4neo DataBase.
        """
        self.db = db
        if self.db.connection is None:
            raise RuntimeError("k4neo.Queries needs k4neo.DataBase with established connection")

    def get_project_id(self, sample_name: str) -> str:
        """
        For a given sample_name search the matching study_id. Studies
        are organized in tables requiring to query the matching study table for each sample

        Args:
            sample_name (str): _description_

        Returns:
            str: _description_
        """

        project_id = pd.read_sql(
            f"SELCT study_id FROM samples WHERE sample_name == ?",
            self.db.connection,
            params=[
                sample_name,
            ],
        )
        if project_id.empty:
            return
        return project_id.get("study_id", None)

    def get_sample_study(self) -> pd.DataFrame:
        table = pd.read_sql("SELECT sample_name, study_id FROM samples", self.db.connection)
        return table

    def annotate_samples_of_project(self, samples: pd.DataFrame) -> pd.DataFrame:
        """
        Given a dataframe containing all sample hist of a single study, query the database
        for sample documents and annotate each sample with tissue, developmental_stage and disease
        :param samples:
        :return:
        """
        study = samples["study_id"].unique()
        assert (
            len(study) == 1
        ), "More than one project in samples dataframe. Function supports only query for one study"
        study = study.item()
        sample_query = samples.sample_name.unique().tolist()
        # annotation = table.search(self.query.sample_name.one_of(sample_query))
        placeholders = ",".join("?" for _ in sample_query)
        annotation = pd.read_sql_query(
            f"""SELECT s.sample_name,
                s.developmental_stage,
                s.disease, 
                tm.tissue
                FROM samples s
                LEFT JOIN tissue_map tm ON s.tissue = tm.tissue_public
                WHERE sample_name IN ({placeholders});
            """,
            self.db.connection,
            params=sample_query,
        )
        # In case there is noting matching the query
        if annotation.empty:
            samples["tissue"] = np.nan
            samples["developmental_stage"] = np.nan
            samples["disease"] = np.nan
            return samples
        # Merge group df with annotation features
        samples = pd.merge(samples, annotation, on="sample_name", how="left")
        return samples

    def document_to_pd(self, document: dict) -> pd.DataFrame:
        """
        Convert document to pandas dataframe
        :param document:
        :return:
        """
        df = pd.DataFrame.from_dict(document)
        return df

    def annotate_tissue_counts(self, samples: pd.DataFrame) -> pd.DataFrame:
        """
        Given a dataframe containing all tissue hits in a single study, query the database
        for precomputed tissue counts and annotate
        :param samples:
        :return:
        """

        study = samples["study_id"].unique()
        assert len(study) == 1, "Dataframe can only contain one study to search for "
        study_id = study.item()
        # Get corresponding study table
        tissue_counts = pd.read_sql_query(
            "SELECT * FROM tissue_counts WHERE study_id = ?",
            self.db.connection,
            params=[
                study_id,
            ],
        )
        # If query doesn't return a match initialize empty df as return
        if tissue_counts.empty:
            return pd.DataFrame(columns=samples.columns.tolist() + ["total"])

        tissue_counts.rename(columns={"count": "total"}, inplace=True)
        samples = pd.merge(samples, tissue_counts, how="left")
        breakpoint()
        return samples

    def get_tissue_counts(self) -> pd.DataFrame:
        """
        Return precomputed total tissue counts per developmental state.
        """
        tissue_counts = pd.read_sql_query("SELECT * FROM tissue_counts;", self.db.connection)
        tissue_counts.rename(columns={"count": "total"}, inplace=True)
        return tissue_counts
