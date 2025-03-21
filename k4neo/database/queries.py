import pandas as pd
from tinydb import TinyDB, Query
from k4neo.database.database import DataBase
import numpy as np
from loguru import logger


class Queries:
    """
    Generic query class that provides methods to query the document database and
    add document information to annotation table
    """

    def __init__(self, db: DataBase):
        self.db = db.database
        self.query = Query()

    def get_project_id(self, sample_name):
        """
        For a given sample_name search the matching study_id. Studies
        are organized in tables requiring to query the matching study table for each sample
        :param sample_name:
        :return study_id
        """
        table = self.db.table("sample_study_table")
        project_id = table.get(self.query.sample_name == sample_name)
        if project_id is None:
            return
        return project_id.get("study_id")

    def get_sample_study(self):
        table = self.db.table("sample_study_table")
        return pd.DataFrame.from_dict(table)

    def annotate_samples_of_project(self, samples: pd.DataFrame):
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
        # Get corresponding study table
        table = self.db.table(study)
        sample_query = samples.sample_name.unique().tolist()
        # Query project table for study ids and convert returned document to DataFrame
        annotation = table.search(self.query.sample_name.one_of(sample_query))
        # In case there is noting matching the query
        if not annotation:
            samples["tissue"] = np.nan
            samples["developmental_stage"] = np.nan
            samples["disease"] = np.nan
            return samples
        # Merge group df with annotation features
        annotation = pd.DataFrame.from_dict(annotation)
        samples = pd.merge(samples, annotation, on="sample_name", how="left")
        return samples

    def document_to_pd(self, document: dict):
        """
        Convert document to pandas dataframe
        :param document:
        :return:
        """
        df = pd.DataFrame.from_dict(document)
        return df

    def annotate_tissue_counts(self, samples: pd.DataFrame):
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
        table = self.db.table("tissue_counts")
        tissue_counts = table.search(self.query.study_id == study_id)
        # If query doesn't return a match initialize empty df as return
        if not tissue_counts:
            return pd.DataFrame(columns=samples.columns.tolist() + ["total"])

        tissue_counts = pd.DataFrame.from_dict(tissue_counts)
        tissue_counts.rename(columns={"count": "total"}, inplace=True)
        samples = pd.merge(samples, tissue_counts, how="left")
        return samples

    def get_tissue_counts(self):
        """
        Return precomputed total tissue counts per developmental state.
        """
        table = self.db.table("tissue_counts")
        tissue_counts = pd.DataFrame.from_dict(table)
        tissue_counts.rename(columns={"count": "total"}, inplace=True)
        return tissue_counts
