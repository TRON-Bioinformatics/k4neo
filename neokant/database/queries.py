import sys
import os
import pandas as pd
from tinydb import TinyDB, Query
from neokant.database.database import DataBase


class Queries:
    def __init__(self, db: DataBase):
        self.db = db.database
        self.query = Query()

    def get_project_id(self, sample_name):
        table = self.db.table("sample_study_table")
        project_id = table.get(self.query.sample_name == sample_name)
        if project_id is None:
            return
        return project_id.get("study_id")

    def annotate_samples_of_project(self, samples: pd.DataFrame):
        study = samples['study_id'].unique()
        assert len(study) == 1, "More than one project in samples dataframe. Function supports only query for one sample"
        study = study.item()
        table = self.db.table(study)
        sample_query = samples.sample_name.unique().tolist()
        # Query project table for study ids and convert returned document to DataFrame
        annotation = table.search(self.query.sample_name.one_of(sample_query))
        annotation = pd.DataFrame.from_dict(annotation)
        # Merge group df with annotation features
        samples = pd.merge(samples, annotation, on="sample_name", how='left')
        return samples

    def document_to_pd(self, document: dict):
        df = pd.DataFrame.from_dict(document)
        return df


## Deprecated code
'''
    def get_sample_study_table(self):
        table = self.db.table("sample_study_table")
        return pd.DataFrame(table)

    def get_cts_index_hits(self, samples_per_study):
        cts_result = []
        for study, sample in samples_per_study:
            table = self.db.table(study)
            result = table.get(self.query.sample_name == sample)
            if result:
                cts_result.append(result)
        return pd.DataFrame(result)

    def tissue_count(self):
        table = self.db.table("tissue_counts")
        return pd.DataFrame(table)

    def get_sample_from_study(self, study, sample):
        pass

    def get_study_sample_mapping(self, sample_names: list):
        """
        Generate a mapping of sample tables to sample names
        {study_id : [id1, id2, id3]}
        """
        study_matching = {}
        sample_studies = [x
            for x in self.db.table("sample_study_table").get(self.query.sample_name.any(sample_names))]
        for sample in sample_studies:
            study = sample["project_id"]
            sample_name = sample["sample_name"]
            if study not in study_matching:
                study_matching[study] = [sample_name]
            else:
                study_matching[study].append(sample_name)
        return study_matching

    def get_index_hit(self, index_ids: list):
        matches = []
        sample_study_mapping = self.get_study_sample_mapping(index_ids)

        for study, samples in sample_study_mapping.items():
            table = self.db.table(study)
            for sample in samples:
                matches.append(table.search(self.query.sample_name == sample))
        return matches
'''