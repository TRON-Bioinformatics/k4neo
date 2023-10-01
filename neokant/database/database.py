from neokant.parser.parser import Parser
from tinydb import TinyDB, Query
from logzero import logger as loggy


class DataBase:
    def __init__(self, db_file, data_set_file, tissue_map):
        """
        Database initialization
        """
        self.db_file = db_file
        self.database = TinyDB(self.db_file)
        self.data_set_file = data_set_file
        self.tissue_map = Parser.parse_tissuemap_into_document(tissue_map)

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
            exists = sample_study_table.contains((sample_study_query.sample_name == this_sample["sample_name"]) &
                                                 (sample_study_query.project_id == this_sample["project_id"]))
            if not exists:
                sample_study_query.insert({x: this_sample[x] for x in ("sample_name", "project_id")})

    def _add_samples(self, study_id: str, study_annot: dict, sample_count: str):
        """
        Check consistency of metadata database. Compare loaded database to sample tables shipped
        with neokant. Add or update documents in database.
        """
        study_query = Query()
        study_table = self.database.table(study_id)
        # Check if sampe table is already initialiazed
        if len(study_table) == sample_count:
            loggy.info(f"Sample table {study_id} already in database")
            return
        if len(study_table) == 0:
            loggy.info(f"Initializing table for study {study_id}")
            study_table.insert_multiple(study_annot)
            loggy.info(f"Loaded {len(study_annot)} documents into study table")
        # Query DB for missing entries
        elif len(study_table) != sample_count:
            counter = 0
            for element in study_annot:
                exists = study_table.contains(study_query.index_name == element["index_name"])
                if not exists:
                    study_table.insert(element)
                    counter += 1
            loggy.info(f"Loaded {counter} missing documents into {study_id} table")
        else:
            loggy.error("Don't know what to do here")

    def _add_tissues(self):
        counter = 0
        tissue_table = self.database.table('tissue_map')
        tissue_query = Query()
        for element in self.tissue_map:
            exists = tissue_table.contains(
                (tissue_query.tissue_public == element['tissue_public']) &
                (tissue_query.tissue_public == element['tissue']) &
                (tissue_query.tissue_public == element['subtissue'])
            )
            if not exists:
                tissue_table.insert(element)
                counter += 1
        loggy.info(f"Loaded {counter} tissue map documents into tissue_map table")

    @staticmethod
    def _update_sample_document_with_tissue(sample: dict, tissue_map: list) -> dict:
        tissue_public = sample['tissue']
        tissue_match = {}
        found = False
        i = 0
        # Find the first hit in the list
        while not found:
            if tissue_map[i]['tissue_public'] == tissue_public:
                found = True
                tissue_match = tissue_map[i]
            i += 1
        if not tissue_match:
            loggy.error(f"Could not find for sample {sample['index_name']} a tissue match. Ignoring for annotation")
            return

        sample['tissue'] = tissue_match['tissue']
        sample['subtissue'] = tissue_match['subtissue']
        return sample

    def _init_database(self):
        """
        Initialize when establishing database handle
        """
        loggy.info("Adding tissue mapping into database")
        self._add_tissues()
        loggy.info("Adding samples into database")
        for study_id, study_annot, sample_count in self._parse_study_table():
             # If study is not in database parse table into document format
            study_elements = Parser.parse_into_document(study_annot)
            # Update samples with tissue mapping and add subtissue section
            study_elements = [DataBase._update_sample_document_with_tissue(x, self.tissue_map) for x in study_elements]
            self._sample_study_table(study_id, study_elements)
            self._add_samples(study_id, study_elements, sample_count)

    def setup_db(self):
        """
        Prepare document db to work on it
        :return:
        """
        self._init_database()


class Query:
    def __init__(self, db: DataBase):
        self.db = db
        self.query = Query()

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

    def get_tissue_per_sample(self, sample_name,):
        matches = self.db.get(self.query.sample_name == sample_name)
        try:
            tissue = matches["tissue"]
        except KeyError as error:
            loggy.error("Tissue key is missing")
