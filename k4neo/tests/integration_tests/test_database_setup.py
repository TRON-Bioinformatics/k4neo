from k4neo.database.database import CreateDataBase
import pathlib
from tinydb import TinyDB, Query


class TestDatabaseIntegration:
    data_set_file = (
        pathlib.Path(__file__).parent.parent / "resources" / "index_data" / "study_annotation.txt"
    )
    tissue_map = (
        pathlib.Path(__file__).parent.parent
        / "resources"
        / "index_data"
        / "tissue_mapping"
        / "mapping.tsv"
    )
    db_file = pathlib.Path(__file__).parent.parent / "resources" / "index_data" / "example.db"
    expected_db = TinyDB(db_file)

    def test_database_setup(self):
        self.db = CreateDataBase(
            db_file=None, data_set_file=self.data_set_file, tissue_map=self.tissue_map, test=True
        )
        self.db.setup_db()
        self.db.precomputations()
        # Check that the four tables created during setup are exactly the same.
        assert (
            self.expected_db.table("PRJNA1234").all() == self.db.database.table("PRJNA1234").all()
        )
        assert (
            self.expected_db.table("sample_study_table").all()
            == self.db.database.table("sample_study_table").all()
        )
        assert (
            self.expected_db.table("tissue_counts").all()
            == self.db.database.table("tissue_counts").all()
        )
        assert (
            self.expected_db.table("tissue_map").all() == self.db.database.table("tissue_map").all()
        )
