import pytest
import pandas as pd
from k4neo.database_sqlite.database import CreateDataBase
from k4neo.database_sqlite.queries import Queries

@pytest.fixture
def test_queries_integration(tmp_path):
    # Create temporary database
    # Two studies with two samples each
    # Two tissues in tissue map

    tissue_map = tmp_path / "tissue.tsv"
    study_file = tmp_path / "study.txt"
    tissue_map.write_text(
        (
            "tissue_description_found_in_public_data\ttissue\tsubtissue\tOntology_UBERON\tHigher_grouping\tstudies\n"
            "A\tB\tC\tUBERON:0001013\tGroup1\tSTUDY1\n"
            "D\tE\tF\tUBERON:0001014\tGroup2\tSTUDY1\n"
        )
    )

    study_file.write_text(
        "".join([
            "STUDY1\t{}\t2\n".format(tmp_path / "samples.tsv"),
            "STUDY2\t{}\t2\n".format(tmp_path / "samples2.tsv")
        ])
    )
    
    samples_file = tmp_path / "samples.tsv"
    samples_file.write_text(
        "sample_name\truns\ttissue\tdevelopmental_stage\tdisease\tsex\n"
        "S1\tR1\tA\tadult\thealthy\tM\n"
        "S2\tR2\tD\tadult\ttumor\tF\n"
    )

    samples_file_2 = tmp_path / "samples2.tsv"
    samples_file_2.write_text(
        "sample_name\truns\ttissue\tdevelopmental_stage\tdisease\tsex\n"
        "S3\tR3\tD\tadult\thealthy\tM\n"
        "S4\tR4\tD\tadult\thealthy\tF\n"
    )

    db = CreateDataBase(
        db_file=None,
        data_set_file=study_file,
        tissue_map=tissue_map,
        test=True,
    )
    db.connect()
    db.create_static_tables()
    db.setup_db()
    db.precomputations()
    yield db
    db.close()

def test_queries_runtime_error(tmp_path):
    
    tissue_map = tmp_path / "tissue.tsv"
    study_file = tmp_path / "study.txt"
    
    db = CreateDataBase(
        db_file=None,
        data_set_file=study_file,
        tissue_map=tissue_map,
        test=True,
    )
    with pytest.raises(RuntimeError):
        Queries(db=db)

def test_queries_functions(test_queries_integration):
    
    queries = Queries(test_queries_integration)

    # Test get_sample_study
    sample_study = queries.get_sample_study()
    assert sample_study.shape[0] == 4

    # Test annotate_samples_of_project
    annotated = queries.annotate_samples_of_project(sample_study[sample_study["study_id"] == "STUDY1"])
    assert "tissue" in annotated.columns
    assert set(annotated["tissue"]) == {"B", "E"}

    annotated = queries.annotate_samples_of_project(sample_study[sample_study["study_id"] == "STUDY2"])
    assert "tissue" in annotated.columns
    assert set(annotated["tissue"]) == {"E", "E"}

    # Test get_tissue_counts
    tc = queries.get_tissue_counts()
    assert "total" in tc.columns
    assert tc["total"].sum() == 4

    # Test annotate_tissue_counts
    #annotated_tc = queries.annotate_tissue_counts(annotated)
    #assert "total" in annotated_tc.columns
    #assert annotated_tc.shape[0] == 
