import pytest
import pandas as pd
from k4neo.database_sqlite.database import CreateDataBase
from k4neo.database_sqlite.queries import Queries

@pytest.fixture
def setup_db(tmp_path):
    db = CreateDataBase(
        db_file=None,
        data_set_file=tmp_path / "study.txt",
        tissue_map=tmp_path / "tissue.tsv",
        test=True,
    )
    db.connect()
    db.create_static_tables()
    db.insert_tissue_table([
        {"tissue_public": "A", "tissue": "B", "subtissue": "C"},
        {"tissue_public": "D", "tissue": "E", "subtissue": "F"},
    ])
    db.insert_sample_table([
        {
            "sample_name": "S1",
            "study_id": "STUDY1",
            "runs": "R1",
            "tissue": "A",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": "M",
        },
        {
            "sample_name": "S2",
            "study_id": "STUDY1",
            "runs": "R2",
            "tissue": "D",
            "developmental_stage": "adult",
            "disease": "tumor",
            "sex": "F",
        },
    ])
    db.precomputations()
    yield db
    db.close()

def test_get_project_id(setup_db):
    queries = Queries(setup_db)
    result = queries.get_project_id("S1")
    assert result.iloc[0] == "STUDY1"
    result = queries.get_project_id("S2")
    assert result.iloc[0] == "STUDY1"

    result = queries.get_project_id("S3")
    assert result is None

def test_get_sample_study(setup_db):
    queries = Queries(setup_db)
    df = queries.get_sample_study()
    assert set(df.columns) == {"sample_name", "study_id"}
    assert df.shape[0] == 2

def test_annotate_samples_of_project(setup_db):
    queries = Queries(setup_db)
    df = pd.DataFrame([{"sample_name": "S1", "study_id": "STUDY1"}])
    annotated = queries.annotate_samples_of_project(df)
    
    assert "tissue" in annotated.columns
    assert "developmental_stage" in annotated.columns
    assert "disease" in annotated.columns
    
    assert annotated.loc[0, "tissue"] == "B"
    assert annotated.loc[0, "developmental_stage"] == "adult"
    assert annotated.loc[0, "disease"] == "healthy"

    # Check for a sample not in the database returns NaN annotations
    df = pd.DataFrame([{"sample_name": "S3", "study_id": "STUDY3"}])
    annotated = queries.annotate_samples_of_project(df)
    
    assert "tissue" in annotated.columns
    assert "developmental_stage" in annotated.columns
    assert "disease" in annotated.columns
    
    assert pd.isna(annotated.loc[0, "tissue"]) is True
    assert pd.isna(annotated.loc[0, "developmental_stage"]) is True
    assert pd.isna(annotated.loc[0, "disease"]) is True

def test_document_to_pd():
    from k4neo.database_sqlite.queries import Queries
    doc = {"a": [1, 2], "b": [3, 4]}
    df = Queries.document_to_pd(None, doc)
    assert isinstance(df, pd.DataFrame)
    assert set(df.columns) == {"a", "b"}

def test_annotate_tissue_counts(setup_db):
    queries = Queries(setup_db)
    df = pd.DataFrame([
        { 
            "sample_name": "S1",
            "study_id": "STUDY1",
            "runs": "R1",
            "tissue": "A",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": "M",
        }
    ])
    annotated = queries.annotate_tissue_counts(df)
    assert "total" in annotated.columns
    assert annotated.shape[0] == 1

    queries = Queries(setup_db)
    df = pd.DataFrame([
        { 
            "sample_name": "S3",
            "study_id": "STUDY3",
            "runs": "R3",
            "tissue": "A",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": "M",
        }
    ])
    annotated = queries.annotate_tissue_counts(df)
    assert "total" in annotated.columns
    assert annotated.shape[0] == 0
    assert pd.isna(annotated["total"]).all()

def test_get_tissue_counts(setup_db):
    queries = Queries(setup_db)
    df = queries.get_tissue_counts()
    
    assert df.shape[0] >= 1
    
    assert set(df.columns).issuperset(
        {"tissue", "developmental_stage", "disease", "total", "study_id"}
    )

    assert df.loc[
        (df["study_id"] == "STUDY1") & (df["tissue"] == "B"), "total"
    ].iloc[0] == 1

    assert df.loc[
        (df["study_id"] == "STUDY1") & (df["tissue"] == "E"), "total"
    ].iloc[0] == 1