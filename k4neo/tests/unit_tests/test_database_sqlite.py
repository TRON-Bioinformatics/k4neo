import pytest
import pandas as pd
from k4neo.database_sqlite.database import CreateDataBase


@pytest.fixture
def tissue_records():
    return [
        {"tissue_public": "A", "tissue": "B", "subtissue": "C"},
        {"tissue_public": "D", "tissue": "E", "subtissue": "F"},
    ]


@pytest.fixture
def sample_records():
    return [
        # Valid sample
        {
            "sample_name": "S1",
            "study_id": "STUDY1",
            "runs": "R1",
            "tissue": "A",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": "M",
        },
        # Valid sample (missing optional 'sex')
        {
            "sample_name": "S2",
            "study_id": "STUDY1",
            "runs": "R2",
            "tissue": "D",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": None,
        },
        # Invalid sample (missing required 'tissue')
        {
            "sample_name": "S3",
            "study_id": "STUDY2",
            "runs": "R3",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": "F",
            "tissue": "A",
        },
    ]


@pytest.fixture
def sample_records_different_disease_and_developmental_stage():
    # Expected counts:
    # Aggregated tissue counts: 2 A (tissue B) healthy, 1 A (tissue B) tumor, 1 D (tissue E) healthy, 1 D (tissue E) fetal
    return [
        # Valid sample
        {
            "sample_name": "S1",
            "study_id": "STUDY1",
            "runs": "R1",
            "tissue": "A",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": "M",
        },
        # Valid sample (missing optional 'sex')
        {
            "sample_name": "S2",
            "study_id": "STUDY1",
            "runs": "R2",
            "tissue": "A",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": None,
        },
        {
            "sample_name": "S3",
            "study_id": "STUDY1",
            "runs": "R3",
            "tissue": "A",
            "developmental_stage": "adult",
            "disease": "tumor",
            "sex": None,
        },
        {
            "sample_name": "S12",
            "study_id": "STUDY2",
            "runs": "R12",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": "F",
            "tissue": "D",
        },
        {
            "sample_name": "S13",
            "study_id": "STUDY2",
            "runs": "R13",
            "developmental_stage": "fetal",
            "disease": "healthy",
            "sex": "F",
            "tissue": "D",
        },
    ]


@pytest.fixture
def sample_record_single():
    return [
        # Valid sample
        {
            "sample_name": "S1",
            "study_id": "STUDY1",
            "runs": "R1",
            "tissue": "A",
            "developmental_stage": "adult",
            "disease": "healthy",
            "sex": "M",
        },
    ]


@pytest.fixture
def in_memory_db(tmp_path):

    # Set up dummy database with no data
    db = CreateDataBase(
        db_file=None,
        data_set_file=tmp_path / "study.txt",
        tissue_map=tmp_path / "tissue.tsv",
        test=True,
    )
    db.connect()
    db.create_static_tables()
    yield db
    db.close()


def test_create_tables(in_memory_db):
    """Test that the necessary tables are created in the database."""
    cursor = in_memory_db.connection.cursor()
    tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
    table_names = {t[0] for t in tables}
    expected = {
        "sample_study_mapping",
        "tissue_map",
        "samples",
        "tissue_counts",
        "aggregated_tissue_counts",
    }
    assert expected.issubset(table_names)


def test_insert_and_query_tissue_map(in_memory_db, tissue_records):
    """Test inserting and querying tissue_map table."""
    in_memory_db.insert_tissue_table(tissue_records)
    df = pd.read_sql("SELECT * FROM tissue_map", in_memory_db.connection)
    assert df.shape[0] == 2
    assert set(df["tissue_public"]) == {"A", "D"}
    assert set(df["tissue"]) == {"B", "E"}
    assert set(df["subtissue"]) == {"C", "F"}


def test_insert_and_query_samples(in_memory_db, tissue_records, sample_record_single):
    # Required for FOREIGN KEY constraint
    in_memory_db.insert_tissue_table(tissue_records)

    in_memory_db.insert_sample_table(sample_record_single)
    df = pd.read_sql("SELECT * FROM samples", in_memory_db.connection)
    assert df.shape[0] == 1
    assert df.iloc[0]["sample_name"] == "S1"
    assert len(df.columns) == 7


def test_insert_and_query_valid_samples(in_memory_db, tissue_records, sample_records):
    # Required for FOREIGN KEY constraint
    in_memory_db.insert_tissue_table(tissue_records)
    in_memory_db.insert_sample_table(sample_records)
    df = pd.read_sql("SELECT * FROM samples", in_memory_db.connection)
    assert df.shape[0] == 3
    assert len(df.columns) == 7

    assert df.iloc[0]["sample_name"] == "S1"
    assert df.iloc[0]["sex"] == "M"

    assert df.iloc[1]["sample_name"] == "S2"
    assert df.iloc[1]["sex"] is None

    assert df.iloc[2]["sample_name"] == "S3"
    assert df.iloc[2]["sex"] == "F"


def test_precomputation(in_memory_db, tissue_records, sample_records):
    in_memory_db.insert_tissue_table(tissue_records)
    in_memory_db.insert_sample_table(sample_records)

    # Runs precomputations and stores results in the database
    in_memory_db.precomputations()

    df = pd.read_sql("SELECT * FROM tissue_counts", in_memory_db.connection)
    assert df.shape[0] == 3
    # Study specific counts
    assert df.loc[(df["study_id"] == "STUDY1") & (df["tissue"] == "B"), "count"].iloc[0] == 1
    assert df.loc[(df["study_id"] == "STUDY1") & (df["tissue"] == "E"), "count"].iloc[0] == 1
    assert df.loc[(df["study_id"] == "STUDY2") & (df["tissue"] == "B"), "count"].iloc[0] == 1

    df = pd.read_sql("SELECT * FROM aggregated_tissue_counts", in_memory_db.connection)
    assert df.shape[0] == 2
    # Aggregated counts
    assert (
        df.loc[(df["developmental_stage"] == "adult") & (df["tissue"] == "B"), "count"].iloc[0] == 2
    )
    assert (
        df.loc[(df["developmental_stage"] == "adult") & (df["tissue"] == "E"), "count"].iloc[0] == 1
    )


def test_precomputation_with_different_developmental_stages_and_disease_states(
    in_memory_db, tissue_records, sample_records_different_disease_and_developmental_stage
):
    in_memory_db.insert_tissue_table(tissue_records)
    in_memory_db.insert_sample_table(sample_records_different_disease_and_developmental_stage)

    # Runs precomputations and stores results in the database
    in_memory_db.precomputations()

    df = pd.read_sql("SELECT * FROM tissue_counts", in_memory_db.connection)
    assert df.shape[0] == 4
    # Study specific counts
    assert (
        df.loc[
            (df["study_id"] == "STUDY1") & (df["tissue"] == "B") & (df["disease"] == "healthy"),
            "count",
        ].iloc[0]
        == 2
    )
    assert (
        df.loc[
            (df["study_id"] == "STUDY1") & (df["tissue"] == "B") & (df["disease"] == "tumor"),
            "count",
        ].iloc[0]
        == 1
    )
    assert (
        df.loc[
            (df["study_id"] == "STUDY2")
            & (df["tissue"] == "E")
            & (df["developmental_stage"] == "adult"),
            "count",
        ].iloc[0]
        == 1
    )
    assert (
        df.loc[
            (df["study_id"] == "STUDY2")
            & (df["tissue"] == "E")
            & (df["developmental_stage"] == "fetal"),
            "count",
        ].iloc[0]
        == 1
    )

    df = pd.read_sql("SELECT * FROM aggregated_tissue_counts", in_memory_db.connection)
    assert df.shape[0] == 4
    # Aggregated counts
    assert (
        df.loc[
            (df["developmental_stage"] == "adult")
            & (df["tissue"] == "B")
            & (df["disease"] == "healthy"),
            "count",
        ].iloc[0]
        == 2
    )
    assert (
        df.loc[
            (df["developmental_stage"] == "adult")
            & (df["tissue"] == "B")
            & (df["disease"] == "tumor"),
            "count",
        ].iloc[0]
        == 1
    )
