import pathlib
import pandas as pd
import pytest
from k4neo.annotator.annotator import Annotator
from k4neo.prepare.prepare import Prepare
from k4neo.database_sqlite.database import CreateDataBase
from k4neo.database_sqlite.queries import Queries


@pytest.fixture
def setup_test_environment(tmp_path):
    """
    Set up a complete test environment with database, prepared sequences, and test data.
    """
    # Create test sequence table
    test_seq_table = tmp_path / "test_sequences.tsv"
    test_seq_table.write_text(
        "junc_id\tcts_id\tcts_seq\tpos\tquery_length\n"
        "test1\tseq1\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t30\t40\n"
        "test2\tseq2\tTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\t30\t40\n"
    )

    # Prepare sequences using Prepare class
    working_dir = tmp_path / "working"
    prepare = Prepare(
        working_dir=working_dir,
        sequence_table=test_seq_table,
        index_kmer_size=21,
    )
    prepare_output = prepare.do_prepare()

    # Create test database
    tissue_map = tmp_path / "tissue.tsv"
    tissue_map.write_text(
        "tissue_description_found_in_public_data\ttissue\tsubtissue\tOntology_UBERON\tHigher_grouping\tstudies\n"
        "brain\tbrain\twhole\tUBERON:0001013\tCNS\tSTUDY1\n"
        "liver\tliver\twhole\tUBERON:0001014\tDigestive\tSTUDY1\n"
        "lung\tlung\twhole\tUBERON:0001015\tRespiratory\tSTUDY2\n"
    )

    study_file = tmp_path / "study.txt"
    study_file.write_text(
        "STUDY1\t{}\t2\n".format(tmp_path / "samples1.tsv")
        + "STUDY2\t{}\t2\n".format(tmp_path / "samples2.tsv")
    )

    samples_file_1 = tmp_path / "samples1.tsv"
    samples_file_1.write_text(
        "sample_name\truns\ttissue\tdevelopmental_stage\tdisease\tsex\n"
        "S1\tR1\tbrain\tadult\thealthy\tM\n"
        "S2\tR2\tliver\tadult\thealthy\tF\n"
    )

    samples_file_2 = tmp_path / "samples2.tsv"
    samples_file_2.write_text(
        "sample_name\truns\ttissue\tdevelopmental_stage\tdisease\tsex\n"
        "S3\tR3\tlung\tadult\thealthy\tM\n"
        "S4\tR4\tlung\tfetal\thealthy\tF\n"
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

    return {
        "prepare_output": prepare_output,
        "working_dir": working_dir,
        "annotation_yaml": working_dir / "annotation_input.yaml",
        "db": db,
        "tmp_path": tmp_path,
    }


def test_annotator_initialization(setup_test_environment):
    """
    Test that Annotator initializes correctly with config YAML.
    """
    env = setup_test_environment
    
    # Initialize Annotator with config YAML
    annotator = Annotator(
        config_yaml=env["annotation_yaml"],
        index_kmer_size=21,
        working_dir=env["working_dir"],
    )

    # Verify initialization
    assert annotator.config is not None
    assert annotator.sequence_table is not None
    assert isinstance(annotator.sequence_table, pd.DataFrame)
    assert len(annotator.sequence_table) > 0
    assert annotator.index_kmer_size == 21
    assert annotator.working_dir == env["working_dir"]

    # Verify sequence table has required columns
    required_columns = ["cts_id", "cts_seq", "pos", "query_length", "query_sequence", "query_cts_id"]
    for col in required_columns:
        assert col in annotator.sequence_table.columns

    env["db"].close()


def test_annotator_count_aggregation(setup_test_environment):
    """
    Test the _count_aggregation method.
    """
    env = setup_test_environment
    
    annotator = Annotator(
        config_yaml=env["annotation_yaml"],
        index_kmer_size=21,
    )

    # Create mock parsed results
    parsed_results = pd.DataFrame({
        "cts_id": ["seq1", "seq1", "seq1", "seq2", "seq2"],
        "study_id": ["STUDY1", "STUDY1", "STUDY2", "STUDY1", "STUDY2"],
        "sample_name": ["S1", "S2", "S3", "S1", "S4"],
        "tissue": ["brain", "liver", "lung", "brain", "lung"],
        "developmental_stage": ["adult", "adult", "adult", "adult", "fetal"],
        "disease": ["healthy", "healthy", "healthy", "healthy", "healthy"],
    })

    result = annotator._count_aggregation(parsed_results)

    # Verify result structure
    assert "count" in result.columns
    assert "cts_id" in result.columns
    assert "study_id" in result.columns
    assert "tissue" in result.columns
    assert "developmental_stage" in result.columns
    assert "disease" in result.columns

    # Verify counts are aggregated correctly
    assert result["count"].dtype == int
    assert result["count"].min() >= 0

    env["db"].close()


def test_annotator_split_found(setup_test_environment):
    """
    Test the _split_found method that separates non-detected sequences.
    """
    env = setup_test_environment
    
    # Create mock parsed results with some sequences not found
    parsed_results = pd.DataFrame({
        "cts_id": ["seq1", "seq2", "seq3"],
        "sample_name": ["S1", None, "S3"],
    })

    not_expressed = Annotator._split_found(parsed_results)

    # Verify structure
    assert "cts_id" in not_expressed.columns
    assert "count" in not_expressed.columns
    assert "total" in not_expressed.columns
    assert "disease" in not_expressed.columns
    assert "developmental_stage" in not_expressed.columns
    assert "tissue" in not_expressed.columns
    assert "study_id" in not_expressed.columns

    # Verify only non-found sequences are included
    assert len(not_expressed) == 1
    assert not_expressed.iloc[0]["cts_id"] == "seq2"
    assert not_expressed.iloc[0]["count"] == 0
    assert not_expressed.iloc[0]["total"] == 0

    env["db"].close()


def test_annotator_calculate_healthy_sample_rate(setup_test_environment):
    """
    Test the _calculate_healthy_sample_rate method.
    """
    env = setup_test_environment
    
    # Create mock data
    parsed_results = pd.DataFrame({
        "cts_id": ["seq1", "seq1", "seq2"],
        "developmental_stage": ["adult", "adult", "adult"],
        "tissue": ["brain", "liver", "brain"],
        "disease": ["healthy", "healthy", "healthy"],
        "count": [2, 3, 1],
    })

    tissue_counts = pd.DataFrame({
        "study_id": ["STUDY1", "STUDY1", "STUDY2"],
        "disease": ["healthy", "healthy", "healthy"],
        "developmental_stage": ["adult", "adult", "adult"],
        "tissue": ["brain", "liver", "lung"],
        "total": [10, 20, 15],
    })

    result = Annotator._calculate_healthy_sample_rate(parsed_results, tissue_counts)

    # Verify result structure
    assert "sample_rate" in result.columns
    assert "cts_id" in result.columns
    assert "tissue" in result.columns
    assert "developmental_stage" in result.columns
    assert "samples_per_tissue" in result.columns

    # Verify sample rates are calculated
    assert result["sample_rate"].notna().all()
    assert (result["sample_rate"] >= 0).all()
    assert (result["sample_rate"] <= 1).all()

    env["db"].close()


def test_annotator_calculate_tumor_sample_rate(setup_test_environment):
    """
    Test the _calculate_tumor_sample_rate method.
    """
    env = setup_test_environment
    
    # Create mock data
    parsed_results = pd.DataFrame({
        "cts_id": ["seq1", "seq1", "seq2"],
        "disease": ["primary solid tumor", "primary solid tumor", "metastatic"],
        "tissue": ["LUAD", "BRCA", "SKCM"],
        "count": [5, 10, 3],
    })

    tissue_counts = pd.DataFrame({
        "study_id": ["TCGA", "TCGA", "TCGA"],
        "disease": ["primary solid tumor", "primary solid tumor", "metastatic"],
        "tissue": ["LUAD", "BRCA", "SKCM"],
        "total": [100, 200, 50],
    })

    result = Annotator._calculate_tumor_sample_rate(parsed_results, tissue_counts)

    # Verify result structure
    assert "sample_rate" in result.columns
    assert "cts_id" in result.columns
    assert "tissue" in result.columns
    assert "disease" in result.columns
    assert "index_count" in result.columns

    # Verify sample rates are calculated
    assert result["sample_rate"].notna().all()
    assert (result["sample_rate"] >= 0).all()
    assert (result["sample_rate"] <= 1).all()

    env["db"].close()


def test_annotator_annotate_cts(setup_test_environment):
    """
    Test the annotate_cts method with mock data.
    """
    env = setup_test_environment
    
    annotator = Annotator(
        config_yaml=env["annotation_yaml"],
        index_kmer_size=21,
    )

    queries = Queries(env["db"])

    # Create mock parsed results
    parsed_results = pd.DataFrame({
        "cts_id": ["seq1", "seq1", "seq2", "seq3"],
        "sample_name": ["S1", "S2", None, "S3"],
    })

    result = annotator.annotate_cts(parsed_results, queries)

    # Verify result structure
    assert isinstance(result, pd.DataFrame)
    assert "cts_id" in result.columns
    assert "count" in result.columns
    assert "total" in result.columns
    assert "tissue" in result.columns
    assert "developmental_stage" in result.columns
    assert "disease" in result.columns
    assert "study_id" in result.columns

    # Verify all sequences are included (found and not found)
    unique_seqs = result["cts_id"].unique()
    assert "seq1" in unique_seqs
    assert "seq2" in unique_seqs
    assert "seq3" in unique_seqs

    env["db"].close()


def test_annotator_annotate_sequences(setup_test_environment):
    """
    Test the annotate_sequences method that merges annotations with original sequences.
    """
    env = setup_test_environment
    
    annotator = Annotator(
        config_yaml=env["annotation_yaml"],
        index_kmer_size=21,
    )

    # Create mock annotated results
    annotated_cts = pd.DataFrame({
        "cts_id": [annotator.sequence_table.iloc[0]["query_cts_id"]],
        "count": [5],
        "total": [10],
        "tissue": ["brain"],
        "developmental_stage": ["adult"],
        "disease": ["healthy"],
        "study_id": ["STUDY1"],
    })

    result = annotator.annotate_sequences(annotated_cts)

    # Verify result structure
    assert isinstance(result, pd.DataFrame)
    assert "cts_id" in result.columns
    assert "cts_seq" in result.columns
    assert "query_sequence" in result.columns
    assert "count" in result.columns
    assert "total" in result.columns
    assert "tissue" in result.columns

    # Verify data was merged correctly
    assert len(result) > 0
    assert result["total"].dtype in [int, 'int64']

    env["db"].close()


def test_annotator_annotate_sample_rate2(setup_test_environment):
    """
    Test the annotate_sample_rate2 method.
    """
    env = setup_test_environment
    
    annotator = Annotator(
        config_yaml=env["annotation_yaml"],
        index_kmer_size=21,
    )

    queries = Queries(env["db"])

    # Create mock annotated results
    annotated_cts = pd.DataFrame({
        "cts_id": [annotator.sequence_table.iloc[0]["query_cts_id"]],
        "count": [1],
        "total": [1],
        "tissue": ["brain"],
        "developmental_stage": ["adult"],
        "disease": ["healthy"],
        "study_id": ["STUDY1"],
    })

    healthy_rate, tumor_rate = annotator.annotate_sample_rate2(annotated_cts, queries, min_total=1)

    # Verify healthy sample rate result
    assert isinstance(healthy_rate, pd.DataFrame)
    assert "sample_rate" in healthy_rate.columns
    assert "cts_id" in healthy_rate.columns
    assert "tissue" in healthy_rate.columns
    assert "developmental_stage" in healthy_rate.columns
    assert len(healthy_rate) > 0

    # Verify tumor sample rate result
    assert isinstance(tumor_rate, pd.DataFrame)
    assert "sample_rate" in tumor_rate.columns
    assert "cts_id" in tumor_rate.columns
    assert "tissue" in tumor_rate.columns
    assert "disease" in tumor_rate.columns

    # Verify sample rates are within valid range
    if len(healthy_rate) > 0:
        print(healthy_rate)
        assert (healthy_rate["sample_rate"] >= 0).all()
        assert (healthy_rate["sample_rate"] <= 1).all()

    if len(tumor_rate) > 0:
        assert (tumor_rate["sample_rate"] >= 0).all()
        assert (tumor_rate["sample_rate"] <= 1).all()

    env["db"].close()
