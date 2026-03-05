import pathlib
import pandas as pd
import yaml
from k4neo.prepare.prepare import Prepare, PrepareOutput


def test_prepare(tmp_path):
    """
    Test that Prepare class correctly processes sequences and generates required output files.
    """
    test_sequence_table = (
        pathlib.Path(__file__).parent.parent / "resources" / "example_seq_table.tsv"
    )

    # Initialize Prepare with tmp_path as working directory
    prepare = Prepare(working_dir=tmp_path, sequence_table=test_sequence_table, index_kmer_size=21)
    result = prepare.do_prepare()

    assert isinstance(result, PrepareOutput)

    # Check that all expected output files exist
    assert result.query_fasta.exists(), "Query FASTA file should exist"
    assert result.seq_to_short_output.exists(), "Seq-to-short output file should exist"
    assert result.sequence_table_output.exists(), "Sequence table output file should exist"
    assert result.working_dir.exists(), "Working directory should exist"

    yaml_output = tmp_path / "annotation_input.yaml"
    assert yaml_output.exists(), "annotation_input.yaml should be created"

    # Verify content of annotation_input.yaml
    with open(yaml_output, "r") as file_handle:
        yaml_content = yaml.safe_load(file_handle)

    assert "query_fasta" in yaml_content
    assert "seq_to_short_output" in yaml_content
    assert "sequence_table_output" in yaml_content
    assert "working_dir" in yaml_content

    # Verfiy FASTA file for pipeline was created
    with open(result.query_fasta, "r") as file_handle:
        fasta_content = file_handle.read()

    assert fasta_content.count(">") > 0, "FASTA file should contain sequences"

    # Verfiy sequence table
    seq_table = pd.read_csv(result.sequence_table_output, sep="\t")

    required_columns = [
        "cts_id",
        "cts_seq",
        "pos",
        "query_length",
        "query_sequence",
        "query_cts_id",
    ]
    for col in required_columns:
        assert col in seq_table.columns, f"Column {col} should be in sequence table"

    assert seq_table["query_sequence"].notna().all(), "All sequences should have query_sequence"
    assert seq_table["query_cts_id"].notna().all(), "All sequences should have query_cts_id"

    # Check that sequences are not too short
    min_length = 21 + 4
    assert (
        seq_table["query_sequence"].str.len() >= min_length
    ).all(), f"All query sequences should be at least {min_length} bp"

    # Verfiy seq_to_short file is empty and
    seq_to_short = pd.read_csv(result.seq_to_short_output, sep="\t")
    for col in required_columns:
        assert col in seq_to_short.columns, f"Column {col} should be in seq_to_short"

    # Verify cts_to_query_cts mapping for debugging
    cts_mapping = tmp_path / "cts_to_query_cts.tsv"
    assert cts_mapping.exists(), "cts to query_cts mapping file should exist"

    mapping_df = pd.read_csv(cts_mapping, sep="\t")
    assert "cts_id" in mapping_df.columns
    assert "query_cts_id" in mapping_df.columns
    assert len(mapping_df) == len(
        seq_table
    ), "Mapping should have same number of entries as sequence table"

    work_dir = tmp_path / "workDir"
    assert work_dir.exists(), "workDir subdirectory should be created"
    assert work_dir.is_dir(), "workDir should be a directory"


def test_prepare_with_short_sequences(tmp_path):
    """
    Test that Prepare correctly handles sequences that are too short to query.
    """
    # Create a test data with one short sequence and one good sequence
    short_seq_table = tmp_path / "short_sequences.tsv"
    short_seq_table.write_text(
        "junc_id\tcts_id\tcts_seq\tpos\tquery_length\n"
        "test1\tshort1\tACGTACGTACGTACGTACGTACGT\t12\t10\n"  # Short sequence
        "test2\tgood1\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t20\t40\n"  # Good length
    )

    prepare = Prepare(
        working_dir=tmp_path / "output", sequence_table=short_seq_table, index_kmer_size=21
    )

    result = prepare.do_prepare()

    # Read the generated files
    seq_table = pd.read_csv(result.sequence_table_output, sep="\t")
    seq_to_short = pd.read_csv(result.seq_to_short_output, sep="\t")

    # test2 should be in sequence_table
    assert len(seq_table) == 1, "Should have exactly one queryable sequence"

    # Sequence in seq_table should be long enough
    min_length = 21 + 4
    for _, row in seq_table.iterrows():
        assert (
            len(row["query_sequence"]) >= min_length
        ), f"Sequence {row['cts_id']} should be at least {min_length} bp"

    # test1 should be in seq_to_short table
    assert len(seq_to_short) == 1, "Should have exactly one short sequence"

    for _, row in seq_to_short.iterrows():
        assert (
            len(row["query_sequence"]) < min_length
        ), f"Short sequence {row['cts_id']} should be less than {min_length} bp"
