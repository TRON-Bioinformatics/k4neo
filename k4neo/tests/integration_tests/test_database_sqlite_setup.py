import pandas as pd
from k4neo.database_sqlite.database import CreateDataBase


def test_full_setup_and_precomputations(tmp_path):
    # Prepare mock files
    tissue_map = tmp_path / "tissue.tsv"
    study_file = tmp_path / "study.txt"
    tissue_map.write_text(
        (
            "tissue_description_found_in_public_data\ttissue\tsubtissue\tOntology_UBERON\tHigher_grouping\tstudies\n"
            "subcutaneous adipose tissue\tadipose tissue\tsubcutaneous\tUBERON:0001013\tConnective & Soft tissue\tENCODE-ENTEX\n"
        )
    )

    study_file.write_text("STUDY1\t{}\t1\n".format(tmp_path / "samples.tsv"))
    samples_file = tmp_path / "samples.tsv"
    samples_file.write_text(
        "sample_name\truns\ttissue\tdevelopmental_stage\tdisease\tsex\n"
        "S1\tR1\tsubcutaneous adipose tissue\tadult\thealthy\tM\n"
    )

    with CreateDataBase(
        db_file=None,
        data_set_file=study_file,
        tissue_map=tissue_map,
        test=True,
    ) as db:
        db.setup_db()
        db.precomputations()

        # Check tissue_map
        tissue_df = pd.read_sql("SELECT * FROM tissue_map", db.connection)
        assert tissue_df.shape[0] == 1

        # Check samples
        samples_df = pd.read_sql("SELECT * FROM samples", db.connection)
        assert samples_df.shape[0] == 1

        # Check precomputed tissue_counts
        tc_df = pd.read_sql("SELECT * FROM tissue_counts", db.connection)
        assert tc_df.shape[0] >= 1

        # Check aggregated_tissue_counts
        agg_df = pd.read_sql("SELECT * FROM aggregated_tissue_counts", db.connection)
        assert agg_df.shape[0] >= 1
