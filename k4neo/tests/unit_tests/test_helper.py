from k4neo.helper.helper import InputValidation, SequenceOperation
import pandas as pd


class TestInputValidation:
    def test_columns_missing_all_present(self):
        df = pd.DataFrame(columns=["name", "age", "email"])
        expected = (False, [])
        result = InputValidation.columns_missing(df, ["name", "age", "email"])
        assert result == expected

    def test_columns_missing_some_missing(self):
        df = pd.DataFrame(columns=["name", "age"])
        expected = (True, ["email"])
        result = InputValidation.columns_missing(df, ["name", "age", "email"])
        assert result == expected

    def test_columns_missing_all_missing(self):
        df = pd.DataFrame(columns=["id", "status"])
        expected = (True, ["name", "email"])
        result = InputValidation.columns_missing(df, ["name", "email"])
        assert result == expected

    def test_columns_missing_empty_check_list(self):
        df = pd.DataFrame(columns=["name", "email"])
        expected = (False, [])
        result = InputValidation.columns_missing(df, [])
        assert result == expected

    def test_columns_missing_empty_dataframe(self):
        df = pd.DataFrame()
        expected = (True, ["name", "email"])
        result = InputValidation.columns_missing(df, ["name", "email"])
        assert result == expected


class TestSequenceOperation:

    SEQ = "ATGGTAGATA"  # Length = 10,

    def test_standard_middle_slice(self):
        result = SequenceOperation.subset_cts(self.SEQ, pos=5, length=5)
        # Index at 5 ("A"), expected indices to be rturned 3–8 → "GTAGA"
        assert result == "GTAGA"

    def test_near_start(self):
        result = SequenceOperation.subset_cts(self.SEQ, pos=1, length=4)
        # start = max(0, 1 - 2) = 0, stop = 0 + 4 = 4 → "ATGG"
        assert result == "ATGG"

    def test_near_end(self):
        result = SequenceOperation.subset_cts(self.SEQ, pos=9, length=4)
        # start = 9 - 2 = 7, stop = min(7+4, 10) = 10 → "ATA"
        assert result == "ATA"

    def test_full_sequence(self):
        result = SequenceOperation.subset_cts(self.SEQ, pos=5, length=20)
        #
        assert result == self.SEQ

    def test_none_pos_or_length_returns_original(self):
        assert SequenceOperation.subset_cts(self.SEQ, pos=None, length=5) == self.SEQ
        assert SequenceOperation.subset_cts(self.SEQ, pos=5, length=None) == self.SEQ
        assert SequenceOperation.subset_cts(self.SEQ, pos=pd.NA, length=5) == self.SEQ
        assert SequenceOperation.subset_cts(self.SEQ, pos=5, length=pd.NA) == self.SEQ

    def test_empty_sequence(self):
        assert SequenceOperation.subset_cts("", pos=3, length=5) == ""

    def test_start_stop_exact_edges(self):
        # pos = 1, length = 3 → start = 0, stop = 3 → "ATG"
        result = SequenceOperation.subset_cts(self.SEQ, pos=1, length=3)
        assert result == "ATG"

        # pos = 9, length = 3 → start = 8, stop = 10 → "ATA"
        result = SequenceOperation.subset_cts(self.SEQ, pos=9, length=3)
        assert result == "ATA"

    def test_kmer(self):
        result = SequenceOperation.get_kmers(self.SEQ, 5)
        #ATGGTAGATA
        assert result == ["ATGGT", "TGGTA", "GGTAG", "GTAGA", "TAGAT", "AGATA"]
    
    def test_canonicalize(self):
        # Test rev complemenent is returned
        result = SequenceOperation.canonicalize("ATGGT")
        assert result == "ACCAT"
        # Test k-mer is returned
        result = SequenceOperation.canonicalize("AGATA")
        assert result == "AGATA"
