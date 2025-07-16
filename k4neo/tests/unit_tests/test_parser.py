import pytest
import re
from six import PY2
from collections import ChainMap
from k4neo.parser.parser import Parser, IndexResultParser
from k4neo.parser.index_parser import IndexResultParser2, BinaryKmerIndexResultParser

# Tests for Parser class

@pytest.fixture
def mocker_tissue_map(mocker):
    data_lines = [
            'tissue_description_found_in_public_data\ttissue\tsubtissue',
            'subcutaneous adipose tissue\tadipose tissue\tsubcutaneous',
            'mesenteric fat pad\tadipose tissue\tmesenteric fat pad',
            'adrenal tissue\tadrenal gland\tNA']
    mocked_data = mocker.mock_open(read_data='\n'.join(data_lines))
    builtin_open = "__builtin__.open" if PY2 else "builtins.open"
    mocker.patch(builtin_open, mocked_data)

def test_read_tissue_map(mocker_tissue_map):
    tissue_map = Parser._read_tissue_map(tissue_map='fakefile')
    assert isinstance(tissue_map, list)
    assert len(tissue_map) == 3
    assert list(tissue_map[0].keys()) == ['tissue_public', 'tissue', 'subtissue']
    assert list(tissue_map[0].values()) == ['subcutaneous adipose tissue', 'adipose tissue', 'subcutaneous']


@pytest.fixture
def mocker_sample_file(mocker):
    data_lines = [
         'sample_name\truns\ttissue\tdevelopmental_stage\tdisease',
          'SAMEA962332\tERR030872,ERR030903\tthyroid\tadult\tcancer',
          'SAMEA962347\tERR030876,ERR030899\tskeletal muscle\tadult\thealthy',
          'IT_N_101_sRL\t\tSAMN21521836\tTestis\tadult\thealthy']

    mocked_data = mocker.mock_open(read_data='\n'.join(data_lines))
    builtin_open = "__builtin__.open" if PY2 else "builtins.open"
    mocker.patch(builtin_open, mocked_data)

def test_read_sample_file(mocker_sample_file):
    s1_expected = {'sample_name': 'SAMEA962332', 'runs': 'ERR030872,ERR030903', 'tissue': 'thyroid', 'developmental_stage': 'adult', 'disease': 'cancer'}
    samples = Parser._read_sample_file(data_table='fakefile')
    assert isinstance(samples, list)
    assert len(samples) == 3
    assert samples[0] == s1_expected


# Test for IndexResultParser class.

@pytest.fixture
def example_table_rows():
   """
   Pytest fixture representing a batch of rows returns by parrow.iter_batches
   """
   rows_raptor = [
       {"__index_level_0__": "cts_1", "s1": 0.7, "s2": 0.7, "s3": 0.7},
       {"__index_level_0__": "cts_2", "s1": None, "s2": None, "s3": 0.7},
       {"__index_level_0__": "cts_3", "s1": None, "s2": None, "s3": None},]
   rows_kmindex = [
       {"__index_level_0__": "cts_1", "s1": 0.9123, "s2": 0.895, "s3": 0.9997},
       {"__index_level_0__": "cts_2", "s1": 0.000, "s2": 0.734, "s3": 0.111},
       {"__index_level_0__": "cts_3", "s1": 0.000, "s2": 0.000, "s3": 0.000},]
   return rows_raptor, rows_kmindex


def test_parse_table_rows_raptor(example_table_rows):
    parsed_raptor_results = list(map(lambda p: IndexResultParser._parse_table_row(table_row=p, method="raptor", kmer_ratio=0.7), example_table_rows[0]))
    assert len(parsed_raptor_results) == 3
    # Merge returned dicts together for easy testing of cases
    query_results = ChainMap(*parsed_raptor_results)
    # Test all samples contain the sequence
    assert query_results['cts_1'] == set(["s1", "s2", "s3"])
    # Test only a single sample detected
    assert query_results['cts_2'] == set(["s3",])
    # Test detected in None of the samples
    assert query_results['cts_3'] == set([None,])


def test_parse_table_rows_kmindex(example_table_rows):
    parsed_kmindex_results = list(map(lambda p: IndexResultParser._parse_table_row(table_row=p, method="kmindex", kmer_ratio=0.7), example_table_rows[1]))
    assert len(parsed_kmindex_results) == 3
    # Merge returned dicts together for easy testing of cases
    query_results = ChainMap(*parsed_kmindex_results)
    # Test all samples contain the sequence
    assert query_results['cts_1'] == set(["s1", "s2", "s3"])
    # Test only a single sample detected
    assert query_results['cts_2'] == set(["s2",])
    # Test detected in None of the samples
    assert query_results['cts_3'] == set([None,])

# test for IndexResultParser2 and BinaryKmerIndexResultParser
