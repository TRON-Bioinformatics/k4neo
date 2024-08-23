import pytest
import re
from six import PY2
from collections import ChainMap
from k4neo.database.database import CreateDataBase
from k4neo.parser.parser import Parser
import pandas as pd



@pytest.fixture
def example_database(mocker):
    fake_data = [
        {'tissue_public' : 'subcutaneous adipose tissue', 'tissue':'adipose tissue', 'subtissue':'subcutaneous'},
        {'tissue_public' : 'mesenteric fat pad', 'tissue':'adipose tissue', 'subtissue':'mesenteric fat pad'},
        {'tissue_public' : 'adrenal tissue', 'tissue':'adrenal gland', 'subtissue':'NA'},
    ]
    database = CreateDataBase(db_file=None, data_set_file='', tissue_map='', test=True)
    mocker.patch("k4neo.parser.parser.Parser.parse_tissuemap_into_document", return_value=fake_data)
    database._add_tissues()
    return database

@pytest.fixture
def mocker_study_file(mocker):
    data_lines = [
          'E-MTAB-2836\thealthy_tissue_studies/E-MTAB-2836/EMTAB2836_kmer_annotation.tsv\t122',
          'E-MTAB-513\thealthy_tissue_studies/E-MTAB-513/EMTAB513_kmer_annotation.tsv\t16',
          'E-MTAB-6814\thealthy_tissue_studies/E-MTAB-6814/EMTAB6814_kmer_annotation.tsv\t313']

    mocked_data = mocker.mock_open(read_data='\n'.join(data_lines))
    builtin_open = "__builtin__.open" if PY2 else "builtins.open"
    mocker.patch(builtin_open, mocked_data)


def test_parse_study_table(mocker_study_file):
    database = CreateDataBase(db_file=None, data_set_file='fakefile', tissue_map='', test=True)
    assert len(list(database._parse_study_table())) == 3

def test_tissue_table_initialisation(mocker):
    fake_data = [
        {'tissue_public' : 'subcutaneous adipose tissue', 'tissue':'adipose tissue', 'subtissue':'subcutaneous'},
        {'tissue_public' : 'mesenteric fat pad', 'tissue':'adipose tissue', 'subtissue':'mesenteric fat pad'},
        {'tissue_public' : 'adrenal tissue', 'tissue':'adrenal gland', 'subtissue':'NA'},
    ]
    database = CreateDataBase(db_file=None, data_set_file='', tissue_map='', test=True)
    mocker.patch("k4neo.parser.parser.Parser.parse_tissuemap_into_document", return_value=fake_data)
    database._add_tissues()
    inserted_elements = pd.DataFrame(database.database.table('tissue_map'))
    assert inserted_elements.shape[0] == 3
    assert inserted_elements.get("tissue", pd.Series()).unique().shape[0] == 2
    assert inserted_elements.get("tissue", pd.Series()).unique().tolist() == ['adipose tissue', 'adrenal gland']

def test_sample_study_table_initialisation(mocker):
    pass

def test_sample_table_initialisation(mocker):
    pass

def test_update_sample_document_with_tissue(example_database):
    db = example_database
    fake_data = [
        {"sample_name": "S1", "tissue": "subcutaneous adipose tissue"},
        {"sample_name": "S2", "tissue": "Adrenal_Gland"},
    ]
    expected_output = [
        ({"sample_name": "S1", "tissue": "adipose tissue", "subtissue": "subcutaneous"}, True),
        ({"sample_name": "S2", "tissue": "Adrenal_Gland"}, False)
    ]
    updated_data = list(map(db._update_sample_document_with_tissue, fake_data))
    assert updated_data == expected_output