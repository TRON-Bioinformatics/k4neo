import pytest
import re
from six import PY2
from collections import ChainMap
from k4neo.parser.parser import Parser
from k4neo.parser.index_parser import (
    IndexResultParser2,
    BinaryKmerIndexResultParser,
    QuantitativeKmerIndexParser,
)

# Tests for Parser class


@pytest.fixture
def mocker_tissue_map(mocker):
    data_lines = [
        "tissue_description_found_in_public_data\ttissue\tsubtissue",
        "subcutaneous adipose tissue\tadipose tissue\tsubcutaneous",
        "mesenteric fat pad\tadipose tissue\tmesenteric fat pad",
        "adrenal tissue\tadrenal gland\tNA",
    ]
    mocked_data = mocker.mock_open(read_data="\n".join(data_lines))
    builtin_open = "__builtin__.open" if PY2 else "builtins.open"
    mocker.patch(builtin_open, mocked_data)


def test_read_tissue_map(mocker_tissue_map):
    tissue_map = Parser._read_tissue_map(tissue_map="fakefile")
    assert isinstance(tissue_map, list)
    assert len(tissue_map) == 3
    assert list(tissue_map[0].keys()) == ["tissue_public", "tissue", "subtissue"]
    assert list(tissue_map[0].values()) == [
        "subcutaneous adipose tissue",
        "adipose tissue",
        "subcutaneous",
    ]


@pytest.fixture
def mocker_sample_file(mocker):
    data_lines = [
        "sample_name\truns\ttissue\tdevelopmental_stage\tdisease",
        "SAMEA962332\tERR030872,ERR030903\tthyroid\tadult\tcancer",
        "SAMEA962347\tERR030876,ERR030899\tskeletal muscle\tadult\thealthy",
        "IT_N_101_sRL\t\tSAMN21521836\tTestis\tadult\thealthy",
    ]

    mocked_data = mocker.mock_open(read_data="\n".join(data_lines))
    builtin_open = "__builtin__.open" if PY2 else "builtins.open"
    mocker.patch(builtin_open, mocked_data)


def test_read_sample_file(mocker_sample_file):
    s1_expected = {
        "sample_name": "SAMEA962332",
        "runs": "ERR030872,ERR030903",
        "tissue": "thyroid",
        "developmental_stage": "adult",
        "disease": "cancer",
    }
    samples = Parser._read_sample_file(data_table="fakefile")
    assert isinstance(samples, list)
    assert len(samples) == 3
    assert samples[0] == s1_expected


# test for IndexResultParser2 and BinaryKmerIndexResultParser


@pytest.mark.parametrize(
    "kmer_ratio,expected",
    [
        (
            0.1,
            {
                "d2e9c05e86abc70c804fc48a7f25aad8": {"MCF7", "SKBR3"},
                "b410c5db97c66dc7a4d200ec44ce4a1d": {"MCF7", "SKBR3"},
            },
        ),
        (
            0.25,
            {
                "d2e9c05e86abc70c804fc48a7f25aad8": {"MCF7", "SKBR3"},
                "b410c5db97c66dc7a4d200ec44ce4a1d": {"MCF7"},
            },
        ),
        (
            0.5,
            {
                "d2e9c05e86abc70c804fc48a7f25aad8": {"MCF7"},
                "b410c5db97c66dc7a4d200ec44ce4a1d": {"MCF7"},
            },
        ),
        (
            0.8,
            {
                "d2e9c05e86abc70c804fc48a7f25aad8": {"MCF7"},
                "b410c5db97c66dc7a4d200ec44ce4a1d": {None},
            },
        ),
        (
            1.0,
            {
                "d2e9c05e86abc70c804fc48a7f25aad8": {"MCF7"},
                "b410c5db97c66dc7a4d200ec44ce4a1d": {None},
            },
        ),
    ],
)
def test_parse_kmindex(mocker, kmer_ratio, expected):
    mock_data = [
        {"samples": "samples_2:d2e9c05e86abc70c804fc48a7f25aad8", "MCF7": "1", "SKBR3": "0.25"},
        {"samples": "samples_2:b410c5db97c66dc7a4d200ec44ce4a1d", "MCF7": "0.7", "SKBR3": "0.1"},
    ]

    # mock csv.DictReader
    mocker.patch("k4neo.parser.index_parser.DictReader", return_value=mock_data)
    mocker.patch("builtins.open", mocker.mock_open(read_data="ignored"))

    parser = BinaryKmerIndexResultParser("dummy.tsv", method="kmindex", kmer_ratio=kmer_ratio)

    results = parser._parse_kmindex()

    assert results == expected


@pytest.fixture
def example_raptor_output():
    """
    1 = sample to minimiser mapping
    2 = Raptor example output
    3 = Expected output
    """
    return (
        "\n".join(
            [
                "minimiser_id\tsample_name",
                "/minimiser/0115_128_CE_TMA_1TR1_SL1_S8_L005_R1_001.minimiser\tMCF7",
                "/minimiser/0115_043_CE_TMA_1TR4_SL2_S5_L004_R1_001.minimiser\tSKBR3",
            ]
        ),
        "\n".join(
            [
                "## Index is compressed = true",
                "## Index is HIBF = true",
                "#0\t/minimiser/0115_128_CE_TMA_1TR1_SL1_S8_L005_R1_001.minimiser",
                "#1\t/minimiser/0115_043_CE_TMA_1TR4_SL2_S5_L004_R1_001.minimiser",
                "#QUERY_NAME\tUSER_BINS",
                "d2e9c05e86abc70c804fc48a7f25aad8\t0",
                "b410c5db97c66dc7a4d200ec44ce4a1d",
                "61c9762eae1072735d7ff13834da8942",
                "9221002eafd64f81d677be49d5b52885\t0,1",
            ]
        ),
        {
            "d2e9c05e86abc70c804fc48a7f25aad8": {"MCF7"},
            "b410c5db97c66dc7a4d200ec44ce4a1d": {None},
            "61c9762eae1072735d7ff13834da8942": {None},
            "9221002eafd64f81d677be49d5b52885": {"MCF7", "SKBR3"},
        },
    )


def test_parse_raptor(mocker, example_raptor_output):
    from io import StringIO

    # Mock minimiser file
    mock_mapping_file = StringIO(example_raptor_output[0])
    # Mock open for raptor resultsw
    mock_search_file = StringIO(example_raptor_output[1])
    mocker.patch("builtins.open", side_effect=[mock_mapping_file, mock_search_file])

    parser = BinaryKmerIndexResultParser("dummy.tsv", "raptor", "dummy2.txt", kmer_ratio=0.7)
    results = parser._parse_raptor()

    assert results == example_raptor_output[2]


# Parsing of subindex results is happening in parallel.
# Here we test that method results are correctly put together for further processing.


def test_add_real_sample_to_placeholder_set():
    """
    Test that adding samples for default datastructure works
    """
    target = {None}
    new = {"sample1", "sample2"}
    IndexResultParser2.update_sample_set(target, new)
    assert target == {"sample1", "sample2"}


def test_add_placeholder_to_placeholder_set():
    """
    Test that placeholder stays if not detected in any sample
    """
    target = {None}
    new = {None}
    IndexResultParser2.update_sample_set(target, new)
    assert target == {None}


def test_add_placeholder_to_real_sample():
    """
    Test that placeholder of next subindex (if cts was not detected in any sample) result is ignored
    """
    target = {"sample1"}
    new = {None}
    IndexResultParser2.update_sample_set(target, new)
    assert target == {"sample1"}


def test_add_mix_of_placeholder_and_real_sample_to_real_sample():
    """
    Test that placeholder is also removed when mixed with real sample name
    """
    target = {"sample1"}
    new = {"sample2", None}
    IndexResultParser2.update_sample_set(target, new)
    assert target == {"sample1", "sample2"}


def test_parse_jellyfish_basic(mocker):
    # Simulate a jellyfish output file with two cts_ids and multiple kmer counts
    file_content = (
        "cts1\tAAA\t5\n" "cts1\tAAC\t7\n" "cts2\tAAG\t2\n" "cts2\tAAT\t3\n" "cts2\tACC\t4\n"
    )
    mock_open = mocker.patch("builtins.open", mocker.mock_open(read_data=file_content))
    parser = QuantitativeKmerIndexParser("dummy.txt", "jellyfish")
    result = parser.parse_jellyfish()
    assert result == {
        "cts1": [5, 7],
        "cts2": [2, 3, 4],
    }
    mock_open.assert_called_once_with("dummy.txt", "r")


def test_parse_jellyfish_empty_file(mocker):
    # Simulate an empty jellyfish output file
    mock_open = mocker.patch("builtins.open", mocker.mock_open(read_data=""))
    parser = QuantitativeKmerIndexParser("dummy.txt", "jellyfish")
    result = parser.parse_jellyfish()
    assert result == {}
    mock_open.assert_called_once_with("dummy.txt", "r")


def test_parse_jellyfish_single_entry(mocker):
    # Simulate a file with a single entry
    file_content = "ctsX\tAAA\t42\n"
    mock_open = mocker.patch("builtins.open", mocker.mock_open(read_data=file_content))
    parser = QuantitativeKmerIndexParser("dummy.txt", "jellyfish")
    result = parser.parse_jellyfish()
    assert result == {"ctsX": [42]}
    mock_open.assert_called_once_with("dummy.txt", "r")
