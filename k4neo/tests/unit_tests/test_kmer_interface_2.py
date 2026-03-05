from k4neo.index.kmer_index import KmerIndex, KmerMetaIndex
from k4neo.index.index_loader import load_metaindex_from_manifest
import pathlib
import yaml
import pytest
from pydantic import ValidationError


@pytest.fixture
def valid_yaml_manifest_single_index():
    index_manifest = pathlib.Path(__file__).parent.parent / "resources" / "index_manifest.yaml"
    with open(index_manifest, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data


@pytest.fixture
def valid_yaml_manifest_multi_index():
    index_manifest = (
        pathlib.Path(__file__).parent.parent / "resources" / "index_multi_manifest.yaml"
    )
    with open(index_manifest, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data


@pytest.fixture
def partial_valid_yaml_file(tmp_path: pathlib.Path):

    content = {
        # Invalid entry. Neither index or sample mapping exist
        "test_index": {
            "samples": 47,
            "path": "tests/resources/bla.idx",
            "sample_mapping": "tests/resources/index_mapping.txt",
            "method": "toolA",
        },
        # Valid k-mer index
        "test_index2": {
            "samples": 20,
            "path": "k4neo/tests/resources/index/raptor.index",
            "sample_mapping": "k4neo/tests/resources/index/index_mapping.txt",
            "method": "raptor",
        },
        # Invalid index. Method not supported
        "test_index3": {
            "samples": 20,
            "path": "k4neo/tests/resources/index/raptor.index",
            "sample_mapping": "k4neo/tests/resources/index/index_mapping.txt",
            "method": "reindeer",
        },
    }

    file = tmp_path / "partial_valid.yaml"

    with file.open("w") as f:
        yaml.safe_dump(content, f)
    return file


@pytest.fixture
def invalid_yaml_file(tmp_path: pathlib.Path):

    content = {
        "test_index": {
            "samples": 47,
            "path": "tests/resources/bla.idx",
            "sample_mapping": "tests/resources/index_mapping.txt",
            "method": "toolA",
        },
        "test_index2": {
            "samples": 0,
            "path": "k4neo/tests/resources/index/raptor.index",
            "sample_mapping": "k4neo/tests/resources/index/index_mapping.txt",
            "method": "raptor",
        },
        "test_index3": {
            "samples": 20,
            "path": "k4neo/tests/resources/index/raptor.index",
            "sample_mapping": "k4neo/tests/resources/index/index_mapping.txt",
            "method": "reindeer",
        },
    }

    file = tmp_path / "invalid.yaml"

    with file.open("w") as f:
        yaml.safe_dump(content, f)
    return file


def test_KmerIndex_valid_config(mocker):
    with mocker.patch("pathlib.Path.exists", return_value=True):
        km = KmerIndex(
            samples=10,
            path="some.idx",
            sample_mapping="map.txt",
            method="raptor",
        )

        assert km.samples == 10


def test_KmerIndex_invalid_config_method(mocker):
    with mocker.patch("pathlib.Path.exists", return_value=True), pytest.raises(ValidationError):
        KmerIndex(
            samples=10,
            path="some.idx",
            sample_mapping="map.txt",
            method="bla",
        )


def test_KmerIndex_invalid_config_path_not_exists():
    with pytest.raises(ValidationError):
        KmerIndex(
            samples=10,
            path="some.idx",
            sample_mapping="map.txt",
            method="raptor",
        )


def test_KmerMetaIndex_works_with_single_index(valid_yaml_manifest_single_index):

    km = KmerMetaIndex.model_validate(valid_yaml_manifest_single_index)

    # One raptor index
    assert "test_index" in km.root

    assert len(km.root.keys()) == 1
    assert km.total_samples() == 2

    assert km.get_all_index_methods() == set(
        [
            "raptor",
        ]
    )

    # Test get methods and if
    tindex = km.get_single_index("test_index")
    assert isinstance(tindex, KmerIndex)

    assert tindex.path == pathlib.Path("k4neo/tests/resources/index/raptor.index")
    assert tindex.sample_mapping == pathlib.Path("k4neo/tests/resources/index/index_mapping.txt")
    assert tindex.method == "raptor"


def test_KmerMetaIndex_works_with_meta_index(valid_yaml_manifest_multi_index):

    km = KmerMetaIndex.model_validate(valid_yaml_manifest_multi_index)

    # Test multiple indices of different types
    assert km.get_all() == km.root

    assert len(km.root.keys()) == 4
    assert km.total_samples() == 427

    assert km.check_existence("test_index")
    assert km.check_existence("test_index2")
    assert km.check_existence("test_index3")
    assert km.check_existence("test_index4")

    # Test get'er methods and helper methods
    assert km.get_all_index_methods() == set(["raptor", "kmindex"])

    assert km.get_index_method("test_index") == "raptor"
    assert km.get_index_method("test_index4") == "kmindex"

    assert km.get_index_to_method_mapping() == {
        "test_index": "raptor",
        "test_index2": "raptor",
        "test_index3": "raptor",
        "test_index4": "kmindex",
    }

    path_exp = pathlib.Path("k4neo/tests/resources/index/index_mapping.txt")
    assert km.get_index_to_sample_mapping() == {
        "test_index": path_exp,
        "test_index2": path_exp,
        "test_index3": path_exp,
        "test_index4": path_exp,
    }


def test_load_metaindex_from_manifest():

    index_manifest = (
        pathlib.Path(__file__).parent.parent / "resources" / "index_multi_manifest.yaml"
    )
    km = load_metaindex_from_manifest(index_manifest)
    assert len(km.root.keys()) == 4


def test_load_metaindex_from_manifest_works_with_partial_invalid_manifest(partial_valid_yaml_file):
    # One index is correct. Other should be excluded from MetaIndex
    km = load_metaindex_from_manifest(partial_valid_yaml_file)
    assert len(km.root.keys()) == 1
    assert "test_index2" in km.root


def test_load_metaindex_from_manifest_with_invalid_manifest(invalid_yaml_file):
    # Raise ValueError as no index passed validation
    with pytest.raises(ValueError):
        load_metaindex_from_manifest(invalid_yaml_file)
