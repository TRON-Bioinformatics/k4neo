import yaml
import pathlib
from k4neo.index.kmer_index import KmerIndex, KmerMetaIndex
from loguru import logger
import pydantic_core


def load_metaindex_from_manifest(manifest_path: pathlib.Path) -> KmerMetaIndex:
    """Load k-mer metaindex from yaml

    Load k4neo k-mer based metaindex yaml and create instance of `KmerMetaIndex`.
    Individual indices are validated and only indices passing pydantic validation
    (see `KmerIndex`) are registered into k4neo metaindex.

    Args:
        manifest_path (pathlib.Path): A path to k4neo metaindex yaml file

    Raises:
        ValueError: If none of the indices passes validation and nothing is left for query

    Returns:
        KmerMetaIndex: Instance of k4neo metaindex
    """
    keep_indices = {}
    with open(manifest_path, "r", encoding="utf-8") as file_handle:
        manifest_content = yaml.safe_load(file_handle)

    for this_index, this_index_properties in manifest_content.items():
        try:
            KmerIndex(**this_index_properties)
        except pydantic_core._pydantic_core.ValidationError as e:
            logger.warning(f"Excluding {this_index} as it failed validation")
            continue
        keep_indices[this_index] = this_index_properties

    if not keep_indices:
        raise ValueError(f"{manifest_path} contains no valid entries.")

    return KmerMetaIndex.model_validate(keep_indices)
