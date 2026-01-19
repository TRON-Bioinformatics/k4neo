"""k4neo k-mer index

This module aims to replace the old KmerIndex interface with a more modular
alternative. The reason for reimplementation is the very unflexible structure
of the old class. In k4neo we provide a meta k-mer index represented in a single object
instance. However, this does not allow to use and query indices independently.
Moreover the old method was to dependent on the query pipeline attributes.
With the new structure the pipeline and the parser could be replaced with another
implementation yielding the same return types.



"""

import pathlib
from k4neo.pipeline import  K4NEO_SUPPORTED_TOOLS
from typing import Literal, Dict, List
from pydantic import BaseModel, RootModel, Field, ValidationError, field_validator, field_serializer


class KmerIndex(BaseModel):
    """Representation of a single k-mer index

    Attributes:
        samples (int | None): Number of samples in k-mer index.
            Must be greater than 0.

        path (pathlib.Path): Path to the directory/file containing the k-mer index.

        sample_mapping (pathlib.Path | None): Path to the file mapping sample to their index id.
            Only required when using Raptor.

        method: (str): K-mer index method.
            One of raptor, kmindex, jellyfish

        kmer_depth (int): K-mer depth of indexed samples.
            Only required for quantitative k-mer indices. Defaults to None.

    Raises:
        ValueError: If `samples` id less than or equal to 0 or
            if `sample_mapping` and `path` do not exist.

    """

    samples: int | None = Field(default=None, gt=0)
    path: pathlib.Path
    sample_mapping: pathlib.Path | None = None
    method: Literal[*K4NEO_SUPPORTED_TOOLS]
    kmer_depth: int | None = None

    @field_validator("path", "sample_mapping")
    def file_must_exist(cls, path: pathlib.Path | None, info):
        # sample_mapping is optional and not required for jellyfish.
        if path is None:
            return path
        if not path.exists():
            raise ValueError(f"File does not exist: {path}")
        return path

    @field_serializer("path", "sample_mapping", mode="plain")
    def serialize_path(self, value: pathlib.Path) -> str:
        return str(value)


class KmerMetaIndex(RootModel):
    """Representation of a k-mer meta index

    `KmerMetaIndex` represents a collection of `KmerIndex` instances that
    can be queried together or independently. It allows to manage multiple
    k-mer indices of different types in a single object as a metaindex.

    Attributes:
        root (Dict[str, KmerIndex]): A collection of `KmerIndex` intances
            that form the metaindex.
    """

    root: Dict[str, KmerIndex]

    def get_single_index(self, name: str) -> KmerIndex:
        """Get single KmerIndex instance from KmerMetaIndex

        Args:
            name (str): Identifier of KmerIndex

        Returns:
            KmerIndex: KmerIndex instance of `name`
        """
        return self.root[name]

    def get_all(self) -> Dict[str, KmerIndex]:
        """Get all KmerIndex instances in KmerMetaIndex

        Returns:
            Dict[str, KmerIndex]: All KmerIndex instances in KmerMetaIndex
        """
        return self.root

    def check_existence(self, name: str) -> bool:
        """Check if KmerIndex exists in KmerMetaIndex

        Args:
            name (str): Identifier of KmerIndex

        Returns:
            bool: True if KmerIndex exists in KmerMetaIndex
        """
        return name in self.root.keys()

    def get_index_method(self, name: str) -> str:
        """Get the index method for a specific KmerIndex

        Args:
            name (str): Identifier of KmerIndex

        Returns:
            str: Index method of the specified KmerIndex
        """
        index = self.root[name]
        return index.method

    def get_all_index_methods(self):
        """Get all unique index methods in the meta index

        Returns:
            set: A set of unique index methods
        """
        return set(index.method for index in self.root.values())

    def get_index_to_method_mapping(self) -> Dict[str, str]:
        """Index to method mapping for all indices in meta index

        Returns:
            Dict[str, str]: Mapping of index name to method
        """
        return {key: index.method for key, index in self.root.items()}

    def get_index_to_sample_mapping(self) -> Dict[str, pathlib.Path]:
        """Index to sample mapping for all indices in meta index

        Returns:
            Dict[str, pathlib.Path]: Mapping of index name to sample mapping file
        """
        return {key: index.sample_mapping for key, index in self.root.items()}

    def get_kmer_depth_mapping(self) -> Dict[str, int]:
        """Obtain k-mer depth mapping for all indices in meta index

        Returns:
            Dict[str, int]: Mapping of index name to k-mer depth
        """
        return {
            key: index.kmer_depth
            for key, index in self.root.items()
            if index.kmer_depth is not None
        }

    def total_samples(self) -> int:
        """Number of indexed samples

        Returns:
            int:
        """
        return sum(i.samples for i in self.root.values())
