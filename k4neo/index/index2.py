"""k4neo k-mer index

This module aims to replace the old KmerIndex interface with a more modular
alternative. The reason for reimplementation is the very unflexible structure 
of the old class. In k4neo we provide a metaindex represented in a single object
instance. However, this does not allow to use and query indices independently.
Moreover the old method was to dependent on the query pipeline attributes. 
With the new structure the pipeline and the parser could be replaced with another 
implementation yielding the same return types.


Todo:
    Better error and value handling
"""

import pathlib
import yaml
from k4neo.pipeline.query_pipeline import QueryPipeline, QueryPipelineConfig
from k4neo.pipeline import TARGET_RULES_OF_METHODS, K4NEO_SUPPORTED_TOOLS
from k4neo.parser.parser import IndexResultParser
from k4neo.parser.index_parser import IndexResultParser2
from loguru import logger
from typing import Literal, Dict, List

from pydantic import  BaseModel, RootModel, Field, ValidationError, field_validator
import tempfile


class KmerIndex(BaseModel):
    """Representation of a single k-mer index

    Attributes:
        samples (int): Number of samples in k-mer index.
            Must be greater than 0.
        path (pathlib.Path): Path to the directory/file containing the k-mer index.
        sample_mapping (pathlib.Path): Path to the file mapping sample to their index id.
            Only required when using Raptor.
        method: (str): K-mer index method.
            One of raptor, kmindex, jellyfish

    Raises:
        ValueError: If `samples` id less than or equal to 0 or
            if `sample_mapping` and `path` do not exist.

    """
    samples: int = Field(gt=0)
    path: pathlib.Path
    sample_mapping: pathlib.Path
    method: Literal[*K4NEO_SUPPORTED_TOOLS]

    @field_validator("path", "sample_mapping")
    def file_must_exist(cls, path: pathlib.Path):
        if not path.exists():
            raise ValueError(f"File does not exist: {path}")
        return path


class KmerMetaIndex(RootModel):
    """Representation of a k-mer meta index

    `KmerMetaIndex` represents a k-mer meta index composed of
    multiple `KmerIndex` instances. Allows the combintation of
    multiple indices and index types e.g. raptor and kmindex

    Attributes:
        root (Dict[str, KmerIndex]): A collection of k-mer indices of type `KmerIndex`
            that form the metaindex.
    """
    root: Dict[str, KmerIndex]

    def get(self, name: str) -> KmerIndex:
        """Get KmerIndex from KmerMetaIndex

        Args:
            name (str): Identifier of KmerIndex

        Returns:
            KmerIndex: KmerIndex instance of `name`
        """
        return self.root[name]

    def get_all(self) -> Dict[str, KmerIndex]:
        """Get root field of KmerMetaIndex

        Returns:
            Dict[str, KmerIndex]: Root object
        """
        return self.root
    
    def check_existence(self, name: str) -> bool:
        """Check if KmerIndex is part of KmerMetaIndex

        Args:
            name (str): Identifier of KmerIndex

        Returns:
            bool: A boolean indicating presence/absence
        """
        return name in self.root.keys()
    
    def get_index_method(self, name:str) -> str:
        """_summary_

        Args:
            name (str): _description_

        Returns:
            str: _description_
        """
        index = self.root[name]
        return index.method

    def get_all_index_methods(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        return set(
            index.method for index in self.root.values() 
        )
    
    def get_index_to_method_mapping(self) -> Dict[str, pathlib.Path]:
        """_summary_

        Returns:
            Dict[str, pathlib.Path]: _description_
        """
        return { key: index.method for key, index in self.root.items() }
    
    def get_index_to_sample_mapping(self) -> Dict[str, pathlib.Path]:
        """_summary_

        Returns:
            Dict[str, pathlib.Path]: _description_
        """
        return { key: index.sample_mapping for key, index in self.root.items() }

    def get_kmer_depth_mapping(self) -> Dict[str, int]:
        """Not implemented

        Returns:
            Dict[str, int]: _description_
        """
        pass
    
    def total_samples(self) -> int:
        return sum(i.samples for i in self.root.values())

class KmerIndexProcessor:
    def __init__(self, meta_index: KmerMetaIndex, pipeline: pathlib.Path, workflow_profile: pathlib.Path, quantitative: bool = False):
        self.meta_index = meta_index
        self.pipeline = pipeline
        self.workflow_profile = workflow_profile
        self.quantitative = quantitative

        self.index_to_method_mapping = self.meta_index.get_index_to_method_mapping()
        
        if not quantitative:
            # This is only required for raptor
            self.index_to_sample_mapping = self.meta_index.get_index_to_sample_mapping()

        if quantitative:
            self.kmer_depth_mapping = self.meta_index.get_kmer_depth_mapping()

    def _get_pipeline_target_rules(self) -> str:
        """
        """
        target_rules = set()
        methods = set(self.index_to_method_mapping.values())
        for this_method in methods:
            target_rules.add(TARGET_RULES_OF_METHODS[this_method])
        return " ".join(target_rules)
    
    def call_pipeline(self, pipeline_config: QueryPipelineConfig, target_rules: str, working_dir, slurm, cores):
        
        logger.info("-> Searching index for context sequences")
        pipeline = QueryPipeline(self.pipeline,
                                 self.workflow_profile,
                                 pipeline_config,
                                 working_dir, target_rule=target_rules)
        result = pipeline.run_pipeline(slurm=slurm, cores=cores)
        return result


    def search_index(
        self,
        subindex_id: str,
        query_sequences: pathlib.Path,
        working_dir: pathlib.Path,
        slurm: bool = True,
        kmer_ratio: float = 0.7,
        cores: int = 8
    ):
        """Execute k-mer query pipeline

        Method calls the QueryPipeline to execute a search of sequences in fasta format
        against sub-indices provided in the index manifest

        Args:
            query_sequences: File path of sequences in FASTA format
            working_dir: Path of pipeline workdir
            slurm: If True, QueryPipeline will submit jobs to slurm scheduler.
            cores: Number of threads

        Returns:
            result: A dictionary mapping search sequences to detected samples

        """
        with tempfile.NamedTemporaryFile(mode="w", delete=True, delete_on_close=False, dir=self.working_dir ) as temp_manifest:
            yaml.dump(self.meta_index, temp_manifest)
            temp_manifest.close()
            
            # Configure pipeline for search
            pipeline_config = QueryPipelineConfig(
                index=self.temp_manifest,
                kmer_ratio=kmer_ratio,
                index_to_method_mapping=self.index_to_method_mapping,
            )
            pipeline_config["query"].update({"query_fasta": str(query_sequences)})
            result = self.call_pipeline(pipeline_config, self._get_pipeline_target_rules(), working_dir, slurm, cores)
        
        return result
    












#    def result_parser2(self, query_pipeline_results, cores=8):
#        """
#        Parse results returned by k-mer index
#        """
#        # query_pipeline_results [("method", "subindex_name", "result_path")]
#        parser_compatible_structure = []
#        for this_result in query_pipeline_results.query_path:
#            parser_compatible_structure.append(
#                (
#                    this_result[0],
#                    this_result[1],
#                    self.index_to_sample_mapping[this_result[1]],
#                    this_result[2],
#                )
#            )
#        parser = IndexResultParser2(query_pipeline_results=parser_compatible_structure, cores=cores)
#        query_hits = parser.parse_result2(kmer_ratio=self.kmer_ratio)
#        return query_hits
