import os
import sys
import pathlib
import yaml
from k4neo.pipeline.query_pipeline import QueryPipeline, QueryPipelineConfig
from k4neo.parser.parser import IndexResultParser
from logzero import logger


class KmerIndex(object):
    """
    K-mer (pipeline) index interface
    """
    def __init__(self,
                 pipeline: str,
                 index_manifest: str,
                 kmer_ratio: float = 0.7):

        # Generate config representation that can be passed directly to the
        # snakemake call
        self.pipeline = pipeline
        self.index_manifest = index_manifest,
        self.kmer_ratio = kmer_ratio

        self.pipeline_config = QueryPipelineConfig(index=self.index_manifest, kmer_ratio=self.kmer_ratio)
        self.index_struct = self.read_index_struct(self.index_manifest)

    def read_index_struct(file):
        """
        Read meta index of different raptor indices
        """
        with open(file, 'r') as file_handle:
            index_struct = yaml.safe_load(file_handle)
        # ToDo add sanity checks if meta index file is correctly formatted
        return index_struct
    
    def get_index_methods(self):
        """
        Results of indices of the same method are aggregated into one query file.
        """
        index_methods = set()
        for index_id, index_properties in self.index_struct.items():
            method = index_properties.get('method', None)
            if method is None:
                logger.warning(f"k-mer method not specified for index: {index_id}. Pipeline will probably fail")
            index_methods.add(method)
        return index_methods

    def search_index(self, query_sequences: pathlib.Path, working_dir: pathlib.Path, slurm: bool = True, cores: int = 8):
        """
        Execute k-mer query pipeline on query sequences
        """
        #  To make this class reusable for multiple queries, we extend in this function
        # the config with the search sequence. Here execution specific modifications can be applied
        pipeline_config = self.pipeline_config.config.copy()
        pipeline_config["query"].update({'query_fasta': query_sequences})
        pipeline = QueryPipeline(self.pipeline, pipeline_config, working_dir, self.methods)
        logger.info("-> Searching index for context sequences")
        result = pipeline.run_pipeline(slurm=slurm, cores=cores)
        return result

    def result_parser(self, result):
        """
        Parse results returned by k-mer index
        """
        parser = IndexResultParser(result.query_path,
                                   self.method,
                                   kmindex_cutoff=self.kmer_ratio
                                   )
        query_hits = parser.parse_results()
        return query_hits

