import os
import sys
import pathlib
from k4neo.pipeline.query_pipeline import QueryPipeline, QueryPipelineConfig
from k4neo.parser.parser import IndexResultParser
from logzero import logger


class KmerIndex(object):
    """
    K-mer (pipeline) index interface
    """
    def __init__(self,
                 pipeline: str,
                 index: str,
                 method: str,
                 raptor_sample_mapping: str = None,
                 kmer_ratio: float = 0.7):

        # Generate config representation that can be passed directly to the
        # snakemake call
        self.pipeline = pipeline
        self.method = method
        self.index = index
        self.raptor_sample_mapping = raptor_sample_mapping
        self.kmer_ratio = kmer_ratio

        self.pipeline_config = QueryPipelineConfig(index=self.index, method=self.method, kmer_ratio=self.kmer_ratio)

    def search_index(self, query_sequences: pathlib.Path, working_dir: pathlib.Path, slurm: bool = True, cores: int = 8):
        """
        Execute k-mer query pipeline on query sequences
        """
        #  To make this class reusable for multiple queries, we extend in this function
        # the config with the search sequence. Here execution specific modifications can be applied
        pipeline_config = self.pipeline_config.config.copy()
        pipeline_config["query"].update({'query_fasta': query_sequences})
        pipeline = QueryPipeline(self.pipeline, pipeline_config, working_dir)
        logger.info("Searching index for context sequences")
        result = pipeline.run_pipeline(slurm=slurm, cores=cores)
        return result

    def result_parser(self, result):
        """
        Parse results returned by k-mer index
        """
        sample_mapping = None
        if self.method == 'raptor':
            sample_mapping = self.raptor_sample_mapping
        parser = IndexResultParser(result.query_path,
                                   self.method,
                                   raptor_sample_mapping=sample_mapping,
                                   kmindex_cutoff=self.kmer_ratio
                                   )
        query_hits = parser.parse_results()
        return query_hits

