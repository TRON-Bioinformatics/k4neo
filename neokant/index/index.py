import os
import sys
import pathlib
from neokant.pipeline.query_pipeline import QueryPipeline, PipelineConfig
from neokant.parser.parser import IndexResultParser
from logzero import logger


class KmerIndex(object):
    def __init__(self,
                 pipeline: str,
                 index: str,
                 method: str,
                 reindeer_sample_mapping: str = None,
                 raptor_sample_mapping: str = None,
                 kmindex_cutoff: float = 0.7,
                 kmer_ratio: float = 0.45):

        # Generate config representation that can be passed directly to the
        # snakemake call
        self.pipeline = pipeline
        self.method = method
        self.index = index
        self.reindeer_sample_mapping = reindeer_sample_mapping
        self.raptor_sample_mapping = raptor_sample_mapping
        self.kmindex_cutoff = kmindex_cutoff

        self.pipeline_config = PipelineConfig(index=index,
                                              method=method,
                                              kmer_ratio=kmer_ratio)

    def search_index(self, query_sequences, working_dir, slurm=True):
        # To make this class reusable for multiple queries, we extend in this function
        # the config with the search sequence. Here execution specific modifications can be applied
        pipeline_config = self.pipeline_config.config.copy()
        pipeline_config["query"].update({'query_fasta': query_sequences})
        pipeline = QueryPipeline(self.pipeline, pipeline_config, working_dir)
        logger.info("Searching index for context sequences")
        result = pipeline.run_pipeline(slurm=slurm)
        return result

    def result_parser(self, result):
        sample_mapping = None
        if self.method == 'reindeer':
            sample_mapping = self.reindeer_sample_mapping
        if self.method == 'raptor':
            sample_mapping = self.raptor_sample_mapping
        parser = IndexResultParser(result.query_path,
                                   self.method,
                                   reindeer_sample_mapping=sample_mapping,
                                   raptor_sample_mapping=sample_mapping,
                                   kmindex_cutoff=self.kmindex_cutoff
                                   )
        query_hits = parser.parse_results()
        return query_hits

