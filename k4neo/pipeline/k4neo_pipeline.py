import os
import sys
from k4neo.database.database import DataBase, CreateDataBase
from k4neo.annotator.annotator import Annotator
from k4neo.parser.parser import IndexResultParser
from k4neo.pipeline.query_pipeline import IndexPipeline, IndexPipelineConfig
from k4neo.plotter.plotter import Plotter
from k4neo.setup_logging import setup_logging
from k4neo.helper.helper import DiskIO

class K4neoPipeline:
    """
    Generic class to represent an entire k4neo annotation run.
        1. Preparing the context sequences, 
        2. k-mer pipeline
        3. Parse results
        4. Annotate the results
    """
    def __init__(self,
                 pipeline,
                 workflow_profile, 
                 index_manifest,
                 kmer_ratio,
                 cores,
                 slurm: bool = False,
                 working_dir,
                 queries,
                 metadata_db,

    ):
        self.annotator = Annotator(working_dir, queries, metadata_db)
        self.kmer_index = KmerIndex(pipeline=pipeline, 
                                    workflow_profile=workflow_profile,
                                    index_manifest=index_manifest,
                                    kmer_ratio=kmer_ratio)
        self.cores = cores
        self.working_dir = working_dir
        self.slurm = slurm

    def prepare(self):
        pass

    def search_index():
        
        query_pipeline_results = self.kmer_index.search_index(
            self.queries, self.working_dir, slurm=self.slurm, cores=self.cores
        )
        return query_pipeline_results

    def parse_results(query_pipeline_results):
        
        parser_compatible_structure = []
        for this_result in query_pipeline_results.query_path:
            parser_compatible_structure.append(
                (
                    this_result[0],
                    this_result[1],
                    self.kmer_index.index_to_sample_mapping[this_result[1]],
                    this_result[2],
                )
            )
        parser = IndexResultParser2(query_pipeline_results=parser_compatible_structure, cores=self.cores)
        query_hits = parser.parse_results(kmer_ratio=self.kmer_ratio)
        return query_hits
    
    def annotate_results():
        pass
    