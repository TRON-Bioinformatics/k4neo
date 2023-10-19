import os
import pathlib
import sys
from logzero import logger
from snakemake import snakemake
import time
from neokant.exceptions import NeoKantPipelineException
from dataclasses import dataclass


@dataclass
class QueryPipelineResult:
    """
    Data class to hold results from query pipeline
    """
    query_path: str

@dataclass
class IndexPipelineResult:
    """
    Data class to hold results from indexing pipeline
    """
    index_path: str


class Pipeline:
    def __init__(self, worklow: str, config: dict, working_dir: str, target_rule: str = "all"):
        self.workflow = pathlib.Path(worklow).resolve()
        self.working_dir = working_dir
        self.config = config
        self.target_rule = target_rule

    def run(self, dryrun: bool = False, slurm: bool = True) -> bool:
        return_code = snakemake(
                                         self.workflow,
                                         workdir=self.working_dir,
                                         config=self.config,
                                         dryrun=dryrun,
                                         slurm=slurm
                                         )
        return return_code

    def determine_final_query(self):
        result = ''
        method = self.config.get('query', dict()).get('method', '')
        if not method:
            raise ValueError("Missing attribute method... Can not determine final output file")
        match method:
            case 'cobs':
                result = os.path.join(self.working_dir, 'query/cobs/cobs_search.txt')
            case 'raptor':
                result = os.path.join(self.working_dir, 'query/raptor/raptor_search.txt')
            case 'reindeer':
                result = os.path.join(self.working_dir, 'query/reindeer/reindeer_search.txt')
            case'kmindex':
                result = os.path.join(self.working_dir, 'query/kmindex/kmindex_search.txt')
            case _:
                logger.error("Tool not supported by neoKant pipeline")
        return result



class QueryPipeline(Pipeline):
    def __init__(self, workflow: str, config: dict, working_dir: str):
        logger.info("Initialising neokant query pipeline...")
        super().__init__(workflow, config, working_dir, target_rule="query")

    def run_pipeline(self, slurm: bool = True) -> QueryPipelineResult:
        final_query = self.determine_final_query()
        logger.info("Submitting query pipeline...")
        return_code = self.run(dryrun=False, slurm=slurm)
        if not return_code:
            raise NeoKantPipelineException("Pipeline failed to execute")
        logger.info("Finished query pipeline")
        return QueryPipelineResult(query_path=final_query)

    def _test_pipeline(self):
        """
        execute dryrun of  pipeline with specified config file for testing purposes
        :return:
        """
        logger.info("Starting dry run..")
        return_code = self.run(dryrun=True, slurm=False)
        if not return_code:
            raise NeoKantPipelineException("Pipeline failed to execute")

class PipelineConfig:
    """
    Generate generic object to hold pipeline config
    Snakemake accepts a dictionary representing the content of the yaml config file
    """
    def __init__(self,
                 index: str,
                 method: str,
                 kmer_ratio: float,
                 verbose = False):
        # Generate config object to be used with
        self.config = {"query": {
                                "index": index,
                                "method": method,
                                "kmer_ratio": kmer_ratio}}
        if verbose:
            self.log_configuration

    def log_configuration(self):
        logger.info("Query pipeline configuration")
        for k, v in self.config.items():
            logger.info("{}={}".format(k, v))