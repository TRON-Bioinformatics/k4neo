import os
import pathlib
import sys
from logzero import logger
from snakemake import snakemake
import time
from neokant.exceptions import NeoKantPipelineException


@dataclass
class QueryPipelineResult:
    query_path: str


class Pipeline():
    def __init__(self, worklow: str, config: dict, target_rule: str = "all"):
        self.workflow = pathlib.Path(worklow).resolve()
        self.config = config
        self.target_rule = target_rule

    def run(self) -> bool:
        return True

    def determine_final_query(self):
        method = self.config['method']
        if method == 'cobs':
            return ''
        elif method == 'raptor':
            return ''
        elif method == 'reindeer':
            return ''
        elif method == 'kmindex':
            return ''
        else:
            return


class QueryPipeline(Pipeline):
    def __init__(self, workflow: str, config: dict):
        logger.info("Initialising neokant query pipeline...")
        super().__init(workflow, config, target_rule="query")


    def run_pipeline(self) -> QueryPipelineResult:
        final_query = self.determine_final_query()
        logger.info("Submitting query pipeline...")
        return_code = self.run()
        if not return_code:
            raise NeoKantPipelineException("Pipeline failed to execute")
        logger.info("Finished query pipeline")
        return QueryPipelineResult(query_path=final_query)



class IndexPipeline(Pipeline):
    pass