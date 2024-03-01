import pathlib
from logzero import logger
from snakemake import snakemake
from k4neo.exceptions import K4neoPipelineException
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
    def __init__(self, worklow: str, config: dict, working_dir: pathlib.Path, target_rule: str = "all"):
        self.workflow = pathlib.Path(worklow).resolve()
        self.working_dir = working_dir
        self.config = config
        self.target_rule = target_rule

    def run(self, dryrun: bool = False, slurm: bool = True, cores: int = 8) -> bool:
        """
        Run snakemake pipeline using the selected executor
        """
        return_code = snakemake(
                                self.workflow,
                                workdir=self.working_dir,
                                config=self.config,
                                dryrun=dryrun,
                                slurm=slurm,
                                rerun_triggers="mtime",
                                cores=cores)
        return return_code

    def determine_final_query(self):
        """
        Find query output of selected query method
        """
        result = ''
        method = self.config.get('query', dict()).get('method', '')
        if not method:
            raise ValueError("Pipeline config is missing attribute 'method'. Can not determine final output file")
        match method:
            case 'cobs':
                result = self.working_dir / 'query' / 'cobs' / 'cobs_search.txt'
            case 'raptor':
                result = self.working_dir / 'query' / 'raptor' / 'raptor_search.txt'
            case 'reindeer':
                result = self.working_dir / 'query' / 'reindeer' / 'reindeer_search.txt'
            case'kmindex':
                result = self.working_dir / 'query' / 'kmindex' / 'kmindex_search.txt'
            case _:
                logger.error("Tool not supported by k4neo pipeline")
        return result

    def determine_final_index(self):
        """
        Find index output of selected indexing method
        """
        result = ''
        method = self.config.get('indexing', dict()).get('method', '')
        if not method:
            raise ValueError("Pipeline config is missing attribute 'method'. Can not determine final index file")
        match method:
            case 'raptor':
                result = self.working_dir / 'index' / 'raptor' / 'raptor.index'
            case'kmindex':
                result = self.working_dir / 'index' / 'kmindex' / 'global_index'
            case _:
                logger.error("Tool not supported by k4neo indexing pipeline")
        return result


class QueryPipeline(Pipeline):
    def __init__(self, workflow: str, config: dict, working_dir: pathlib.Path):
        logger.info("Initialising k4neo query pipeline...")
        super().__init__(workflow, config, working_dir, target_rule="query")

    def run_pipeline(self, slurm: bool = True, cores: int = 8) -> QueryPipelineResult:
        """
        Execute query pipeline
        """
        final_query = self.determine_final_query()
        logger.info("Submitting query pipeline...")
        return_code = self.run(dryrun=False, slurm=slurm, cores=cores)
        if not return_code:
            raise K4neoPipelineException("Pipeline failed to execute")
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
            raise K4neoPipelineException("Pipeline failed to execute")

class IndexPipeline(Pipeline):
    """
    Wrapper method to call k4neo k-mer indexing pipeline
    """
    def __init__(self, worklow: str, config: dict, working_dir: pathlib.Path):
        super().__init__(worklow, config, working_dir, target_rule="index")

    def run_pipeline(self, slurm: bool = True) -> IndexPipelineResult:
        """
        Execute indexing pipeline
        """
        final_index = self.determine_final_index()
        logger.info("Submitting indexing pipeline")
        return_code = self.run(dryrun=False, slurm=slurm)
        if not return_code:
            raise K4neoPipelineException("Indexing pipeline failed to execute")
        logger.info("Finished indexing pipeline")
        return IndexPipelineResult(index_path=final_index)

    def _test_pipeline(self):
        """
        execute dryrun of  pipeline with specified config file for testing purposes
        :return:
        """
        logger.info("Starting dry run..")
        return_code = self.run(dryrun=True, slurm=False)
        if not return_code:
            raise K4neoPipelineException("Pipeline failed to execute")

class PipelineConfig:
    """
    Generate generic object to hold pipeline config
    Snakemake accepts a dictionary representing the content of the traditional yaml config file
    """
    def __init__(self,
                 query: bool,
                 indexing: bool,
                 verbose = True):
        # Generate config object to be used with
        self.config = {"modus": {
                                 "query": query,
                                 "indexing": indexing}}
        if verbose:
            self.log_configuration()

    def log_configuration(self):
        logger.info("Query pipeline configuration")
        for item in self.config:
            for k, v in self.config[item].items():
                logger.info("{}={}".format(k, v))

class QueryPipelineConfig(PipelineConfig):
    """
    Generate config for query pipeline
    """
    def __init__(self,
                 index: str,
                 method: str,
                 kmer_ratio: float,
                 verbose = True):
        super().__init__(query=True, indexing=False)
        # Generate config object to be used with
        self.config = {"query": {
                                "index": index,
                                "method": method,
                                "kmer_ratio": kmer_ratio}}
        if verbose:
            self.log_configuration()

class IndexPipelineConfig(PipelineConfig):
    """
    Generate config for indexing pipeline
    """
    def __init__(self,
                 samples: str,
                 method: str,
                 kmer_size: float = 21,
                 cutoff: int = 2,
                 fpr: float = 0.05,
                 verbose = True):
        super().__init__(query=False, indexing=True)
        # Generate config object to be used with
        self.config = {"indexing": {
                                "samples": samples,
                                "method": method,
                                "kmer_size": kmer_size,
                                "cutoff": cutoff,
                                "fpr": fpr}}
        if verbose:
            self.log_configuration()
