import pathlib
import tempfile
import subprocess
import yaml
from logzero import logger
from k4neo.exceptions import K4neoPipelineException
from dataclasses import dataclass


@dataclass
class QueryPipelineResult:
    """
    Data class to hold results from query pipeline
    """
    query_path: dict

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
    
    @staticmethod
    def execute_cmd(cmd, working_dir = "."):
        """This function runs a command into a subprocess."""
        logger.info("-> Executing CMD: {}".format(" ".join(cmd)))
        p = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = working_dir, shell=False)
        if p.returncode != 0:
            logger.error(p.stderr)
        return p.returncode

    def run(self, dryrun: bool = False, slurm: bool = True, cores: int = 8) -> int:
        """
        Run snakemake pipeline using the selected executor to orchestrate search arcoss subindices
        """
        with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=self.workdir) as temp_config:
            yaml.dump(self.config, temp_config)
            temp_config.close()
        cmd = ['snakemake',
               '--snakefile', self.workflow,
               '--local-cores', str(self.cores),
               '--jobs', str(self.cores),
               '--configfile', str(temp_config.name),
               '--use-conda',
               '--directory', str(self.self.working_dir),
               '--rerun-triggers', 'mtime']
        if slurm:
            cmd.extend(['--executor', 'slurm'])
        return_code = Pipeline.execute_cmd(cmd)
        return return_code

    def determine_final_query(self):
        """
        Based on selected methods in index manifest, find query output that qould be created by pipeline
        """
        results = {}
        methods = self.config.get('query', dict()).get('methods', [])
        for this_method in methods:
            match this_method:
                case 'raptor':
                    results['raptor'] = self.working_dir / 'query' / 'raptor' / 'search.parquet'
                case'kmindex':
                    results['kmindex'] = self.working_dir / 'query' / 'kmindex' / 'search.parquet'
                case _:
                    raise ValueError(f"Tool {this_method} not supported by k4neo query pipeline")
        return results

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
        if return_code != 0:
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
        if return_code != 0:
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
                 kmer_ratio: float,
                 methods: set,
                 verbose = True):
        super().__init__(query=True, indexing=False)
        # Generate config object to be used with
        self.config = {"query": {
                                "index": index,
                                "kmer_ratio": kmer_ratio,
                                "methods": list(methods)}}
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
