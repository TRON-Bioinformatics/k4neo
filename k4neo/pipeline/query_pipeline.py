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
    """Generic representation of SnakeMake pipeline"""
    def __init__(self, worklow: pathlib.Path, config: dict, working_dir: pathlib.Path, target_rule: str = "all"):
        """Parameter initialization

        Args:
            worklow (pathlib): The path to the tronmake-kmer-pipeline Snakefile.
            config (dict): Pipeline configuration object.
            working_dir (pathlib.Path): Working directory of pipeline.
            target_rule (str, optional): Target rule to request from workflow. Defaults to "all".
        """
        self.workflow = worklow.resolve()
        self.working_dir = working_dir.resolve()
        self.config = config
        self.target_rule = target_rule
    
    @staticmethod
    def execute_cmd(cmd, working_dir = "."):
        """Run shell command in a subprocess and return the exit-code."""
        logger.info("-> Executing CMD: {}".format(" ".join(cmd)))
        p = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = working_dir, shell=False)
        if p.returncode != 0:
            logger.error(p.stderr)
        return p.returncode

    def run(self, dryrun: bool = False, slurm: bool = True, cores: int = 8) -> int:
        """Run snakemake pipeline using configuration object"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=self.working_dir) as temp_config:
            yaml.dump(self.config, temp_config)
            temp_config.close()
        cmd = ['snakemake',
               '--snakefile', str(self.workflow),
               '--local-cores', str(cores),
               '--jobs', str(cores),
               '--configfile', str(temp_config.name),
               '--directory', str(self.working_dir),
               '--rerun-triggers', 'mtime', '--workflow-profile', 'none']
        if slurm:
            cmd.extend(['--executor', 'slurm'])
        return_code = Pipeline.execute_cmd(cmd)
        return return_code

    def determine_final_query(self):
        """Based on selected methods in index manifest, find query output that qould be created by pipeline"""
        results = {}
        methods = self.config.get('query', dict()).get('methods', [])
        for this_method in methods:
            match this_method:
                case 'raptor':
                    results['raptor'] = self.working_dir / 'query' / 'raptor' / 'search.parquet'
                case'kmindex':
                    results['kmindex'] = self.working_dir / 'query' / 'kmindex' / 'search.parquet'
                case _:
                    raise ValueError(f"-> Tool {this_method} not supported by k4neo query pipeline")
        return results

    def determine_final_index(self):
        """Find index output of user-specified k-mer indexing method"""
        result = ''
        method = self.config.get('indexing', dict()).get('method', '')
        if not method:
            raise ValueError("-> Pipeline config is missing attribute 'method'. Can not determine final index file")
        match method:
            case 'raptor':
                result = self.working_dir / 'index' / 'raptor' / 'raptor.index'
            case'kmindex':
                result = self.working_dir / 'index' / 'kmindex' / 'global_index'
            case _:
                logger.error("-> Tool not supported by k4neo indexing pipeline")
        return result


class QueryPipeline(Pipeline):
    """Query pipeline wrapper"""
    def __init__(self, workflow: str, config: dict, working_dir: pathlib.Path):
        """Parameter initialization"""
        logger.info("-> Initialising k4neo query pipeline.")
        super().__init__(workflow, config, working_dir, target_rule="query")

    def run_pipeline(self, slurm: bool = True, cores: int = 8) -> QueryPipelineResult:
        """Execute query pipeline"""
        final_query = self.determine_final_query()
        logger.info("-> Submitting query pipeline.")
        return_code = self.run(dryrun=False, slurm=slurm, cores=cores)
        if return_code != 0:
            raise K4neoPipelineException("-> Pipeline failed to execute")
        logger.info("-> Finished query pipeline")
        return QueryPipelineResult(query_path=final_query)

    def _test_pipeline(self):
        """Execute dryrun of  pipeline with specified config file for testing purposes"""
        logger.info("-> Starting dry run")
        return_code = self.run(dryrun=True, slurm=False)
        if return_code != 0:
            raise K4neoPipelineException("Pipeline failed to execute")

class IndexPipeline(Pipeline):
    """Indexing pipeline wrapper"""
    def __init__(self, worklow: str, config: dict, working_dir: pathlib.Path):
        """Parameter initialization"""
        super().__init__(worklow, config, working_dir, target_rule="index")

    def run_pipeline(self, slurm: bool = True) -> IndexPipelineResult:
        """Execute indexing pipeline"""
        final_index = self.determine_final_index()
        logger.info("Submitting indexing pipeline")
        return_code = self.run(dryrun=False, slurm=slurm)
        if not return_code:
            raise K4neoPipelineException("Indexing pipeline failed to execute")
        logger.info("Finished indexing pipeline")
        return IndexPipelineResult(index_path=final_index)

    def _test_pipeline(self):
        """Execute dryrun of  pipeline with specified config file for testing purposes"""
        logger.info("Starting dry run..")
        return_code = self.run(dryrun=True, slurm=False)
        if not return_code:
            raise K4neoPipelineException("Pipeline failed to execute")

class PipelineConfig:
    """
    Class to represent SnakeMake pipeline configuration
    """
    def __init__(self,
                 query: bool,
                 indexing: bool,
                 verbose = True):
        """Parameter initialization

        Args:
            query (bool): Pipeline should run in query mode.
            indexing (bool): Pipeline should run in indexing mode.
            verbose (bool, optional): Print configuation details. Defaults to True.
        """
        # Generate config object to be used with
        self.config = {"modus": {
                                 "query": query,
                                 "indexing": indexing}}
        if verbose:
            self.log_configuration()

    def log_configuration(self):
        """Log configuration dictionary on command line"""
        logger.info("Query pipeline configuration")
        for item in self.config:
            for k, v in self.config[item].items():
                logger.info("{}={}".format(k, v))

class QueryPipelineConfig(PipelineConfig):
    """Generate config for query modus"""
    def __init__(self,
                 index: pathlib.Path,
                 kmer_ratio: float,
                 methods: set,
                 verbose = True):
        """Parameter initialization

        Args:
            index (pathlib.Path): Path to yaml based k-mer manifest file.
            kmer_ratio (float): K-mer ratio used by raptor to determine presence/absence in indexed samples.
            methods (set): k-mer methods of indices described in manifest.
            verbose (bool, optional): Print configuation details. Defaults to True.
        """
        super().__init__(query=True, indexing=False, verbose=False)
        # Generate config object to be used with
        self.config["query"] = {
                                "index": str(index),
                                "kmer_ratio": kmer_ratio,
                                "methods": list(methods)}
        if verbose:
            self.log_configuration()

class IndexPipelineConfig(PipelineConfig):
    """Generate config for indexing modus"""
    def __init__(self,
                 samples: pathlib.Path,
                 method: str,
                 kmer_size: float = 21,
                 cutoff: int = 2,
                 fpr: float = 0.05,
                 verbose = True):
        """Parameter initialization

        Args:
            samples (str): Path to TSV files with samples to include in k-mer index.
            method (str): k-mer method used for indexing.
            kmer_size (float, optional): k-mer size of index. Defaults to 21.
            cutoff (int, optional): A value to define solid and weak k-mers. k-mers with lower occurence are omitted from index. Defaults to 2.
            fpr (float, optional): The theoretical false positive rate of bloom filters. Defaults to 0.05.
            verbose (bool, optional): Print configuation details. Defaults to True.
        """
        super().__init__(query=False, indexing=True, verbose=False)
        # Generate config object to be used with
        self.config["indexing"] = {
                                "samples": samples,
                                "method": method,
                                "kmer_size": kmer_size,
                                "cutoff": cutoff,
                                "fpr": fpr}
        if verbose:
            self.log_configuration()
