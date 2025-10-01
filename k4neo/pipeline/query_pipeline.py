import pathlib
import tempfile
import yaml
from dataclasses import dataclass
from k4neo.exceptions import K4neoPipelineException
from k4neo.helper.helper import ShellExec
from loguru import logger


@dataclass
class QueryPipelineResult:
    """
    Data class to hold results from query pipeline
    """

    query_path: list[
        tuple[str, str, pathlib.Path]
    ]  # List of tuples with method and path to query result


@dataclass
class IndexPipelineResult:
    """
    Data class to hold results from indexing pipeline
    """

    index_path: str


class Pipeline:
    """
    Generic representation of a snakemake pipeline. The pipeline class
    holds the path to workflow (snakefile), a configuration object, the working directory
    and the requested target rule. Upon exeuction the config object is written into a
    yaml file and the pipeline executed. The class `pipeline` is intended to be a reusable
    interface for specific pipelines.
    """

    def __init__(
        self,
        worklow: pathlib.Path,
        workflow_profile: pathlib.Path,
        config: dict,
        working_dir: pathlib.Path,
        target_rule: str = "all",
    ):
        """Parameter initialization

        Args:
            worklow (pathlib): The path to the tronmake-kmer-pipeline Snakefile.
            config (dict): Pipeline configuration object.
            working_dir (pathlib.Path): Working directory of pipeline.
            target_rule (str, optional): Target rule to request from workflow. Defaults to "all".
        """
        self.workflow = worklow.resolve()
        self.workflow_profile = workflow_profile.resolve()
        self.working_dir = working_dir.resolve()
        self.config = config
        self.target_rule = target_rule

    def run(self, dryrun: bool = False, slurm: bool = True, cores: int = 8) -> int:
        """Build and execute pipeline

        The run method implements the execution of the pipeline. It takes care
        of writing the config object into yaml file and constructing
        the snakemake shell command for execution in a subprocess.

        Args:
            dryrun (bool, optional): Execute pipeline without performing any operation. Defaults to False.
            slurm (bool, optional): Execute pipeline with slurm support. Defaults to True.
            cores (int, optional): Number of cores. Defaults to 8.

        Returns:
            int: _description_
        """
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, dir=self.working_dir
        ) as temp_config:
            yaml.dump(self.config, temp_config)
            temp_config.close()
        cmd = [
            "snakemake",
            "--snakefile",
            str(self.workflow),
            "--local-cores",
            str(cores),
            "--jobs",
            str(cores),
            "--configfile",
            str(temp_config.name),
            "--directory",
            str(self.working_dir),
            "--rerun-triggers",
            "mtime",
            "--workflow-profile",
            "none",
            "--profile",
            str(self.workflow_profile.parent),
            "--until",
            self.target_rule,
        ]
        if slurm:
            cmd.extend(["--executor", "slurm"])
        if dryrun:
            cmd.extend(["--dry-run"])
        return_code = ShellExec.execute_cmd(cmd)
        return return_code


class QueryPipeline(Pipeline):
    """
    The QueryPipeline class implements the search mode of the TronMake k-mer pipeline.
    """

    def __init__(
        self,
        workflow: pathlib.Path,
        workflow_profile: pathlib.Path,
        config: dict,
        working_dir: pathlib.Path,
        target_rule: str
    ):
        """Parameter initialization"""
        logger.debug("Initialising of k4neo query pipeline.")
        super().__init__(workflow, workflow_profile, config, working_dir, target_rule=target_rule)

    def determine_final_query(self):
        """Based on selected methods in index manifest, find query output that could be created by pipeline"""
        results = []
        for this_index_name, this_method in self.config.get(
            "index_to_method_mapping", dict
        ).items():
            match this_method:
                case "raptor":
                    # query/raptor/{subindex}/search.tsv"
                    results.append(
                        (
                            this_method,
                            this_index_name,
                            self.working_dir / "query" / "raptor" / this_index_name / "search.tsv",
                        )
                    )
                case "kmindex":
                    # "query/kmindex/{subindex}/search.tsv",
                    results.append(
                        (
                            this_method,
                            this_index_name,
                            self.working_dir / "query" / "kmindex" / this_index_name / "search.tsv",
                        )
                    )
                case "jellyfish":
                    # query/jellyfish/{subindex}/search.tsv"
                    results.append(
                        (
                            this_method,
                            this_index_name,
                            self.working_dir / "query" / "jellyfish" / this_index_name / "quantitative_search.tsv",
                        )
                    )
                case _:
                    raise ValueError(f"-> Tool {this_method} not supported by k4neo query pipeline")
        return results

    def run_pipeline(self, slurm: bool = True, cores: int = 8) -> QueryPipelineResult:
        """Execute query pipeline"""
        final_query = self.determine_final_query()
        logger.info("Submitting query pipeline.")
        return_code = self.run(dryrun=False, slurm=slurm, cores=cores)
        if return_code != 0:
            logger.error("Pipeline command returned non zero exit code.")
            raise K4neoPipelineException("Pipeline failed to execute")
        logger.info("Finished query pipeline.")
        return QueryPipelineResult(query_path=final_query)

    def _test_pipeline(self):
        """Execute dryrun of  pipeline with specified config file for testing purposes"""
        logger.info("Starting dry run")
        return_code = self.run(dryrun=True, slurm=False)
        if return_code != 0:
            raise K4neoPipelineException("Pipeline failed to execute")


class IndexPipeline(Pipeline):
    """
    The IndexPipeline class implements the indexing of the TronMake k-mer pipeline.
    """

    def __init__(
        self,
        worklow: pathlib.Path,
        workflow_profile: pathlib.Path,
        config: dict,
        working_dir: pathlib.Path,
    ):
        """Parameter initialization"""
        super().__init__(worklow, workflow_profile, config, working_dir, target_rule="index")

    def determine_final_index(self):
        """Find index output of user-specified k-mer indexing method"""
        result = ""
        method = self.config.get("indexing", dict()).get("method", "")
        if not method:
            raise ValueError(
                "-> Pipeline config is missing attribute 'method'. Can not determine final index file"
            )
        match method:
            case "raptor":
                result = self.working_dir / "index" / "raptor" / "raptor.index"
            case "kmindex":
                result = self.working_dir / "index" / "kmindex" / "global_index"
            case _:
                logger.error("-> Tool not supported by k4neo indexing pipeline")
        return result

    def run_pipeline(self, slurm: bool = True) -> IndexPipelineResult:
        """Execute indexing pipeline"""
        final_index = self.determine_final_index()
        logger.info("-> Submitting indexing pipeline")
        return_code = self.run(dryrun=False, slurm=slurm)
        if not return_code:
            raise K4neoPipelineException("-> Indexing pipeline failed to execute")
        logger.info("-> Finished indexing pipeline")
        return IndexPipelineResult(index_path=final_index)

    def _test_pipeline(self):
        """Execute dryrun of  pipeline with specified config file for testing purposes"""
        logger.info("Starting dry run..")
        return_code = self.run(dryrun=True, slurm=False)
        if not return_code:
            raise K4neoPipelineException("Pipeline failed to execute")


class PipelineConfig:
    """
    Class to represent a generic pipeline configuration.
    """

    def __init__(self, verbose=True):
        """Parameter initialization

        Args:
            verbose (bool, optional): Print configuation details. Defaults to True.
        """

        if verbose:
            self.log_configuration()

    def log_configuration(self):
        """Log configuration dictionary on command line"""
        logger.info("Pipeline configuration:")
        for item in self.config:
            for k, v in self.config[item].items():
                logger.info("  {}={}".format(k, v))


class KmerPipelineConfig(PipelineConfig):
    """
    Class to represent a base TronMake k-mer pipeline configuration.
    """

    def __init__(self, query: bool, indexing: bool, verbose=True):
        """Parameter initialization

        Args:
            query (bool): Pipeline should run in query mode.
            indexing (bool): Pipeline should run in indexing mode.
            verbose (bool, optional): Print configuation details. Defaults to True.
        """
        super().__init__(verbose=verbose)
        self.config = {"modus": {"query": query, "indexing": indexing}}


class QueryPipelineConfig(KmerPipelineConfig):
    """Generate config for query modus"""

    def __init__(
        self, index: pathlib.Path, kmer_ratio: float, index_to_method_mapping: dict, verbose=True
    ):
        """Parameter initialization

        Args:
            index (pathlib.Path): Path to yaml based k-mer manifest file.
            kmer_ratio (float): K-mer ratio used by raptor to determine presence/absence in indexed samples.
            methods (set): k-mer methods of indices described in manifest.
            verbose (bool, optional): Print configuation details. Defaults to True.
        """
        super().__init__(query=True, indexing=False, verbose=False)
        # SnakeMake pipeline options - specific for TronMake k-mer pipeline
        self.config["query"] = {
            "index": str(index),
            "kmer_ratio": kmer_ratio,
        }
        # Not for pipeline. This is used by the QueryPipeline class to determine the files
        # and associated methods.
        self.config["index_to_method_mapping"] = index_to_method_mapping
        if verbose:
            self.log_configuration()


class IndexPipelineConfig(KmerPipelineConfig):
    """Generate config for indexing modus"""

    def __init__(
        self,
        samples: pathlib.Path,
        method: str,
        kmer_size: float = 21,
        cutoff: int = 2,
        fpr: float = 0.05,
        verbose=True,
    ):
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
            "fpr": fpr,
        }
        if verbose:
            self.log_configuration()
