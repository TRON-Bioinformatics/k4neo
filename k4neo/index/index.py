import pathlib
import yaml
from k4neo.pipeline.query_pipeline import QueryPipeline, QueryPipelineConfig
from k4neo.parser.parser import IndexResultParser
from loguru import logger


class KmerIndex(object):
    """
    K-mer (pipeline) index interface
    """

    def __init__(
            self, pipeline: pathlib.Path, workflow_profile: pathlib.Path, index_manifest: pathlib.Path, kmer_ratio: float = 0.7
    ):

        # Generate config representation that can be passed directly to the
        # snakemake call
        self.pipeline = pipeline
        self.workflow_profile = workflow_profile
        self.index_manifest = index_manifest
        self.kmer_ratio = kmer_ratio

        self.index_struct = self.read_index_struct()
        self.index_methods = self._get_index_methods()
        self.index_names = self.index_struct.keys()
        logger.info(
            f"-> Executing queries against {' & '.join(self.index_methods)} k-mer indices"
        )
        self.pipeline_config = QueryPipelineConfig(
            index=self.index_manifest,
            kmer_ratio=self.kmer_ratio,
            methods=self.index_methods,
            index_names=self.index_names,
        )

    def read_index_struct(self):
        """
        Read k4neo meta index definition file.
        """
        with open(self.index_manifest, "r") as file_handle:
            index_struct = yaml.safe_load(file_handle)
        # ToDo add sanity checks if meta index file is correctly formatted
        return index_struct

    def _get_index_methods(self):
        """Collect k-mer index methods

        Results of sub-indices of the same method are aggregated into one query file.

        """
        index_methods = set()
        for index_id, index_properties in self.index_struct.items():
            method = index_properties.get("method", None)
            if method is None:
                logger.warning(
                    f"-> k-mer method is missing for index: {index_id}. Ignoring"
                )
                continue
            index_methods.add(method)
        return index_methods

    def search_index(
        self,
        query_sequences: pathlib.Path,
        working_dir: pathlib.Path,
        slurm: bool = True,
        cores: int = 8,
    ):
        """Execute k-mer query pipeline

        Method calls the QueryPipeline to execute a search of sequences in fasta format
        against sub-indices provided in the index manifest

        Args:
            query_sequences: File path of sequences in FASTA format
            working_dir: Path of pipeline workdir
            slurm: If True, QueryPipeline will submit jobs to slurm scheduler.
            cores: Number of threads

        Returns:
            result: A dictionary mapping search sequences to detected samples

        """
        #  To make this class reusable for multiple queries, we extend in this function
        # the config with the search sequence. Here execution specific modifications can be applied
        pipeline_config = self.pipeline_config.config.copy()
        pipeline_config["query"].update({"query_fasta": str(query_sequences)})
        pipeline = QueryPipeline(self.pipeline, self.workflow_profile, pipeline_config, working_dir)
        logger.info("-> Searching index for context sequences")
        result = pipeline.run_pipeline(slurm=slurm, cores=cores)
        return result

    def result_parser(self, query_pipeline_results, cores=8):
        """
        Parse results returned by k-mer index
        """
        parser = IndexResultParser(
            query_pipeline_results=query_pipeline_results, cores=cores
        )
        query_hits = parser.parse_results(kmer_ratio=self.kmer_ratio)
        return query_hits
