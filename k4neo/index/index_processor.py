"""k4neo k-mer index processor module

This module provides the KmerIndexProcessor class, which facilitates querying
a k-mer meta index using a specified pipeline and workflow profile. It supports
both binary and quantitative k-mer searches. The class acts as an interface
between the meta index and the query pipeline, coordinating the different
components. It serves as a reference implementation for users who want to
build their own k-mer index querying workflows using k4neo components.

"""

import pathlib
import yaml
import tempfile

from k4neo.pipeline.query_pipeline import QueryPipeline, QueryPipelineConfig
from k4neo.parser.index_parser import IndexResultParser2
from k4neo.pipeline import TARGET_RULES_OF_METHODS
from k4neo.index.kmer_index import KmerMetaIndex

from loguru import logger


class KmerIndexProcessor:
    """K-mer index processor

    This class manages the querying of a k-mer meta index using a specified
    pipeline and workflow profile. It supports both binary and quantitative
    k-mer searches. It acts as an interface between the meta index and the
    query pipeline and acts like service class to coordinate the different
    components. It acts as a reference implementation for users who want to
    build their own k-mer index querying workflows using k4neo components.

    """

    def __init__(
        self,
        meta_index: KmerMetaIndex,
        pipeline: pathlib.Path,
        workflow_profile: pathlib.Path,
        quantitative: bool = False,
        pipeline_cls=QueryPipeline,
        pipeline_config_cls=QueryPipelineConfig,
    ):
        """Initialize KmerIndexProcessor

        Args:
            meta_index (KmerMetaIndex): KmerMetaIndex instance
            pipeline (pathlib.Path): Path to the pipeline
            workflow_profile (pathlib.Path): Path to the workflow profile
            quantitative (bool, optional): Whether the search is quantitative. Defaults to False.
            pipeline_cls (QueryPipeline, optional): Pipeline class. Defaults to QueryPipeline.
            pipeline_config_cls (QueryPipelineConfig, optional): Pipeline config class. Defaults to QueryPipelineConfig.
        """
        self.meta_index = meta_index
        self.pipeline = pipeline
        self.workflow_profile = workflow_profile
        self.quantitative = quantitative
        self.pipeline_cls = pipeline_cls
        self.pipeline_config_cls = pipeline_config_cls

        self.index_to_method_mapping = self.meta_index.get_index_to_method_mapping()

        if not quantitative:
            # This is only required for raptor
            self.index_to_sample_mapping = self.meta_index.get_index_to_sample_mapping()

        if quantitative:
            self.kmer_depth_mapping = self.meta_index.get_kmer_depth_mapping()

            # assert (
            #    self.kmer_depth_mapping
            # ), "At least one index requires kmer_depth for quantitative search"

    def _get_pipeline_target_rules(self) -> str:
        """Obtain target rules for pipeline

        Returns:
            str: Target rules passed to snakemake pipeline (--until)
        """
        methods = set(self.index_to_method_mapping.values())
        return " ".join(TARGET_RULES_OF_METHODS[this_method] for this_method in methods)

    def _serialize_manifest(self, working_dir: pathlib.Path) -> str:
        """Create a temporary manifest file for the meta index."""
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, delete_on_close=False, dir=working_dir
        ) as temp_manifest:
            yaml.dump(self.meta_index.model_dump(), temp_manifest)
            temp_manifest.close()
            return temp_manifest.name

    def _build_pipeline_config(
        self, manifest_path: str, query_sequences: pathlib.Path, kmer_ratio: float
    ) -> QueryPipelineConfig:
        """Build pipeline configuration."""
        config = self.pipeline_config_cls(
            index=manifest_path,
            kmer_ratio=kmer_ratio,
            index_to_method_mapping=self.index_to_method_mapping,
        )
        config.config["query"].update({"query_fasta": str(query_sequences)})
        return config

    def search_index(
        self,
        query_sequences: pathlib.Path,
        working_dir: pathlib.Path,
        slurm: bool = True,
        kmer_ratio: float = 0.7,
        cores: int = 8,
    ) -> dict:
        """Execute k-mer query pipeline

        Method calls the QueryPipeline to execute a search of sequences in fasta format
        against sub-index instances in meta_index

        Args:
            query_sequences: File path of sequences in FASTA format
            working_dir: Path of pipeline workdir
            slurm: If True, QueryPipeline will submit jobs to slurm scheduler.
            cores: Number of threads

        Returns:
            result: A dictionary mapping search sequences to detected samples

        """
        manifest_path = self._serialize_manifest(working_dir)
        pipeline_config = self._build_pipeline_config(manifest_path, query_sequences, kmer_ratio)
        target_rule = self._get_pipeline_target_rules()
        pipeline = self.pipeline_cls(
            self.pipeline,
            self.workflow_profile,
            pipeline_config,
            working_dir,
            target_rule=target_rule,
        )
        logger.info("-> Searching index for context sequences")
        return pipeline.run_pipeline(slurm=slurm, cores=cores)

    def result_parser2(self, query_pipeline_results, cores, kmer_ratio):
        """
        Parse results returned by k-mer index
        """
        # query_pipeline_results [("method", "subindex_name", "result_path")]
        if self.quantitative:
            raise NotImplementedError(
                "Quantitative result parsing not implemented yet. Please use QuantitativeKmerIndexParser directly."
            )
        parser_compatible_structure = []
        for this_result in query_pipeline_results.query_path:
            parser_compatible_structure.append(
                (
                    this_result[0],
                    this_result[1],
                    self.index_to_sample_mapping[this_result[1]],
                    this_result[2],
                )
            )
        parser = IndexResultParser2(query_pipeline_results=parser_compatible_structure, cores=cores)
        query_hits = parser.parse_result2(kmer_ratio=kmer_ratio)
        return query_hits
