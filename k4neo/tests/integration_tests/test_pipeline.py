from k4neo.pipeline.query_pipeline import QueryPipeline, QueryPipelineConfig, QueryPipelineResult
from k4neo.index.index import KmerIndex
import pathlib


class TestPipeline:
    index_manifest = pathlib.Path(__file__).parent.parent / "resources" / "index_manifest.yaml"
    kmer_ratio = 0.7
    pipeline = (
        pathlib.Path(__file__).parent.parent.parent
        / "pipeline"
        / "tronmake-kmer-pipeline"
        / "workflow"
        / "Snakefile"
    )
    profile = pathlib.Path(__file__).parent.parent.parent / "pipeline" / "config.yaml"
    fasta = pathlib.Path(__file__).parent.parent / "resources" / "example.fasta"

    def test_raptor_query_pipeline(self):
        qconfig = QueryPipelineConfig(
            index=self.index_manifest,
            kmer_ratio=self.kmer_ratio,
            index_to_method_mapping={"test_index": "raptor"},
        )
        qconfig.config["query"].update({"query_fasta": str(self.fasta)})
        pipeline = QueryPipeline(
            self.pipeline,
            self.profile,
            qconfig,
            pathlib.Path("./"),
            target_rule="query_raptor",
        )
        result = pipeline.run_pipeline(slurm=False, cores=1)
        assert isinstance(result, QueryPipelineResult)
        assert isinstance(result.query_path, list)
        assert result.query_path == [
            (
                "raptor",
                "test_index",
                pathlib.Path("./").resolve() / "query/raptor/test_index/search.tsv",
            )
        ]

    def test_jellyfish_query_pipeline(self):
        quant_manifest = (
            pathlib.Path(__file__).parent.parent / "resources" / "quant_index_manifest.yaml"
        )
        qconfig = QueryPipelineConfig(
            index=quant_manifest,
            kmer_ratio=self.kmer_ratio,
            index_to_method_mapping={"IT_N_103": "jellyfish"},
        )
        qconfig.config["query"].update({"query_fasta": str(self.fasta)})

        pipeline = QueryPipeline(
            self.pipeline, self.profile, qconfig, pathlib.Path("./"), target_rule="all"
        )
        result = pipeline.run_pipeline(slurm=False, cores=1)
        assert isinstance(result, QueryPipelineResult)
        assert isinstance(result.query_path, list)
        assert result.query_path == [
            (
                "jellyfish",
                "IT_N_103",
                pathlib.Path("./").resolve() / "query/jellyfish/IT_N_103/quantitative_search.tsv",
            )
        ]
