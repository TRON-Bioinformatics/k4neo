from k4neo.pipeline.query_pipeline import QueryPipelineResult
from k4neo.index.index_loader import load_metaindex_from_manifest
from k4neo.index.index_processor import KmerIndexProcessor
import pathlib


class TestKmerIndexProcessor:
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

    def test_kmer_index_processor(self):
        index_manifest = pathlib.Path(__file__).parent.parent / "resources" / "index_manifest.yaml"
        # Load meta index from manifest using helper function. Typical use case
        meta_index = load_metaindex_from_manifest(index_manifest)

        processor = KmerIndexProcessor(
            meta_index=meta_index,
            pipeline=self.pipeline,
            workflow_profile=self.profile,
            quantitative=False,
        )

        assert processor.index_to_method_mapping == {"test_index": "raptor"}
        assert processor.index_to_sample_mapping == {
            "test_index": pathlib.Path("k4neo/tests/resources/index/index_mapping.txt")
        }

        result = processor.search_index(
            query_sequences=self.fasta,
            working_dir=pathlib.Path("./"),
            slurm=False,
            kmer_ratio=self.kmer_ratio,
            cores=1,
        )
        assert isinstance(result, QueryPipelineResult)
        assert isinstance(result.query_path, list)

        parsed_results = processor.result_parser2(result, cores=1, kmer_ratio=self.kmer_ratio)

        assert isinstance(parsed_results, dict)
        assert "raptor" in parsed_results
        assert set(parsed_results.keys()) == {
            "raptor",
        }
        assert set(parsed_results["raptor"].keys()) == {
            "dde52f9b214f3e15",
            "b855999df5ed524d",
            "845644f2045d1b52",
            "09f4de7a7f434f4d",
            "8b0ee74a7458e72c",
            "3eb0221fe5e90c3c",
        }

    def test_kmer_index_processor_quantitative(self):
        index_manifest = (
            pathlib.Path(__file__).parent.parent / "resources" / "quant_index_manifest.yaml"
        )
        # Load meta index from manifest using helper function. Typical use case
        meta_index = load_metaindex_from_manifest(index_manifest)

        processor = KmerIndexProcessor(
            meta_index=meta_index,
            pipeline=self.pipeline,
            workflow_profile=self.profile,
            quantitative=True,
        )

        assert processor.index_to_method_mapping == {"IT_N_103": "jellyfish"}
        assert processor.kmer_depth_mapping == {"IT_N_103": 22925084}

        result = processor.search_index(
            query_sequences=self.fasta,
            working_dir=pathlib.Path("./"),
            slurm=False,
            kmer_ratio=self.kmer_ratio,
            cores=1,
        )
        assert isinstance(result, QueryPipelineResult)
        assert isinstance(result.query_path, list)
