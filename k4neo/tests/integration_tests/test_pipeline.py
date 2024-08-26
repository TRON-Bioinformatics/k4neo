from k4neo.pipeline.query_pipeline import QueryPipeline, QueryPipelineConfig, QueryPipelineResult
import tempfile
import pathlib

class TestPipeline():
    index_manifest = pathlib.Path(__file__).parent.parent /'resources' / 'index_manifest.yaml'
    kmer_ratio = 0.7
    pipeline = pathlib.Path(__file__).parent.parent.parent /'pipeline' / 'tronmake-kmer-pipeline' / 'workflow' / 'Snakefile'
    fasta = pathlib.Path(__file__).parent.parent /'resources' / 'example.fasta'

    def test_raptor_query_pipeline(self):
        qconfig = QueryPipelineConfig(index=self.index_manifest, kmer_ratio=self.kmer_ratio, methods=['raptor'])
        qconfig.config["query"].update({'query_fasta': str(self.fasta)})
        pipeline = QueryPipeline(self.pipeline, qconfig.config, pathlib.Path("./"))
        result = pipeline.run_pipeline(slurm=False, cores=1)
        assert isinstance(result, QueryPipelineResult)
        assert isinstance(result.query_path, dict)
        assert "raptor" in result.query_path.keys()
        assert isinstance(result.query_path['raptor'], pathlib.PosixPath)