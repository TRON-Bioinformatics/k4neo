from k4neo.index.index import KmerIndex
import pathlib


class TestKmerIndexInterface():
    index_manifest = pathlib.Path(__file__).parent.parent /'resources' / 'index_manifest.yaml'
    kmer_ratio = 0.7
    pipeline = pathlib.Path(__file__).parent.parent.parent /'pipeline' / 'tronmake-kmer-pipeline' / 'workflow' / 'Snakefile'
    profile = pathlib.Path(__file__).parent.parent.parent /'pipeline' / 'config.yaml'
    fasta = pathlib.Path(__file__).parent.parent /'resources' / 'example.fasta'
    kmer_index = KmerIndex(pipeline=pipeline, workflow_profile=profile, index_manifest=index_manifest, kmer_ratio=kmer_ratio)

    def test_kmer_interface(self):
        # Test that index manifest is read correctly    
        struct_parsed = self.kmer_index.read_index_struct()
        assert struct_parsed == {'test_index': {'samples': 2, 'path': 'k4neo/tests/resources/index/raptor.index', 'sample_mapping': 'k4neo/tests/resources/index/index_mapping.txt', 'method': 'raptor'}}
        # Test that getting methods from manifest works
        methods_parsed_from_struct = self.kmer_index._get_index_methods()
        assert methods_parsed_from_struct == set(['raptor',])
        # Test that indices with missing index method are ignored
        example_struct = {'test_index': 
            {'samples': 2, 'path': 'k4neo/tests/resources/index/raptor.index', 'sample_mapping': 'k4neo/tests/resources/index/index_mapping.txt', 'method': 'raptor'},
            'test_index2': 
            {'samples': 2, 'path': 'k4neo/tests/resources/index/raptor.index', 'sample_mapping': 'k4neo/tests/resources/index/index_mapping.txt', 'method': 'kmindex'},
            'test_index3': 
            {'samples': 2, 'path': 'k4neo/tests/resources/index/raptor.index', 'sample_mapping': 'k4neo/tests/resources/index/index_mapping.txt'}
            }
        self.kmer_index.index_struct = example_struct
        methods_parsed_from_struct = self.kmer_index._get_index_methods()
        assert methods_parsed_from_struct == set(['raptor', 'kmindex'])

    
