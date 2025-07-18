import pathlib
from k4neo.parser.index_parser import IndexResultParser2, BinaryKmerIndexResultParser

#class TestIndexResultParser():
#    pass



class TestBinaryKmerIndexResultPaser:
    
    def test_parsing_raptor(self):
        parser = BinaryKmerIndexResultParser(
            pathlib.Path(__file__).parent.parent / "resources" / "example_index_output" / "search_index_raptor.tsv",
            "raptor",
            pathlib.Path(__file__).parent.parent / "resources" / "example_index_output" / "sample2minimiser.tsv",
            kmer_ratio = 0.7
        )
        
        results = parser.parse_results()

        assert results == {
            "d2e9c05e86abc70c804fc48a7f25aad8": {"MCF7"},
            "b410c5db97c66dc7a4d200ec44ce4a1d": {None},
            "61c9762eae1072735d7ff13834da8942": {None},
            "9221002eafd64f81d677be49d5b52885": {"MCF7", "SKBR3"},
            "4d2420d9a64523b789dee5a2ff0e1a3a": {"MCF7"},
            "ae618273613b012a6008ade29038e65c": {"MCF7"},
            "dc7611c0669d2afb1961b4e6c21f453d": {"MCF7"},
            "fe4058182debdb142dc93b3c0223640a": {"MCF7"},
            "8aebbe56abfaba55fcc2821dbfc2e9c9": {"MCF7"},
            "3b853efc1d5fdea096ac98f6c783cdf5": {"MCF7", "SKBR3"},
            "7bcda9de85306099b5e5dd15c5824249": {"MCF7", "SKBR3"},
            "eeb983a06674aaa92b850e0c94076030": {"MCF7", "SKBR3"},
            "a9588b5674a5d5d1c679e34a88a17d7a": {"MCF7", "SKBR3"},
            "96dbc2491a9573ec1ae08d6a271b913d": {"MCF7"}
        }
    
    def test_parsing_kmindex(self):
        
        parser = BinaryKmerIndexResultParser(
            pathlib.Path(__file__).parent.parent / "resources" / "example_index_output" / "search_index_kmindex.tsv",
            "kmindex",
            "",
            kmer_ratio = 0.7
        )
        results = parser.parse_results()

        assert results == {
            "d2e9c05e86abc70c804fc48a7f25aad8": {"MCF7"},
            "b410c5db97c66dc7a4d200ec44ce4a1d": {None},
            "61c9762eae1072735d7ff13834da8942": {None},
            "9221002eafd64f81d677be49d5b52885": {"MCF7", "SKBR3"},
            "4d2420d9a64523b789dee5a2ff0e1a3a": {"MCF7"},
            "ae618273613b012a6008ade29038e65c": {"MCF7"},
            "dc7611c0669d2afb1961b4e6c21f453d": {"MCF7"},
            "fe4058182debdb142dc93b3c0223640a": {"MCF7"},
            "8aebbe56abfaba55fcc2821dbfc2e9c9": {"MCF7"},
            "3b853efc1d5fdea096ac98f6c783cdf5": {"MCF7", "SKBR3"},
            "7bcda9de85306099b5e5dd15c5824249": {"MCF7", "SKBR3"},
            "eeb983a06674aaa92b850e0c94076030": {"MCF7", "SKBR3"},
            "a9588b5674a5d5d1c679e34a88a17d7a": {"MCF7", "SKBR3"},
            "96dbc2491a9573ec1ae08d6a271b913d": {"MCF7"}
        }
