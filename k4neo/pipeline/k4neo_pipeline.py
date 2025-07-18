import os
import sys

class K4neoPipeline:
    """
    Genereic class to represent entire k4neo run.
    Starting with preparing the context sequences, 
    submitting the k-mer pipeline and annotating the results
    """
    def __init__(self, pipeline, workflow_profile, index_manifest, kmer_ratio, cores, slurm, output):
        self.annotator = ""
        pass

    def _prepare(self):
        self.annotator = ""
        pass

    def _run_kmer_search(self):
        pass

    def _annotate(self, chunk_size):
    