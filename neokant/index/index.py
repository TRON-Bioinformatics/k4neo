import os
import sys
import pathlib
from neokant.pipeline.query_pipeline import QueryPipeline
from logzero import logger





class KmerIndex(object):
    def __init__(self, index, index_options):
        self.index = pathlib.Path(index)
        self.index_options = pathlib.Path(index_options)

    def search_index(self):


