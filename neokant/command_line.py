import os
import sys
import pathlib
from argparse import ArgumentParser
sys.path.append(XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
from logzero import logger
import neokant
from neokant.database.database import DataBase, CreateDataBase
from neokant.annotator.annotator import Annotator
from neokant.parser.parser import IndexResultParser

sys.path.append(XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX/neokant")


def build_database():
    parser = ArgumentParser(
        description=f"neoKant {neokant.VERSION} database builder"
    )
    parser.add_argument(
        '--sample-tables',
        dest='sample_table',
        help='Archive with standardized sample metadata tables to include in annotation db',
    )
    parser.add_argument(
        '--database',
        dest='database',
        help='Database file to create'
    )
    parser.add_argument(
        '--tissue-map',
        dest='tissue_map',
        help='Mapping of different tissue identifiers found across public data to their corresponding GTEx identifier'

    )
    args = parser.parse_args()
    logger.info("Starting database creation...")
    db = CreateDataBase(args.database, args.sample_table, args.tissue_map)
    db.setup_db()
    db.precomputations()
    logger.info("Created database")


def register_index():
    parser = ArgumentParser(
        description=f"neoKant {neokant.VERSION} index register"
    )
    parser.add_argument(
        '--database',
        dest='database',
        help='Annotation database file.',
        required=True
    )
    parser.add_argument(
        '--index',
        dest='index',
        help='kmer index to query.',
        required=True
    )
    parser.add_argument(
        '--method',
        dest='method',
        help='K-mer indexing method',
        required=True
    )
    parser.add_argument(
        '--global-index',
        dest='global_index',
        help="Directory where index should be registered"
    )
    args = parser.parse_args()
    logger.info("Index registration...")


def parse_output():
    parser = ArgumentParser(
        description=f"neoKant {neokant.VERSION} result parser"
    )
    parser.add_argument(
        '--index-results',
        dest='index_table',
        help='Query results from k-mer index',
    )
    parser.add_argument(
        '--tool',
        dest='tool',
        help='indexing method'
    )
    parser.add_argument(
        '--sample-table',
        dest='sample_table',
        help='Sample table mapping index ids to samples. Only required for Reindeer and Raptor HIBF'

    )
    parser.add_argument(
        '--output-table',
        dest='output_table',
        help='Output table'

    )
    parser.add_argument(
        '--kmindex-cutoff',
        dest='kmindex_cutoff',
        help='Kmindex cutoff to filter results',
        default=0.7,
        type=float
    )
    args = parser.parse_args()
    reindeer_sample_mapping = None
    raptor_sample_mapping = None
    if args.tool == 'raptor':
        raptor_sample_mapping = args.sample_table
    if args.tool == 'reindeer':
        reinder_sample_mapping = args.sample_table
    logger.info("Starting query result parsing...")
    parser = IndexResultParser(args.index_table, args.tool,
                               reindeer_sample_mapping=reindeer_sample_mapping,
                               raptor_sample_mapping=raptor_sample_mapping,
                               kmindex_cutoff=args.kmindex_cutoff)
    results = parser.parse_results()
    parser.write_result(results, args.output_table)
    logger.info("Parsed query results into table")


def build_index():
    parser = ArgumentParser(
        description=f"neoKant {neokant.VERSION} indexing pipeline"
    )
    parser.add_argument(
        '--samples',
        dest='samples',
        help='Sample sheet describing fastqs to index',
    )
    parser.add_argument(
        '--index',
        dest='index',
        help='Path to index directory'
    )
    parser.add_argument(
        '--method',
        dest='method',
        help='Indexing method to use',
        default="raptor"
    )
    args = parser.parse_args()
    if args.method not in ['raptor', 'cobs', 'reindeer', 'kmindex']:
        raise ValueError()

    logger.info("Starting indexing pipeline...")
    logger.info(f"Indexing method: {args.method}")
    return


def annotate():
    parser = ArgumentParser(
        description=f"neoKant {neokant.VERSION} annotator"
    )
    parser.add_argument(
        '--database',
        dest='database',
        help='Annotation database file.',
        required=True
    )
    parser.add_argument(
        '--index',
        dest='index',
        help='kmer index to query.',
        required=True
    )
    parser.add_argument(
        '--queries',
        dest='queries',
        help='Tabular format with context sequence and position of interest',
        required=True
    )
    parser.add_argument(
        '--output',
        dest='output',
        help='Annotated output sequences'
    )
    parser.add_argument(
        '--method',
        dest='method',
        help='K-mer indexing method',
        required=True
    )
    parser.add_argument(
        '--ratio',
        dest='kmer_ratio',
        help='Number of shared ki-mers between query and sample to report as hit',
        default=0.45
    )
    parser.add_argument(
        '--working-dir',
        dest='working_dir',
        help='Working directory of neoKant pipeline',
        default="./neokant_query"
    )
    parser.add_argument(
        '--workflow',
        dest='workflow',
        help='path to neokant pipeline',
        default=pathlib.Path(__file__).parent / 'pipeline' / 'workflow' / 'Snakefile'
    )
    parser.add_argument(
        '--sample-table',
        dest='sample_table',
        help='Sample table mapping index ids to samples. Only required for Reindeer and Raptor HIBF'
    )
    parser.add_argument(
        '--kmindex-cutoff',
        dest='kmindex_cutoff',
        help='Kmindex cutoff to filter results',
        default=0.7,
        type=float
    )
    args = parser.parse_args()
    reindeer_sample_mapping = None
    raptor_sample_mapping = None
    if args.method == 'raptor':
        raptor_sample_mapping = args.sample_table
    if args.method == 'reindeer':
        reinder_sample_mapping = args.sample_table
    working_dir = pathlib.Path(args.working_dir).resolve()
    pipeline = pathlib.Path(args.workflow).resolve()
    annotator = Annotator(working_dir,
                          args.queries,
                          args.database)
    annotator.prepare_cts()
    result = annotator.search_cts(pipeline=pipeline,
                                  index=pathlib.Path(args.index),
                                  method=args.method,
                                  reindeer_sample_mapping=reindeer_sample_mapping,
                                  raptor_sample_mapping=raptor_sample_mapping,
                                  kmer_ratio=args.kmer_ratio,
                                  kmindex_cutoff=args.kmindex_cutoff)
    results = annotator.annotate_cts(result)
    results = annotator.annotate_sequences(results)
    with open(args.output, "w") as file_handle:
        results.to_csv(file_handle, sep="\t", index=False)


def query_pipeline():
    parser = ArgumentParser(
        description=f"Run neoKant {neokant.VERSION} query pipeline for testing"
    )
    parser.add_argument(
        '--query-fasta',
        dest='query_fasta',
        help='FASTA file with query sequences',
    )
    parser.add_argument(
        '--index',
        dest='index',
        help='Index/Indexing directory'
    )
    parser.add_argument(
        '--method',
        dest='method',
        help='Method used to create index',
        default="raptor"
    )
    args = parser.parse_args()
    if args.method not in ['raptor', 'cobs', 'reindeer', 'kmindex']:
        raise ValueError()

    logger.info("Starting query pipeline...")

    return

## To test the functionality

def main():
    #build_database()
    #parse_output()
    annotate()

if __name__ == '__main__':
    main()
