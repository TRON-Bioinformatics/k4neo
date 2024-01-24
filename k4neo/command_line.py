
import pathlib
from argparse import ArgumentParser
from logzero import logger
import k4neo
from k4neo.database.database import DataBase, CreateDataBase
from k4neo.annotator.annotator import Annotator
from k4neo.parser.parser import IndexResultParser
from k4neo.pipeline.query_pipeline import IndexPipeline, IndexPipelineConfig


def build_database():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} database builder"
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
        description=f"k4neo {k4neo.VERSION} index register"
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
        description=f"k4neo {k4neo.VERSION} result parser"
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
        description=f"k4neo {k4neo.VERSION} indexing pipeline"
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
        '--kmer_size',
        dest='kmer_size',
        help='Size of k-mers',
        default=21
    )
    parser.add_argument(
        '--cutoff',
        dest='cutoff',
        help='Cutoff to define solid and weak k-mers',
        default=2
    )
    parser.add_argument(
        '--fpr',
        dest='fpr',
        help='Theoretic False Positive Rate of individual bloom filters',
        default=0.05
    )
    parser.add_argument(
        '--workflow',
        dest='workflow',
        help='path to k4neo pipeline',
        default=pathlib.Path(__file__).parent / 'pipeline' / 'workflow' / 'Snakefile'
    )
    parser.add_argument(
        '--slurm',
        dest='slurm',
        help='Submit pipeline jobs to slurm',
        action='store_true'
    )
    args = parser.parse_args()
    samples = args.samples
    working_dir = pathlib.Path(args.index).resolve()
    pipeline = pathlib.Path(args.workflow).resolve()
    pipeline_config = IndexPipelineConfig(samples=samples,
                                          method="raptor",
                                          kmer_size=args.kmer_size,
                                          cutoff=args.cutoff,
                                          fpr=args.fpr)
    logger.info("Starting indexing pipeline...")
    index_pipeline = IndexPipeline(pipeline,
                                   working_dir=working_dir,
                                   config=pipeline_config.config)
    index_pipeline.run_pipeline(slurm=args.slurm)
    logger.info("Finished indexing pipeline...")

    return

def annotate():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} annotator"
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
        help='Number of shared k-mers between query and sample to report as hit',
        default=0.45
    )
    parser.add_argument(
        '--working-dir',
        dest='working_dir',
        help='Working directory of k4neo pipeline',
        default="./k4neo_query"
    )
    parser.add_argument(
        '--workflow',
        dest='workflow',
        help='path to k4neo pipeline',
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
        description=f"Run k4neo {k4neo.VERSION} query pipeline for testing"
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
    #annotate()
    build_index()

if __name__ == '__main__':
    main()
