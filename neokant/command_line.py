import os
import sys
import neoKant
from argparse import ArgumentParser
from logzero import logger


def build_database():
    parser = ArgumentParser(
        description=f"neoKant {neoKant.version} database builder"
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
    return


def build_index():
    parser = ArgumentParser(
        description=f"neoKant {neoKant.version} indexing pipeline"
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
        description=f"neoKant {neoKant.version} annotator"
    )
    parser.add_argument(
        '--database',
        dest='database',
        help='Annotation database file.'
    )
    parser.add_argument(
        '--index',
        dest='index',
        help='kmer index to query.'
    )
    parser.add_argument(
        '--search',
        dest='search',
        help='Tabular format with context sequence and position of interest'
    )

    args = parser.parse_args()
    logger.info("Starting annotation ...")
    return


def query_pipeline():
    parser = ArgumentParser(
        description=f"Run neoKant {neoKant.version} query pipeline for testing"
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

