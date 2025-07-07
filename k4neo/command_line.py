import sys
import pathlib
import gzip
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import k4neo
from k4neo.database.database import DataBase, CreateDataBase
from k4neo.annotator.annotator import Annotator
from k4neo.parser.parser import IndexResultParser
from k4neo.pipeline.query_pipeline import IndexPipeline, IndexPipelineConfig
from k4neo.plotter.plotter import Plotter
#from loguru import logger
from k4neo.setup_logging import setup_logging
from rich.console import Console
from tqdm import tqdm
import gc
import pandas as pd
import time

console = Console()
#logger.remove()
#logger.add(lambda msg: tqdm.write(msg, end=""), level="INFO", colorize=True)

epilog = "Copyright (c) 2025 TRON gGmbH (See LICENSE for licensing details)"


def build_database():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} database builder",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--sample-tables",
        dest="sample_table",
        help="Archive with standardized sample metadata tables to include in annotation db",
        required=True,
    )
    parser.add_argument(
        "--database", dest="database", help="Database file to create", required=True
    )
    parser.add_argument(
        "--tissue-map",
        dest="tissue_map",
        help="Mapping of different tissue identifiers found across public data to their corresponding GTEx identifier",
        required=True,
    )
    args = parser.parse_args()
    logger.info("-> Starting metadata database creation.")
    db = CreateDataBase(args.database, args.sample_table, args.tissue_map)
    db.setup_db()
    db.precomputations()
    logger.info("-> Finished metadata database creation.")


def parse_output():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} parser",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--index-results",
        dest="index_table",
        help="Query results from k-mer index",
    )
    parser.add_argument("--tool", dest="tool", help="indexing method")
    parser.add_argument("--output-table", dest="output_table", help="Output table")
    parser.add_argument(
        "--kmindex-cutoff",
        dest="kmindex_cutoff",
        help="Kmindex cutoff to filter results",
        default=0.7,
        type=float,
    )
    args = parser.parse_args()
    logger.info("Starting query result parsing...")
    parser = IndexResultParser(
        args.index_table,
        args.tool,
        raptor_sample_mapping=args.raptor_sample_mapping,
        kmindex_cutoff=args.kmindex_cutoff,
    )
    results = parser.parse_results()
    parser.write_result(results, args.output_table)
    logger.info("Parsed query results into table")

def annotate():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} annotator",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--database", dest="database", help="Annotation database file.", required=True
    )
    parser.add_argument(
        "--index", dest="index_manifest", help="k-mer index to query.", required=True
    )
    parser.add_argument(
        "--queries",
        dest="queries",
        help="Tabular format with context sequence and position of interest",
        required=True,
    )
    parser.add_argument("--output", dest="output", help="Output prefix for annotated sequences")
    parser.add_argument(
        "--ratio",
        dest="kmer_ratio",
        help="Number of shared k-mers between query and sample to report as hit",
        default=0.7,
        type=float,
    )
    parser.add_argument(
        "--working-dir",
        dest="working_dir",
        help="Working directory of k4neo pipeline",
        default="./k4neo_query",
    )
    parser.add_argument(
        "--workflow",
        dest="workflow",
        help="path to tronmake k-mer pipeline",
        default=pathlib.Path(__file__).parent
        / "pipeline"
        / "tronmake-kmer-pipeline"
        / "workflow"
        / "Snakefile",
    )
    parser.add_argument(
        "--profile",
        dest="workflow_profile",
        help="A yaml file containing snakemake options for execution",
        default=pathlib.Path(__file__).parent / "pipeline" / "default_profile.yaml",
    )
    parser.add_argument(
        "--kmer",
        dest="kmer_size",
        help="K-mer size of search index",
        default=21,
        type=int,
    )
    parser.add_argument(
        "--cpu",
        dest="cpu",
        help="Number of cpus for local execution",
        default=16,
        type=int,
    )
    parser.add_argument(
        "--slurm", dest="slurm", help="Submit query job to slurm", action="store_true"
    )
    parser.add_argument(
        "--chunk-size",
        dest="chunk_size",
        help="Chunk size for processing input sequences",
        type=int,
        default=10000,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Verbose logs (default: False)",
    )
    args = parser.parse_args()
    console.print(k4neo.logo, style="bold red")


    # Create pipeline workdir if not existent
    working_dir = pathlib.Path(args.working_dir).resolve()
    working_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logger
    output_directory = pathlib.Path(args.output).parent
    log_file_name = pathlib.Path(output_directory) / "k4neo.log"

    logger = setup_logging(log_file_name, args.verbose)

    pipeline = pathlib.Path(args.workflow).resolve()
    workflow_profile = pathlib.Path(args.workflow_profile).resolve()
    index_manifest = pathlib.Path(args.index_manifest).resolve()
    
    annotator = Annotator(working_dir, args.queries, args.database)
    annotator.prepare_cts()
    result_dict = annotator.search_cts(
        pipeline=pipeline,
        workflow_profile=workflow_profile,
        index_manifest=index_manifest,
        kmer_ratio=args.kmer_ratio,
        cores=args.cpu,
        slurm=args.slurm,
    )

    # Write non-queryable sequences to disk. These are the sequences that shorter than the k-mer size
    output_non_queryable = pathlib.Path(args.output + "_non_querable.tsv.gz")
    if len(annotator.non_queryable.index > 0):
        logger.debug("-> Writing non-queryable sequences to disk")
        with gzip.open(output_non_queryable, "wb") as file_handle:
            annotator.non_queryable.to_csv(file_handle, sep="\t", index=False)

    # Result could come from different indexing methods such as binary and quantitative indices.
    for this_method, this_df in result_dict.items():
        logger.info(f"Annotating {this_method} query results")

        # Group df by cts_id and determine number of chunks
        grouped_df = this_df.groupby("cts_id")
        group_keys = list(grouped_df.groups.keys())
        num_chunks = (len(group_keys) + args.chunk_size - 1) // args.chunk_size

        # Define final output files
        output_annotated = pathlib.Path(args.output + f"_annotated_{this_method}.tsv.gz")
        output_healthy_rate = pathlib.Path(
            args.output + f"_healthy_sample_rate_{this_method}.tsv.gz"
        )
        output_tumor_rate = pathlib.Path(args.output + f"_tumor_sample_rate_{this_method}.tsv.gz")
        
        first_chunk = True
        with (
            gzip.open(output_annotated, "wb") as f_annotated,
            gzip.open(output_healthy_rate, "wb") as f_healthy,
            gzip.open(output_tumor_rate, "wb") as f_tumor,
        ):

            for i in tqdm(
                range(0, len(group_keys), args.chunk_size), total=num_chunks, desc="Processing chunks"
            ):
                chunk_t0 = time.perf_counter()
                # Collect cts_ids belonging to currently processed chunk and get rows
                chunk_keys = group_keys[i : i + args.chunk_size]
                chunk = pd.concat([grouped_df.get_group(k) for k in chunk_keys])

                # Run k4neo annotator functions on chunk
                results = annotator.annotate_cts(chunk)
                sample_hits = annotator.annotate_sequences(results)
                healthy_sample_rate, tumor_sample_rate = annotator.annotate_sample_rate2(results)
                
                chunk_t1 = time.perf_counter()
                logger.debug(f"[chunk {i}] required: {(chunk_t1 - chunk_t0):.2f} s")

                # Write to disk
                sample_hits.to_csv(f_annotated, sep="\t", index=False, header=first_chunk)
                healthy_sample_rate.to_csv(f_healthy, sep="\t", index=False, header=first_chunk)
                tumor_sample_rate.to_csv(f_tumor, sep="\t", index=False, header=first_chunk)
                
                first_chunk = False  # Stop writing header after first chunk

                del results
                del sample_hits
                del healthy_sample_rate
                del tumor_sample_rate
                gc.collect()
    
        logger.info("Finished k4neo analysis.")


def plot():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} plotter",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--input", dest="input", help="Input prefix for annotated sequences", required=True
    )
    parser.add_argument("--output", dest="output", help="Output file", required=True)
    
    args = parser.parse_args()
    console.print(k4neo.logo, style="bold red")

    pl = Plotter(
        args.input + "_healthy_sample_rate_raptor.tsv.gz",
        args.input + "_tumor_sample_rate_raptor.tsv.gz",
    )
    pl.plot(args.output)
