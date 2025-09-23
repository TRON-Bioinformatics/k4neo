import sys
import pathlib
import gzip
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import k4neo
from k4neo.database.database import DataBase, CreateDataBase
from k4neo.annotator.annotator import Annotator
from k4neo.annotator.reference_annotation import ReferenceIndexer, KmerUniquenessAnnotator
from k4neo.parser.parser import IndexResultParser
from k4neo.parser.index_parser import IndexResultParser2
from k4neo.pipeline.query_pipeline import IndexPipeline, IndexPipelineConfig
from k4neo.plotter.plotter import Plotter
from k4neo.setup_logging import setup_logging
from k4neo.helper.helper import DiskIO
from k4neo.helper.async_writer import AsyncDFWriter
from rich.console import Console
from tqdm import tqdm
import gc
import pandas as pd
import time

console = Console()

epilog = "Copyright (c) 2025 TRON gGmbH (See LICENSE for licensing details)"


def process_chunk(chunk, annotator):
    results = annotator.annotate_cts(chunk)
    sample_hits = annotator.annotate_sequences(results)
    healthy_sample_rate, tumor_sample_rate = annotator.annotate_sample_rate2(results)
    return sample_hits, healthy_sample_rate, tumor_sample_rate


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

    log_file_name = pathlib.Path(args.database).parent / "k4neo_db_build.log"

    logger = setup_logging(log_file_name, verbose=True)

    logger.info("-> Starting metadata database creation.")
    db = CreateDataBase(args.database, args.sample_table, args.tissue_map)
    db.setup_db()
    db.precomputations()
    logger.info("-> Finished metadata database creation.")


def build_ref_index():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} reference index",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--genome",
        dest="genome",
        help="Genome fasta file",
        required=True,
    )
    parser.add_argument(
        "--transcriptome",
        dest="transcriptome",
        help="Transcriptome fasta file",
        required=True,
    )
    parser.add_argument("--kmer", dest="kmer_size", help="K-mer size", default=21)
    parser.add_argument(
        "--output",
        dest="output",
        help="Output directory",
        required=True,
    )
    args = parser.parse_args()

    log_file_name = pathlib.Path(args.output).parent / "k4neo_ref_index.log"

    logger = setup_logging(log_file_name, verbose=True)

    ReferenceIndexer.prepare_data(args.genome, args.transcriptome, args.kmer_size, args.output)


def annotate_uniqueness():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} uniqueness annotation",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--fasta",
        dest="queries",
        help="FASTA file",
        required=True,
    )
    parser.add_argument(
        "--reference_indices",
        dest="ref_index",
        help="JellyFish indices of genome and transcriptome",
        required=True,
    )
    parser.add_argument(
        "--output",
        dest="output",
        help="Tabular output with uniqueness annotation",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Verbose logs (default: False)",
    )
    args = parser.parse_args()

    log_file_name = pathlib.Path(args.output).parent / "k4neo_uniqueness_annot.log"

    logger = setup_logging(log_file_name, args.verbose)

    logger.info("Calculating uniqueness of sequences")

    uniq = KmerUniquenessAnnotator(args.ref_index)
    results = uniq.annotate_fasta(args.queries)

    counter = 0
    with open(args.output, "w") as file_handle:
        file_handle.write(
            "cts_id\tcts_unique_rate\tcts_ref_rate\tcts_ref_single_gene_locus_rate\tcts_ref_multi_gene_locus_rate\tcts_ref_single_transcript_rate\tcts_ref_multi_transcript_rate\n"
        )
        for seq_id, rates in results.items():
            file_handle.write(
                f"{seq_id}\t{rates['cts_unique_rate']:.4f}\t{rates['cts_ref_rate']:.4f}\t{rates['cts_ref_single_gene_locus_rate']:.4f}\t{rates['cts_ref_multi_gene_locus_rate']:.4f}\t{rates['cts_ref_single_transcript_rate']:.4f}\t{rates['cts_ref_multi_transcript_rate']:.4f}\n"
            )
            counter += 1
    logger.info(f"Annotated uniqueness of {counter} sequences")


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
    parser.add_argument(
        "--compress",
        dest="compression",
        help="Compress final output files with gzip",
        action="store_true",
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
    # Write non-queryable sequences to disk
    if len(annotator.non_queryable.index > 0):
        logger.info("-> Writing non-queryable sequences to disk")
        output_non_queryable = (
            pathlib.Path(args.output + "_non_querable.tsv.gz")
            if args.compression
            else pathlib.Path(args.output + "_non_querable.tsv")
        )
        DiskIO.write_df(annotator.non_queryable, output_non_queryable, args.compression)

    # Write a debug table that maps cts_ids to query_ids
    annotator.sequence_table[["cts_id", "query_cts_id"]].to_csv(
        pathlib.Path(args.output + "_cts_to_query_cts.tsv"), sep="\t", index=False
    )

    for method_name in result_dict:
        total_cts = len(result_dict[method_name])
        logger.info(f"-> Annotating query results of method: {method_name}")
        pbar = tqdm(total=total_cts, desc=f"Processing {method_name}")

        if args.compression:
            output_annotated = pathlib.Path(args.output + f"_annotated_{method_name}.tsv.gz")
            output_healthy_rate = pathlib.Path(
                args.output + f"_healthy_sample_rate_{method_name}.tsv.gz"
            )
            output_tumor_rate = pathlib.Path(
                args.output + f"_tumor_sample_rate_{method_name}.tsv.gz"
            )
        else:
            output_annotated = pathlib.Path(args.output + f"_annotated_{method_name}.tsv")
            output_healthy_rate = pathlib.Path(
                args.output + f"_healthy_sample_rate_{method_name}.tsv"
            )
            output_tumor_rate = pathlib.Path(args.output + f"_tumor_sample_rate_{method_name}.tsv")

        first_chunk = True

        # Start writer threads
        healthy_writer = AsyncDFWriter(output_healthy_rate, compression=args.compression)
        healthy_writer.start()

        tumor_writer = AsyncDFWriter(output_tumor_rate, compression=args.compression)
        tumor_writer.start()

        annot_writer = AsyncDFWriter(output_annotated, compression=args.compression)
        annot_writer.start()

        for _, batch_len, chunk in IndexResultParser2.generate_dataframe_in_batches(
            {method_name: result_dict[method_name]}, batch_size=args.chunk_size
        ):

            # k4neo Annotation functions
            results = annotator.annotate_cts(chunk)
            sample_hits = annotator.annotate_sequences(results)
            healthy_sample_rate, tumor_sample_rate = annotator.annotate_sample_rate2(results)

            # Send metrics to background writer threads
            annot_writer.write(
                sample_hits,
                [
                    "cts_id",
                    "count",
                    "total",
                    "disease",
                    "developmental_stage",
                    "tissue",
                    "study_id",
                ],
                append=not first_chunk,
                header=first_chunk,
            )

            healthy_writer.write(
                healthy_sample_rate,
                ["cts_id", "developmental_stage", "tissue", "sample_rate"],
                append=not first_chunk,
                header=first_chunk,
            )

            tumor_writer.write(
                tumor_sample_rate,
                ["cts_id", "disease", "tissue", "sample_rate"],
                append=not first_chunk,
                header=first_chunk,
            )

            first_chunk = False  # turn off headers after first write
            pbar.update(batch_len)

        pbar.close()
        logger.info("Waiting for writer threads to finish")
        # Wait for writer threads to finish and close
        annot_writer.wait_until_done()
        healthy_writer.wait_until_done()
        tumor_writer.wait_until_done()

        annot_writer.stop()
        healthy_writer.stop()
        tumor_writer.stop()


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
