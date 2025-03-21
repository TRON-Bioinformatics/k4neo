import sys
import pathlib
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from logzero import logger
import k4neo
from k4neo.database.database import DataBase, CreateDataBase
from k4neo.annotator.annotator import Annotator
from k4neo.parser.parser import IndexResultParser
from k4neo.pipeline.query_pipeline import IndexPipeline, IndexPipelineConfig

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
        description=f"k4neo {k4neo.VERSION} result parser",
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


def build_index():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} indexing pipeline",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--samples",
        dest="samples",
        help="Sample sheet describing fastqs to index",
    )
    parser.add_argument("--index", dest="index", help="Path to index directory")
    parser.add_argument(
        "--kmer_size", dest="kmer_size", help="Size of k-mers", default=21
    )
    parser.add_argument(
        "--cutoff",
        dest="cutoff",
        help="Cutoff to define solid and weak k-mers",
        default=2,
    )
    parser.add_argument(
        "--fpr",
        dest="fpr",
        help="Theoretic False Positive Rate of individual bloom filters",
        default=0.05,
    )
    parser.add_argument(
        "--workflow",
        dest="workflow",
        help="path to k4neo pipeline",
        default=pathlib.Path(__file__).parent / "pipeline" / "workflow" / "Snakefile",
    )
    parser.add_argument(
        "--slurm",
        dest="slurm",
        help="Submit pipeline jobs to slurm",
        action="store_true",
    )
    args = parser.parse_args()
    samples = args.samples
    working_dir = pathlib.Path(args.index).resolve()
    pipeline = pathlib.Path(args.workflow).resolve()
    pipeline_config = IndexPipelineConfig(
        samples=samples,
        method="raptor",
        kmer_size=args.kmer_size,
        cutoff=args.cutoff,
        fpr=args.fpr,
    )
    logger.info("-> Starting indexing pipeline.")
    index_pipeline = IndexPipeline(
        pipeline, working_dir=working_dir, config=pipeline_config.config
    )
    index_pipeline.run_pipeline(slurm=args.slurm)
    logger.info("-> Finished indexing pipeline.")

    return


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
    parser.add_argument(
        "--output", dest="output", help="Output prefix for annotated sequences"
    )
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
    args = parser.parse_args()
    working_dir = pathlib.Path(args.working_dir).resolve()
    pipeline = pathlib.Path(args.workflow).resolve()
    index_manifest = pathlib.Path(args.index_manifest).resolve()
    annotator = Annotator(working_dir, args.queries, args.database)
    annotator.prepare_cts()
    result_dict = annotator.search_cts(
        pipeline=pipeline,
        index_manifest=index_manifest,
        kmer_ratio=args.kmer_ratio,
        cores=args.cpu,
        slurm=args.slurm,
    )
    # Write non-queryable sequences to disk
    output_non_queryable = pathlib.Path(args.output + "_non_querable.tsv.gz")
    if len(annotator.non_queryable.index > 0):
        logger.info("-> Writing non-queryable sequences to disk")
        with open(output_non_queryable, "w") as file_handle:
            annotator.non_queryable.to_csv(file_handle, sep="\t", index=False)

    for this_method, this_df in result_dict.items():
        logger.info(f"-> Annotating query results of method: {this_method}")
        results = annotator.annotate_cts(this_df)
        sample_hits = annotator.annotate_sequences(results)
        healthy_sample_rate, tumor_sample_rate = annotator.annotate_sample_rate2(
            results
        )

        # Write index hits to output
        output_annotated = pathlib.Path(
            args.output + f"_annotated_{this_method}.tsv.gz"
        )
        output_healthy_rate = pathlib.Path(
            args.output + f"_healthy_sample_rate_{this_method}.tsv.gz"
        )
        output_tumor_rate = pathlib.Path(
            args.output + f"_tumor_sample_rate_{this_method}.tsv.gz"
        )
        with open(output_annotated, "w") as file_handle:
            sample_hits.to_csv(file_handle, sep="\t", index=False)
        # Write sample rate to output
        with open(output_healthy_rate, "w") as file_handle:
            healthy_sample_rate.to_csv(file_handle, sep="\t", index=False)
        with open(output_tumor_rate, "w") as file_handle:
            tumor_sample_rate.to_csv(file_handle, sep="\t", index=False)
