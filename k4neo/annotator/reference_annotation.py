#!/usr/bin/env python3

import argparse
import subprocess
import os
import json
from Bio import SeqIO
from typing import List, Dict
from k4neo.helper.helper import JellyFishHelper
from loguru import logger


class KmerUniquenessAnnotator:

    def __init__(self, manifest: str, k: int = None):

        with open(os.path.join(data_dir, "metadata.json")) as f:
            meta = json.load(f)

        self.k = k or meta["kmer"]
        if self.k != meta["kmer"]:
            raise ValueError(f"K-mer size mismatch! Expected {meta['kmer']}, got {self.k}")

        self.genome_index = os.path.join(meta["data_dir"], "genome.jf")
        self.transcriptome_index = os.path.join(meta["data_dir"], "transcriptome.jf")

    def _get_kmers(self, sequence: str) -> List[str]:
        return [sequence[i : i + self.k] for i in range(len(sequence) - self.k + 1)]

    def annotate_sequence(self, sequence: str) -> Dict[str, float]:
        kmers = self._get_kmers(sequence)
        genome_counts = [JellyFishHelper.query_index(kmer, self.genome_index) for kmer in kmers]
        transcriptome_counts = [
            JellyFishHelper.query_index(kmer, self.transcriptome_index) for kmer in kmers
        ]

        gene_specific = sum(
            1 for g, t in zip(genome_counts, transcriptome_counts) if g <= 1 and t >= 1
        )
        transcript_specific = sum(
            1 for g, t in zip(genome_counts, transcriptome_counts) if g <= 1 and t == 1
        )
        chimera_specific = sum(
            1 for g, t in zip(genome_counts, transcriptome_counts) if g == 0 and t == 0
        )

        total_kmers = len(kmers)
        if total_kmers == 0:
            return {
                "gene_specific_rate": 0.0,
                "transcript_specific_rate": 0.0,
                "chimera_specific_rate": 0.0,
            }

        return {
            "gene_specific_rate": gene_specific / total_kmers,
            "transcript_specific_rate": transcript_specific / total_kmers,
            "chimera_specific_rate": chimera_specific / total_kmers,
        }

    def annotate_fasta(self, fasta_file: str) -> Dict[str, Dict[str, float]]:
        results = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            rates = self.annotate_sequence(str(record.seq))
            results[record.id] = rates
        return results


def prepare_data(genome_fasta: str, transcriptome_fasta: str, k: int, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    genome_index = os.path.join(outdir, "genome.jf")
    transcriptome_index = os.path.join(outdir, "transcriptome.jf")
    db_path = os.path.join(outdir, "annotation.db")

    logger.info("🔄 Creating genomic reference index with jellyfish...")
    ShellExec.execute_cmd(
        ["jellyfish", "count", "-m", str(k), "-C", "-s", "3G", "-o", genome_index, genome_fasta],
        check=True,
    )

    logger.info("🔄 Creating transcriptomic reference index with jellyfish...")
    ShellExec.execute_cmd(
        [
            "jellyfish",
            "count",
            "-m",
            str(k),
            "-s",
            "100M",
            "-o",
            transcriptome_index,
            transcriptome_fasta,
        ],
        check=True,
    )

    metadata = {
        "kmer": k,
        "genome_fasta": os.path.basename(genome_fasta),
        "transcriptome_fasta": os.path.basename(transcriptome_fasta),
        "data_dir": os.path.abspath(outdir)
    }

    with open(os.path.join(outdir, "metadata.json"), "w") as f:
        json.dump(metadata, f, indent=2)


def run_annotation(args):
    annotator = KmerAnnotator(manifest=args.manifest, k=args.kmer)
    results = annotator.annotate_fasta(args.fasta)

    if args.format == "json":
        output = json.dumps(results, indent=2)
    else:
        lines = ["ID\tgene_specific_rate\ttranscript_specific_rate"]
        for seq_id, rates in results.items():
            lines.append(
                f"{seq_id}\t{rates['gene_specific_rate']:.4f}\t{rates['transcript_specific_rate']:.4f}"
            )
        output = "\n".join(lines)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
    else:
        print(output)


def main():
    parser = argparse.ArgumentParser(description="Annotator CLI Tool")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # PREPARE
    prepare = subparsers.add_parser("prepare", help="Erzeuge Jellyfish-Indizes und GTF-Datenbank")
    prepare.add_argument("--genome", required=True, help="Genom FASTA-Datei")
    prepare.add_argument("--transcriptome", required=True, help="Transkriptom FASTA-Datei")
    prepare.add_argument("-k", "--kmer", type=int, default=31, help="k-mer Größe (default: 31)")
    prepare.add_argument(
        "-o", "--outdir", required=True, help="Zielverzeichnis für vorbereitete Dateien"
    )

    # ANNOTATE
    annotate = subparsers.add_parser(
        "annotate", help="Annotiere Sequenzen aus FASTA mit vorbereiteten Daten"
    )
    annotate.add_argument(
        "-m", "--manifest", required=True, help="Manifest"
    )
    annotate.add_argument("-f", "--fasta", required=True, help="FASTA-Datei mit Sequenzen")
    annotate.add_argument(
        "-k", "--kmer", type=int, required=False, help="k-mer Größe (muss zur Vorbereitung passen)"
    )
    annotate.add_argument("-o", "--output", help="Ausgabedatei")
    annotate.add_argument("--format", choices=["json", "tsv"], default="json", help="Ausgabeformat")

    args = parser.parse_args()

    if args.command == "prepare":
        prepare_data(args.genome, args.transcriptome, args.kmer, args.outdir)
    elif args.command == "annotate":
        run_annotation(args)


if __name__ == "__main__":
    main()
