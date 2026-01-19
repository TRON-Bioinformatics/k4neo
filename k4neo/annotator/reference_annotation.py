#!/usr/bin/env python3

import argparse
import subprocess
import os
import json
from Bio import SeqIO
from typing import List, Dict
from k4neo.helper.helper import JellyFishHelper, SequenceOperation
from loguru import logger
from probables import CountingBloomFilter


class KmerUniquenessAnnotator:

    def __init__(self, manifest: str, k: int = None):

        with open(manifest, "r") as file_handle:
            self.meta = json.load(file_handle)

        self.k = k or self.meta["kmer_size"]
        if self.k != self.meta["kmer_size"]:
            raise ValueError(
                f"K-mer size mismatch! Expected {self.meta['kmer_size']}, got {self.k}"
            )

        self.canonical = False

        self.genome_index = os.path.join(self.meta["data_dir"], "genome.jf")
        self.transcriptome_index = os.path.join(self.meta["data_dir"], "transcriptome.jf")

        self.gene_index = CountingBloomFilter(
            filepath=os.path.join(self.meta["data_dir"], "gene.bf")
        )
        self.transcript_index = CountingBloomFilter(
            filepath=os.path.join(self.meta["data_dir"], "transcriptome.bf")
        )

    def annotate_sequence(self, sequence: str) -> Dict[str, float]:
        """Uniqueness annotation of sequence

        Returns three metrics describing different levels of uniqueness.
        K-mer uniqueness is calculated with respect to genome and transcriptome.

        Genome specific rate: Rate of k-mers that are specific to one position in the genome but can occur several times in transcriptome
        Transcript specific rate: Rate of k-mers that are specific to one transcript at one genomic locus.
        Chimera specific: Rate of unique k-mers.


        Args:
            sequence (str): A nucleotide sequence

        Returns:
            Dict[str, float]: A dict with uniqueness rates.
        """
        kmers = SequenceOperation.get_kmers(sequence, self.k)
        total_kmers = len(kmers)

        # For genome we search the canoncial representation
        genome_counts = [
            JellyFishHelper.query_index(SequenceOperation.canonicalize(kmer), self.genome_index)
            for kmer in kmers
        ]
        # For transcriptome annotation we require the strand of the k-mer to match as this matters for uniqueness of antisense transcripts etc.
        transcriptome_counts = [
            JellyFishHelper.query_index(kmer, self.transcriptome_index) for kmer in kmers
        ]
        # Check if k-mer is unique to one genic location or multiple genes have a transcript containing this k-mer
        gene_counts = [self.gene_index.check(this_kmer) for this_kmer in kmers]
        transcript_counts = [self.transcript_index.check(this_kmer) for this_kmer in kmers]

        # gene_specific = sum(
        #    1 for g, t in zip(genome_counts, transcriptome_counts) if g <= 1 and t >= 1
        # )
        # transcript_specific = sum(
        #    1 for g, t in zip(genome_counts, transcriptome_counts) if g <= 1 and t == 1
        # )
        # Unique k-mers of the sequence
        cts_unique = sum(
            1 for g, t in zip(genome_counts, transcriptome_counts) if g == 0 and t == 0
        )
        # k-mers occuring in the reference
        cts_ref = total_kmers - cts_unique
        # k-mers in reference from a single genic loci
        cts_ref_single_gene_locus_rate = sum(1 for count in gene_counts if count == 1)
        # k-mers in reference from multiple genic loci
        cts_ref_multi_gene_locus_rate = total_kmers - cts_ref_single_gene_locus_rate
        # k-mers specific to a single transcript from the reference
        cts_ref_single_transcript_rate = sum(
            1
            for gx_count, tx_count in zip(gene_counts, transcript_counts)
            if gx_count == 1 and tx_count == 1
        )
        # k-mers from multiple transcripts
        cts_ref_multi_transcript_rate = sum(
            1
            for gx_count, tx_count in zip(gene_counts, transcript_counts)
            if gx_count == 1 and tx_count > 1
        )

        if total_kmers == 0:
            return {
                "cts_unique_rate": 0.0,
                "cts_ref_rate": 0.0,
                "cts_ref_single_gene_locus_rate": 0.0,
                "cts_ref_multi_gene_locus_rate": 0.0,
                "cts_ref_single_transcript_rate": 0.0,
                "cts_ref_multi_transcript_rate": 0.0,
            }

        return {
            "cts_unique_rate": cts_unique / total_kmers,
            "cts_ref_rate": cts_ref / total_kmers,
            "cts_ref_single_gene_locus_rate": cts_ref_single_gene_locus_rate / total_kmers,
            "cts_ref_multi_gene_locus_rate": cts_ref_multi_gene_locus_rate / total_kmers,
            "cts_ref_single_transcript_rate": cts_ref_single_transcript_rate / total_kmers,
            "cts_ref_multi_transcript_rate": cts_ref_multi_transcript_rate / total_kmers,
        }

    def annotate_fasta(self, fasta_file: str) -> Dict[str, Dict[str, float]]:
        """Annotate a fasta file

        Args:
            fasta_file (str): Path to nucleotide fasta file

        Returns:
            Dict[str, Dict[str, float]]: A nested dict mapping sequences to uniqueness rates.
        """
        results = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            rates = self.annotate_sequence(str(record.seq))
            results[record.id] = rates
        return results


class ReferenceIndexer:

    @staticmethod
    def prepare_data(genome_fasta: str, transcriptome_fasta: str, kmer_size: int, outdir: str):
        os.makedirs(outdir, exist_ok=True)
        genome_index = os.path.join(outdir, "genome.jf")
        transcriptome_index = os.path.join(outdir, "transcriptome.jf")

        if not os.path.exists(genome_index):
            logger.info("Creating genomic reference index with jellyfish...")
            JellyFishHelper.generate_index(
                genome_fasta, genome_index, bf_size="3G", canonical=True, kmer_size=kmer_size
            )
        if not os.path.exists(transcriptome_index):
            logger.info("Creating transcriptomic reference index with jellyfish...")
            JellyFishHelper.generate_index(
                transcriptome_fasta,
                transcriptome_index,
                bf_size="100M",
                canonical=False,
                kmer_size=kmer_size,
            )

        metadata = {
            "kmer_size": kmer_size,
            "genome_fasta": os.path.basename(genome_fasta),
            "transcriptome_fasta": os.path.basename(transcriptome_fasta),
            "data_dir": os.path.abspath(outdir),
        }

        with open(os.path.join(outdir, "metadata.json"), "w") as f:
            json.dump(metadata, f, indent=2)
