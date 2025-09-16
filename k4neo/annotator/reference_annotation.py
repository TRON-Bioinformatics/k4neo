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

        with open(manifest, "r") as file_handle:
            self.meta = json.load(file_handle)

        self.k = k or self.meta["kmer_size"]
        if self.k != self.meta["kmer_size"]:
            raise ValueError(f"K-mer size mismatch! Expected {self.meta['kmer_size']}, got {self.k}")
        
        self.canonical = False
        
        self.genome_index = os.path.join(self.meta["data_dir"], "genome.jf")
        self.transcriptome_index = os.path.join(self.meta["data_dir"], "transcriptome.jf")

    def _get_kmers(self, sequence: str) -> List[str]:
        """Extract k-mers from sequence

        Args:
            sequence (str): A nucleotide sequence

        Returns:
            List[str]: A list of len(sequence) - self.k + 1 k-mers 
        """
        return [sequence[i : i + self.k] for i in range(len(sequence) - self.k + 1)]

    @staticmethod
    def canonicalize(kmer):
        """Canonical representation of k-mer

        A canonical k-mer is defined as the lexicographic smaller string of
        a k-mer and it's reverse complement

        Args:
            kmer (str): The k-mer to bring into canonical form

        Returns:
            str: The canonical k-mer
        """
        reverse_complement = {'A':'T',
                              'T':'A',
                              'C':'G',
                              'G':'C'}
        kmer_rev = kmer[::-1]
        kmer_rev = [reverse_complement[char] for char in kmer_rev]
        kmer_rev = ''.join(kmer_rev)
        if kmer < kmer_rev:
            return kmer
        return kmer_rev


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
        kmers = self._get_kmers(sequence)
        # For genome we search the canoncial representation
        genome_counts = [JellyFishHelper.query_index(KmerUniquenessAnnotator.canonicalize(kmer), self.genome_index) for kmer in kmers]
        # For transcriptome annotation we require the strand of the k-mer to match as this matters for uniqueness of antisense transcripts etc.
        transcriptome_counts = [
            JellyFishHelper.query_index(kmer, self.transcriptome_index) for kmer in kmers
        ]

        gene_specific = sum(
            1 for g, t in zip(genome_counts, transcriptome_counts) if g <= 1 and t >= 1
        )
        transcript_specific = sum(
            1 for g, t in zip(genome_counts, transcriptome_counts) if g <= 1 and t == 1
        )
        total_specific = sum(
            1 for g, t in zip(genome_counts, transcriptome_counts) if g == 0 and t == 0
        )
        total_unspecific = sum(
            1 for g, t in zip(genome_counts, transcriptome_counts) if g >= 1 or t >= 1
        )

        total_kmers = len(kmers)
        if total_kmers == 0:
            return {
                "gene_specific_rate": 0.0,
                "transcript_specific_rate": 0.0,
                "total_specific_rate": 0.0,
                "total_unspecific_rate": 0.0
            }

        return {
            "gene_specific_rate": gene_specific / total_kmers,
            "transcript_specific_rate": transcript_specific / total_kmers,
            "total_specific_rate": total_specific / total_kmers,
            "total_unspecific_rate": total_unspecific / total_kmers,
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
            JellyFishHelper.generate_index(genome_fasta, genome_index, bf_size="3G", canonical=True, kmer_size=kmer_size)
        if not os.path.exists(transcriptome_index):
            logger.info("Creating transcriptomic reference index with jellyfish...")
            JellyFishHelper.generate_index(transcriptome_fasta, transcriptome_index, bf_size="100M", canonical=False, kmer_size=kmer_size)

        metadata = {
            "kmer_size": kmer_size,
            "genome_fasta": os.path.basename(genome_fasta),
            "transcriptome_fasta": os.path.basename(transcriptome_fasta),
            "data_dir": os.path.abspath(outdir)
        }

        with open(os.path.join(outdir, "metadata.json"), "w") as f:
            json.dump(metadata, f, indent=2)
