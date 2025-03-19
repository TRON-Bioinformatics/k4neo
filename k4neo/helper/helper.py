#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class FastaHandler:
    """
    Generic class to handle FASTA operations
    """
    @staticmethod
    def write_fasta(entries: pd.DataFrame, fasta_file: pathlib.Path):
        """
        Write query context sequences into fasta format
        """
        assert all([x in entries.columns for x in ['query_cts_id', 'query_sequence']]), "Columns are missing"
        sequences = entries[['query_cts_id', 'query_sequence']].drop_duplicates()
        fasta_entries = []
        with open(fasta_file, 'w') as file_handle:
            for row in sequences.itertuples(index=False):
                record = SeqRecord(
                    Seq(row.query_sequence),
                    id = row.query_cts_id,
                    name = row.query_cts_id,
                    description='')
                fasta_entries.append(record)
            SeqIO.write(fasta_entries, file_handle, "fasta")

class InputValidation:
    @staticmethod
    def columns_exist:
        pass

class SequenceOperation:
    """
    Generic methods to shorten and filter context sequences
    """
    @staticmethod
    def subset_cts (cts_seq: str, pos: int = None, length: int = None) -> str:
        """Generate query sequence

        Generate sequence of interest to query in kmer index e.g. the fusion breakpoint, splice junction. If length
        is null the original sequence is returned for search in index e.g. a full length transcript or SNP sequence.

        Args:
            cts_seq (str): The full-length context sequence.
            pos (int, optional): Position in context sequence to extract query seqeunce from. Defaults to None.
            length (int, optional): Length of query sequence. +/- length/2 around pos. Defaults to None.

        Returns:
            str: Context sequence string
        """
        if pd.isnull(pos) or pd.isnull(length):
            return cts_seq

        start = max(0, pos - round(length / 2))
        stop = min(start + length, len(cts_seq) - 1)
        sequence_of_interest = cts_seq[start:stop]
        return sequence_of_interest

    def filter_short_sequence
        pass
        
    def calculate_cts_id(cts_seq):
        return xxhash.xxh64(cts_seq).hexdigest()
    