#!/usr/bin/env python3

class FastaHandler:
    @staticmethod
    def write_fasta(entries: dict, fasta_file: str):
        with open(fasta_file, 'w') as file_handle:
            for header, sequence in entries.items():
                file_handle.write(f"{header}\n")
                file_handle.write(f"{sequence}\n")
                
    @staticmethod
    def read_fasta(fasta_file):
        pass




class Annotator:
    def __init__(self, working_dir: str, sequence_table) -> None:
        self.working_dir = working_dir
        self.sequence_table = self.read_context_seq(sequence_table)

    def read_context_seq(self):
        pass

    def read_vcf_file(self):
        pass

    def annotator(self):
        pass

    def call_pipeline(self):
        """
        Implement a call to a pipeline that searches the context sequences in the index
        """
        pass