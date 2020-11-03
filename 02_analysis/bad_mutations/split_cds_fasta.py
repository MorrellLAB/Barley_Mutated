#!/usr/bin/env python3
"""This script takes in a CDS fasta file and splits the sequences up
into multiple files (one sequence per file). This allows for parallelizing
when running BAD_Mutations.

Usage: ./split_cds_fasta.py [fasta_file] [out_dir]

Where:
1) [fasta_file] is the full filepath to the fasta file
2) [out_dir] is the full filepath to where we want to store the split fasta files.
"""

import sys
import os
from Bio import SeqIO

def split_by_sequence(fasta_records, out_dir):
    """Split fasta by one sequence per file."""
    for record in fasta_records:
        out_prefix = record.id
        # Strip trailing slash to prevent path related errors
        out_fp = os.path.normpath(os.path.expanduser(out_dir)) + "/" + out_prefix + ".fa"
        # Write one record to file
        SeqIO.write(record, out_fp, "fasta")


def main(fasta_fp, out_dir_fp):
    """Driver function."""
    # Load the fasta file
    print("Loading fasta file...")
    fasta_dat = list(SeqIO.parse(os.path.expanduser(fasta_fp), "fasta"))
    # Split fasta file
    print("Splitting fasta file...")
    split_by_sequence(fasta_dat, out_dir_fp)
    print("Done.")


main(sys.argv[1], sys.argv[2]) # Run the program
