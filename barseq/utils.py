#!/usr/bin/env python3

"""
Module that provides helper functions for barseq. Basically for methods that have no home...

"""

import os
import re
import sys
import csv
import logging
import subprocess
import pandas as pd
from pathlib import Path

# Module imports
from barseq.exceptions import SampleMapError, ExperimentDirectoryExistsError

__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'

DNA_COMPLEMENT_BASES = {'A': 't', 'G': 'c', 'T': 'a', 'C': 'g'}


class Cd:
    """ Context manager for moving between directories. """

    def __init__(self, new_path):
        self.new_path = new_path

    def __enter__(self):
        self.old_path = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.old_path)


class BarSeqLogger:
    # Get logger
    logger = logging.getLogger('Barseq')
    logger.setLevel(logging.INFO)
    # Set up stream handler for stdout
    FORMATTER = logging.Formatter(datefmt="%H:%M:%S",
                                  fmt="[%(asctime)s] %(levelname)s - %(module)s - %(message)s")

    # Set up console logger
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(FORMATTER)

    logger.addHandler(sh)

    logger.debug('Logger initiated.')

    @staticmethod
    def add_file_handler(logfile_path: Path):
        if not logfile_path.is_file():
            logfile_path.touch()
        fh = logging.FileHandler(logfile_path, mode='w')
        fh.setFormatter(BarSeqLogger.FORMATTER)
        BarSeqLogger.logger.addHandler(fh)


def format_filename(name: str) -> str:
    """
    TAKEN FROM PYINSEQ
    Converts input into a valid name for file and recording results
    :param name: Input name
    :return output: Formatted name
    """
    # Strip extension
    name = "".join(name.split(".")[0])
    return re.sub(r"(?u)[^-\w]", "", name.strip().replace(" ", "_"))


def write_output(sample_dict: dict, barcode_dict: dict, output_path: Path) -> None:
    """
    Convert results file into a pandas dataframe with following
    structure.

    :param sample_dict: Dictionary with samples and barcode counts
    :param barcode_dict: Dictionary with barcode and gene names
    :param output_path: Relative path for new csv file
    :return None

    |    Gene    | Barcode |   Sample 1  |  Sample 2   | ... |
    |------------|---------|-------------|-------------|-----|
    |   Gene 1   | ATCGCGT |     500     |     10      | ... |

    """
    # Barcode index wil be used, gene names will be treated as additional column
    barcode_index = [barcode for barcode in barcode_dict.keys()]
    gene_d = {barcode: d['gene'] for barcode, d in barcode_dict.items()}

    df = pd.DataFrame(index=barcode_index).join(pd.DataFrame.from_dict(gene_d, orient='index', columns=['Gene']))
    df.index.name = 'Barcodes'
    for (sample, t), count_d in sample_dict.items():
        sample_count_dict = dict()
        for barcode, bar_count_d in count_d.items():
            sample_count_dict[barcode] = bar_count_d['count']
        df = df.join(pd.DataFrame.from_dict(data=sample_count_dict, orient='index', columns=[sample]))
    # Write to output
    df.to_csv(output_path)
    return


def create_directories(dir_path: Path, force_it: bool):
    """
    Create barseq directories to store results
    :param dir_path: path to dir
    :param force_it: bool that forces program to overwrite dir
    :return:
    """
    if dir_path.is_dir() and not force_it:
        raise ExperimentDirectoryExistsError(dir_path)
    dir_path.mkdir(parents=True, exist_ok=True)
    return


def uncompress_fastq_gz(filename: Path) -> None:
    """
    Uncrompresses fastq.gz files using screed module
    :param filename: path to fastq.gz file
    :return:
    """
    # Check if gz compressed
    if filename.suffix in ["gz"]:
        # Build pipe for uncompressed file
        subprocess.Popen(("gzip", "-dk", filename))
    return


def complement_dna(sequence: str) -> str:
    """
    Gives complement sequence of dna string input (5`-> 3`)
    :param sequence: Sequence to complement
    :return: complement sequence
    """
    for bp, c in DNA_COMPLEMENT_BASES.items():
        sequence = sequence.replace(bp, c)
    return sequence.upper()


def reverse_barcodes(barcode) -> str:
    """ Given barcode, returns the reverse"""
    complement_barcode = complement_dna(barcode)
    return complement_barcode[::-1]


def convert_args_to_paths(args) -> dict:
    """
    Converts argparse.Namespace object into a dict of {arg: Path(value)}
    :param args:
    :return:
    """
    return {arg: Path(value) for arg, value in vars(args).items() if value}


def read_sample_map(filename: str) -> dict():
    """
    Reads sample map file into a dictionary to rename sample names in barseq run.
    :param filename:
    :return sample_map: dict {Old_name: {new_name: Value, Type: Value}}
    """
    sample_map = dict()
    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter=find_del(filename))
        for line in reader:
            sample_map[line['Sample']] = {'Name': line['Name'], 'Type': line['Type']}
    return sample_map


def find_del(filename: str) -> str:
    """
    Helper for finding delimiter in table file.
    :param filename:
    :return:
    """
    with open(filename, 'r') as f:
        line = repr(f.readline())
        if ',' in line:
            return ','
        if '\\t' in line:
            return '\\t'

def total_reads(seq_file: Path) -> int:
    return
if __name__ == '__main__':
    pass
