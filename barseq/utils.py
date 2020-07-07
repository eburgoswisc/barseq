#!/usr/bin/env python3

"""
Script that provides helper functions for package.

"""

import os
import re
import sys
import csv
import logging
import subprocess
import pandas as pd
from pathlib import Path

__author__ = "Emanuel Burgos"
__email__ = "eburgos@wisc.edu"

# Stdout logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(module)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M",
)
# Get logger
logger = logging.getLogger("barseq")

class Cd:
    """ Context manager for moving between directories. """
    def __init__(self, new_path):
        self.new_path = new_path

    def __enter__(self):
        self.old_path = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.old_path)


def read_barcodes(barcodes_file: Path) -> dict:
    """
    Read in barcodes from file

    :param barcodes_file: path to csv file with barcodes and gene name
    :return barcode_dict:

    barcode_dict = {
        barcode_1 : {"gene": Gene_1, "count": 0}
        barcode_2 : {"gene": Gene_2, "count": 0}
    }
    """
    # Store barcodes
    barcode_dict = dict()
    with open(barcodes_file, "r") as csv_barcode:
        # Skip Header
        next(csv_barcode)
        for line in csv.reader(csv_barcode):
            # Ignore comments
            if not line[0].startswith("#"):
                gene = line[1]
                barcode = line[0].upper()
                # Check for duplicate barcode
                if barcode not in barcode_dict:
                    barcode_dict[barcode] = {"gene": gene, "count": 0}
                else:
                    logger.error(f"Barcode {barcode} already in dictionary.")
                    raise IOError(f"Duplicate error: {barcode} already in dictionary")
    # Add _other for barcode that do not match
    barcode_dict["_other"] = {"gene": "_other", "count": 0}
    return barcode_dict


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


def format_filename(name: str) -> str:
    """
    Converts input into a valid name for file and recording results

    :param name: Input name
    :return output: Formatted name
    """
    # Strip extension
    name = "".join(name.split(".")[0])
    # Taken from Pyinseq software
    return re.sub(r"(?u)[^-\w]", "", name.strip().replace(" ", "_"))


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


if __name__ == '__main__':
    pass
