#!/usr/bin/env python3

import csv
import logging
from pathlib import Path

# Module import
from barseq.exceptions import DuplicateBarcodeError

__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'

# Get logger
logger = logging.getLogger()

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
    with open(barcodes_file, "r") as f:
        for line in csv.DictReader(f):
            # Grab needed items
            gene = line['Gene']
            barcode = line['Barcode'].upper()
            strain_type = line['Type'].upper()
            # Check for duplicate barcode
            if barcode not in barcode_dict:
                barcode_dict[barcode] = {"gene": gene, "count": 0, 'Type': strain_type}
            else:
                logger.error(f"Barcode {barcode} already in dictionary.")
                raise DuplicateBarcodeError(f"Duplicate barcode error: {barcode} already in dictionary")
    # Add _other for barcode that do not match
    barcode_dict["_other"] = {"gene": "_other", "count": 0, "Type": None}
    return barcode_dict


if __name__ == '__main__':
    pass
