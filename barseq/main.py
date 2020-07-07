#!/usr/bin/env python3

"""
Main pipeline for barseq software

"""

import logging
from copy import deepcopy
from pathlib import Path

# Module import
from barseq.process_reads import count_barcodes
from barseq.read_barcodes import read_barcodes
from barseq.utils import *

__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'

# Get logger
logger = logging.getLogger()

class BarSeqRun():
    def __init__(self, input: str, barcodes: str, output_dir: str, experiment: str, sample_map: str=None, force_it: bool=False):

        self.program = 'BarSeq'
        self.experiment = experiment
        self.sequence_dir = Path(input)
        self.barcodes_path = Path(barcodes)
        self.output_dir = Path(output_dir).joinpath(f'results/{self.experiment}')
        self.logfile = self.output_dir.joinpath(f'{self.experiment}.log')
        self.counts_results_path = self.output_dir.joinpath('barcode_counts_table.csv')

        # Create the output directory
        create_directories(self.output_dir, force_it)

        # Add logfile handler to logger
        #BarSeqLogger.add_file_handler(self.logfile)

        if sample_map:
            self.sample_map_path = Path(sample_map)



    def get_sequence_files(self):
        """ Gets sequence files relative path and returns list with Path objects
        :return: list with relative path objects from cwd for sequences
        """
        # Relative to cwd Path objects (cwd/sequences/filename.fastq[.gz])
        return sorted([Path(n) for n in self.sequence_dir.glob('*fastq*')])

    def rename_sample_columns(self):
        new_barcode_count = dict()
        for sample, d in self.sample_map.items():
            if sample in self.barcode_count_samples_dict.keys():
                new_barcode_count[(d['Name'], d['Type'])] = self.barcode_count_samples_dict[sample]
        self.barcode_count_samples_dict = new_barcode_count
        return

    def run(self) -> None:
        """
        Main pipeline for analyzing barseq data.
        """
        logger.info("***** Starting barseq *****")

        logger.info(f"Collecting sequence files from {self.sequence_dir}")
        self.sequences_list = self.get_sequence_files()

        logger.info(f"Reading in barcodes from {self.barcodes_path.name}")
        self.barcode_dict = read_barcodes(self.barcodes_path)
        self.barcode_count_samples_dict = dict()

        if self.sample_map_path:
            self.sample_map = read_sample_map(self.sample_map_path)

        # Process each sequencing file
        for seq_file in self.sequences_list:
            # Get sample name without suffixes
            suffixes = ''.join(seq_file.suffixes)
            sample = seq_file.name.replace(suffixes, '')

            if 'R2' in seq_file.name.split('_'):
                logger.warning(f'Skipping {sample} since software does not support paired-end reads yet.')
                continue

            logger.info(f"Mapping barcoded reads in {sample}")
            self.barcode_count_samples_dict[sample] = deepcopy(self.barcode_dict)
            # Seq file needs to be just name
            count_barcodes(seq_file, self.barcode_count_samples_dict[sample])

        # Change names
        if self.sample_map_path:
            logger.info(f"Renaming columns using sample mapping found in {self.sample_map_path.name}")
            self.rename_sample_columns()

        # Write to output
        logger.info(f"Storing results in {self.output_dir}")
        logger.info(f"Writing count data into {self.counts_results_path}")
        write_output(self.barcode_count_samples_dict, self.barcode_dict, self.counts_results_path)

        # Confirm completion of barseq
        logger.info("***** barseq is complete! *****")
        return


if __name__ == "__main__":
    pass
