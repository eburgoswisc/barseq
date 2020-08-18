#!/usr/bin/env python3

"""
Main pipeline for barseq software

"""

from copy import deepcopy

# Module import
from barseq.process_reads import count_barcodes
from barseq.read_barcodes import read_barcodes
from barseq.utils import *

__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'

# Get logger
logger = BarSeqLogger.logger

class BarSeqRun():
    def __init__(self, input: str, barcodes: str, output_dir: str, experiment: str, sample_map: str=None, force_it: bool=False):

        logger.info("***** Starting barseq *****")

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
        BarSeqLogger.add_file_handler(self.logfile)

        if sample_map:
            self.sample_map_path = Path(sample_map)
            self.sample_map = read_sample_map(self.sample_map_path)

        logger.info(f"Reading in barcodes from {self.barcodes_path.name}")
        self.barcode_dict = read_barcodes(self.barcodes_path)
        self.barcode_count_samples_dict = dict()

        logger.info(f"Collecting sequence files from {self.sequence_dir}")
        self.get_sequence_files()

    def get_sequence_files(self):
        """ Gets sequence files relative path and returns list with Path objects """
        # Check that sequence dir exists
        if not self.sequence_dir.exists():
            raise NotADirectoryError("Sequence directory does not exist.")
        # Relative to cwd Path objects (cwd/sequences/filename.fastq[.gz])
        seq_files = list(self.sequence_dir.glob('*fastq*'))
        for sample in self.sample_map:
            f_path = [f for f in seq_files if sample in f.name][0]
            self.sample_map[sample]['filepath'] = f_path
        return

    def rename_sample_columns(self):
        """ Rename samples in map using names from user provided file """
        new_barcode_count = dict()
        for sample, d in self.sample_map.items():
            if sample in self.barcode_count_samples_dict.keys():
                new_barcode_count[(d['Name'], d['Type'])] = self.barcode_count_samples_dict[sample]
        self.barcode_count_samples_dict = new_barcode_count
        return

    def run(self) -> None:
        """ Main pipeline for analyzing barseq data """
        # Process each sequencing file
        for sample, metadata in self.sample_map.items():
            # Get sample name without suffixes
            seq_file = metadata['filepath']
            suffixes = ''.join(seq_file.suffixes)
            sample = seq_file.name.replace(suffixes, '')

            if 'R2' in seq_file.name.split('_'):
                logger.warning(f'Skipping {sample} since software does not support paired-end reads yet.')
                continue
            # Need a newline
            print()
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
