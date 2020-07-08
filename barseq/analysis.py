#!/usr/bin/env python3

"""
Script for performing mathematical operations on the results count table

"""

import logging
import numpy as np
import pandas as pd
from pathlib import Path

# Module imports
from barseq.read_barcodes import read_barcodes
from barseq.utils import *

__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'

# Get logger
logger = BarSeqLogger.logger


class BarSeqAnalysis(object):
    def __init__(self, cts_table: str, sample_map: str, barcodes: str, output_dir: str):

        self.count_file = Path(cts_table)
        self.output_dir = Path(output_dir)

        # Drop Barcode column and _other row
        self.raw_counts_df = pd.read_csv(cts_table, index_col=[1])
        self.raw_counts_df.drop(labels=['_other'], inplace=True)
        self.raw_counts_df.drop(columns=['Barcodes'], inplace=True)

        self.barcode_map = barcodes
        self.sample_map = sample_map

        return

    def run(self) -> pd.DataFrame:
        # Find ids for Wild type barcodes and input samples
        self.n_samples = len(self.raw_counts_df.columns)
        self.WT_id = [d['gene'] for d in self.barcode_map.values() if d['Type'] == 'WT']
        self.input_id = [d['Name'] for d in self.sample_map.values() if d['Type'] == 'Input'][0]

        # Add pseudo count of 1 to zeros
        self.raw_counts_df.replace(to_replace=0, value=1, inplace=True)
        self.sum_samples = self.raw_counts_df.sum(axis=0)

        logger.info("Normalize raw counts to input barcode counts.")
        self.normalize_counts_to_input()

        logger.info("Calculate relative frequency for each barcode.")
        self.calculate_rf()

        logger.info("Calculate the competitive index.")
        self.calculate_competitive_index()

        # Save final df into results file
        self.final_df.to_csv(self.output_dir.joinpath('barcode-counts-analysis.txt'), sep='\t')

        return

    def normalize_counts_to_input(self):
        self.rf_input = dict()
        for strain, counts in self.raw_counts_df[self.input_id].iteritems():
            self.rf_input[strain] = counts / self.sum_samples[self.input_id]

        # Loop through each sample column, return series with strain as index
        self.normalized_counts = pd.DataFrame(index=self.raw_counts_df.index)
        for col in self.raw_counts_df:
            normalized_sample = dict()
            # Skip gene column

            for strain, count in self.raw_counts_df[col].iteritems():
                # Place in dictionary
                normalized_sample[strain] = (count / self.rf_input[strain]) / self.n_samples
            self.normalized_counts = self.normalized_counts.join(pd.Series(normalized_sample,
                                                                           index=self.raw_counts_df.index,
                                                                           name=col))

    def calculate_rf(self):
        # Loop through each col, get relative frequency
        self.rf_normalized_output = pd.DataFrame(index=self.raw_counts_df.index)
        for col in self.normalized_counts:
            # Get sum of all counts in sample
            sum_counts_sample = self.normalized_counts[col].sum()
            rf_sample = dict()
            for strain, count in self.normalized_counts[col].iteritems():
                rf_sample[strain] = count / sum_counts_sample
            self.rf_normalized_output = self.rf_normalized_output.join(pd.Series(rf_sample, name=col))

    def calculate_competitive_index(self):
        # Grab only WT rows and take mean
        self.rf_norm_mean_WT = {col: wt_rf for col, wt_rf in self.rf_normalized_output.loc[self.WT_id, :].mean().iteritems()}
        self.final_df = pd.DataFrame(index=self.rf_normalized_output.index)
        wt_input = self.rf_norm_mean_WT[self.input_id]
        for sample in self.rf_normalized_output:
            wt_sample = self.rf_norm_mean_WT[sample]
            sample_ci = dict()
            for strain, rf in self.rf_normalized_output[sample].iteritems():
                numerator_ci = rf / wt_sample
                denominator_ci = self.rf_normalized_output.loc[strain, self.input_id] / wt_input
                sample_ci[strain] = np.log10(numerator_ci / denominator_ci)
            self.final_df = self.final_df.join(pd.Series(sample_ci, name=sample))

