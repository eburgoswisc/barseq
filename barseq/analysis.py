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
from barseq.utils import read_sample_map

__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'

# Get logger
logger = logging.getLogger()


class BarSeqAnalysis(object):
    def __init__(self, cts_table: str, sample_map: str, barcodes: str, output_dir: str, rename_sample=False,
                 read_bar=False):

        self.count_file = Path(cts_table)
        self.output_dir = Path(output_dir)

        # Drop Barcode column and _other row
        self.raw_counts_df = pd.read_csv(cts_table, index_col=[1])
        self.raw_counts_df.drop(labels=['_other'], inplace=True)
        self.raw_counts_df.drop(columns=['Barcodes'], inplace=True)

        if read_bar:
            self.barcode_map = read_barcodes(barcodes)
        else:
            self.barcode_map = barcodes

        if rename_sample:
            self.sample_map = read_sample_map(sample_map)
        else:
            self.sample_map = sample_map

        # Find identifiers
        self.n_samples = len(self.raw_counts_df.columns)
        self.WT_id = [d['gene'] for d in self.barcode_map.values() if d['Type'] == 'WT']

        self.input_id = []
        for values in self.sample_map.values():
            if values['Type'] == 'Input':
                self.input_id = values['Name']

        self.sum_sample_strains = None
        self.rf_input = {}
        self.rf_norm_mean_WT = {}

        return

    def run(self) -> pd.DataFrame:
        # Mapping tables
        self.sum_sample_strains = self.raw_counts_df.sum(axis=0)

        self.get_rf_input()

        self.normalize_counts()

        self.calculate_rf()

        self.average_wt()
        self.get_ci_values()

        # Save final df into results file
        self.final_df.to_csv(self.output_dir.joinpath('barcode-counts-analysis.txt'), sep='\t')
        return self.final_df

    def calculate_rf(self):
        # Loop through each col, get relative frequency
        self.rf_normalized_output = pd.DataFrame(index=self.raw_counts_df.index)
        for col in self.normalized_df:
            # Get sum of all counts in sample
            sum_counts_sample = self.normalized_df[col].sum()
            rf_sample = dict()
            for strain, count in self.normalized_df[col].iteritems():
                rf_sample[strain] = count / sum_counts_sample
            self.rf_normalized_output = self.rf_normalized_output.join(pd.Series(rf_sample, name=col))

    def get_rf_input(self):
        self.rf_input = dict()
        for strain, counts in self.raw_counts_df[self.input_id].iteritems():
            self.rf_input[strain] = counts / self.sum_sample_strains[self.input_id]
        return

    def normalize_counts(self):
        # Loop through each sample column, return series with strain as index
        self.normalized_df = pd.DataFrame(index=self.raw_counts_df.index)
        for col in self.raw_counts_df:
            normalized_sample = dict()
            # Skip gene column
            if col != 'Gene':
                for strain, count in self.raw_counts_df[col].iteritems():
                    # Place in dictionary
                    normalized_sample[strain] = (count / self.rf_input[strain]) / self.n_samples
                self.normalized_df = self.normalized_df.join(pd.Series(normalized_sample,
                                                                       index=self.raw_counts_df.index, name=col))

    def average_wt(self):
        # Grab only WT rows and take mean
        rf_norm_out_mean_WT = self.rf_normalized_output.loc[self.WT_id, :].mean()
        self.rf_norm_mean_WT = {col: wt_rf for col, wt_rf in rf_norm_out_mean_WT.iteritems()}

    def get_ci_values(self):
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

