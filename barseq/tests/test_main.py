#!/usr/bin/env python3

"""
Test module for testing main barseq pipeline

"""

import filecmp
import subprocess

# Module import
from barseq.tests.test_setup import temp_dir

__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'



def test_barseq(temp_dir):
    # Create testing paths
    input_sequence = temp_dir(['data', 'input', 'sequences', 'no-index'])
    input_barcodes = temp_dir(['data', 'input', 'barcodes-samples.csv'])
    input_sample_map = temp_dir(['data', 'input', 'sample_map_test.csv'])
    expected_output = temp_dir(['data', 'output'])
    experiment = 'barseq_pytest'
    test_output = temp_dir(['results', experiment])

    # Call main with test data
    subprocess.call(['barseq',
                     '-i', input_sequence,
                     '-b', input_barcodes,
                     '-e', experiment,
                     '-s', input_sample_map,
                     '-o', temp_dir()],
                    cwd=temp_dir())
    # Check output
    assert filecmp.dircmp(expected_output, test_output)
    assert filecmp.cmp(expected_output.joinpath('barcode_counts_table.csv'),
                       test_output.joinpath('barcode_counts_table.csv'))

    # Check log files
    cmp_log_file = expected_output.joinpath('barseq_pytest_log.txt')
    test_log_file = test_output.joinpath(f'{experiment}_log.txt')

    with open(cmp_log_file) as cmp_log, open(test_log_file) as test_log:
        for cmp, test in zip(cmp_log, test_log):
            test_line = '-'.join(test.strip().split('-')[3:])
            cmp_line = '-'.join(cmp.strip().split('-')[3:])

            assert cmp_line == test_line
