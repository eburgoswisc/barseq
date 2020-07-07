#!/usr/bin/env python3

"""
Test module for testing functions in barseq.

"""

# Module import
from barseq.exceptions import DuplicateBarcodeError
from barseq.read_barcodes import read_barcodes
from barseq.utils import *

__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'


def test_read_tab_delimited_barcodes():
    barcode_file = 'barseq/tests/data/input/barcodes-samples.csv'
    output_d = {'ATGAAGACTGTTGCCGTA': {'Type': 'WT', 'count': 0, 'gene': 'bar1'},
                'CACGACGCCCTCCGCGGA': {'Type': 'WT', 'count': 0, 'gene': 'bar2'},
                'ACTATTACGCAAAATAAT': {'Type': 'WT', 'count': 0, 'gene': 'bar3'},
                'ATGGAAGATATTATTATT': {'Type': 'WT', 'count': 0, 'gene': 'bar4'},
                'CCTCTCCAACCGGGTCTG': {'Type': 'MUT', 'count': 0, 'gene': 'bar5'},
                'CCCGGTCGCCTAGCCCCG': {'Type': 'MUT', 'count': 0, 'gene': 'bar6'},
                'GGCCCCCCGCCCGTCCCC': {'Type': 'MUT', 'count': 0, 'gene': 'bar7'},
                'GGATCACTGCTAGCGTAT': {'Type': 'MUT', 'count': 0, 'gene': 'bar8'},
                'CCTGCAGCAGCGGCCCGC': {'Type': 'MUT', 'count': 0, 'gene': 'bar9'},
                'ACACATGCAGACATAGAG': {'Type': 'MUT', 'count': 0, 'gene': 'bar10'},
                'CGCGCCATCCGCCGCCCA': {'Type': 'MUT', 'count': 0, 'gene': 'bar11'},
                'AATATTCAGATGGGACGT': {'Type': 'MUT', 'count': 0, 'gene': 'bar12'},
                '_other': {'Type': None, "gene": "_other", "count": 0}
                }
    assert output_d == read_barcodes(barcode_file)


def test_format_filename():
    assert "weird_file-name_" == format_filename("weird:_file-%$name_='][")
    assert "weird_file-name_" == format_filename("      weird:_file-%$name_='][")
    assert "weird_file-name_" == format_filename("weird:_file-%$name_='][     ")


def test_duplicate_barcode_error():
    try:
        barcode_duplicate_file = Path('barseq/tests/data/input/barcode-samples-duplicate.csv')
        read_barcodes(barcodes_file=barcode_duplicate_file)
        # Should not reach this...
        assert False
    except Exception as e:
        assert isinstance(e, DuplicateBarcodeError)


def test_uncompress_fastq_gz():
    compressed_files = os.listdir('barseq/tests/data/input/sequences/compressed')
    for f in compressed_files:
        uncompress_fastq_gz(Path(f'barseq/tests/data/input/sequences/compressed/{f}'))


def test_reverse_barcodes():
    barcode1 = 'ATGAAGACTGTTGCCGTA'
    barcode2 = 'CACGACGCCCTCCGCGGA'
    assert reverse_barcodes(barcode1) == 'TACGGCAACAGTCTTCAT'
    assert reverse_barcodes(barcode2) == 'TCCGCGGAGGGCGTCGTG'


def test_complement_dna():
    assert complement_dna('ATCC') == 'TAGG'
    assert complement_dna('TGCGTGCGTGTGAAAAATACATAT') == 'ACGCACGCACACTTTTTATGTATA'


def test_read_sample_map():
    sample_map = {
        'data1': {'New_name': 'Input-1', 'Type': 'Input'},
        'data2': {'New_name': 'Output-1', 'Type': 'Output'},
        'data3': {'New_name': 'Output-2', 'Type': 'Output'},
    }
    assert read_sample_map('barseq/tests/data/input/sample_map_test.csv') == sample_map