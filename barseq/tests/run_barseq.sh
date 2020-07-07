#!/usr/bin/env bash

## Meant to be run from project root directory
# Setup stuff...
source barseq/tests/setup.sh
echo '****** INITIATE TESTING OF BarSeq PIPELINE ******'

echo
# CREATE A DUMP DIR FOR TESTING
CREATE_DUMP

# Run Barseq main pipeline...
which barseq
# PARAMS for barseq test
INPUT="../barseq/tests/data/input/sequences/no-index"
BARCODES="../barseq/tests/data/input/barcodes-samples.csv"
SAMPLEMAP="../barseq/tests/data/input/sample_map_test.csv"
pwd
if barseq -i $INPUT -b $BARCODES -e "test" --sample-map $SAMPLEMAP; then
:
else
:
fi
# CLEAN UP DUMP
CLEAN_UP