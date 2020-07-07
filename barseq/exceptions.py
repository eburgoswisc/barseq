#!/usr/bin/env python3

"""

Custom exception classes for barseq software

"""


__author__ = 'Emanuel Burgos'
__email__ = 'eburgos@wisc.edu'


class BarseqException(Exception):
    """
    Custom exception class for BarSeq.
    """
    def __init__(self):
        super().__init__()


class DuplicateBarcodeError(BarseqException):
    def __init__(self, msg):
        super().__init__()
        self.msg = msg


class ExperimentDirectoryExistsError(Exception):
    def __init__(self, dir_path):
        self.msg = f" {dir_path} directory already exists. " \
                   f"Delete or rename {dir_path}, or provide new name for the experiment.."
        super().__init__(self.msg)


class SampleMapError(BarseqException):
    def __init__(self, msg):
        super().__init__()
        self.msg = msg

