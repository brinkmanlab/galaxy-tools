#!/usr/bin/env python
import sys
from biopython_convert import get_args, convert
from Bio.SeqIO.InsdcIO import _InsdcWriter

# Quote anticodon qualifiers
_InsdcWriter.FTQUAL_NO_QUOTE = tuple(v for v in _InsdcWriter.FTQUAL_NO_QUOTE if v not in ['anticodon', 'transl_except'])

if __name__ == '__main__':
    convert(*get_args(sys.argv[1:]))

