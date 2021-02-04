#!/usr/bin/env python
import sys
from biopython_convert import get_args, convert
from Bio.SeqIO.InsdcIO import _InsdcWriter, EmblWriter
import Bio

# Quote anticodon qualifiers
_InsdcWriter.FTQUAL_NO_QUOTE = tuple(v for v in _InsdcWriter.FTQUAL_NO_QUOTE if v not in ['anticodon', 'transl_except'])

if Bio.__version__ == "1.78":
    # TODO monkeypatch until https://github.com/biopython/biopython/pull/3476
    _write_the_first_lines_orig = EmblWriter._write_the_first_lines

    def _write_the_first_lines(self, record):
        orig_type = record.annotations.get("molecule_type")
        record.annotations["molecule_type"] = orig_type.upper()
        ret = _write_the_first_lines_orig(self, record)
        record.annotations["molecule_type"] = orig_type
        return ret

    EmblWriter._write_the_first_lines = _write_the_first_lines

if __name__ == '__main__':
    convert(*get_args(sys.argv[1:]))

