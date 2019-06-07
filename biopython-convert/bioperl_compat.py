import sys
from biopython_convert import get_args, convert
from Bio.SeqIO.InsdcIO import _InsdcWriter

# Quote anticodon qualifiers
_InsdcWriter.FTQUAL_NO_QUOTE = tuple(v for v in _InsdcWriter.FTQUAL_NO_QUOTE if v not in ['anticodon', 'transl_except'])

if __name__ == '__main__':
    input_path, input_type, output_path, output_type, jpath, split, stats = get_args(sys.argv[1:])

    with open(input_path, "r") as input:
        convert(input, input_type, output_path, output_type, jpath, split, stats)
