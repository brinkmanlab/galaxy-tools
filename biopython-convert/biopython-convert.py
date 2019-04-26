#!/usr/bin/env python
from Bio import SeqIO
import itertools
import getopt
import gffutils
from gffutils import biopython_integration
import sys

usage = "Use: biopython-convert.py [-s] [-v] input_file input_type output_file output_type\n"\
        "\t-s Split records into seperate files\n"\
        "\t-v Print version and exit\n"

def append_filename(path, s):
    """
    Append a string to a file name before the extension
    :param path: file path
    :param s: string to append
    :return: appended file name
    """
    s = output_path.rsplit(".", 1)
    if len(s) > 1:
        return f"{s[0]}{s}.{s[1]}"
    else:
        return path + s

if __name__ == '__main__':
    split = False
    # Parse arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'vs')
        for opt, val in opts:
            if opt == '-v':
                import __version

                print(__version.__versionstr__)
                exit(0)
            elif opt == '-s':
                split = True

    except getopt.GetoptError as err:
        print("Argument error(", err.opt, "): ", err.msg, file=sys.stderr)
        args = []

    if len(args) < 4:
        print(usage, file=sys.stderr)
        exit(1)

    gff_types = ["gff", "gff3"]
    input_path = args[0]
    input_type = args[1]
    output_path = args[2]
    output_type = args[3]

    with open(input_path, "rU") as input:
        if input_type not in gff_types and output_type not in gff_types:
            SeqIO.convert(sys.stdin, args[1], sys.stdout, args[1])
            exit(0)

        if input_type in gff_types:
            db = gffutils.create_db(input, ":memory:", merge_strategy="create_unique")
            input_records = map(
                lambda x: SeqIO.SeqRecord("", features=list(x)),
                itertools.groupby(
                    map(biopython_integration.to_seqfeature, db.all_features(order_by="seqid")),
                    lambda x: x.id
                )
            )
        else:
            input_records = SeqIO.parse(input_path, input_type)

        output = open(append_filename(output_path, ".0") if split else output_path, "w")

        for i, record in enumerate(input_records):
            for feature in record.features:
                if output_type in gff_types:
                    print(biopython_integration.from_seqfeature(feature), file=output)
                else:
                    SeqIO.write(record, output, args[3])
            if split:
                output.close()
                output = open(append_filename(output_path, "." + (i+1)))

