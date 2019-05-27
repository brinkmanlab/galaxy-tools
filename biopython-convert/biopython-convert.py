#!/usr/bin/env python
"""
Biopython Convert
Convert between any formats that Biopython supports or gffutils.
Provides a means of querying/filtering documents using JMESPath query language.
"""
from Bio import SeqIO
import itertools
import getopt
import gffutils
from gffutils import biopython_integration
import JMESPathGen
import sys

usage = "Use: biopython-convert.py [-s] [-v] [-f JMESPath] input_file input_type output_file output_type\n" \
        "\t-s Split records into seperate files\n" \
        "\t-q JMESPath to select records. Must return list of SeqIO records. Root is list of input SeqIO records.\n" \
        "\t-v Print version and exit\n"

def append_filename(path, s):
    """
    Append a string to a file name before the extension
    :param path: file path
    :param s: string to append
    :return: appended file name
    """
    seg = path.rsplit(".", 1)
    if len(seg) > 1:
        return seg[0] + s + "." + seg[1]
    else:
        return path + s

if __name__ == '__main__':
    split = False
    jpath = None
    # Parse arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'vsq:')
        for opt, val in opts:
            if opt == '-v':
                import __version

                print(__version.__versionstr__)
                exit(0)
            elif opt == '-s':
                split = True
            elif opt == '-q':
                jpath = val

    except getopt.GetoptError as err:
        print("Argument error(" + str(err.opt) + "): " + err.msg, file=sys.stderr)
        args = []

    # Check for minimum number of arguments
    if len(args) < 4:
        print(usage, file=sys.stderr)
        exit(1)

    gff_types = ["gff", "gff3"]
    input_path = args[0]
    input_type = args[1]
    output_path = args[2]
    output_type = args[3]

    with open(input_path, "r") as input:
        if input_type in gff_types:
            # If input is GFF load with gffutils library
            db = gffutils.create_db(input, ":memory:", merge_strategy="create_unique")
            # Wrap features in generator that converts to BioPython SeqRecords
            input_records = map(
                lambda x: SeqIO.SeqRecord("", features=list(x)),
                itertools.groupby(
                    map(biopython_integration.to_seqfeature, db.all_features(order_by="seqid")),
                    lambda x: x.id
                )
            )
        else:
            input_records = SeqIO.parse(input_path, input_type)

        # Open output file with file name suffix if splitting
        output = open(append_filename(output_path, ".0") if split else output_path, "w")

        # Wrap input in JMESPath selector if provided
        if jpath:
            input_records = JMESPathGen.search(jpath, input_records)

        for i, record in enumerate(input_records):
            # TODO allow objects other than SeqRecord, transform to SeqRecord or handle special output (like if output format == txt|json, pretty print object)
            if output_type in gff_types:
                # If output type is GFF, use gffutils library
                for feature in record.features:
                    print(biopython_integration.from_seqfeature(feature), file=output)
            else:
                SeqIO.write(record, output, args[3])
            if split:
                # If splitting, open next file
                output.close()
                output = open(append_filename(output_path, "." + str(i+1)))

