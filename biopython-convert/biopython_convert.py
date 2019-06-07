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

usage = "Use: biopython_convert.py [-s] [-v] [-i] [-q JMESPath] input_file input_type output_file output_type\n" \
        "\t-s Split records into seperate files\n" \
        "\t-q JMESPath to select records. Must return list of SeqIO records. Root is list of input SeqIO records.\n" \
        "\t-i Print out details of records during conversion\n" \
        "\t-v Print version and exit\n"

gff_types = ["gff", "gff3"]

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

def get_args(sysargs):
    split = False
    jpath = None
    stats = False
    # Parse arguments
    try:
        opts, args = getopt.gnu_getopt(sysargs, 'vsiq:')
        for opt, val in opts:
            if opt == '-v':
                import __version

                print(__version.__versionstr__)
                exit(0)
            elif opt == '-s':
                split = True
            elif opt == '-q':
                jpath = val
            elif opt == '-i':
                stats = True

    except getopt.GetoptError as err:
        print("Argument error(" + str(err.opt) + "): " + err.msg, file=sys.stderr)
        args = []

    # Check for minimum number of arguments
    if len(args) < 4:
        print(usage, file=sys.stderr)
        exit(1)

    input_path = args[0]
    input_type = args[1]
    output_path = args[2]
    output_type = args[3]

    return input_path, input_type, output_path, output_type, jpath, split, stats

def convert(input, input_type, output_path, output_type, jpath, split, stats):
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
        input_records = SeqIO.parse(input, input_type)

    # Open output file with file name suffix if splitting
    output = open(append_filename(output_path, ".0") if split else output_path, "w")

    # Wrap input in JMESPath selector if provided
    if jpath:
        input_records = JMESPathGen.search(jpath, input_records)

    if isinstance(input_records, dict):
        # Support generating a new record in JMESPath
        input_records = SeqIO.SeqRecord(**input_records)

    if isinstance(input_records, SeqIO.SeqRecord):
        # Support returning single record from JMESPath
        input_records = (input_records,)

    for i, record in enumerate(input_records):  # type: (Integer, SeqIO.SeqRecord)
        # TODO allow objects other than SeqRecord, transform to SeqRecord or handle special output (like if output format == txt|json, pretty print object)
        if isinstance(input_records, dict):
            # Support generating new records in JMESPath
            record = SeqIO.SeqRecord(**record)

        if stats:
            print(f"${i}\t${record.id}\t${record.name}\t${record.description}\t${len(record)}")

        if output_type in gff_types:
            # If output type is GFF, use gffutils library
            for feature in record.features:
                print(biopython_integration.from_seqfeature(feature), file=output)
        else:
            SeqIO.write(record, output, output_type)
        if split:
            # If splitting, open next file
            output.close()
            output = open(append_filename(output_path, "." + str(i + 1)))

if __name__ == '__main__':
    input_path, input_type, output_path, output_type, jpath, split, stats = get_args(sys.argv[1:])

    with open(input_path, "r") as input:
        convert(input, input_type, output_path, output_type, jpath, split, stats)
