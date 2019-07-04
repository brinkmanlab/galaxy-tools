#!/usr/bin/env python
import sys
from Bio import SeqIO
from collections import defaultdict

usage = """
make_unique_id
Makes all record ids unique across all input data.
All input data must be the same format.

Use: make_unique_id.py [-v] <format> <input1> <output1> [<input2> <output2> ... <inputn> <outputn>]
\t-v Print version and exit 

Valid formats: clustal, embl, fasta, fasta-2line, fastq-sanger, fastq, fastq-solexa, fastq-illumina, genbank, gb, imgt,
nexus, phd, phylip, pir, seqxml, sff, stockholm, tab, qual
"""

if __name__ == '__main__':
    if '-v' in sys.argv:
        print('1.0')
        exit(0)

    if len(sys.argv) < 5:
        print("Missing arguments", file=sys.stderr)
        print(help, file=sys.stderr)
        exit(1)

    format = sys.argv[1]
    ids = defaultdict(int)
    
    def makeUnique(seq):
        global ids
        count = ids[seq.id]
        if count:
            count += 1
            suffix = "_" + str(count)
            print(f"{seq.id}\t{seq.id}{suffix}")
            seq.id += suffix
            seq.name += suffix
        else:
            count += 1

    paths = iter(sys.argv[2:])

    for input, output in zip(paths, paths):
        SeqIO.write(
            map(makeUnique, SeqIO.parse(input, format)),
            format
        )


