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

    if len(sys.argv) < 4:
        print("Missing arguments", file=sys.stderr)
        print(usage, file=sys.stderr)
        exit(1)

    format = sys.argv[1]
    ids = defaultdict(int)
    
    def makeUnique(seq):
        seqlenlen = len(str(len(seq)))
        newid = seq.id[:26 - seqlenlen]  # Genbank has a max length for the id and sequence length number, truncate the sequence id if too long
        count = ids[newid]
        ids[newid] += 1
        if count:
            suffix = "_" + str(count)
            if len(newid) + len(suffix) + 1 + seqlenlen > 26:
               newid = newid[:25 - seqlenlen - len(suffix)]

            newid += suffix
            seq.name += suffix

        if seq.id != newid:
            print(f"{seq.id}\t{newid}")

        seq.id = newid
        return seq


    paths = iter(sys.argv[2:])

    for input, output in zip(paths, paths):
        SeqIO.write(
            map(makeUnique, SeqIO.parse(input, format)),
            output,
            format
        )


