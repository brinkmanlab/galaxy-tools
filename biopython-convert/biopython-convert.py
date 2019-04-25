from Bio import SeqIO
import gffutils
from gffutils import biopython_integration
import sys

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) < 2:
        print("Use: biopython-convert.py [input type] [output type] < [input file] > [output file]", sys.stderr)
        exit(1)

    if args[1] in ["gff", "gff3"]:
        if args[0] in ["gff", "gff3"]:
            print("Just copy it silly..", sys.stderr)
            exit(1)
        else:
            for feature in SeqIO.parse(sys.stdin, args[0]):
                print(biopython_integration.from_seqfeature(feature))
    else:
        if args[0] in ["gff", "gff3"]:
            db = gffutils.create_db(input, ":memory:", merge_strategy="create_unique")
            for feature in db.all_features():
                SeqIO.write(biopython_integration.to_seqfeature(feature), sys.stdout, args[1])
        else:
            SeqIO.convert(sys.stdin, args[0], sys.stdout, args[1])
