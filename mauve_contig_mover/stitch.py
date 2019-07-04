#!/usr/bin/env python
import sys

from Bio import SeqIO, Alphabet
from Bio.Seq import Seq
import csv

usage = """
Mauve Contig Mover - Stitch
Stitch contigs into a single contig.
Compliments reversed sequences and rewrites all feature coordinates.

Use: stitch.py [-v] <padding length> <MauveCM contigs.tab path> <draft file path> <draft file format> [final sequence id]
\t-v Print version and exit
Valid draft file formats:
abi, abi-trim, ace, cif-atom, cif-seqres, clustal, embl, fasta, fasta-2line, fastq-sanger, fastq, fastq-solexa, fastq-illumina,
genbank, gb, ig, imgt, nexus, pdb-seqres, pdb-atom, phd, phylip, pir, seqxml, sff, sff-trim, stockholm, swiss, tab, qual, uniprot-xml, gff3
"""


def getOrder(path):
    """
    Parse MCM contig order file and iterate rows after "Ordered Contigs"
    :param path: path to MCM *_contig.tab
    :return: tuple(type, label, contig_type, strand, left_end, right_end)
    """
    with open(path, "r") as alignment_file:
        alignments = iter(csv.reader(alignment_file, delimiter="\t"))
        try:
            alignment = next(alignments)

            # Jog to beginning of ordered alignments
            while not (len(alignment) and "Ordered Contigs" == alignment[0]):
                alignment = next(alignments)

            # Skip column header
            next(alignments)

            while True:
                yield next(alignments)
        except StopIteration:
            return


def stitch(pad, contigs, order):
    """
    Reduce contigs to single contig by concatenation.
    Compliments reversed sequences and rewrites all feature coordinates.
    :param pad: Seq or SeqRecord instance to insert between contigs 
    :param contigs: dict of SeqRecords keyed on the record name
    :param order: iterable of tuples containing sequence names and orientation 
    :return: concatentated SeqRecord
    """
    result = None
    # Concat in order with padding
    for alignment in order:
        if len(alignment) < 4: continue
        try:
            contig = contigs.pop(alignment[1])  # type: SeqIO.SeqRecord
            if alignment[3] == "complement":
                contig = contig.reverse_complement()
            if result:
                # A lot is happening in the background here. Biopython handles the feature coordinates implicitly.
                result += pad + contig
            else:
                result = contig
                pad.alphabet = result.seq.alphabet
        except KeyError:
            pass

    # Concat remaining in arbitrary order
    for unordered in contigs.values():
        if result:
            result += pad + unordered
        else:
            result = unordered

    return result


if __name__ == '__main__':
    if '-v' in sys.argv:
        print('1.0')
        exit(0)

    if len(sys.argv) < 5:
        print("Missing arguments", file=sys.stderr)
        print(help, file=sys.stderr)
        exit(1)

    pad_len = int(sys.argv[1])
    if pad_len < 0:
        print("Padding length must be >= 0", file=sys.stderr)
        print(help, file=sys.stderr)
        exit(1)

    contig_path = sys.argv[2]
    draft_path = sys.argv[3]
    draft_format = sys.argv[4]

    order = getOrder(contig_path)
    pad = Seq('N'*pad_len)
    contigs = {seq.name: seq for seq in SeqIO.parse(draft_path, draft_format)}

    result = stitch(pad, contigs, order)

    if result:
        # Ensure there is only one 'source' feature
        # TODO

    if result and len(sys.argv) > 5:
        result.id = sys.argv[5]
        result.description = ""

    SeqIO.write(result, sys.stdout, draft_format)
