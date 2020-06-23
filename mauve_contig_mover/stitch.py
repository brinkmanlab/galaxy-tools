#!/usr/bin/env python
import sys
import csv
import getopt

from Bio import SeqIO, Alphabet
from Bio.Seq import Seq

usage = """
Mauve Contig Mover - Stitch
Stitch contigs into a single contig.
Compliments reversed sequences and rewrites all feature coordinates.

Use: stitch.py [-v] [-s 'final sequence id'] <padding length> <draft file path> <draft file format> [MauveCM contigs.tab path]
\t-v Print version and exit
\t-s Provide an ID for the final sequence, the first sequence ID will be used otherwise
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
    seqid = None
    # Parse arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'vs:iq:')
        for opt, val in opts:
            if opt == '-v':
                print('1.0')
                exit(0)
            elif opt == '-s':
                seqid = val
    except getopt.GetoptError as err:
        print("Argument error(" + str(err.opt) + "): " + err.msg, file=sys.stderr)
        args = []
    
    # Check for minimum number of arguments
    if len(args) < 3:
        print(usage, file=sys.stderr)
        exit(1)    

    pad_len = int(args[0])
    if pad_len < 0:
        print("Padding length must be >= 0", file=sys.stderr)
        print(help, file=sys.stderr)
        exit(1)

    draft_path = args[1]
    draft_format = args[2]
    
    if len(args) < 4:
        order = ()
    else:
        order = getOrder(args[3])
        
    pad = Seq('N'*pad_len)
    contigs = {seq.name: seq for seq in SeqIO.parse(draft_path, draft_format)}

    result = stitch(pad, contigs, order)

    if result:
        # Ensure there is only one 'source' feature
        # TODO
        pass

    if result and seqid:
        result.id = seqid
        result.description = ""

    result.seq.alphabet = Alphabet.generic_dna # TODO Investigate why this is required for some datasets
    SeqIO.write(result, sys.stdout, draft_format)
