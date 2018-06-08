#!/usr/bin/env python
import os
import argparse
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline


def refineposition(refasta, queryseq, position):
    # blast 100 bp of query sequence around the start and end or each inversion against reference, to refine position
    blast_command = NcbiblastnCommandline(cmd='blastn', task='megablast', subject=refasta, outfmt=6)
    joinseq = ''.join(queryseq[(position-150):(position+150)])
    stdout, stderr = blast_command(stdin=joinseq)
    qstop = stdout.rstrip('\n').split('\t')[6]
    return int(qstop) + (position-150)


def main():
    #python ../../inversionfix.py -r 4-0010.rdiff -v 4-0010-temp.vcf -q 4-0010.qdiff -f $GROUPHOME/data/genomes/4-0010.fasta -p $GROUPHOME/resources/H37Rv.fasta -t 4-0010-temp.fasta
    parser = argparse.ArgumentParser(description='finds inversions and adds them to vcf output of biodiff',
                                     epilog='Runs nucmer and show-diff on fasta files to find inversions')
    parser.add_argument('-r', '--rdiff', help='dnadiff output, reference coordinates', default='out.rdiff')
    parser.add_argument('-q', '--qdiff', help='dnadiff output, query coordinates', default='out.qdiff')
    parser.add_argument('-v', '--vcf', help='temp vcf file', default='out-temp.vcf')
    parser.add_argument('-f', '--qfasta', help='query fasta file')
    parser.add_argument('-p', '--rfasta', help='reference fasta file')
    parser.add_argument('-t', '--tfasta', help='temp fasta file', default='out-temp.fasta')
    args = parser.parse_args()

    # Find inversions in terms of Query position, and output temp.fasta with those inversions undone
    # parse args.qdiff to get inversions in query coordinates
    inversions = []  # stores locations of inversions in query sequence
    with open(args.qdiff) as qfile:
        line2inv = False  # dnadiff records inversions in 2 lines. first occurrence is start, 2nd is end.
        for line in qfile:
            cleanline = line.strip('\n')
            qdiff = cleanline.split('\t')
            if qdiff[1] == 'INV' and not line2inv:  # first line of INV pair, append new inversion with Start key
                vardict = dict(Start=int(qdiff[2]))
                inversions.append(vardict)
                line2inv = True
            elif qdiff[1] == 'INV':  # 2nd linst of INV pair, set End key of last inversion dict in list
                inversions[-1]['End'] = int(qdiff[2])
                line2inv = False

    # parse query and reference sequence
    qrecord = SeqIO.parse(args.qfasta, 'fasta').next()  # read in the fasta file, should read only the first contig
    qseq = qrecord.seq  # get the sequence
    refrecord = SeqIO.parse(args.rfasta, 'fasta').next()
    rseq = refrecord.seq

    # refine inversion coordinates using blast
    for inv in inversions:
        inv['Start'] = refineposition(refasta=args.rfasta, queryseq=qseq, position=inv['Start']) - 1
        # -1 is to convert from blast's 1 based counting to python's 0 based counting
        inv['End'] = refineposition(refasta=args.rfasta, queryseq=qseq, position=inv['End']) - 1

    # undo the inversions found in args.qdiff
    fixedrecord = qrecord
    for invdict in inversions:
        before = qseq[0:invdict['Start']]
        fixedinv = qseq[invdict['Start']:invdict['End']].reverse_complement()
        after = qseq[invdict['End']:]
        qseq = before + fixedinv + after
    # output noninverted fasta to args.tfasta, which biodiff will then use
    fixedrecord.seq = qseq
    SeqIO.write(fixedrecord, args.tfasta, 'fasta')

    # Find inversions in terms of Reference position, and output them in VCF format
    chromname = qrecord.id
    # parse args.rdiff and make args.vcf with just the inversions, translocations, etc. in vcf format (no header)
    # which biodiff will eventually combine with its own vcf from the other differences it finds
    vcflines = []
    with open(args.rdiff) as rfile:
        line2inv = False  # dnadiff records inversions in 2 lines. first occurrence is start, 2nd is end.
        invstart = 0
        for line in rfile:
            cleanline = line.strip('\n')
            rdiff = cleanline.split('\t')
            # rdiff[2] is genome position, rdiff[1] is variation type (GAP, INV, etc)
            if rdiff[1] == 'INV' and not line2inv:
                invstart = refineposition(args.qfasta, rseq, int(rdiff[2]))  # blast position in reference to query sequence
                vcfstring = chromname + '\t' + str(invstart) + '\t.\tN\tINV\t.\tPASS\t'
                line2inv = True
            elif rdiff[1] == 'INV':
                invend = refineposition(args.qfasta, rseq, int(rdiff[2]))
                vcfstring += 'SVLEN=' + str(int(invend - invstart)) + ';SVTYPE=INV\n'
                vcflines.append(vcfstring)
                line2inv = False
    with open(args.vcf, 'w') as vfile:
        for line in vcflines:
            vfile.write(line)

if __name__ == '__main__':
    main()

