#!/usr/bin/env python
import os
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description='creates new fasta with inversions for testing')
    parser.add_argument('-f', '--fasta', help='fasta file', default='bit-H37Rv.fasta')
    parser.add_argument('-t', '--ifasta', help='inversion fasta file', default='ibit-H37Rv.fasta')
    parser.add_argument('-c', '--cfasta', help='normal fasta file with same line output', default='cbit-H37Rv.fasta')
    args = parser.parse_args()

    record = SeqIO.parse(args.fasta, 'fasta').next()  # read in the fasta file, should read only the first contig
    SeqIO.write(record, args.cfasta, 'fasta') # write out same record in fasta with same bp per line for quick compare
    seq = record.seq  # get the sequence of the genome

    inversion1 = seq[200:570].reverse_complement()  # so inversion should be from 201 to 571 in 1 based counting
    iseq = seq[0:200] + inversion1 + seq[570:]
    inversion2 = seq[1140:1590].reverse_complement()  # 1141 to 1591
    iseq2 = iseq[0:1140] + inversion2 + iseq[1590:]

    irecord = record
    irecord.seq = iseq2
    SeqIO.write(irecord, args.ifasta, 'fasta')


if __name__ == '__main__':
    main()