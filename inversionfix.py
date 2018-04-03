#!/usr/bin/env python
import os
import argparse

GROUPHOME = os.environ['GROUPHOME']

def main():
    parser = argparse.ArgumentParser(description='finds inversions and adds them to vcf output of biodiff',
                                     epilog='Runs nucmer and show-diff on fasta files to find inversions')
    parser.add_argument('-r', '--rdiff', help='dnadiff output, reference coordinates', default='out.rdiff')
    parser.add_argument('-q', '--qdiff', help='dnadiff output, query coordinates', default='out.qdiff')
    parser.add_argument('-v', '--vcf', help='temp vcf file', default='out-temp.vcf')
    parser.add_argument('-t', '--tfasta', help='temp fasta file', default='out-temp.fasta')
    parser.add_argument('-f', '--qfasta', help='query fasta file')
    args = parser.parse_args()

    # parse args.rdiff and make args.vcf with just the inversions, translocations, etc. in vcf format (no header)
    # which biodiff will eventually combine with its own vcf from the other differences it finds

    # parse args.qdiff to get inversions in querry coordinates
    # parse args.qfasta to get fasta of querry sequence,
    # undo the inversions found in args.qdiff
    # output noninverted fasta to args.tfasta, which biodiff will then use



if __name__ == '__main__':
    main()

