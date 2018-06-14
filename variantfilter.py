#!/usr/bin/env python
import argparse


def main():
    parser = argparse.ArgumentParser(description='Adjusts output of udiff2vcf',
                                     epilog='Report only position and length for variants over 10bp long.')
    parser.add_argument('-i', '--invcf', help='input vcf file from udiff2vcf')
    parser.add_argument('-o', '--ovcf', help='output vcf file (no header)')
    args = parser.parse_args()

    outlist = list()  # list to store modified lines for output as new vcf file
    inversionlist = list()  # list to store inversion positions, to check if variants are inside or near them
    with open(args.invcf) as infile:
        for line in infile:
            cleanline = line.strip('\n')
            vcfbits = cleanline.split('\t')
            if vcfbits[4] == '<INV>':  # parse inversion entries to create list of inversion positions
                infobits = vcfbits[7].split(';')
                invlength = infobits[0].split('=')
                inversiondict = dict(Start=int(vcfbits[1]), End=int(vcfbits[1])+int(invlength[1]))
                inversionlist.append(inversiondict)
                outlist.append(line)
            # Note: All SV entrees are at the start of the input file
            else:  # parse other entries and adjust based on length of REF/ALT and position relative to inversions
                if len(vcfbits[3]) > 10 or len(vcfbits[4]) > 10:
                    vcfbits[3] = '.'
                    vcfbits[4] = '.'
                    vcfbits[7] = 'SVLEN=' + str(len(vcfbits[4]))
                intrainv = False  # flag to track if variant is inside an inversion
                nearinv = False  # flag to track if variant is near the edge of an inversion
                for inv in inversionlist:
                    if inv['Start'] < int(vcfbits[1]) < inv['End']:
                        intrainv = True
                    if inv['Start'] - 10 < int(vcfbits[1]) < inv['Start'] + 1:
                        nearinv = True
                    if inv['End'] - 1 < int(vcfbits[1]) < inv['End'] + 10:
                        nearinv = True
                if intrainv and vcfbits[7] == '.':  # replace the . if the INFO column is empty, append otherwise
                    vcfbits[7] = 'Ambiguity=Inside_Inversion'
                elif intrainv:
                    vcfbits[7] += ';Ambiguity=Inside_Inversion'
                if nearinv and vcfbits[7] == '.':  # replace the . if the INFO column is empty, append otherwise
                    vcfbits[7] = 'Ambiguity=Near_Inversion'
                elif nearinv:
                    vcfbits[7] += ';Ambiguity=Near_Inversion'
                newline = '\t'.join(vcfbits) + '\n'  # convert entry back to string for writing
                outlist.append(newline)

    # write new temp vcf
    with open(args.ovcf, 'w') as vfile:
        for oline in outlist:
            vfile.write(oline)


if __name__ == '__main__':
    main()
