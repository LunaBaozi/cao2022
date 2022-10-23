#!/usr/bin/env python

import sys
import os
import argparse
import re


def parse_arguments( argv ):
    argv_tmp = sys.argv
    sys.argv = argv
    description = 'Collecting all the dna sequences from the DNAwork outputs.'
    parser = argparse.ArgumentParser( description = description )
    parser.add_argument('-rename_fname', type=str, default='rename.list', help='the list file with keys from the first stript')
    parser.add_argument('-output', type=str, default='DNA_sequence.list', help="final output file to order.")
    args = parser.parse_args()
    sys.argv = argv_tmp
    return args

def get_dnaseq( fname ):
    s = ''
    with open( fname ) as f:
        lines = [line.strip() for line in f]
    count = 0
    for line in lines:
        if -1 != line.find( "The DNA sequence" ):
            count += 1
            continue
        if count == 0:
            continue
        if -1 != line.find('-------------------------------------'):
            count += 1
            continue
        if count == 3 : break
        s += line
    return re.sub('[\s\r\n\t0-9]', '', s)

def main( argv ):
    args = parse_arguments( argv )

    if not os.path.exists( args.rename_fname ):
        print( "Can not find the rename file which contains the key values I need!!")
        exit(0)

    final_seqs = []
    with open( args.rename_fname ) as f:
        # key name protein_seq
        lines = [ line.strip() for line in f ]

    for line in lines:
        temp = line.split()
        dnawork_out_fname = temp[0] + '.dwo'
        dna_seq = get_dnaseq( dnawork_out_fname )
        final_seqs.append( "%s %s" % (line, dna_seq) )

    with open( args.output, 'w' ) as f:
        for line in final_seqs:
            f.write( line + '\n' )

if __name__ == '__main__':
        main( sys.argv )
