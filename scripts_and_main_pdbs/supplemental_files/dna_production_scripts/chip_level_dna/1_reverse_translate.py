#!/usr/bin/env python

# Written by Longxing Cao
#   touched up by Brian Coventry

import sys
import argparse
import hashlib
from collections import namedtuple
from string import Template
import os

job_template = Template('''solutions 1
repeat 8
LOGFILE "${fname}"
pattern
  BamHI GGATCC
  NdeI CATATG
  XhoI CTCGAG
  NheI GCTAGC
  BsaI GGTCTC
  BsaI GAGACC
  PolyA AAAAAAAA
  PolyG GGGGG
  PolyT TTTTTTTT
  PolyC CCCCCCCC
  Aarl  CACCTGC
//
                      
codon ${codon}

protein
${sequence}
//'''
)

def parse_arguments( argv ):
    argv_tmp = sys.argv
    sys.argv = argv
    description = 'Reverse Translate accepts a protein sequences and uses DNAWorks to generate a DNA sequence representing the most likely non-degenerate coding sequence.'
    parser = argparse.ArgumentParser( description = description )
    parser.add_argument('-seqs', type=str, nargs='*', help='one or more sequences')
    parser.add_argument('-seq_list', type=str, help='a list file of all sequences')
    parser.add_argument('-organism', type=str, default='yeast', choices=['yeast', 'ecoli'], help="yeast display or ecoli expression")
    parser.add_argument('-padding', type=int, default=-1, help='pad all the sequences to the same length')
    parser.add_argument('-padding_seq', type=str, default='GGS', help='linker sequence to add')
    parser.add_argument('-rename_fname', type=str, default="rename.list", help="the name of the output file")
    parser.add_argument('-jobs_fname', type=str, default="dnaworks_commands.list", help="the name of the output file")
    parser.add_argument('-dnaworks', type=str, default="", help="Path to the dnaworks executable")
    args = parser.parse_args()
    sys.argv = argv_tmp
    return args

def get_key( md5_str, existing_keys, key_len = 25):
    index = 0
    while True:
        key = md5_str[ index : index+key_len ]
        if key not in existing_keys:
            break
        else:
            index += 1
            if index >= 8:
                print( "Key conflict, you are so lucky. There must be lots of duplicates" )
                exit(0)

    return key

def main( argv ):
    args = parse_arguments( argv )

    if ( args.dnaworks == "" or not os.path.exists(args.dnaworks) ):
        print("You must pass -dnaworks with a path to the dnaworks executable.")

    SeqInfo = namedtuple('SeqInfo', 'key name sequence')
    seqs_dict = { } # use a dict to store all sequences
    md5 = hashlib.md5()
    commands_all = list()
    if args.seqs == None:
        assert (args.seq_list != None)
        with open(args.seq_list) as f:
            all_seqs = [line.strip() for line in f]
        for seq in all_seqs:
            temp = seq.split()
            prot_name = temp[1]
            prot_seq = temp[0]
            if args.padding != -1:
                prot_seq = (prot_seq+args.padding_seq*100)[:args.padding]
            # currently the first is name and the second is the sequence
            md5.update(prot_name.encode('ascii'))
            key = get_key( md5.hexdigest(), seqs_dict.keys() )
            seqs_dict[key] = SeqInfo(key=key, name = prot_name, sequence = prot_seq)
    else:
        all_seqs = args.seqs
        for seq in all_seqs:
            prot_seq = seq
            if args.padding != -1:
                prot_seq = (prot_seq+args.padding*100)[:args.padding]
            md5.update(prot_seq.encode('ascii'))
            key = get_key( md5.hexdigest(), seqs_dict.keys() )
            seqs_dict[key] = SeqInfo(key=key, name=seq, sequence=prot_seq)

    pwd = os.getcwd()
    for key, val in seqs_dict.items():
        job_fname = val.key + '.dwj'
        output_fname = val.key + '.dwo'
        job_str = job_template.substitute( fname = output_fname, sequence = val.sequence , codon = 'S. cerevesiae' if args.organism == 'yeast' else 'E. coli')
        with open(job_fname, 'w') as f:
            f.write( job_str )
        commands_all.append( 'cd %s; %s %s' %(pwd, args.dnaworks, job_fname ) )

    with open(args.jobs_fname, 'w') as f:
        for com in commands_all:
            f.write(com + '\n')
            
    with open( args.rename_fname, 'w') as f:
        for key, val in seqs_dict.items():
            f.write( "%s %s %s\n" % (val.key, val.name, val.sequence) )

    print("")
    print( "Job file generation done!!!" )
    print( args.jobs_fname + " contains all of the command to run!!!" )
    print( "\nGood luck!!!\n")

if __name__ == '__main__':
        main( sys.argv )
