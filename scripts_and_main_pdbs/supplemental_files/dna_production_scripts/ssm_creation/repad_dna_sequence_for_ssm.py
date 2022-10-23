#!/usr/bin/env python

import os
import sys
from collections import defaultdict



dna_name_file = sys.argv[1]
seq_name_file = sys.argv[2]
ssm_file = sys.argv[3]




name_to_dna = {}
with open(dna_name_file) as f:
    for line in f:
        line = line.strip()
        if (len(line) == 0):
            continue
        sp = line.split()

        if ( len(sp) != 2 ):
            print("File not in SEQUENCE NAME format: ", dna_name_file)
            assert(False)

        name = sp[1]
        dna = sp[0]

        name_to_dna[name] = dna

seq_name_pairs = []
with open(seq_name_file) as f:
    for line in f:
        line = line.strip()
        if (len(line) == 0):
            continue
        sp = line.split()
        seq_name_pairs.append(sp)


seq_name_pairs_by_parent = defaultdict(list)
with open(ssm_file) as f:
    for line in f:
        line = line.strip()
        if (len(line) == 0):
            continue
        sp = line.split()
        parent = "__".join(sp[1].split("__")[:-2])
        seq_name_pairs_by_parent[parent].append(sp)



table = { 
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
}

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def translate_dna(dna, force_partial=None):
    translated = ""
    for i in range(0, len(dna), 3):
        part = dna[i:i+3]

        if ( len(part) < 3 ):
            if ( force_partial is None ):
                eprint("Partial sequence!!!")
            else:
                translated += force_partial
        else:
            part = part.upper()
            if ( part not in table ):
                eprint("Unrecognized codon! " + part)
                assert(False)
            translated += table[part]
    return translated



for prot_seq, name in seq_name_pairs:
    full_dna = name_to_dna[name]

    translated = translate_dna(full_dna)

    assert(prot_seq in translated)

    prot_offset = translated.find(prot_seq)
    dna_offset = prot_offset*3
    dna_end = dna_offset + len(prot_seq)*3

    prot_dna = full_dna[dna_offset:dna_end]
    assert(translate_dna(prot_dna) == prot_seq)

    pad_start = full_dna[:dna_offset]
    pad_end = full_dna[dna_end:]

    assert(name in seq_name_pairs_by_parent)
    child_pairs = seq_name_pairs_by_parent[name]

    print(full_dna, name + "__native")

    for child_dna_seq, child_name in child_pairs:
        print(pad_start + child_dna_seq + pad_end, child_name)








