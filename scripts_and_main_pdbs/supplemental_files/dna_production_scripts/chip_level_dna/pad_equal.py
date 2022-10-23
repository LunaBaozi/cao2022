#!/usr/bin/env python

# A script to pad protein sequences out to a specific length
#  using a repetitive pattern of amino acids

# to_order.seq should contain 2 space separated columns
#    SEQUENCE NAME
# with no header.

# Typical usage:

# ./pad_equal.py GGS 65 to_order.seq > padded_to_order.seq


import os
import sys
import numpy as np


pattern = sys.argv[1]
total_length = int(sys.argv[2])
sequence_file = sys.argv[3]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

sequences = []
names = []
with open(sequence_file) as f:
    for line in f:
        line = line.strip()
        if (len(line) == 0 ):
            continue
        sp = line.split()
        if ( len(sp) != 2 ):
            eprint("File not in SEQUENCE NAME format")
            assert(False)
        sequences.append(sp[0])
        names.append(sp[1])

mul = int(total_length / len(pattern) + 1)
full_pad = pattern * mul


for sequence, name in zip(sequences, names):

    missing = total_length - len(sequence)
    if ( missing < 0 ):
        eprint("Warning! A sequence was too long! Truncating to %i"%total_length)
        sequence = sequence[:total_length]
        missing = total_length - len(sequence)
    assert( missing >= 0 )

    missing_left = int(np.floor(missing/2))
    missing_right = int(np.ceil(missing/2))

    left_pad = ""
    if ( missing_left > 0 ):
        left_pad = full_pad[-missing_left:]
    right_pad = ""
    if ( missing_right > 0 ):
        right_pad = full_pad[:missing_right]

    final = left_pad + sequence + right_pad
    assert(len(final) == total_length)
    print(final, name)



