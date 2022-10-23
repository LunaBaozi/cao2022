#!/usr/bin/env python

import sys
from Bio.Seq import Seq

if len( sys.argv ) != 2:
    print( "Usage: " + sys.argv[0] + " DNA_sequence.list")
    exit(0)

fname = sys.argv[1]

with open(fname) as f:
    lines = [ line.strip() for line in f ]

for line in lines:
    temp = line.split()
    dna = temp[-1]
    prot = temp[-2]

    if len(dna) % 3 != 0:
        print ("This sequence is not a multiple of three: %s" % line )
        sys.exit(0)

    coding_dna = Seq(dna)
    tr = str(coding_dna.translate())
    if not tr == prot:
        print( "Protein seq and DNA sequence don't match: %s" % line )
        sys.exit(0)
    else:
        pass

print( "Everything is OK, just order!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
