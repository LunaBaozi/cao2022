#!/usr/bin/env python

import sys
from Bio.Seq import Seq

yeast_petcon3_primer5 = "GGGTCGGCTTCGCATATG"
yeast_petcon3_primer3 = "CTCGAGGGTGGAGGTTCC"

if len(sys.argv) != 3:
    print( "Usage: " + sys.argv[0] + " input_file output_file")
    exit(0)

assert( len(sys.argv) == 3 )
fname = sys.argv[1]
fout  = sys.argv[2]

with open(fname) as f:
    lines = [ line.strip() for line in f ]

new_lines = []
for line in lines:
    temp = line.split()
    temp[-1] = (yeast_petcon3_primer5 + temp[-1] + yeast_petcon3_primer3)[:230]
    
    new_lines.append( ' '.join(temp) )

for index in range( len(lines) ):
    if not lines[index].split()[-1] in new_lines[index].split()[-1]:
        print( "Something is wrong!!!!!!!" )
        exit(0)
    if len( new_lines[index].split()[-1] ) > 230:
        print( "DNA length exceeds the maximun length!!!!")
        exit(0)

with open(fout, 'w') as f:
    for line in new_lines:
        f.write(line+'\n')


print( "Everything is OK, adaptors added!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
