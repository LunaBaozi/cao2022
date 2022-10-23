import sys, os
import pymol
from pymol import cmd
import gzip, bz2
import re
import string


pdbname = cmd.get_object_list()[0]

def read_pdb( fname ):
    if os.path.exists( fname + ".pdb.gz" ):
        return gzip.open(fname + ".pdb.gz" )
    elif os.path.exists( fname + ".pdb" ):
        return open( fname + ".pdb" )
    elif os.path.exist( fname + ".pdb.bz2" ):
        return open( fname + ".pdb.bz2" )
    else:
        return False

f = read_pdb( pdbname )
if not f:
    print ("No pdb files found, exit!")
    exit(0)

pdb_infos = [ line.split()[2] for line in f if -1 != line.find("RIFRES") ]


sele_str = "not resi " + '+'.join(pdb_infos) + " and chainA"

cmd.select("sele", sele_str)
cmd.hide("sticks","sele")

sele_str = "resi " + '+'.join(pdb_infos) + " and chainA"


cmd.select("sele", sele_str)
cmd.show("sticks","sele")
util.cbay("sele")
cmd.zoom("sele", 20)
cmd.delete("sele")

print ("  ".join(pdb_infos))
