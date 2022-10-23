#!/usr/bin/env python

import os
import sys
import argparse

from pyrosetta import *
from pyrosetta.rosetta import *

init("-beta_nov16")




parser = argparse.ArgumentParser(description="Open this file with a text editor for very important information!")
parser.add_argument("original_pdb", type=str, help="The PDB file to maken an SSM for. Binder = Chain A")
parser.add_argument("--output_pdbs", action="store_true", help="If you'd rather have pdbs than a silent file, use this.")

args = parser.parse_args(sys.argv[1:])


og_pdb = args.original_pdb
name = os.path.basename(og_pdb).replace(".gz", "").replace(".pdb", "")


def to_silent(pose, tag):
    silent_name = "%s.silent"%name
    sfd_out = core.io.silent.SilentFileData( silent_name, False, False, "binary", core.io.silent.SilentFileOptions())
    struct = sfd_out.create_SilentStructOP()
    struct.fill_struct(pose, tag)
    sfd_out.add_structure(struct)
    sfd_out.write_all(silent_name, False)


def break_disulfide(pose, seqpos):
    if ( pose.residue(seqpos).name1() != "C" ): return
    try:
        other = core.conformation.get_disulf_partner(pose.conformation(), seqpos)
    except:
        return

    try:
        core.conformation.break_disulfide(pose.conformation(), seqpos, other)
    except:
        pass


print(og_pdb)
assert(os.path.exists(og_pdb))

og_pose = pose_from_file(og_pdb)
og_pose.pdb_info(core.pose.PDBInfo(og_pose))


sequence = og_pose.sequence()[:og_pose.conformation().chain_end(1)]

scorefxn = get_fa_scorefxn()

natives_written = set()

for pos in range(1, og_pose.conformation().chain_end(1)+1):

    # Skip disulfides because they make the next few steps crash. 
    # there is certainly an easy work around
    for name1 in "ACDEFGHIKLMNPQRSTVWY":


        pose = og_pose.clone()
        if ( pose.residue(pos).name1() == "C" ):
            break_disulfide( pose, pos )
            

        protocols.toolbox.pose_manipulation.repack_this_residue(pos, pose, scorefxn, False, name1)

        pdb_info = core.pose.PDBInfo(pose)
        pdb_info.add_reslabel(pos, "MUT")
        pose.pdb_info(pdb_info)

        mut = name + "__%i__%s"%(pos, name1)

        if ( args.output_pdbs ):
            pose.dump_pdb(mut + ".pdb")
        else:
            to_silent(pose, mut)











