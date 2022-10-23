#!/usr/bin/env python

# Color your PDB by shannon entropy from counts data

# This is the super-simple no-assumptions method

# It is expected that for all the pdbs you pass, the pooled_counts.list
#  has a description column with PDB-NAME__SEQPOS__LETTER
# The native should be listed as PDB-NAME__native


# By default, the color-scheme is the same for all pdbs you run here.
#   i.e. --high_entropy and --low_entropy are determined and used across the board
# If instead you want per-pdb entropy, use --per_pdb_entropy


# From Cao 2021

# Brian Coventry 2021



import os
import sys
import pandas as pd
import xarray as xr
import numpy as np
import argparse
import scipy.stats
import warnings


parser = argparse.ArgumentParser(description="")
parser.add_argument("pooled_counts.list", type=str, help="A file with a header row (including description) and counts. ")
parser.add_argument("which_column", type=str, help="Which column in pooled_counts.list should we use?")
parser.add_argument("parent_pdbs", type=str, nargs="+", help="The pdb files you wish to perform this analysis on. Name must coordinate with"
                                                            +" pooled_counts.list")
parser.add_argument("--high_entropy", type=float, default=-1, help="Assign a max-entropy value so that colors are consistent across"
                                                                +" experiments. Disabled by default.")
parser.add_argument("--low_entropy", type=float, default=-1, help="Assign a low-entropy value so that colors are consistent across"
                                                                +" experiments. Disabled by default")
parser.add_argument("--per_pdb_entropy", action="store_true", help="Colors are relative for each PDB. Each output will have a different"
                                                                +" min and max entropy.")
parser.add_argument("--no_data_red", action="store_true", help="By default, no-data is colored blue. This will color it red.")
parser.add_argument("--output_folder", type=str, default="./", help="Where to output pdbs?")

args = parser.parse_args(sys.argv[1:])


print("Loading pyrosetta")

from pyrosetta import *
from pyrosetta.rosetta import *
init("-mute all -beta_nov16")


counts = pd.read_csv(args.__getattribute__("pooled_counts.list"), sep="\s+")

counts['description'] = counts['description'].str.replace("([^_])_([0-9]+)_([A-Z])$", r"\1__\2__\3")

assert(args.which_column in list(counts))
assert("description" in list(counts))
counts['my_counts'] = counts[args.which_column]
counts = counts[['my_counts', "description"]].copy()


def load_pose_data(pdbs):
    sequences = {}
    poses = {}

    for pdb in pdbs:
        print("    " + pdb)

        name = os.path.basename(pdb)
        if ( name.endswith(".gz") ):
            name = name[:-len(".gz")]
        if ( name.endswith(".pdb") ):
            name = name[:-len(".pdb")]

        pose = pose_from_file(pdb)
        
        full_seq = pose.sequence()
        sequence = full_seq[:pose.conformation().chain_end(1)]


        sequences[name] = sequence
        poses[name] = pose

    return sequences, poses



print("Loading pdbs")

sequences, poses = load_pose_data(args.parent_pdbs)


# Make the KD rows for the parent
parent_rows = []
for idx, row in counts[counts['description'].str.endswith("__native")].iterrows():
    name = row['description'].replace("__native", "")
    if ( name not in sequences ):
        continue
    seq = sequences[name]
    
    for il, let in enumerate(seq):
        seqpos = il + 1
        
        new_row = row.copy()
        new_row['native'] = True
        new_row['description'] = name + "__%i__%s"%(seqpos, let)
        parent_rows.append(new_row)

        
extra_counts = pd.DataFrame(parent_rows)

# Add parent rows to the counts and then remove the originals
counts['native'] = False
counts = pd.concat((counts, extra_counts))
old_size = len(counts)
counts = counts.drop_duplicates(subset=['description'], keep='last')
if ( len(counts) != old_size ):
    print("!!!!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!")
    print("  There were duplicate names in the dataframe")
    print("  Either you truly do have duplicate names in your experimental data")
    print("  Or you have experimental datapoints for mutations that were native")
    print("  This can also trigger if you have the wrong sequence")
    print("")
counts = counts[~counts['description'].str.endswith("__native")].copy()

parts = counts['description'].str.extract(r"(.*)__([0-9]+)__([A-Z])")
counts['ssm_parent'] = parts[0]
counts['ssm_seqpos'] = parts[1].astype(str)
counts['ssm_letter'] = parts[2]


# Convert to xarray because it's way easier to do calcualtions
ds = pd.pivot_table(counts, index=['ssm_parent', 'ssm_letter', 'ssm_seqpos']).to_xarray()

# Find probability of each aa at each seqpos
p_of_letter = ds['my_counts'].fillna(0) / ds['my_counts'].fillna(0).sum(dim='ssm_letter')

# Definition of Shannon Entropy
per_pos_entropy = -( p_of_letter * np.log( p_of_letter ) ).sum(dim='ssm_letter')

ds['entropy'] = per_pos_entropy

# Convert back to dataframe and remove the extra rows that xarray added
df = ds.to_dataframe().reset_index()
df['description'] = df['ssm_parent'] + "__" + df['ssm_seqpos'] + "__" + df['ssm_letter']
df = df.merge(counts[['description']], 'inner', 'description')


min_entropy = np.nanmin(df['entropy'])
max_entropy = np.nanmax(df['entropy'])

print("")
print("Minimum entropy: %.2f"%min_entropy)
print("Maximum entropy: %.2f"%max_entropy)

if ( args.low_entropy >= 0 ):
    print("")
    print("Reassigning minimum entropy to %.2f"%(args.low_entropy))
    min_entropy = args.low_entropy

if ( args.high_entropy >= 0 ):
    print("")
    print("Reassigning maximum entropy to %.2f"%(args.high_entropy))
    max_entropy = args.max_entropy


print("")
print("Outputing colored pdbs to %s:"%(args.output_folder))


for name, subdf in df.groupby('ssm_parent'):

    print("    " + name)

    local_min = min_entropy
    local_max = max_entropy
    if ( args.per_pdb_entropy ):
        local_min = np.nanmin(subdf['entropy'])
        local_max = np.nanmax(subdf['entropy'])
        print("      --per_pdb_entropy: Using min = %.2f, max=%.2f"%(local_min, local_max))


    pose = poses[name]
    pdb_info = pose.pdb_info()

    for seqpos in range(1, pose.size()+1):

        rows = subdf[subdf['ssm_seqpos'] == str(seqpos)]
        if ( len(rows) == 0 ):
            entropy = local_max if args.no_data_red else local_min
        else:
            row = rows.iloc[0]
            if ( np.isnan(row['entropy'])):
                entropy = local_max if args.no_data_red else local_min
            else:
                entropy = np.clip(row['entropy'], local_min, local_max)

        res = pose.residue(seqpos)
        for iatom in range(1, res.natoms()+1):
            pdb_info.bfactor(seqpos, iatom, entropy)

        # The last residue gets the min and max set so we ensure the final pdb has the correct range
        pdb_info.bfactor(pose.size(), 1, local_min)
        pdb_info.bfactor(pose.size(), 2, local_max)


    pose.pdb_info(pdb_info)

    path = os.path.join(args.output_folder, name + "_bfactor.pdb")
    pose.dump_pdb(path)





print("")
print("Done! Entropies inserted into b-factors")
print(" Visualize in pymol with:")
print("    spectrum b, blue_white_red ")
print("")
print(" Red is variable. Blue is conserved. No-data is %s"%("Red" if args.no_data_red else "Blue"))


