#!/usr/bin/env python

# Color your PDB by shannon entropy from counts data

# This script should give more accurate entropy estimates as it will pick the best pool
#  to calculate entropy from (See --entropy_pool_conc_factor).
#  Alternatively, you can calculate pseudo-entropy based on
#  the estimated KD values. The pseudo-entropies have the advantage of looking at multiple
#  sorts and negates differences in DNA synthesis.

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
import re


parser = argparse.ArgumentParser(description="")
parser.add_argument("pooled_counts.list", type=str, help="Same as estimate_affinity_from_ngs.py")
parser.add_argument("sorts.csv", type=str, help="Same as estimate_affinity_from_ngs.py")
parser.add_argument("affinity_estimate_out", type=str, help="The output from estimate_affinity_from_ngs.py. Your SSM mutations must be labeled"
                                                            +" such that they are <parent_name>__<seqpos>__<letter>. And the parent design should"
                                                            +" be labeled <parent_name>__native. See sorting_ngs_data/*_ssm/pooled_counts.list"
                                                            +" for examples.")
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

parser.add_argument("--entropy_pool_conc_factor", type=float, default=0.1, help="Which pool should we use for sequence entropy?"
                                                                        +" The default setting of 0.1 will try to pick the pool closest to 10-fold"
                                                                        +" lower than the parent Yeast-KD. 1.0 would pick the pool closest to the"
                                                                        +" Yeast-KD, and 10 10-fold above.")

parser.add_argument("--use_pseudo_entropies", action="store_true", help="Instead of true sequence entropy. Use the KD-estimates to produce"
                                                                        +" the expected number of counts.")


args = parser.parse_args(sys.argv[1:])


print("Loading pyrosetta")

from pyrosetta import *
from pyrosetta.rosetta import *
init("-mute all -beta_nov16")



# definition of KD
def kd_to_frac(kd, conc):
    return conc / (kd + conc)

def assign_rounds(sorts):
    # Assing round by starting from expression and moving forwards
    sorts['expression_parent'] = ""
    sorts.loc[sorts['special'] == 'naive', 'round'] = -1
    expression_mask = sorts['special'] == 'expression'
    sorts.loc[expression_mask, 'round'] = 0
    sorts.loc[expression_mask, 'expression_parent'] = sorts[expression_mask]['pool_name']
    missing = True
    iterr = 0
    while ( missing ):
        missing = False
        iterr +=1
        for idx, row in sorts.iterrows():
            parent = row['parent_pool']
            if ( parent == "" ):
                continue
            if ( np.isnan(sorts.loc[parent]['round']) ):
                missing = True
                continue
            sorts.loc[idx, 'round'] = sorts.loc[parent]['round'] + 1
            if ( row['expression_parent'] == "" ):
                sorts.loc[idx, 'expression_parent'] = sorts.loc[parent]['expression_parent']
            
        if ( iterr > 100 ):
            sys.exit("Could not assign rounds. Parent structure does not make sense. All pools must be derived from expression"
                    +" (With the exception of a naive pool)")


kd_df = pd.read_csv(args.affinity_estimate_out, sep="\s+")

if ( "target" in list(kd_df) ):
    target = kd_df['target'].iloc[0]
else:
    target = "none"
assert("kd_lb" in list(kd_df))
assert("kd_ub" in list(kd_df))
assert("lowest_conc" in list(kd_df))
assert("highest_conc" in list(kd_df))
assert("low_conf" in list(kd_df))
assert("avid_doesnt_agree" in list(kd_df))
assert("description" in list(kd_df))
kd_df = kd_df[['kd_lb', 'kd_ub', 'lowest_conc', 'highest_conc', 'low_conf', 'avid_doesnt_agree', 'description' ]].copy()

kd_df['description'] = kd_df['description'].str.replace("([^_])_([0-9]+)_([A-Z])$", r"\1__\2__\3")

sorts = pd.read_csv(args.__getattribute__("sorts.csv"), sep=",")
counts = pd.read_csv(args.__getattribute__("pooled_counts.list"), sep="\s+")
sorts.index = sorts['pool_name']

counts['description'] = counts['description'].str.replace("([^_])_([0-9]+)_([A-Z])$", r"\1__\2__\3")

# change nan to empty string
sorts['parent_pool'] = sorts['parent_pool'].fillna("")
sorts['special'] = sorts['special'].fillna("")
sorts['avidity'] = sorts['avidity'].fillna("")
sorts['notes'] = sorts['notes'].fillna("").astype(str)
sorts.loc[sorts['avidity'] == 'avi', 'avidity'] = 'avid'
sorts['round'] = np.nan
if ( "num_cells" in list(sorts) and "collected_cells" not in list(sorts)):
    sorts['collected_cells'] = sorts['num_cells']

if ( (sorts['special'] == 'expression').sum() == 0 ):
    sys.exit("You must have a value in sorts.csv where the special column == \"expression\"")

if ( ~np.all((sorts['avidity'] == 'avid') | (sorts['avidity'] == "") )):
    sys.exit("Valid choices for the avidity column are \"avid\" and \"\"")


# Give each sort a name and assign round 1, 2, 3, etc
assign_rounds(sorts)

# Only keep expression and non-avdity standard sorts
useful_mask = (sorts['avidity'] == "") & ((sorts['special'] == "") | (sorts['special'] == "expression"))

# If two sorts are at the same concentration. We want the later one.
sorts = sorts[useful_mask].sort_values('round')
sorts = sorts.drop_duplicates('concentration', keep='last')


# Remake counts but in terms of concentrations
new_df = pd.DataFrame(index=range(len(counts)))
for idx, row in sorts.iterrows():
    if ( row['special'] == 'expression'):
        continue
    
    name = "%.3f"%(row['concentration'])
    
    new_df[name + "_counts"] = counts[row['pool_name']]
    new_df[name + "_enrich"] = counts[row['pool_name']] / counts[row['expression_parent']]

new_df['description'] = counts['description']
counts = new_df

# Merge new counts into kd_df
to_print = "We have %i points with KD. %i with counts. And %%i with both."%(len(kd_df), len(counts))
kd_df = kd_df.merge(counts, 'inner', 'description')
print(to_print%len(kd_df))


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
for idx, row in kd_df[kd_df['description'].str.contains("_+native$")].iterrows():
    name = re.sub( "_+native$", "", row['description'])
    if ( name not in sequences ):
        continue
    seq = sequences[name]
    
    for il, let in enumerate(seq):
        seqpos = il + 1
        
        new_row = row.copy()
        new_row['native'] = True
        new_row['description'] = name + "__%i__%s"%(seqpos, let)
        parent_rows.append(new_row)

# del counts
# del sorts

extra_kd_df = pd.DataFrame(parent_rows)
    
# Add parent rows to the kd_df and then remove the originals
kd_df['native'] = False
kd_df = pd.concat((kd_df, extra_kd_df))
old_size = len(kd_df)
kd_df = kd_df.drop_duplicates(subset=['description'], keep='last')
if ( len(kd_df) != old_size ):
    print("!!!!!!!!!!!!!!!!!! Warning !!!!!!!!!!!!!!!!!")
    print("  There were duplicate names in the dataframe")
    print("  Either you truly do have duplicate names in your experimental data")
    print("  Or you have experimental datapoints for mutations that were native")
    print("  This can also trigger if you have the wrong sequence")
    print("")
kd_df = kd_df[~kd_df['description'].str.contains("_+native$")].reset_index().copy()

parts = kd_df['description'].str.extract(r"(.*[^_])_+([0-9]+)_+([A-Z])$")
kd_df['ssm_parent'] = parts[0]
kd_df['ssm_seqpos'] = parts[1].astype(str)
kd_df['ssm_letter'] = parts[2]




print("Concentrations we are calculating entropy from")

# Figure out what pool we're going to use for entropy calculations and then gather data
kd_df['chosen_counts'] = 0

for parent, subdf in kd_df.groupby('ssm_parent'):
    if ( parent not in sequences ):
        continue
    natives = subdf[subdf['native'].astype(bool)]
    if ( len(natives) == 0 ):
        print("  No native design for " + parent)
        continue
    parent_kd = np.clip(np.sqrt(natives['kd_ub'].mean() * natives['kd_lb'].mean()), natives['lowest_conc'].iloc[0]/10, 
                                                                                        natives['highest_conc'].iloc[0]*1e8)


    look_kd = parent_kd * args.entropy_pool_conc_factor


    if ( args.use_pseudo_entropies ):

        print("  %8.3f -- %s"%(look_kd, parent))

        # It doesn't matter what the total number is here, but we'll set it to 100,000
        total_cells = 100000

        # We use the upper-bound on the kd estimate

        expected_frac = kd_to_frac( subdf['kd_ub'].clip(natives['lowest_conc'].iloc[0]/10, 
                                                        natives['highest_conc'].iloc[0]*1e8), look_kd )
        expected_cells = total_cells * expected_frac


        kd_df.loc[subdf.index, 'chosen_counts'] = expected_cells

    else:
        # Find all the pool concentrations we could use to perform the entropy calculation.
        # Originally this script accepted multiple experiments, in practice, all the pools here will be valid
        options = []
        for name in list(subdf):
            if ( not name.endswith("_counts") ):
                continue
            if ( name.startswith("chosen")):
                continue
            valids = (~subdf[name].isnull()).sum()
            if ( valids / len(subdf) > 0.2 ):
                options.append(name)
                
        
        conces = np.array([float(x.split("_")[0]) for x in options])
        
        # Error is the ratio of the concentration to the look_concentration
        errors = conces / look_kd
        mask = errors < 1
        errors[mask] = 1 / errors[mask]
        
        smallest = np.argmin(errors)
        picked = options[smallest]
        my_ratio = conces[smallest] / parent_kd
        
        
        kd_df.loc[subdf.index, 'chosen_counts'] = subdf[picked]
        
        
        print("  %8.3f -- %s"%(conces[smallest], parent))
    
    


# Convert to xarray because it's way easier to do calcualtions
ds = pd.pivot_table(kd_df, index=['ssm_parent', 'ssm_letter', 'ssm_seqpos']).to_xarray()

# Find probability of each aa at each seqpos
with np.errstate(divide='ignore'):
    p_of_letter = ds['chosen_counts'].fillna(0) / ds['chosen_counts'].fillna(0).sum(dim='ssm_letter')

# Definition of Shannon Entropy
with np.errstate(divide='ignore'):
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

    if ( name not in sequences ):
        continue

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

if ( not args.use_pseudo_entropies ):
    print("")
    print("Also consider examining --use_pseudo_entropies . It may give less (or more) noise.")


