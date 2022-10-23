#!/usr/bin/env python

import os
import sys

from pyrosetta import *
from pyrosetta.rosetta import *

init()

input_name = sys.argv[1]
output_name = sys.argv[2]

pose = pose_from_pdb(sys.argv[1])

pose.pdb_info(core.pose.PDBInfo(pose))

pose.dump_pdb(output_name)
