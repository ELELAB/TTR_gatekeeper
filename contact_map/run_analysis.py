#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-



# Standard library
import os
import subprocess
# Third-party package
import MDAnalysis as mda
from MDAnalysis.analysis import contacts

# Python interpreter to be used
interpreter = "python3.7"
# Script to perform the contact analysis
script = "analyze_contacts.py"

# Directory where to store the matrices
matrix_dir = os.path.join(os.getcwd(), "matrices")
os.makedirs(matrix_dir, exist_ok = True)

# Path to the topologies
top_path = "/data/user/shared_projects/ttr_gatekeep/productive/{:s}/9-md/md_prot.tpr"
# Path to the trajectories
traj_path = "/data/user/shared_projects/ttr_gatekeep/productive/{:s}/9-md/traj_centered_150ns.xtc"
# Path to the PDBs
pdb_path = "reference_structures/1F41_altA{:s}_processed.pdb"

# Names of the systems in the PDBs
pdb_systems = ["", "_K35N", "_K35T", "_R34G", "_R34T"]
# Names of the systems in the trajectories
traj_systems = ["wt", "K35N", "K35T", "R34G", "R34T"]

# List of PDBs
pdbs = [pdb_path.format(s) for s in pdb_systems]
# List of topologies
tops = [top_path.format(s) for s in traj_systems]
# List of trajectories
trajs = [traj_path.format(s) for s in traj_systems]

# String for atom selection
sel = "protein_and_not_backbone"
# Distance cut-off for two atoms to be in contact
dist_cut = 4.5
# Sequence cut-off for two residues to be in contact
seq_cut = 1
# Number of atom pairs that need to be in contact for two
# residues to be considered in contact
atom_pairs = 1

# For each system
for top, traj, pdb, traj_sys in zip(tops, trajs, pdbs, traj_systems):

    # Set the output matrices' names
    op = os.path.join(matrix_dir, f"{traj_sys}_pdb.dat")
    ot = os.path.join(matrix_dir, f"{traj_sys}_traj.dat")

    # Perform the contact analysis on the starting structure
    # and write out the resulting matrix
    subprocess.run([interpreter,
                    script,
                    "-f", pdb,
                    "-s", pdb,
                    "-o", op,
                    "--sel", sel,
                    "--dist-cut", str(dist_cut),
                    "--seq-cut", str(seq_cut),
                    "--atom-pairs", str(atom_pairs)])

    # Perform the contact analysis on the trajectory 
    # and write out the resulting matrix
    subprocess.run([interpreter,
                    script,
                    "-f", traj,
                    "-s", top,
                    "-o", ot,
                    "--sel", sel,
                    "--dist-cut", str(dist_cut),
                    "--seq-cut", str(seq_cut),
                    "--atom-pairs", str(atom_pairs)])

    