#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-



# Standard library
import os
import subprocess



# Python interpreter to be used
interpreter = "python3.7"
# Script to perform the contact analysis
script = "plot_contacts.py"

# Directory where to store the plots
plot_dir = os.path.join(os.getcwd(), "plots")
os.makedirs(plot_dir, exist_ok = True)

# Path to the PDBs (= reference structures)
pdb_path = "reference_structures/1F41_altA{:s}_processed.pdb"
# Path to the matrices of contacts for the reference structures
matrix_pdb_path = "matrices/{:s}_pdb.dat"
# Path to the matrices of contacts for the trajectories
matrix_traj_path = "matrices/{:s}_traj.dat"

# Names of the systems in the PDBs
pdb_systems = ["", "_K35N", "_K35T", "_R34G", "_R34T"]

# Names of the systems in the output PDF files
names = ["wt", "K35N", "K35T", "R34G", "R34T"]

# List of PDBs
pdbs = [pdb_path.format(s) for s in pdb_systems]
# List of matrices for the reference structures
matrices_pdb = [matrix_pdb_path.format(s) for s in names]
# List of matrices for the trajectories
matrices_traj = [matrix_traj_path.format(s) for s in names]

# Color map to be used in the heatmaps
cmap = "Greens"

# Font to be used for the text elements of the heatmaps
font = "fonts/Helvetica.ttf"

# Regions of the system whose contacts will be plotted
regions = {"chainA_36-40" : (26, 30), 
           "chainA_55-64" : (45, 54),
           "chainB_36-40" : (142, 146),
           "chainB_55-64" : (161, 170)}

# For each system
for pdb, mpdb, mtraj, name in \
    zip(pdbs, matrices_pdb, matrices_traj, names):

    # For each region
    for region_name, region in regions.items():

        # Convert the region definition into a string
        region = ",".join([str(i) for i in region])

        # Set the output plots' names
        op = os.path.join(plot_dir, f"{name}_pdb_{region_name}.pdf")
        ot = os.path.join(plot_dir, f"{name}_traj_{region_name}.pdf")

        # Plot the heatmap for the contacts in the selected
        # region in the reference structure
        subprocess.run([interpreter,
                        script,
                        "-m", mpdb,
                        "-r", pdb,
                        "-o", op,
                        "--cmap", cmap,
                        "--font", font,
                        "--region", region])

        # Plot the heatmap for the contacts in the selected
        # region in the trajectories
        subprocess.run([interpreter,
                        script,
                        "-m", mtraj,
                        "-r", pdb,
                        "-o", ot,
                        "--cmap", cmap,
                        "--font", font,
                        "--region", region])
