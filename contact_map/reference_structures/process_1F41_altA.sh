#!/bin/bash

in_pdb=../../pdbs/1F41_altA.pdb
out_pdb=1F41_altA_processed.pdb

python remove_heteroatoms.py -i $in_pdb -o $out_pdb
