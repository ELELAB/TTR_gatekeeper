#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-



# Standard library
import sys
# Third-party packages
import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import numpy as np



class ContactAnalysis:

    """Perform a residue-residue contact analysis on a trajectory.
    """

    def __init__(self, sel, dist_cut, seq_cut, atom_pairs):

        # String representing the selection of atoms to
        # be included in the contact analysis
        self.sel = sel

        # Distance cutoff - the distance between two atoms must
        # be under this value for the atoms to be in contact
        self.dist_cut = dist_cut

        # Sequence cutoff - pairs of residues within this sequence
        # distance will not be considered for contact analysis
        self.seq_cut = seq_cut

        # Number of atom pairs that must be in contact for two
        # residues to be considered in contact
        self.atom_pairs = atom_pairs

        # Matrix containing the residue-residue contacts'
        # frequencies over the whole trajectory 
        self.matrix = None


    def _get_nb_atoms_per_res(self, sel):
        """Get the neighboring atoms of each residue in the
        system as an UpdatingAtomGroup.
        """
        
        # Split the atom selection by residue    
        splitted_byres = sel.atoms.split(level = "residue")
        
        # Generate a list of length equal to len(splitted_byres) 
        # with tuples containing each an AtomGroup instance with the
        # atoms of each residue and an UpdatingAtomGroup with the
        # neighoring atoms of each residue.
        # Use no periodic boundary conditions.
        nb_atoms_per_res = []
        for atoms_resi in splitted_byres:
            nb_atoms = \
                sel.select_atoms(f"around {self.dist_cut} group sg ",
                                 sg = atoms_resi,
                                 updating = True,
                                 periodic = False)
            
            # Append the atoms belonging to the residue and the
            # neighboring atoms to the list
            nb_atoms_per_res.append((atoms_resi, nb_atoms))
        
        # Return the list
        return nb_atoms_per_res


    def _compute_contacts(self, nb_atoms_per_res, matrix_shape):
        """Compute contacts for a single frame.
        """

        # Initialize the matrix of contacts
        matrix = np.ones(shape = matrix_shape, dtype = np.float64) * \
                 -(np.inf)
        

        #------------------------- Residue 'i' -----------------------#


        # For each tuple of (atoms of residue 'i', 
        # neighboring atoms or residue 'i')
        for atoms_resi, nb_atoms in nb_atoms_per_res:
            
            # Get the residue and segment (= chain) index
            resindex_i = atoms_resi[0].residue.ix
            segindex_i = atoms_resi[0].segment.ix
     
            # Group the neighboring atoms by residue
            nb_atoms_byres = nb_atoms.split(level = "residue")


            #----------------------- Residue 'i' ---------------------#


            # For each group of neighboring atoms (each group
            # belonging to a different residue 'j')
            for atoms_resj in nb_atoms_byres:
                
                # Get the residue and segment (= chain) index
                resindex_j = atoms_resj[0].residue.ix
                segindex_j = atoms_resj[0].segment.ix


                #------------------ Sequence cut-off ----------------#


                # Ignore residue 'j' if it is within the sequence
                # cutoff with respect to residue 'i' and they belong 
                # to the same segment (= chain)
                is_same_segment = (segindex_i == segindex_j)
                is_within_seq_cut = \
                    resindex_j >= (resindex_i - self.seq_cut) and \
                    resindex_j <= (resindex_i + self.seq_cut)

                # If the two residues are within the sequence cutoff
                # and belong to the same segment, skip them
                if is_same_segment and is_within_seq_cut:
                    continue


                #------------------- Distance array ------------------#


                # Do not compute twice the pairs for each couple
                # of residues 'i' and 'j', therefore compute only
                # if the symmetric cell is still empty
                if matrix[resindex_j, resindex_i] == -(np.inf):
                    d_array = dist.distance_array(atoms_resi.positions,
                                                  atoms_resj.positions)
                    
                    # Compute the number of atom pairs between the two
                    # residues having a distance lower or equal than
                    # the distance cutoff
                    atom_pairs = np.sum(d_array <= self.dist_cut)
                    
                    # If there are sufficient atom pairs in contact
                    # for the residues to be in contact, save the
                    # contact
                    if atom_pairs >= self.atom_pairs:
                        # Update the matrix
                        matrix[resindex_i,resindex_j] = 1

        # Symmetrize the matrix
        matrix = np.maximum(matrix, matrix.transpose())

        # Convert all negative infinites to 0.0
        matrix[matrix == -(np.inf)] = 0.0

        # Return the matrix
        return matrix


    def _iter_frames(self, u, sel, matrix_shape):
        """Return a generator iterating over all frames of the
        trajectory and computing the contacts for each of them.
        """

        # Get the neighboring atoms per residue
        nb_atoms_per_res = self._get_nb_atoms_per_res(sel = sel)
        
        # For each frame
        for fs in u.trajectory:
            
            # Log the progress in analyzing frames
            sys.stdout.write(f"\rAnalyzing frame: {fs.frame}")
            sys.stdout.flush()
            
            # Yield the matrix of contacts
            yield self._compute_contacts(\
                    matrix_shape = matrix_shape,
                    nb_atoms_per_res = nb_atoms_per_res)        


    def run(self, u):
        """Run the contact analysis.
        """

        # Select the atoms based on the selection string
        sel = u.select_atoms(self.sel)
 
        # Set the shape of the final matrix (# residues x # residues)
        matrix_shape = (len(u.residues), len(u.residues))
        
        # Create the iterator for the matrices of contacts for
        # each frame
        matrices = \
            self._iter_frames(u = u,
                              sel = sel,
                              matrix_shape = matrix_shape)  
        
        # Set the average matrix to be the first matrix
        avg_matrix = next(matrices)
        
        # Set a counter for storing the total number of matrices
        tot_matrices = 1
        
        # For each matrix
        for i, matrix in enumerate(matrices):
            # Sum the current matrix to the average one
            avg_matrix = np.add(avg_matrix, matrix)
            # Increament the counter of the total matrices
            tot_matrices += 1
        
        # Save the matrix of contacts' frequencies
        self.matrix = avg_matrix / tot_matrices


if __name__ == "__main__":
    
    import argparse

    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Add the arguments
    f_help = "Input trajectory."
    parser.add_argument("-f", 
                        type = str,
                        required = True,
                        help = f_help)

    s_help = "Input topology."
    parser.add_argument("-s",
                        type = str,
                        required = True,
                        help = s_help)

    o_help = "Output matrix."
    parser.add_argument("-o",
                        type = str,
                        required = True,
                        help = o_help)

    sel_default = "protein_and_not_backbone"
    sel_help = f"String for atom selection (default: {sel_default})."
    parser.add_argument("--sel",
                        type = str,
                        default = sel_default,
                        help = sel_help)

    dist_cut_default = 4.5
    dist_cut_help = f"Distance cut-off (default : {dist_cut_default})."
    parser.add_argument("--dist-cut",
                        type = float,
                        default = dist_cut_default,
                        help = dist_cut_help)

    seq_cut_default = 1
    seq_cut_help = f"Sequence cut-off (default : {seq_cut_default})."
    parser.add_argument("--seq-cut",
                        type = int,
                        default = seq_cut_default,
                        help = seq_cut_help)

    atom_pairs_default = 1
    atom_pairs_help = \
        f"Number of atom pairs in contact for two residues to " \
        f"be considered in contact (default: {atom_pairs_default})."
    parser.add_argument("--atom-pairs",
                        type = int,
                        default = atom_pairs_default,
                        help = atom_pairs_help)

    # Parse the arguments
    args = parser.parse_args()
    sel = args.sel.replace("_", " ")

    # Create the Universe object
    u = mda.Universe(args.s, args.f)
    
    # Initialize the contact analysis
    ca = ContactAnalysis(sel = sel,
                         dist_cut = args.dist_cut,
                         seq_cut = args.seq_cut,
                         atom_pairs = args.atom_pairs)

    # Run the contact analysis
    ca.run(u)

    # Save the matrix of contacts' frequencies
    np.savetxt(args.o, ca.matrix)



