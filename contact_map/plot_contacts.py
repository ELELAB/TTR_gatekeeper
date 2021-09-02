#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-



# Third party packages
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns



# Mapping between the three-letters name and the one-letter
# name for the 20 canonical amino acids
three2one = \
    {"ALA" : "A", "CYS" : "C", "ASP" : "D", "GLU" : "E",
     "PHE" : "F", "GLY" : "G", "HIS" : "H", "ILE" : "I",
     "LYS" : "K", "LEU" : "L", "MET" : "M", "ASN" : "N",
     "PRO" : "P", "GLN" : "Q", "ARG" : "R", "SER" : "S",
     "THR" : "T", "VAL" : "V", "TRP" : "W", "TYR" : "Y"} 


def plot_contacts_heatmap(df, region, out, cmap, font):

    # Get the font properties to use for the various text
    # elements    
    fp_cbarlabel = fm.FontProperties(fname = font, size = 12)
    fp_cbarticklabels = fm.FontProperties(fname = font, size = 10)
    fp_axlabels = fm.FontProperties(fname = font, size = 12)
    fp_axticklabels = fm.FontProperties(fname = font, size = 10)

    # If a specific region needs to be plotted
    if region is not None:
        # Get the starting and ending points of the region
        start, end = region
        # Take only the portion of the dataframe corresponding
        # to the region of interest
        df = df.iloc[start:end+1, start:end+1]

    # Generate the heatmap
    ax = sns.heatmap(df,
                     cbar = False,
                     vmin = 0.0,
                     vmax = 1.0,
                     center = 0.5,
                     cmap = cmap,
                     square = True,
                     linecolor = "black",
                     linewidths = 0)

    # Generate the color bar
    cbar = plt.colorbar(mappable = ax.get_children()[0],
                        extend = "both",
                        extendrect = False,
                        orientation = "vertical",
                        pad = 0.2)

    # Set the color bar label
    cbar.ax.set_ylabel(ylabel = "Contact persistence",
                       fontproperties = fp_cbarlabel)

    # Set the font properties for the color bar tick
    # labels
    for t in cbar.ax.get_yticklabels():
        t.set_font_properties(fp_cbarticklabels)

    # Set the x-axis' label
    ax.set_xlabel(xlabel = "Residues",
                  fontproperties = fp_axlabels)
    
    # Set the x-axis' tick labels
    ax.set_xticklabels(labels = ax.get_xticklabels(),
                       fontproperties = fp_axticklabels)

    # Set the x-axis' tick parameters
    ax.tick_params(axis = "x",
                   which = "both",
                   length = 0)

    # Set the y-axis' label
    ax.set_ylabel(ylabel = "Residues",
                  fontproperties = fp_axlabels)
    
    # Set the y-axis' tick labels
    ax.set_yticklabels(labels = ax.get_yticklabels(),
                       fontproperties = fp_axticklabels)

    # Set the y-axis' tick parameters
    ax.tick_params(axis = "y",
                   which = "both",
                   length = 0)

    # Save the figure to the PDF output file
    plt.savefig(out,
                dpi = 300,
                transparent = True,
                bbox_inches = "tight")
            
    # Clear the figure
    plt.clf()
            
    # Close the current figure window
    plt.close()



if __name__ == "__main__":
    
    import argparse

    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Add the arguments
    m_help = "Input matrix."
    parser.add_argument("-m", 
                        type = str,
                        required = True,
                        help = m_help)

    r_help = "Input reference structure."
    parser.add_argument("-r",
                        type = str,
                        required = True,
                        help = r_help)

    o_help = "Output PDF file."
    parser.add_argument("-o",
                        type = str,
                        required = True,
                        help = o_help)

    cmap_default = "PuOr"
    cmap_help = f"Color map to be used. Default is: {cmap_default}."
    parser.add_argument("--cmap",
                        type = str,
                        default = cmap_default,
                        help = cmap_help)

    font_default = None
    font_help = \
        "Font to be used for text elements. Default is: the default " \
        "matplotlib font."
    parser.add_argument("--font",
                        type = str,
                        default = font_default,
                        help = font_help)

    region_default = None
    region_help = \
        "Residue indexes of the region to be plotted. Default is: " \
        "contacts for the whole system are plotted."
    parser.add_argument("--region",
                        type = str,
                        default = region_default,
                        help = region_help)

    # Parse the arguments
    args = parser.parse_args()

    # Get the region of interest, if provided
    region = \
        [int(i) for i in args.region.split(",")] if args.region \
        is not None else args.region

    # Get the residues in the reference structure
    residues = mda.Universe(args.r).residues

    # Format the residues' names
    resnames = [f"{three2one[r.resname]}{r.resnum}" for r in residues]

    # Load the matrix of contacts
    matrix = np.loadtxt(args.m)

    # Create a dataframe from the matrix of contacts, where
    # rows and columns are named after the residues
    df = pd.DataFrame(data = matrix,
                      index = resnames,
                      columns = resnames)

    # Plot the contacts' heatmap
    plot_contacts_heatmap(df = df,
                          region = region,
                          out = args.o,
                          cmap = args.cmap,
                          font = args.font)