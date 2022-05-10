#!/usr/bin/env python3
"""
========================================================================================
 (C) (or copyright) 2022. Triad National Security, LLC. All rights reserved.

 This program was produced under U.S. Government contract 89233218CNA000001 for Los
 Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
 for the U.S. Department of Energy/National Nuclear Security Administration. All rights
 in the program are reserved by Triad National Security, LLC, and the U.S. Department
 of Energy/National Nuclear Security Administration. The Government is granted for
 itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
 license in this material to reproduce, prepare derivative works, distribute copies to
 the public, perform publicly and display publicly, and to permit others to do so.
========================================================================================

Explores a pio file.
Usage: python3 pioExplore.py <dumpFile>

By default prints a list of variables.
Uncomment different lines to do different things

"""
import numpy as np

def usage():
    import sys

    print(f"""
    
    Usage:  {sys.argv[0]} pio_file 

    Set blocks in code to True / False to see different capabilities

    """)

import sys
from pio import *

if len(sys.argv) < 2:
    # Not enough arguments
    usage()
    sys.exit(1)

#####################################
########## START HERE ###############
#####################################

# Read file (first argument)
try:
    p = pio(sys.argv[1])    
    print(f'numcell={p.numcell}')
except:
    print(sys.argv)
    print("\n\n   Unable to read file")
    usage()
    sys.exit(1)

# Activate blocks depending on what you want to do

# Print list of variables, default action
if True:
    for name in p.names:
        print(name)
    print("")

# Print elements of class PIO
if False: 
    print(dir(p))

# A more complex example that will print the values of
# top level cells whose centers satisfy the conditions:
#    (cx > 1.5) && (cx <= 1.6) && (cy > 1.5) && (cy <= 1.6)
#
# where (cx, cy) are the coordinates of the centers of
# given cell
if True:
    # Steps:
    #  Initialize PIO to handle fractional variables (starting with chunk_)
    #  Print a subset of cells in a piece of space
    #  Generate a mask of ltop cells
    #  Read in the center array
    #  Generate a mask from conditionals on coordinates
    #  Create an index list from mask
    #  index variables by this index list
    
    # Initialize the material variable expansion keys
    # Without this you cannot access fractional variables
    # i.e. the ones that start with "chunk_"
    p.updateCsrIndices("chunk_nummat", "chunk_mat", "vcell", 1)
    
    # Generate ltopMask and ltop array from daughter array
    ltopMask = p.readArray('cell_daughter_0') == 0
    ltop = np.argwhere(ltopMask)[:,0]
    print(f'numtop={len(ltop)} (number of top level cells) ')
    print(f'nummat={p.csrN} (number of materials)')

    # Read in centers of the cells by looping over p.ndim
    # multidimensional arrays in dump files are indexed from 1
    centers = np.array([p.readArray(f'cell_center_{iDim}') for iDim in range(1,p.ndim+1)])

    # define spatial limits and generate index list
    xmin = 0.03; xmax = 0.53
    ymin = 0.03; ymax = 0.53

    # Generate Mask with comparisons
    compareMask = np.logical_and(centers[0][:]>xmin, centers[0][:]<=xmax)
    compareMask = np.logical_and(compareMask, centers[1][:]>ymin)
    compareMask = np.logical_and(compareMask, centers[1][:]<=ymax)
    
    # Mask with only ltop cells
    compareMask = np.logical_and(compareMask, ltopMask)

    # generate the index array for our subset of space
    indices = np.argwhere(compareMask)[:,0]

    # print first ten entries in the index list
    print(f"number of cells that match: {len(indices)}")
    print(f"Indices of match={indices}")

    # Read in volume and pressure and print for our subset of space
    # by indexing with the indices array
    volume = p.readArray('vcell_0')[indices]
    pressure = p.readArray('pres_0')[indices]
    fractional_vol = p.expandCsrVariable('chunk_vol_0', True)
    print("")
    print("-"*83)
    print(f'{"i":>6} {"index":>6} {"cx":>8} {"cy":>8} {"volume":>12} '
          f'{"pressure":>12} {f"fvol_1":>12} {f"fvol_2":>12}')
    print("-"*83)
    for i, index in enumerate(indices):
        # We index centers and fractional volumes with "index" because they are full arrays
        cx = centers[0][index]
        cy = centers[1][index]
        fvol_1 = fractional_vol[0][index]
        fvol_2 = fractional_vol[1][index]
        
        # We index volume and pressure with "i" because they are reduced arrays
        vol = volume[i]
        pres = pressure[i]

        # print the values
        print(f'{i:>6} {index:>6} {cx:>8.4} {cy:>8.4} {vol:>12.4} {pres:>12.4} {fvol_1:>12.4} {fvol_2:>12.4}')
