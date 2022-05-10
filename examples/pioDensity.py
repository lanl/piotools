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

# Read in mass
mass = p.readArray('mass_0')

# read in volume
vcell = p.readArray('vcell_0')

rho = mass/vcell

# print the density for first 10 cells
print(rho[1:10])
