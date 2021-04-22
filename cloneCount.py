#!/usr/bin/env python3
"""
========================================================================================
 (C) (or copyright) 2021. Triad National Security, LLC. All rights reserved.

 This program was produced under U.S. Government contract 89233218CNA000001 for Los
 Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
 for the U.S. Department of Energy/National Nuclear Security Administration. All rights
 in the program are reserved by Triad National Security, LLC, and the U.S. Department
 of Energy/National Nuclear Security Administration. The Government is granted for
 itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
 license in this material to reproduce, prepare derivative works, distribute copies to
 the public, perform publicly and display publicly, and to permit others to do so.
========================================================================================

Takes two arguments: PIO file and numprocs

Use the PIO class to estimate number of
clone cells for a simulation

Uses a slow but reliable algorithm.

** REQUIRES Python3 **

Author: Sriram Swaminarayan (sriram@lanl.gov)
Date: March 24, 2021
Version: 1.0

"""
import sys
import numpy as np
from pio import pio
from math import floor

if len(sys.argv) != 3:
    print(f"\n\n   Usage: {sys.argv[0]} <filename> <nprocs> \n\n")
    exit(1)

p = pio(sys.argv[1])
nprocs = int(sys.argv[2])
ndim = p.ndim

dtr = p.readArray(b"cell_daughter_0")

numcell = p.numcell
nbrs = np.empty([numcell, 2 * ndim])
for i in range(ndim):
    print("processing neighbors in direction", i)
    idx = 2 * i + 1
    nl = p.readArray(b"cell_index_%1d" % (idx))
    nh = p.readArray(b"cell_index_%1d" % (idx + 1))
    for k in range(numcell):
        nbrs[k, idx - 1] = nl[k] - 1
        nbrs[k, idx] = nh[k] - 1

print("counting clones")
block = int(2 ** ndim)
quantum = round(float(numcell) / float(nprocs))
remainder = quantum % block
quantum = quantum - remainder

nEnd = 0
nClones = 0
r = 0
n = [0] * (2 * ndim)
print(f"quantum={quantum}, remainder={remainder}, block={block}")
for i in range(nprocs):
    nStart = nEnd
    if i == nprocs - 1:
        nEnd = numcell
    else:
        r += remainder
        if r >= block:
            offset = block
            r = r - block
        else:
            offset = 0
        nEnd = nStart + quantum + offset
    nEnd = int(nEnd)
    print(f"checking proc={i},s={nStart} e={nEnd} n={nEnd-nStart}")
    for j in range(nStart, nEnd):
        if dtr[j] > 0:
            continue
        n = nbrs[j]
        for k in n:
            if k > nEnd or k < nStart or k == j:
                nClones += 1

print(f"\n  nProcs={nprocs}, nClones={nClones}\n")

#
# comment out next line if you want to check how well
# our partitioning tracks xRage's partitionng.
#
# print(p.readArray(b"global_numcell_0"))
