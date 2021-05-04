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


def cloneCountP(p, nprocs, verbose=False):
    import numpy as np
    from cloneInfo import countClones

    ndim = p.ndim
    numcell = p.numcell
    block = int(2 ** ndim)
    quantum = round(float(numcell) / float(nprocs))
    remainder = quantum % block
    quantum = quantum - remainder

    nEnd = 0
    nClones = 0
    r = 0
    n = [0] * (2 * ndim)
    gns = [0] * nprocs
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
        gns[i] = nEnd - nStart
    return countClones(gns, p, verbose)

def cloneCount(fname, nprocs, verbose=False):
    from pio import pio
    p = pio(fname)
    return cloneCountP(p, nprocs, verbose)

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print(f"\n\n   Usage: {sys.argv[0]} <filename> <nprocs> \n\n")
        exit(1)
    myFile = sys.argv[1]
    myNProcs = int(sys.argv[2])

    print(cloneCount(myFile, myNProcs, 0))
