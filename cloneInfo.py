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
Returns information about an AMHC dump.
Return value is a dictionary with keys:
 numcell: number of cells in the simulation
 nClones: number of clone cells
 nMother: number of mother cells
    nTop: number of top level cells.

Note that the number of real cells in the simulation is nTop and the number
of additional cells is nMother + nClones

Author: Sriram Swaminarayan (sriram@lanl.gov)
Date: April 22, 2021
Version: 1.0

"""


def countClones(cellsPerProc, p, verbose=False):
    """
    Uses the passed global_numcell instead of the one from the pio file
    to generate clones.
      cellsPerProc = array that is nProcs long with number of cells on each processor
      p = PIO file class
    """
    import numpy as np

    ndim = p.ndim

    # read in relevant arrays
    dtr = p.readArray(b"cell_daughter_0")

    if verbose:
        print(f"nProcs={len(cellsPerProc)}")
    # generate procID
    if verbose:
        print("Generating procid")
    id = 0
    procid = []
    last = int(0)
    for x in cellsPerProc:
        x = int(x)
        procid.extend([id] * x)
        print(f"Proc {id:>6}: {last+1:>8} - {last+x:<8}")
        last = len(procid)
        id += 1

    nClones = 0
    nMother = 0
    numcell = p.numcell
    print("numcell=", p.numcell)
    for iDir in range(ndim):
        if verbose:
            print("Counting Clones in direction", iDir)
        idx = 2 * iDir + 1
        nl = p.readArray(b"cell_index_%1d" % (idx)).astype(int) - 1
        nh = p.readArray(b"cell_index_%1d" % (idx + 1)).astype(int) - 1
        for i in range(numcell):
            if dtr[i] > 0:
                if iDir == 0:
                    # don't double count
                    nMother += 1
                continue
            myProc = procid[i]

            xL = nl[i]
            xH = nh[i]
            if xL < 0 or xL == i or procid[xL] != myProc:
                nClones += 1
            if xH < 0 or xH == i or procid[xH] != myProc:
                nClones += 1

    return {
        "numcell": numcell,
        "nClones": nClones,
        "nMother": nMother,
        "nTop": numcell - nMother,
    }


def CloneInfo(f, verbose=False):
    import numpy as np
    from pio import pio
    from math import floor

    # open the file
    p = pio(f)
    ndim = p.ndim

    # read in relevant arrays
    dtr = p.readArray(b"cell_daughter_0")
    gns = p.readArray(b"global_numcell_0")
    if gns is None:
        gns = [p.numcell]
    return countClones(gns, p, verbose=False)


if __name__ == "__main__":
    import sys

    files = sys.argv[1:]
    for f in files:
        print(CloneInfo(f, 0))
