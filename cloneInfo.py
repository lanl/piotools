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

    if verbose:
        print(f'nProcs={len(gns)}')
    # generate procID
    if verbose:
        print("Generating procid")
    id = 0
    procid = []
    last = int(0)
    for x in gns:
        x = int(x)
        procid.extend([id]*x)
        print(f'Proc {id:>6}: {last+1:>8} - {last+x:<8}')
        last = len(procid)
        id += 1
    
    numcell = p.numcell
    nbrs = np.empty([numcell, 2 * ndim],dtype='int')
    for i in range(ndim):
        if verbose:
            print("Generating neighbors in direction",i)
        idx = 2 * i + 1
        nl = p.readArray(b"cell_index_%1d"%(idx))
        nh = p.readArray(b"cell_index_%1d"%(idx+1))
        for k in range(numcell):
            nbrs[k,idx-1] = int(nl[k]-1)
            nbrs[k,idx] = int(nh[k]-1)

    if verbose:
        print("Neighbor generation done.")

    nClones = 0
    nMother = 0
    print('numcell=',p.numcell)
    for i in range(p.numcell):
        myProc = procid[i]
        if dtr[i] > 0:
            nMother += 1
            continue
        for x in nbrs[i]:
            if ( x < 1 or
                 x == i or
                 procid[x] != myProc):
                nClones += 1

    return {'numcell':numcell, 'nClones':nClones, 'nMother':nMother, 'nTop':numcell-nMother}
    

if __name__ == "__main__":
    import sys
    files = sys.argv[1:]
    for f in files:
        print(CloneInfo(f,1))
