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

 Script: addVar.py
 Author: Sriram Swaminarayan (sriram@lanl.gov)
   Date: March 17, 2021
Version: 1.1
  Usage: python3 addProcID.py <input-dmp######> [outbase=tmp]

Input: xRage dump file
Output:
  <outbase>.pio:  file that can be opened by Paraview
  <outbase>-dmp000000: input file with processor_id array added

What it does:
Appends a new cell array called 'processor_id' to a PIO file. that
contains the MPI processor ID on which that cell resided

"""
from __future__ import print_function
from pio import pio
import numpy as np
import os
import sys

def usage():
    print(f"""
    **************************************************************
    Usage: python3 {sys.argv[0]} <input-dmp######> [outbase=tmp]
    **************************************************************

      Input: xRage dump file
      Output:
        <outbase>.pio:  file that can be opened by Paraview
        <outbase>-dmp000000: input file with processor_id array added

    What it does:
      Appends a new cell array called 'processor_id' to a PIO file. 
      thatcontains the MPI processor ID on which that cell resided
    **************************************************************
    """)


if (len(sys.argv) > 3) or (len(sys.argv) <2):
      usage()
else:
    # input file
    filename = sys.argv[1]

    # output file
    try:
        outbase = sys.argv[2]
    except:
        outbase = 'tmp'


    # Read in the dump file metadata
    p = pio(filename)
    
        
    # Create processor ID array
    procID = np.zeros(p.numcell,dtype='int')
    print('numcell=',p.numcell)
    data = p.readArray(b'global_numcell_0')
    dtr = p.readArray(b'cell_daughter_0')
    lastIdx = 0
    ncell = 0
    id = 0
    delta = p.numcell/10;
    print("processing PE info")
    for i in range(len(data)):
        nextIdx = lastIdx + int(data[i])
        # procID[lastIdx:nextIdx] = int(i)
        for j in range(int(data[i])):
            if dtr[id] == 0:
                procID[ncell] = int(i)
                ncell += 1
            id += 1
        lastIdx = nextIdx
        print(f"\rProcessed  %{100.*float(lastIdx)/float(p.numcell):.1f} done", end="")
    print("")

    print(f"  Writing {ncell} entries to disk")
    #with open('/tmp/xxx.bin','wb') as fp:
    np.resize(procID,ncell)
    np.savez_compressed('/tmp/proc.npz', procID=procID)
    xxx = np.load('/tmp/proc.npz')
    print('written cells = ',xxx['procID'].size)
