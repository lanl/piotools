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

    # Ensure that outfile doesn't already exist
    outfile = outbase + '-dmp000000'
    if os.path.exists(outfile):
        print(f"\n\n  *** File {outfile} exists, will not proceed. ***\n\n")
        exit(1)

    # Read in the dump file metadata
    p = pio(filename)
    
    # Write the Paraview file
    with open(outbase+'.pio','w') as ofp:
        ofp.writelines(
            ["DUMP_DIRECTORY .\n",
             f"DUMP_BASE_NAME {outbase}\n"]
            )
        
    # Create processor ID array
    newData = np.zeros(p.numcell,dtype='double')
    data = p.readArray(b'global_numcell_0')
    lastIdx = 0
    for i in range(len(data)):
        nextIdx = lastIdx + int(data[i])
        newData[lastIdx:nextIdx] = i
        lastIdx = nextIdx

    # Write new file with processor id included
    p.writeWithNewCellArray(outfile, "processor_id", newData)

            
