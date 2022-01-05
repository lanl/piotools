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
Uses pio.py to print variables

"""


def usage():
    import sys

    f"""
    
    Usage:  {sys.argv[0]} pio_file var1 [ var_2 var_3 ...]

    """


import sys
from pio import *

if len(sys.argv) < 3:
    usage()
    sys.exit(1)

p = pio(sys.argv[1])
for n in sys.argv[2:]:
    if n == "--list":
        for name in p.names:
            print(name)
    found = 0
    for name in p.names:
        if name.startswith(n):
            c = p.readArray(name)
            print(f"{name} = {c}")
            found += 1
    if found == 0:
        print(f"    *** No instances of {n} were found")
