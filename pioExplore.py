#!/usr/bin/env python3
"""
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
