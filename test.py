import sys
import numpy as np
from pio import pio
from math import floor

# load the PIO file
p = pio(sys.argv[1])


# print numdim and numcell
print(p.ndim)
print(p.numcell)

#print all variables
for x in p.names:
    if x.startswith(b'matdef'): print(x)

#print one variable, note 'b' is REQUIRED before string to make it bytearray
m1 = p.readArray(b"matdef_1")
print(m1)
