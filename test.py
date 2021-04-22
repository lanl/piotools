import sys
import numpy as np
from pio import pio
from math import floor

# load the PIO file
p = pio(sys.argv[1])

# print numdim and numcell
print(p.ndim)
print(p.numcell)

# print all variables
for x in p.names:
    print(x)

# print one variable, note 'b' is REQUIRED before string to make it bytearray
mLen = p.readArrayInt(b"MATIDENT_LEN_0")
print(mLen)

m1 = p.readArray(b"matident_0")
print(m1)
print("")
mdf1 = p.readArray(b"matdef_1")
print(mdf1)
print("")
mdf2 = p.readArray(b"matdef_2")
print(mdf2)
