from pio import pio

import sys

def readOneEntry(p, name, j):
    return p.readArrayRange(name,j,1)[0]
    
p = pio(sys.argv[1])

l1 = p.readArray(b"cell_level_0")
l0 = l1[0]
for i in range(p.numcell-1,0,-1):
    if l1[i] == l0:
        break

c = [0]*p.ndim
d = [0]*p.ndim
e = [0]*p.ndim
r = [0]*p.ndim
offsets = [1, 2, 4]
for dim in range(p.ndim):
    name = b"cell_center_%1d"%(dim+1)
    c[dim] = readOneEntry(p, name, 0)
    d[dim] = readOneEntry(p, name, offsets[dim])
    e[dim] = readOneEntry(p, name, i)
    r[dim] = 1+round((e[dim] - c[dim])/(d[dim] - c[dim]))

print(f"Root mesh size = {r}")
