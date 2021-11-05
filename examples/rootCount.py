def getRootSize(p):
    """
    Returns the root mesh size of pio file p
    """
    l1 = p.readArray("cell_level_0")
    l0 = l1[0]
    for i in range(p.numcell - 1, 0, -1):
        if l1[i] == l0:
            break

    c = [0] * p.ndim
    d = [0] * p.ndim
    e = [0] * p.ndim
    r = [0] * p.ndim
    offsets = [1, 2, 4]
    for dim in range(p.ndim):
        name = "cell_center_%1d" % (dim + 1)
        c[dim] = p.readArrayRange(name, 0, 1)[0]
        d[dim] = p.readArrayRange(name, offsets[dim], offsets[dim] + 1)[0]
        e[dim] = p.readArrayRange(name, i, i + 1)[0]
        r[dim] = 1 + round((e[dim] - c[dim]) / (d[dim] - c[dim]))

    return r


if __name__ == "__main__":
    import sys
    from pio import pio

    p = pio(sys.argv[1])
    print(f"Root mesh size = {getRootSize(p)}")
