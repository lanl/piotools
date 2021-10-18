#!/usr/bin/env python3


def getRiotInfo(fname):
    import h5py

    f = h5py.File(fname)
    info = f["/Info"]
    time = info.attrs["Time"]
    walltime = info.attrs["WallTime"]
    cycle = info.attrs["NCycle"]
    try:
        mesh = f["/Mesh"]
        nb = mesh.attrs["nbtotal"]
        nblock = mesh.attrs["blockSize"]
    except:
        nb = info.attrs["NumMeshBlocks"]
        nblock = info.attrs["MeshBlockSize"]
    cellsPerBlock = nblock[0] * nblock[1] * nblock[2]
    allCellsPerBlock1 = (nblock[0] + 2) * (nblock[1] + 2) * (nblock[2] + 2)
    allCellsPerBlock2 = (nblock[0] + 4) * (nblock[1] + 4) * (nblock[2] + 4)

    numcell = nb * cellsPerBlock
    nclone1 = nb * allCellsPerBlock1
    nclone2 = nb * allCellsPerBlock2

    tdict = {
        "cycle": cycle,
        "walltime": walltime,
        "numcell": numcell,
        "nclone1": nclone1,
        "nclone2": nclone2,
    }
    return time, tdict


if __name__ == "__main__":
    import sys

    timeinfo = {}
    tdict = {}
    for f in sys.argv[1:]:
        try:
            t, tdict = getRiotInfo(f)
            timeinfo[t] = tdict
        except:
            print(f"ERROR reading file {f}")

    print("")
    print(f"t,{','.join(tdict.keys())}")
    for t in sorted(timeinfo.keys()):
        tdict = timeinfo[t]
        print(t, end="")
        for x in tdict:
            print(",", tdict[x], end="")
        print("")
