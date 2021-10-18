#!/usr/projects/env python3
"""
Generate data for number of cells Vs time
"""

import sys, os
from pio import *


def clone(theCycle, theBase):
    from cloneInfo import CloneInfo

    fname = f"{theBase}-dmp{theCycle:>06}"
    try:
        x = CloneInfo(fname)
        return x["nClones"]
    except:
        return None


outputFile = sys.argv[1]
if not outputFile.endswith("-output"):
    print(
        f"""
    
    Usage: {sys.argv[0]} problem-output

    """
    )
    exit(1)
basename = outputFile[:-7]
timeinfo = {}
mpc = None
with open(sys.argv[1], "r") as fp:
    lines = fp.readlines()
    for idx, l in enumerate(lines):
        if l.startswith(
            "    cycle    t             dtnext    timestep  cstb  tpct  epct  ritr  hitr   sumritr    wallhr sumwallhr    sumcpuhr      date:time"
        ):
            data = lines[idx + 1].strip().split()
            c = int(data[0])
            t = float(data[1])
            dt = float(data[2])
            wall = float(data[10])
            wallSum = float(data[11])
            cpuSum = float(data[12])
            tstamp = data[13]
            clones = clone(c, basename)
            timeinfo[t] = {
                "cycle": c,
                "dt": dt,
                "wall": wall,
                "wallSum": wallSum,
                "cpuSum": cpuSum,
                "tstamp": tstamp,
            }
            tdict = timeinfo[t]
            data = lines[idx + 5].strip().split()
            tdict["mxlvl"] = int(float(data[0]))
            tdict["numpe"] = int(float(data[1]))
            tdict["nummat"] = int(float(data[2]))
            tdict["sum_cell"] = int(float(data[3]))
            tdict["sum_top"] = int(float(data[4]))
            if clones is not None:
                tdict["numclone"] = tdict["sum_cell"] + clones
            else:
                tdict["numclone"] = ""
            tdict["avg_mat_per_cell"] = mpc
            mpc = None
        elif l.strip().startswith("materials/cell"):
            mpc = float(l.split()[-1])

print(f"t,{','.join(tdict.keys())}")


for t in sorted(timeinfo.keys()):
    tdict = timeinfo[t]
    print(t, end="")
    for x in tdict:
        print(",", tdict[x], end="")
    print("")

# files = sys.argv[2:]
# for f in files:
#     p = pio(f,verbose=1)
#     #for x in p.names:
#     #    print(x)

#     t = p.readArray(b"timertimes_0")
#     t = p.readArray(b"hist_dt_0")
#     print(t)


print(basename)
