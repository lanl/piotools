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
"""
from __future__ import print_function
import numpy as np
import struct

class pio:
    """
    Author: Sriram Swaminarayan (sriram@lanl.gov)
    Date: March 17, 2021
    Version: 1.0

    This class will read PIO files and allow minimal manipulation.
    It is expected that other scripts will build on the PIO class
    to do heavy-duty lifting.
    """

    def __init__(self, theFile):
        """
        Initializes a PIO class and returns a PIO object.
        Argument is the PIO file to be read
        """

        self.fp = open(theFile, mode="rb")

        self.offset = 0

        # read signature
        s = self.str(8)
        print(f"  File type is: {s}")
        if not s.lower() == b"pio_file":
            raise ValueError("Invalid file type")

        # Ensure second value is two
        d = self.doubles()
        if d != 2.0:
            raise ValueError("Second value is not 2")

        # read version
        self.version = self.ints()
        print(f"version={self.version}")

        # read element lengths
        self.lName = self.ints()
        self.lHeader = self.ints()
        self.lIndex = self.ints()

        # read data/time
        self.date = self.str(16)
        print(self.date)

        # read number of variables and index location
        self.n = self.ints()
        self.position = self.ints()

        # read file signature
        self.signature = self.ints()

         # read the variable index
        self.names = {}
        self.xnames = []
        print('position=',self.position)
        self.seek(self.position)
        for i in range(int(self.n)):
            hdf = self.readArrayHeader()
            idx = hdf["name"].strip() + b"_%d" % hdf["index"]
            self.names[idx] = hdf
            self.xnames.append(hdf)

        cch = self.names[b"cell_center_1"]
        self.numcell = int(cch["length"])
        self.outOffset = -1

    def writeHeader(self, fp):
        """
        Writes PIO header to the file pointer sent.
        """
        self.outOffset = 0
        fp.write(b"pio_file")
        np.array(
            [2.0, self.version, self.lName, self.lHeader, self.lIndex],
            dtype="double",
        ).tofile(fp)
        fp.write(self.date)
        np.array([self.n, self.position, self.signature], dtype="double").tofile(fp)
        np.zeros((self.lHeader-11), dtype='double').tofile(fp)
        self.outOffset = self.lHeader

    def writeWithNewCellArray(self, outName, newName, newData):
        """
        Adds the new cell array with name `newName` and data
        `newData` to the PIO file and writes it to file named
        `outname`.
        """
        
        # now add in a new array
        self.oldPosition = self.position
        
        self.addCellArray(newName)
        
        with open(outName, "wb") as ofp:
            # write the header
            self.writeHeader(ofp)

            # copy rest of file till old Index offset
            self.copyToOffset(ofp, self.lHeader, self.oldPosition)

            # write new data to file
            newData.tofile(ofp)
            self.outOffset += len(newData)

            # write Index
            self.writeIndex(ofp)

            ofp.close()
        self.outOffset = -1
        

    def writeIndex(self, outfp):
        for x in self.xnames:
            outfp.write(x["bytes"])
            self.outOffset += self.lIndex
            
    def copyToOffset(self, outfp, offsetStart, offsetEnd):
        self.seek(offsetStart)
        sz = offsetEnd - offsetStart
        if self.numcell < 1048576:
            bufsize = 1048576
        else:
            bufsize = self.numcell
        written = 0
        left = sz
        while left:
            if left < bufsize:
                bufsize = left
            buf = self.doubles(bufsize)
            buf.tofile(outfp)
            written += bufsize
            left -= bufsize
            print(f"{0.95*(100.*written)/sz:.2f}%  written ")
        self.outOffset += sz

    def addCellArray(self, name):
        cch = self.copyArrayHeader(self.names[b"pres_0"])
        fmt = "%%-%ds" % self.lName
        cch["name"] = bytes(name, "utf8")
        longName = bytes(fmt % name, "utf8")
        cch["offset"] = self.lIndex
        b = cch["bytes"]
        o = self.lName
        cch["bytes"] = longName + bytes(np.array([0, self.numcell, self.position, 0], dtype='double'))
        print(cch["bytes"], len(cch["bytes"]))
        self.position += self.numcell
        self.n += 1
        self.names[name + "_0"] = cch
        self.xnames.append(cch)
        return

    def readArray(self, name):
        if name not in self.names:
            return None
        hdr = self.names[name]
        self.seek(hdr["offset"])
        data = self.doubles(hdr["length"], force=True)
        print("return:", data)
        return data

    def readArrayHeader(self):
        """
        reads array at current header
        """
        x = self.lName
        start = self.offset
        data = self.str(8 * self.lIndex)
        name = data[:x]
        index = int(struct.unpack('d',data[x:x+8])[0])
        length = int(struct.unpack('d',data[x+8:x+16])[0])
        offset = int(struct.unpack('d',data[x+16:x+24])[0])
        return {
            "name": name,
            "index": index,
            "length": length,
            "offset": offset,
            "bytes": data,
        }
    
    def copyArrayHeader(self, src):
        """ Returns a copy of the array header """
        ret = {}
        for x in src:
            ret[x] = src[x]
        return ret

    def seek(self, offset):
        """ Seeks to offset in the input file """
        self.offset = offset
        offset = int(8 * offset)
        self.fp.seek(offset)

    def str(self, count=1, offset=None):
        """ Reads count characters from the input file """
        count = int(count)
        if offset is not None:
            self.seek(offset)
        s = self.fp.read(count)
        self.offset += int(count/8)
        return s

    def doubles(self, count=1, offset=None, force=False):
        """
        Reads count doubles from the input file.
        If count == 1 and force is False, then it will
        return scalars.  
        """
        count = int(count)
        if offset is not None:
            self.seek(offset)
        value = np.fromfile(self.fp, dtype="double", count=count)
        self.offset += count
        if count == 1 and (not force):
            return value[0]
        return value

    def ints(self, count=1, offset=None, force=False):
        """
        Reads count doubles from the input file and returns as ints.
        If count == 1 and force is False, then it will
        return scalars.  
        """
        count = int(count)
        if offset is not None:
            self.seek(offset)
        raw = np.fromfile(self.fp, dtype="double", count=count)
        self.offset += count
        if count == 1 and (not force):
            return int(raw[0])
        value = [int(x) for x in raw]
        return value


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print(f"""
        {sys.argv[0]} takes at least one argument.
        """)
    else:
        filename = sys.argv[1]
        p = pio(filename)
        print('number of cells = ',p.numcell)


    
