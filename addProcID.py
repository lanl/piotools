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

 Script: addProcID.py
 Author: Sriram Swaminarayan (sriram@lanl.gov)
   Date: March 16, 2021
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
import os
import sys
import numpy as np

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
class pio:
    def __init__(self, theFile, outbase='tmp'):
        self.offset = 0
        outfile = outbase+'-dmp000000'
        if os.path.exists(outfile):
            print(f"\n\n  *** File {outfile} exists, will not proceed. ***\n\n")
            return None
        self.fp = open(theFile, mode='rb')
        
        # read signature
        s = self.str(8)
        print(f'  File type is: {s}')
        if not s.lower() == b'pio_file':
            raise ValueError('invalid file type')

        # Ensure second value is two
        d = self.double()
        if d != 2.0:
            raise ValueError('second value is not 2')
        
        # read version
        self.version = self.int64()
        print(f'version={self.version}')

        self.lName = self.int64()
        self.lHeader = self.int64()
        self.lIndex = self.int64()


        self.date = self.str(16)
        print(self.date)

        self.n = self.int64()
        self.position = self.int64()
        self.signature = self.int64()
        
        print(f'numvars={self.n}, offsetindex={8*self.position}')
        
        # seek to end of header

        # read the name
        self.names={}
        self.xnames = []
        offset = int(8*self.position)
        self.seek(offset)
        for i in range(int(self.n)):
            hdf = self.arrayHeader()
            idx = hdf['name'] + b'_%d'%hdf['index']
            self.names[idx] = hdf
            self.xnames.append(hdf)
            
        cch = self.names[b'cell_center_1']
        self.numcell = int(cch['length'])
        
        print(f'__________NUMCELL={self.numcell}')

        # Write the new dump file
        self.writeNewArray(outfile)

        # Write the Paraview file
        with open(outbase+'.pio','w') as ofp:
            ofp.writelines(
                ["DUMP_DIRECTORY .\n",
                 "DUMP_BASE_NAME tmp\n"]
                )
        
    def writeNewArray(self, outName):
        # now add in a new array
        self.oldPosition = self.position
        self.addCellArray('processor_id')
        newData = np.zeros(self.numcell,dtype='double')
        gnh = self.names[b'global_numcell_0']
        data = self.readArray(b'global_numcell_0')
        lastIdx = 0
        for i in range(len(data)):
            nextIdx = lastIdx + int(data[i])
            newData[lastIdx:nextIdx] = i
            lastIdx = nextIdx
            
        with open(outName,'wb') as ofp:
            self.ofp = ofp
            # write the header
            
            ofp.write(b'pio_file')
            np.array([2.0,
                      self.version,
                      self.lName,
                      self.lHeader,
                      self.lIndex],
                     dtype='double').tofile(ofp)
            ofp.write(self.date)
            np.array([self.n, self.position],dtype='double').tofile(ofp)
            
            # copy rest of file till new data position
            self.copyToOldIndex(ofp,10)
            
            # write new data to file
            newData.tofile(ofp)
            
            # write trailer
            self.writeTrailer(ofp)

            ofp.close()

    def writeTrailer(self, outfp):
        for i in range(self.n):
            outfp.write(self.xnames[i]['bytes'])
            

    def doubleBytes(self, value):
        return bytes(np.array(value,dtype='double'))
        
    def copyToOldIndex(self,outfp, offset):
        self.seek(8*offset)
        sz = self.oldPosition - offset
        bufsize = self.numcell
        read = 0
        left = sz
        while left:
            if left < bufsize:
                bufsize = left
            buf = self.double(bufsize)
            buf.tofile(outfp)
            read += bufsize
            left -= bufsize
            print(f'{0.95*(100.*read)/sz:.2f}%  written ')

    def addCellArray(self, name):
        cch = self.copyHeader(self.names[b'pres_0'])
        fmt = '%%-%ds'%self.lName
        cch['name'] = bytes(name,'utf8')
        longName = bytes(fmt%name,'utf8')
        cch['offset'] = self.lIndex
        b = cch['bytes']
        o = self.lName
        cch['bytes'] = longName + self.doubleBytes([0, self.numcell, self.position, 0])
        print(cch['bytes'], len(cch['bytes']))
        self.position += self.numcell
        self.n += 1
        self.names[name+'_0'] = cch
        self.xnames.append(cch)
        return

    def readArray(self, name):
        if name not in self.names:
            return None
        hdr = self.names[name]
        self.seek(int(hdr['offset']))
        data = self.double(hdr['length'],force=True)
        print('return:',data)
        return data
    
    def arrayHeader(self, offset=None):
        """
        reads array at current header
        """
        start = self.offset
        name = self.str(self.lName).strip()
        index = self.int64()
        length = self.int64()
        offset = 8*self.int64()
        self.seek(start)
        data = self.str(8*self.lIndex)
        return {'name':name, 'index':index, 'length':length, 'offset':offset, 'bytes':data}
        # return {'name':name, 'index':index, 'length':length, 'offset':offset}
        
    def copyHeader(self,src):
        ret = {}
        for x in src:
            ret[x] = src[x]
        return ret

    def seek(self, offset):
        offset = int(offset)
        self.fp.seek(offset)
        self.offset = offset

    def str(self, count=1, offset=None):
        count = int(count)
        if offset is not None:
            self.seek(offset)
        s = self.fp.read(count)
        self.offset += count
        return s

    def double(self, count=1, offset=None, force=False):
        count = int(count)
        if offset is not None:
            self.seek(offset)
        value = np.fromfile(self.fp, dtype='double', count=count)
        if count == 1 and (not force):
            return value[0]
        return value
            
    def int64(self, count=1, offset=None, force=False):
        count = int(count)
        if offset is not None:
            self.seek(offset)
        raw = np.fromfile(self.fp, dtype='double', count=count)
        if count == 1 and (not force):
            return int(raw[0])
        
        value = [int(x) for x in raw]
        return value
            

if __name__ == "__main__":
    if (len(sys.argv) > 3) or (len(sys.argv) <2):
          usage()
    else:
        filename = sys.argv[1]
        try:
            outbase = sys.argv[2]
        except:
            outbase = 'tmp'
        p = pio(filename, outbase)
