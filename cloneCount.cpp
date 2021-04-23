//========================================================================================
// (C) (or copyright) 2020. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S. Department of Energy/National Nuclear
// Security Administration. All rights in the program are reserved by Triad
// National Security, LLC, and the U.S. Department of Energy/National Nuclear
// Security Administration. The Government is granted for itself and others
// acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
// in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do
// so.
//========================================================================================

// TO compile:
//   g++ cloneCount.cpp -I. -o cloneCount -O3
//
// Usage: ./cloneCount <filename> [nProcs]
//
// If nProcs is omitted, it will use the partition in the dump file.
// If nProcs is present, it will generate an estimate of the partition
// for that many processors
//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "pio.hpp"

int main(int argc, char **argv) {

  // Check number of arguments
  if (argc < 2 || argc > 3) {
    std::cerr << std::endl
              << "  ERROR: Usage: " << argv[0] << " <filename>  [nProcs]"
              << std::endl
              << std::endl;
    return (1);
  }

  PIO p(argv[1]);
  int ndim = p.ndim();
  int64_t numcell = p.numcell();

  // Generate / read in partitioning
  int64_t nProcs;
  std::vector<int64_t> gns;
  if (argc > 2) {
    // Generate partitioning
    nProcs = (int64_t)std::stoi(argv[2]);
    gns.resize(nProcs);
    int64_t block = (1 << ndim);
    int64_t quantum = round((double)(numcell) / (double)(nProcs));
    int64_t remainder = quantum % block;
    quantum = quantum - remainder;
    int64_t nEnd = 0;
    int64_t nStart = 0;
    int64_t r = 0;
    for (int id = 0; id < nProcs; id++) {
      int64_t offset = 0;
      nStart = nEnd;
      if (id == nProcs - 1) {
        nEnd = numcell;
      } else {
        r += remainder;
      }
      if (r >= block) {
        offset = block;
        r = r - block;
      } else {
        offset = 0;
      }
      nEnd = nStart + quantum + offset;
      gns[id] = nEnd - nStart;
    }
  } else {
    // Read partitioning from file
    nProcs = p.arraySize("global_numcell").l;
    gns = p.variable<int64_t>("global_numcell");
  }
  // Read in daughter array
  std::vector<int64_t> dtr = p.variable<int64_t>("cell_daughter");

  // Read in neighbors, note this *may* create copies
  std::vector<int64_t> nbrs[2][ndim];
  for (int i = 0; i < ndim; i++) {
    nbrs[0][i] = p.variable<int64_t>("cell_index", 2 * i + 1);
    nbrs[1][i] = p.variable<int64_t>("cell_index", 2 * i + 2);
  }

  // Convert cell count to  procLimits array
  int64_t y = 0;
  std::vector<int> procID;
  procID.resize(numcell);
  int64_t j = 0;
  for (int id = 0; id < nProcs; id++) {
    for (int k = 0; k < gns[id]; k++, j++) {
      procID[j] = id;
    }
  }

  // Check neighbor and count clones
  int64_t nClones = 0;
  int64_t nMother = 0;
  int64_t nTop = 0;
  for (int64_t i = 0; i < numcell; i++) {
    for (int idir = 0; idir < ndim; idir++) {
      if (dtr[i] > 0) {
        if (idir == 0) {
          nMother++;
        }
        continue;
      }
      nTop++;
      int64_t xL = nbrs[0][idir][i] - 1;
      int64_t xH = nbrs[1][idir][i] - 1;
      if (xL == i || xL < 0 || procID[xL] != procID[i]) {
        nClones++;
      }
      if (xH == i || xH < 0 || procID[xH] != procID[i]) {
        nClones++;
      }
    }
  }

  std::cout << "numcell = " << numcell << std::endl;
  std::cout << "nClones = " << nClones << std::endl;
  std::cout << "nMother = " << nMother << std::endl;
  std::cout << "nTop = " << nTop << std::endl;
  std::cout << "Ratio = " << (double)(nClones + nMother) / (double)nTop
            << std::endl;
  return 0;
}
