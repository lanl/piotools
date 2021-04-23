#include <iostream>
#include <string>
#include <vector>

#include "pio.hpp"

int main(int argc, char **argv) {
  PIO p(argv[1]);
  int ndim = p.ndim();
  int64_t numcell = p.numcell();

  // Read in partitioning
  int64_t nProcs = p.arraySize("global_numcell").l;
  std::vector<int64_t> gns = p.variable<int64_t>("global_numcell");

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
        nMother++;
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
