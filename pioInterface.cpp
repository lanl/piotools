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

#include <iostream>
#include <string>
#include <vector>

#include "pioInterface.hpp"

void PioInterface::listFields(FILE *fp) { //< lists fields in the file
  std::vector<std::string> names = pd->arrayOrder;
  for ( auto& n : names) {
    fprintf(fp, "%s\n", n.c_str());
  }
  return;
}

void PioInterface::updateDXyz() {
  /* we are going to cheat and use knowledge of the xrage structure to fill
     in the levels.  Basically we know the neighbbors of cell #1 */
  const int inbr[3] = {1, 2,
                       4}; /* the already known level 1 neighbor information */

  dXyz_ = new xyz_t[nLevel_ + 1];
  for (int l = 0; l <= nLevel_; l++) {
    for (int d = 0; d < 3; d++)
      dXyz_[l].xyz[d] = 1.0;
  }

  /* compute level 1 from known neighbors */
  for (int d = 0; d < nDim_; d++)
    dXyz_[1].xyz[d] = 0.5 * (center_[d][inbr[d]] - center_[d][0]);

  /* compute the rest of the levels from level 1 */
  for (int l = 2; l <= nLevel_; l++) {
    for (int d = 0; d < nDim_; d++) {
      dXyz_[l].xyz[d] = 0.5 * dXyz_[l - 1].xyz[d];
    }
  }
  return;
}

void PioInterface::freeField(const char *name) {
  // do nothing
  return;
}

int64_t PioInterface::getFieldLength(const char *field) {
  return pd->arrayDims[field].l;
}

int64_t PioInterface::getFieldWidth(const char *field) {
  return pd->arrayDims[field].w;
}

template <class T> const T *PioInterface::getUniqMap(const T *field) {
  T *fRet = new T[nCell_];

  for (int64_t i = 0; i < nCell_; i++)
    fRet[i] = field[uniqMap_[i]];
  return (const T *)fRet;
}

template <class T>
const T **PioInterface::getUniqMap(const T **field, const int n) {
  T **fRet = new T *[n];
  for (int d = 0; d < n; d++) {
    fRet[d] = new T[nCell_];
    for (int64_t i = 0; i < nCell_; i++) {
      fRet[d][i] = field[d][uniqMap_[i]];
    }
  }
  return (const T **)fRet;
}

template <class T>
void PioInterface::deleteArray(const T **field, const int n) {
  for (int d = 0; d < n; d++)
    delete[] field[d];
  delete[] field;
}


void PioInterface::updateIMap() {
  // returns a two dimensional array with a map of which cells are
  // leaf cells at any given level.  The level with index 0 contains
  // the number of leaf cells at each of the levels.

  int64_t *counts = NULL;
  iMap = (int64_t **)calloc(nLevel_ + 1, sizeof(int64_t *));

  /* compute how many cells at each level */
  counts = iMap[0] = (int64_t *)calloc(nLevel_ + 1, sizeof(int64_t));
  memset(iMap[0], 0, (nLevel_ + 1) * sizeof(int64_t));
  for (int64_t i = 0; i < nCell_; i++) {
    if (daughter_[i] > 0)
      continue;
    counts[(int64_t)level_[i]]++;
  }
  /* allocate space for the map */
  for (int64_t l = 1; l <= nLevel_; l++) {
    iMap[l] = (int64_t *)calloc(counts[l], sizeof(int64_t));
    counts[l] = 0; /* we will use this as index! */
  }

  /* fill in the maps */
  for (int64_t i = 0; i < nCell_; i++) {
    if (daughter_[i] > 0)
      continue;
    int64_t lvl = (int64_t)level_[i];
    iMap[lvl][counts[lvl]] = i;
    counts[lvl]++;
  }

  return;
}

void PioInterface::releaseMapByLevel() {
  for (int64_t i = 0; i <= nLevel_; i++) {
    free(iMap[i]);
  }
  memset(iMap, 0xDEADBEEF, sizeof(void *));
  free(iMap);
  return;
}

void PioInterface::updateNLevel() { // Finds out how many levels we have by
                                     // inspecting
                                     // cell_level variable
  nLevel_ = 0;
  for (int64_t i = 0; i < nCell_; i++) {
    nLevel_ = (nLevel_ < (int64_t)level_[i] ? level_[i] : nLevel_);
  }
}
PioInterface::~PioInterface() { //< our destructor!
  // release all data
  if (pd)
    delete pd;
  if (uniqMap_)
    delete[] uniqMap_;
  delete[] dXyz_;
  if (iMap)
    releaseMapByLevel();
}

static int i2Compare(const void *va, const void *vb) {
  const i2_t *a = (const i2_t *)va;
  const i2_t *b = (const i2_t *)vb;
  if (a->id < b->id)
    return -1;
  else if (a->id == b->id)
    return 0;
  return 1;
}
#define mymax(a, b) ((a) > (b) ? (a) : (b))
void PioInterface::updateUniqMap() {
  /**<
   * generates a unique id for each cell so that we can directly
   * compare two dump files that were generated from runs on different
   * numbers of processors.
   *
   * Here is how it is done:
   * 1: divide each coordinate by the dXyz of the finest level
   * 2: find Nx,Ny,and Nz as the highest index in x,y,and z directions
   * 3: id = (x/dx) + Nx*(y/dy) + Nx*Ny*(z/dz)
   * 4: sort the IDs to get the map so that uniq[i] = id of cell that
   *    is the ith cell in this mapping
   *
   * I'm sure there's a better way to do this, but this will have to
   * do for now
   **/

  /** allocate temporary space for our IDs **/
  i2_t *i2 = new i2_t[nCell_];
  int64_t nMax[nDim_];
  const double *dxyz = dXyz_[nLevel_].xyz;
  const double scale = (double)(1 << nLevel_);
  int64_t(*idxyz)[3] = new int64_t[nCell_][3];
  //  for(int d=0; d<nDim_; d++) idxyz[d] = new int64_t [nCell_];

  /* zero out the Nxyz */
  memset(nMax, 0, sizeof(nMax));

  /* compute coordinates and max of dimensions */
  for (int d = 0; d < nDim_; d++) {
    for (int64_t i = 0; i < nCell_; i++) {
      idxyz[i][d] = (int64_t)floor((scale * center_[d][i] / dxyz[d] + 0.5));
      nMax[d] = mymax(nMax[d], idxyz[i][d]);
    }
  }

  /* compute the temporary ID and fill in index */
  for (int64_t i = 0; i < nCell_; i++) {
    i2[i].index = i;
    i2[i].id = idxyz[i][0];
    if (nDim_ > 1)
      i2[i].id += nMax[0] * idxyz[i][1]; // 2d?
    if (nDim_ > 1)
      i2[i].id += nMax[0] * nMax[1] * idxyz[i][2]; // 3d?
  }

  qsort((void *)i2, nCell_, sizeof(i2_t), i2Compare);

  /** allocate a unique map **/
  uniqMap_ = new int64_t[nCell_];
  for (int64_t i = 0; i < nCell_; i++)
    uniqMap_[i] = i2[i].index;

  /** free data **/
  delete[] i2;
  //  for(int d=0; d<nDim_; d++) delete[] idxyz.xyz[d];
  delete[] idxyz;

  return;
}

std::map<int, std::vector<double>>
PioInterface::getMaterialVariable(const char *field) {
  std::map<int, std::vector<double>> rMap;
  auto data = getField<double>(field);
  if (data.size() == 0) {
    if (verbose_) {
      std::cout << "Unable to find Field: " << field << std::endl;
    }
    if (!strncmp("chunk_", field, 6)) {
      char newField[strlen(field) + 1];
      sprintf(newField, "frac_%s_0", field + 6);
      if (verbose_) {
        std::cout << "Unable to find Field: " << field << " trying " << newField
                  << std::endl;
      }
      return getMaterialVariable(newField);
    } else {
      if (verbose_) {
        std::cout << "Unable to find Field: " << field << std::endl;
      }
      return rMap;
    }
  }
  if (verbose_) {
    std::cout << "Found Field: " << field << std::endl;
  }

  // allocate space for variables
  for (int i = 1; i <= nMat(); i++) {
    rMap[i] = std::vector<double>(nCell_, 0.0);
  }

  if (nMat() == 1 ) {
    memcpy(rMap[1].data(),data.data(), nCell_*sizeof(double));
    return rMap;
  }
		      
  // cycle through cells and fill in fields
  int64_t idx = 0;
  for (int64_t icell = 0; icell < nCell_; icell++) {
    int64_t n = nMatPerCell_[icell];
    for (int iMat = 0; iMat < n; iMat++, idx++) {
      int idMat = static_cast<int>(idMatPerCell_[idx]);
      rMap[idMat][icell] = data[idx];
    }
  }
  if (verbose_) {
    std::cout << "Done generating map" << std::endl;
  }
  freeField(field);
  return rMap;
}

// initializer takes dump file name and request for unique ids
PioInterface::PioInterface(const char *name, const int uniq,
                             const int verbose)
    : uniqMap_(nullptr), dXyz_(nullptr), iMap(nullptr),
      verbose_(verbose) {
  // initializes a class from file name and request for unique map
  try {
    uniq_ = uniq;
    if (verbose) {
      std::cout << "getting piodata\n";
    }
    pd = new PIO(name);

    nDim_ = pd->ndim();
    nCell_ = pd->numcell();

    if (verbose) {
      std::cout << "getting levels\n";
    }
    level_ = getField<int64_t>("cell_level");
    updateNLevel(); // Note this requires the level array to be gathered first

    if (verbose) {
      std::cout << "getting centers\n";
    }
    center_ = getField2D<double>("cell_center");

    if (verbose) {
      std::cout << "getting daughters\n";
    }
    daughter_ = getField<int64_t>("cell_daughter");

    // Get material variable information
    if (verbose) {
      std::cout << "updating material information\n";
    }
    nMat_ = getFieldWidth("matdef");
    nMatPerCell_ = getField<int64_t>("chunk_nummat");
    idMatPerCell_ = getField<int64_t>("chunk_mat");
    if (verbose > 10) {
      std::cout << "number of materials =" << nMat_ << std::endl;
      const int64_t ilen = getFieldLength("chunk_nummat");
      int64_t idx = 0;
      for (int64_t i = 0; i < ilen; i++) {
        int n = nMatPerCell_[i];
        std::cout << "cell=" << i << ":" << n << "::";
        for (int64_t j = 0; j < n; j++, idx++) {
          std::cout << ":" << idMatPerCell_[idx];
        }
        std::cout << std::endl;
      }
    }

    if (verbose) {
      std::cout << "updating Dxyz\n";
    }
    updateDXyz();
    if (uniq)
      updateUniqMap();
    if (verbose) {
      std::cout << "done\n";
    }

  } catch (...) {
    std::cerr << "\n"
              << "Following Error ocured while trying to read " << name << "\n"
              << "\n";
    exit(2);
  }
}

std::vector<std::string> PioInterface::getFieldNames() {
  return pd->arrayOrder;
}

#ifdef DOPIOMAIN
int main(int argc, const char **argv) {
  PioInterface a(argv[1], 0, 0);
  int ndim = a.nDim();
  a.listFields(stdout);
  auto center = a.center();
  std::cout << std::endl << "____________CENTERS________________" << std::endl;
  for ( int i=0; i<10; i++){
    std::cout << i << ":" ;
    for(int j=0; j<ndim; j++) {
      std::cout << " " << center[j][i];
    }
    std::cout << std::endl;
  }
  return 0;
}
#endif
