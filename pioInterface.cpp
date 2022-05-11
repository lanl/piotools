//========================================================================================
// (C) (or copyright) 2022. Triad National Security, LLC. All rights reserved.
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
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <stdlib.h>
#include <string.h>

#include "pioInterface.h"
#include "pioInterface.hpp"

void PioInterface::listFields(FILE *fp) { //< lists fields in the file
  std::vector<std::string> names = pd->arrayOrder;
  for (auto &n : names) {
    fprintf(fp, "%s\n", n.c_str());
  }
  return;
}

void PioInterface::updateDXyz() {
  const int inbr[3] = {1, 2,
                       4}; /* the already known level 1 neighbor information */

  dXyz_ = new double *[nLevel_ + 1];

  for (int l = 0; l <= nLevel_; l++) {
    dXyz_[l] = new double[3];
    for (int d = 0; d < 3; d++) {
      dXyz_[l][d] = 1.0;
    }
  }

  /* compute level 1 from known neighbors */
  for (int d = 0; d < nDim_; d++)
    dXyz_[1][d] = 0.5 * (center_[d][inbr[d]] - center_[d][0]);

  /* compute the rest of the levels from level 1 */
  for (int l = 2; l <= nLevel_; l++) {
    for (int d = 0; d < nDim_; d++) {
      dXyz_[l][d] = 0.5 * dXyz_[l - 1][d];
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
    counts[level_[i]]++;
  }
  /* allocate space for the map */
  for (int l = 1; l <= nLevel_; l++) {
    iMap[l] = (int64_t *)calloc(counts[l], sizeof(int64_t));
    counts[l] = 0; /* we will use this as index! */
  }

  /* fill in the maps */
  for (int64_t i = 0; i < nCell_; i++) {
    if (daughter_[i] > 0)
      continue;
    int lvl = level_[i];
    iMap[lvl][counts[lvl]] = i;
    counts[lvl]++;
  }

  return;
}

void PioInterface::releaseMapByLevel() {
  for (int i = 0; i <= nLevel_; i++) {
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
    nLevel_ = (nLevel_ < level_[i] ? level_[i] : nLevel_);
  }
}
PioInterface::~PioInterface() { //< our destructor!
  // release all data
  if (pd)
    delete pd;
  if (uniqMap_)
    delete[] uniqMap_;
  for (int i = 0; i <= nLevel_; i++) {
    delete[] dXyz_[i];
  }
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
  const double *dxyz = dXyz_[nLevel_];
  const double scale = (double)(((int64_t)1) << nLevel_);
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
  //  for(int d=0; d<nDim_; d++) delete[] idxyz[d];
  delete[] idxyz;

  return;
}

std::map<int, std::vector<double>>
PioInterface::getMaterialVariable(const char *field, int64_t iStart,
                                  int64_t nCount) {
  std::map<int, std::vector<double>> rMap;
  if (nCount < 0) {
    nCount = nCell_;
  }

  // Compute offset and count and get the raw data
  int64_t matOffset = matStartIndex_[iStart];
  int64_t matCount = matStartIndex_[iStart + nCount] - matOffset;
  auto data = getField<double>(field, 0, matOffset, matCount);
  if (data.size() == 0) {
    // Adjust field name in case of error
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
      return getMaterialVariable(newField, iStart, nCount);
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
    rMap[i] = std::vector<double>(nCount, 0.0);
  }

  if (nMat() == 1) {
    memcpy(rMap[1].data(), data.data(), nCount * sizeof(double));
    return rMap;
  }

  // cycle through cells and fill in fields
  for (int64_t icell = iStart; icell < iStart + nCount; icell++) {
    for (int64_t indexMat = matStartIndex_[icell];
         indexMat < matStartIndex_[icell + 1]; indexMat++) {
      int idMat = matIds_[indexMat];
      rMap[idMat][icell - iStart] = data[indexMat - matOffset];
    }
  }
  if (verbose_) {
    std::cout << "Done generating map" << std::endl;
  }
  freeField(field);
  return rMap;
}

// initializer takes dump file name and request for unique ids
PioInterface::PioInterface(const char *name, const int uniq, const int verbose)
    : uniqMap_(nullptr), dXyz_(nullptr), iMap(nullptr), verbose_(verbose) {
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
    level_ = getField<int>("cell_level");
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
    matIds_ = getField<int>("chunk_mat");
    matStartIndex_ = getField<int64_t>("chunk_nummat");
    matStartIndex_.resize(nCell_ + 1);
    {
      // Shift the start indices
      int64_t startIdx;
      startIdx = 0;
      for (int64_t i = 0; i < nCell_; i++) {
        int64_t n = matStartIndex_[i];
        matStartIndex_[i] = startIdx;
        startIdx += n;
      }
      matStartIndex_[nCell_] = startIdx;
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

template <typename T> static T *deepCopy(std::vector<T> &src) {
  T *ret;
  ret = (T *)calloc(src.size(), sizeof(T));
  memcpy(ret, src.data(), src.size() * sizeof(T));
  return ret;
}

/*******************************************
 * C interface to enable fortran integration
 *******************************************/
static std::map<const int, std::shared_ptr<PioInterface>> myFiles;
extern "C" {
void *pio_ftn_index_2d(void *in, int index) {
  return (void *)((void **)in)[index];
}

void pio_release(const int ID) {
  auto it = myFiles.find(ID);
  if (it != myFiles.end()) {
    myFiles.erase(it);
  }
}

void pio_init(const int ID, const char *fname, const int verbose) {
  pio_release(ID);
  myFiles.insert(std::pair<const int, std::shared_ptr<PioInterface>>(
      ID, std::make_shared<PioInterface>(fname, 0, verbose)));
  return;
}

int pio_nCell(const int ID) { return myFiles[ID]->nCell(); }

int pio_nDim(const int ID) { return myFiles[ID]->nDim(); }

int pio_nMat(const int ID) { return myFiles[ID]->nMat(); }

int64_t pio_length(const int ID, const char *field) {
  return myFiles[ID]->getFieldLength(field);
}

int64_t pio_width(const int ID, const char *field) {
  return myFiles[ID]->getFieldWidth(field);
}

int64_t *pio_daughter(const int ID) { return myFiles[ID]->daughter().data(); }

double *pio_center(const int ID, int index) {
  return myFiles[ID]->center()[index].data();
}

double *pio_get_range_d(int ID, const char *var_name, int index, int64_t iStart,
                        int64_t nCount) {
  auto &pd = myFiles[ID];
  auto data = pd->getField<double>(var_name, index, iStart, nCount);
  return deepCopy<double>(data);
}

double *pio_get_d(int ID, const char *var_name, int index) {
  return pio_get_range_d(ID, var_name, index, 0, -1);
}

int64_t *pio_get_range_i64(int ID, const char *var_name, int index,
                           int64_t iStart, int64_t nCount) {
  auto &pd = myFiles[ID];
  auto data = pd->getField<int64_t>(var_name, index, iStart, nCount);
  return deepCopy<int64_t>(data);
}

// PioInterface::getMaterialVariable(const char *field, int64_t iStart, int64_t
// nCount) {
double **pio_get_matvar_range_d(int ID, const char *field, int64_t iStart,
                                int64_t nCount) {
  double **ret;
  auto &pd = myFiles[ID];
  int nMat = pd->nMat();
  auto myVar = pd->getMaterialVariable(field, iStart, nCount);
  ret = (double **)calloc(nMat + 1, myVar[0].size() * sizeof(double *));
  ret[0] = nullptr;
  for (int i = 1; i <= nMat; i++) {
    ret[i] = deepCopy<double>(myVar[i]);
  }
  return ret;
}

double **pio_get_matvar_d(int ID, const char *field) {
  return pio_get_matvar_range_d(ID, field, 0, -1);
}

int64_t *pio_get_i64(int ID, const char *var_name, int index) {
  return pio_get_range_i64(ID, var_name, index, 0, -1);
}

void pio_release_i64(int64_t *ptr) { free(ptr); }

void pio_release_d(double *ptr) { free(ptr); }

void pio_release_matvar_d(double **ptr, int nMat) {
  for (int i = 1; i <= nMat; i++) {
    free(ptr[i]);
  }
  free(ptr);
}
}

#ifdef DOPIOMAIN
/**********************
 * A simple inline test
 **********************/
int main(int argc, const char **argv) {
  {
    std::cout << "FROM C++" << std::endl;
    PioInterface a(argv[1], 0, 0);
    int ndim = a.nDim();
    a.listFields(stdout);

    auto center = a.center();
    std::cout << std::endl
              << "____________CENTERS________________" << std::endl;
    for (int i = 0; i < 10; i++) {
      std::cout << i << ":";
      for (int j = 0; j < ndim; j++) {
        std::cout << " " << center[j][i];
      }
      std::cout << std::endl;
    }
    std::cout << "C++: ndim=" << ndim << std::endl;
    std::cout << std::endl
              << "____________C++ VFRAC________________" << a.nMat()
              << std::endl;
    int64_t iStart = 2927;
    int64_t nCount = 10;
    auto myVar = a.getMaterialVariable("chunk_vol", iStart, nCount);
    for (int i = 0; i < 10; i++) {
      std::cout << i + iStart << ":";
      for (int j = 1; j <= a.nMat(); j++) {
        std::cout << " " << myVar[j][i];
      }
      std::cout << " centers=";
      for (int j = 0; j < a.nDim(); j++) {
        std::cout << center[j][iStart + i] << " ";
      }
      std::cout << std::endl;
    }
  }
  std::cout << std::endl << "______________________" << std::endl;
  { /* C test */
    int id;
    int ndim = -1;
    id = 22;
    std::cout << "FROM C" << std::endl;
    pio_init(id, argv[1], 0);
    auto center = myFiles[id]->center();
    ndim = pio_nDim(id);
    std::cout << "  C: id=" << id << "; ndim=" << ndim << std::endl;
    double *xc1 = pio_get_range_d(id, "cell_center", 1, 4, 6);
    double *xc2 = pio_get_range_d(id, "cell_center", 2, 4, 6);
    for (int i = 0; i < 6; i++) {
      std::cout << i + 4 << ":" << xc1[i] << " " << xc2[i] << std::endl;
    }

    std::cout << std::endl << "_____________Volume chunks 1-10:" << std::endl;
    int64_t iStart = 2927;
    int64_t nCount = 10;
    double **mcv = pio_get_matvar_range_d(id, "chunk_vol", iStart, nCount);
    for (int i = 0; i < 10; i++) {
      std::cout << i + iStart << ": ";
      for (int j = 1; j <= pio_nMat(id); j++) {
        std::cout << mcv[j][i] << " ";
      }
      std::cout << " centers=";
      for (int j = 0; j < myFiles[id]->nDim(); j++) {
        std::cout << center[j][iStart + i] << " ";
      }
      std::cout << std::endl;
    }

    pio_release_matvar_d(mcv, pio_nMat(id));
    pio_release_d(xc1);
    pio_release_d(xc2);
    pio_release(id);

    if (argc > 2) {
      /* C sedond file test */
      id = 22;
      pio_init(id, argv[2], 0);
      ndim = pio_nDim(id);
      std::cout << "  C: id=" << id << "; ndim=" << ndim << std::endl;
      pio_release(id);
    }
  }
  return 0;
}
#endif
