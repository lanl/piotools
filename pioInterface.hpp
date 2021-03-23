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

#ifndef EXAMPLE_AMHC_PIOINTERFACE_HPP_
#define EXAMPLE_AMHC_PIOINTERFACE_HPP_

#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include <math.h>
#include <unistd.h>

#include "pio.hpp"

typedef struct xyz_t {
  double xyz[3];
} xyz_t;

typedef struct i2_t {
  int64_t id;
  int64_t index;
} i2_t;

class PioInterface {
private:
  int uniq_;   //< if set to 1 will provide unique ids across multiple processor
               // runs
  int nDim_;   //< number of dimensions
  int nLevel_; //< number of levels
  int64_t nCell_;    //< number of cells
  int64_t *uniqMap_; //< if uniqId is set, provides mapping from uniq ids to ids
                     // in dump
                     // file [nCell]
  std::vector<int64_t> level_;    //< level of each cell [nCell]
  std::vector<int64_t> daughter_; //< daughter of each cell [nCell]
  std::map<int, std::vector<double>>
      center_;                        //< center of each cell [nDim][nCell]
  xyz_t *dXyz_;                       //< cell size at each level [nLevel]
  int64_t **iMap;                     //< mapping of cells by level
  PIO *pd;                            //< PIO data struct for dmp file
  int64_t nMat_;                      //< Number of materials
  std::vector<int64_t> nMatPerCell_;  //< Number of materials per cell
  std::vector<int64_t> idMatPerCell_; //< Array of material IDs
  int verbose_;                       //< print verbose information

  // private functions
  void updateIMap();
  void updateDXyz();
  void updateNCell();
  void updateNLevel();
  void updateUniqMap();

  void releaseMapByLevel();
  void freeField(const char *name);

public:
  void listFields(FILE *fp); //< prints fields in the dmp file to fp

  /** member access functions **/
  int uniq() { return uniq_; }
  int nDim() { return nDim_; }
  int nMat() { return nMat_; }
  int64_t nCell() { return nCell_; }
  int nLevel() { return nLevel_; }
  const xyz_t *dXyz() { return dXyz_; }
  const int64_t *uniqMap() { return (const int64_t *)uniqMap_; }

  std::map<int, std::vector<double>> &center() { return center_; }
  std::vector<int64_t> &level() { return level_; }
  std::vector<int64_t> &daughter() { return daughter_; }

  std::vector<std::string> getFieldNames();

  int64_t
  getFieldWidth(const char *field); //< Width / Number of instances of a field
  int64_t getFieldLength(const char *field); //< Length of a field

  std::vector<const char *> getVCField(const char *field, int index = 0) {
    /** for any type other than double we need to get double data and then
     * translate **/
    /* given a pio_data field and a field name, returns the data associated with
     * the field
     */
    std::vector<const char *> cVec;
    std::vector<double> origData = getField<double>(field, index);
    int l = getFieldLength(field);
    const char *data =
        strndup((const char *)(origData.data()), l * sizeof(double));
    cVec.push_back(data);
    return cVec;
  }

  std::string getStringField(const char *field, int index = 0) {
    std::vector<double> origData = getField<double>(field, index);
    int l = getFieldLength(field);
    const char *data =
        strndup((const char *)(origData.data()), l * sizeof(double));
    std::string s;
    s = std::string(data);
    free((void *)data);
    return s;
  }

  template <typename T>
  std::vector<T> getVariable(const char *field, int index = 0) {
    // field name *must* exactly match a field in PIO file
    return pd->variable<T>(field, index);
  }

  template <typename T>
  std::vector<T> getField(const char *field, int index = 0) {
    return pd->variable<T>(field, index);
  }

  template <typename T>
  std::map<int, std::vector<T>> getField2D(const char *field) {
    std::map<int, std::vector<T>> data;
    int w = getFieldWidth(field);
    for (int i = 1; i <= w; i++) {
      data[i - 1] = pd->variable<T>(field, i);
    }
    return data;
  }

  std::map<int, std::vector<double>>
  getMaterialVariable(const char *field); //< gets a map of a material variable

  template <class T> const T *getUniqMap(const T *field);
  template <class T> const T **getUniqMap(const T **field, const int n);
  template <class T> void deleteArray(const T **field, const int n);

  std::vector<std::shared_ptr<double>> getDChunkField(const char *field);

  // initializer takes dump file name and request for unique ids
  PioInterface(const char *name, const int uniq = 0, const int verbose = 0);
  ~PioInterface();
};

#endif
