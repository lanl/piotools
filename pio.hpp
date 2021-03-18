//========================================================================================
// (C) (or copyright) 2021. Triad National Security, LLC. All rights reserved.
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

#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#pragma pack(push, 1)
struct PIOHeader {
  char filetype[8];
  double two;
  double version;
  double lengthName;
  double lengthHeader;
  double lengthIndex;
  char date[16];
  double nArrays;
  double position;
  double signature;
};

struct PIOArrayHeader {
  double index;
  double length;
  double position;
};
#pragma pack(pop)

struct arrayDimensions {
  int64_t l;
  int64_t w;
  arrayDimensions() : l(0), w(0) {}
};

class PIO {
public:
  std::vector<std::string> arrayOrder;
  std::map<std::string, PIOArrayHeader> arrays;
  std::map<std::string, arrayDimensions> arrayDims;
  PIO(std::string filename, bool verbose = false) {
    ndim_ = 0;
    numcell_ = 0;
    fp = fopen(filename.c_str(), "rb");
    auto iRead = fread(&header_, sizeof(PIOHeader), 1, fp);
    if (verbose) {
      printf("%8s\n", header_.filetype);
    }
    if (strncmp(header_.filetype, "pio_file", 8)) {
      std::cout << "Unable to open file" << std::endl;
      memset(&header_, -1, sizeof(PIOHeader));
      return;
    }
    if (header_.two != 2.0) {
      std::cout << "File has right signature but wrong format" << std::endl;
      memset(&header_, -1, sizeof(PIOHeader));
      return;
    }
    if (verbose) {
      std::cout << "           two: " << header_.two << std::endl;
      std::cout << "       version: " << header_.version << std::endl;
      std::cout << "   Name Length: " << header_.lengthName << std::endl;
      std::cout << " Header Length: " << header_.lengthHeader << std::endl;
      std::cout << "  Index Length: " << header_.lengthIndex << std::endl;
      std::cout << "          Date: " << header_.date << std::endl;
      std::cout << "      N Arrays: " << header_.nArrays << std::endl;
      std::cout << "  Index Offset: " << header_.position << std::endl;
      std::cout << "File Signature: " << header_.signature << std::endl;
    }
    // Load array headers
    loadArrayHeaders();
  }
  ~PIO() { fclose(fp); }
  PIOHeader header() { return header_; }
  int ndim() { return ndim_; }
  int numcell() { return numcell_; }

  std::vector<double> variableRead(std::string name, int index = 0) {
    return readArray(name + "_" + std::to_string(index));
  }

  std::vector<double> variable(std::string name, int index = 0) {
    return variableRead(name, index);
  }

  template <typename T>
  std::vector<T> variable(std::string name, int index = 0) {
    auto base = variableRead(name, index);
    return std::vector<T>(base.begin(), base.end());
  }

  std::vector<double> readArray(std::string name) {
    std::vector<double> v;
    if (arrays.count(name) > 0) {
      auto h = arrays[name];
      seek(h.position);
      v.resize(h.length);
      auto iRead = fread(v.data(), h.length, sizeof(double), fp);
    }
    return v;
  }
  void seekRaw(size_t offset) {
    auto myOffset = static_cast<size_t>(offset);
    fseek(fp, static_cast<long int>(offset), SEEK_SET);
  }
  void seek(double offset) { seekRaw(8.0 * offset); }

private:
  FILE *fp;
  int ndim_;
  int numcell_;
  PIOHeader header_;

  std::string toString(char *inS, int l) {
    int idx = l;
    while (idx) {
      idx--;
      if (inS[idx] == ' ') {
        inS[idx] = '\0';
      } else {
        break;
      }
    }
    char *s = strndup(inS, idx + 1);
    std::string ret(s);
    free(s);
    return ret;
  }
  void loadArrayHeaders() {
    double pos = header_.position;
    for (int i = 0; i < header_.nArrays; i++) {
      seek(pos);
      PIOArrayHeader a;
      char name[static_cast<int>(header_.lengthName)];
      size_t iRead;
      iRead = fread(name, sizeof(char), header_.lengthName, fp);
      iRead += fread(&a, sizeof(PIOArrayHeader), 1, fp);
      auto baseName = toString(name, header_.lengthName);
      auto nameStr = baseName + "_" + std::to_string(int(a.index));
      arrayOrder.push_back(nameStr);
      arrays[nameStr] = a;
      if (!strncmp("cell_center", baseName.c_str(), 11)) {
        numcell_ = a.length;
        ndim_ = (ndim_ > a.index ? ndim_ : a.index);
      }

      // Set dimensions of array
      // default constructor gives zeros for first time
      auto &x = arrayDims[baseName];
      x.l = a.length;
      x.w++;

      // position set to next header
      pos += header_.lengthIndex;
    }
  }
};

// END
