#include <cstdio>
#include <map>
#include "pio.hpp"

static std::map<std::string,double> pmap;

double writeIndex(FILE *fp, char *name, double nameLength, double index, double length, double position)  {
    // write index
  char s[static_cast<size_t>(nameLength)];
  for (int i=0; i<nameLength; i++) s[i] = ' ';
  for (int i=0; i<strlen(name); i++) s[i] = name[i];
  std::string s2 = std::string(name) + "_" + std::to_string(index);
  std::cout << "Name=" << s << "=" << position << ":" << pmap[s2] << std::endl;
  fwrite(s, nameLength, 1, fp);
  
  double ahdr[4];
  ahdr[0] = index;
  ahdr[1] = length;
  ahdr[2] = position;
  ahdr[3] = 0;
  fwrite(ahdr, sizeof(double), 4, fp);
 
  return length;
}


int main(int argc, char **argv) {
  char *name;
  name = argv[1];
  
  auto p = PIO(name);

  double numcell = p.numcell();

  for  (auto &h : p.arrayDims) {
    auto name = h.first;
    auto dims = h.second;
    printf("%40s: %d x %d\n", name.c_str(),  dims.w, dims.l);
  }

  std::cout << "dimensionality: "<< p.ndim() << std::endl;
  std::cout << "numcell: "<< numcell << std::endl;

  FILE *fp = fopen("tmp-dmp000000","wb");
  PIOHeader hdr = p.header();
  hdr.nArrays = p.ndim()*(3) // cell_center + cell_index
                + 3 + 1;     // cell_modhter/daughter/level + proc_id
  hdr.position = hdr.lengthHeader +  hdr.nArrays * numcell;

  fwrite(&hdr, sizeof(hdr), 1, fp);
  std::cout << "ncell=" << numcell << std::endl;
  std::cout << "ndim=" << p.ndim() << std::endl;
  std::cout << "diff=" << 8*hdr.lengthHeader-sizeof(hdr)/sizeof(double)<<std::endl;
  size_t xlen = hdr.lengthHeader - 11;
  double d[xlen] = {0};
  fwrite(d, sizeof(double), xlen, fp);

  std::cout << "xlen=" << xlen << ":" << hdr.lengthHeader << std::endl;
  double position = hdr.lengthHeader;
  
  std::cout << "hdr done" << std::endl;
  
  for ( int i = 1; i<= p.ndim(); i++) {
    std::vector<double> v = p.variable("cell_center", i);
    fwrite(v.data(), sizeof(double), numcell, fp);
    pmap[std::string("cell_center_") + std::to_string(i)] = position;
    position += numcell;    
    std::cout << "cell_center " << i << " done " << std::endl;
  }
  {
    pmap[std::string("cell_mother_0")] = position;
    std::vector<double> v = p.variable("cell_mother");
    fwrite(v.data(), sizeof(double), numcell, fp);
    position += numcell;
  }
  std::cout << "cell_mother done" << std::endl;
  {
    pmap[std::string("cell_daughter_0")] = position;
    std::vector<double> v = p.variable("cell_daughter");
    fwrite(v.data(), sizeof(double), numcell, fp);
    position += numcell;
  }
  std::cout << "cell_daughter done" << std::endl;
  {
    pmap[std::string("cell_level_0")] = position;
    std::vector<double> v = p.variable("cell_level");
    fwrite(v.data(), sizeof(double), numcell, fp);
    position += numcell;
  }
  std::cout << "cell_level done" << std::endl;
  for ( int i = 1; i<= 2 * p.ndim(); i++) {
    pmap[std::string("cell_index_")] = position;
    std::vector<double> v = p.variable("cell_index", i);
    fwrite(v.data(), sizeof(double), numcell, fp);
    std::cout << "cell_index " << i << " done" << std::endl;
    position += numcell;
  }
  {
    std::vector<double> v = p.variable("cell_index", 1);
    pmap[std::string("cell_index_0")] = position;
    for (int i=0; i<numcell; i++) {
      v[i] = i;
    }
    fwrite(v.data(), sizeof(double), numcell, fp);
    position += numcell;
  }
  double nameLength = hdr.lengthName;
  position = hdr.lengthHeader;
  for( int i=1; i<=p.ndim(); i++) {
    position += writeIndex(fp, "cell_center", nameLength, i, numcell, position);
  }
  position += writeIndex(fp, "cell_mother", nameLength, 0, numcell, position);
  position += writeIndex(fp, "cell_daughter", nameLength, 0, numcell, position);
  position += writeIndex(fp, "cell_level", nameLength, 0, numcell, position);
  for( int i=1; i<= 2*p.ndim(); i++) {
    position += writeIndex(fp, "cell_index", nameLength, i, numcell, position);
  }
  position += writeIndex(fp, "processor_id", nameLength, 0, numcell, position);

  fclose(fp);

  return 0;
}
