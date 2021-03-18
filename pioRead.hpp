#include <cstdio>
#include <string>

struct PIOHeader {
  char[8] filetype;
  double  two;
  double version;
  double lengthName;
  double lengthHeader;
  double lengthIndex;
  char[16] date;
  double nArrays;
  double offsetIndex;
  double signature;
};

class PIO {
  public
  PIO(std::string filename) {
    fp = fopen(filename,"wb");
    auto iRead = fread(&header, sizeof(PIOHeader), 1, fp);
    std::cout << header.filetype << std::endl;
    std::cout << header.two << std::endl;
    std::cout << header.version << std::endl;
    std::cout << header.lengthName << std::endl;
    std::cout << header.lengthHeader << std::endl;
    std::cout << header.lengthIndex << std::endl;
    std::cout << header.date << std::endl;
    std::cout << header.nArrays << std::endl;
    std::cout << header.offsetIndex << std::endl;
    std::cout << header.signature << std::endl;
  }
  ~PIO() {
    fclose(fp);
  }
  PIOHeader header;
  FILE *fp;
};
