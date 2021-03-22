#include <cstring>
#include <iostream>
#include <vector>

#include "pio.hpp"

int main(int argc, char **argv) {
  PIO p(argv[1]);
  for ( auto& x : p.arrayOrder) {
    std::cout << x << std::endl;
  }
  std::vector<int64_t> level=p.variable<int64_t>("cell_level");
  int id=0;
  for ( int i=0; i<10; i++) {
    for (int j=0; j<10; j++,id++) {
      std::cout << level[id] << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
