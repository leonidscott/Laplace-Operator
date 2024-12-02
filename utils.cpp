#include <iostream>
#include <vector>
#include "utils.h"

namespace utils {
  void print_arr(std::vector<double> arr){
    auto prn = [](double el) { std::cout << el <<"," ; } ;
    std::cout << "[";
    std::for_each(arr.begin(), arr.end()-1, prn);
    std::cout << arr[arr.size() -1] << "]\n";
  }

  std::vector<double> range(double start, double end, double step) {
    int out_length =  (end-start)/step ;
    std::vector<double> out(out_length);
    for (int i = 0; i < out.size(); i++) out[i] = (start + i*step);
    return out;
  }
}
