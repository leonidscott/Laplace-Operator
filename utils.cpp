#include <iostream>
#include <vector>
#include "utils.h"

namespace utils {
  void print_arr(std::vector<double> arr){
    auto prn = [](double el) { std::cout << el <<", " ; } ;
    std::cout << "[";
    std::for_each(arr.begin(), arr.end()-1, prn);
    std::cout << arr[arr.size() -1] << "]\n";
  }

  std::vector<double> range(double start, double end, int N, bool st_ex/*=false*/){
    double const step = (end-start)/(st_ex ? N+1 : N);
    double const first_step = (st_ex  ? start+step : start);
    std::vector<double> out(N);
    for (int i = 0; i < out.size(); i++) out[i] = (first_step + i*step);
    return out;
  }
}
