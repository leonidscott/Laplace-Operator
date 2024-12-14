#include <iostream>
#include <vector>
#pragma once

namespace utils {
  void print_arr(std::vector<double>);
  void print_i_arr(std::vector<int>);
  void prn_c_arr(double* arr, int length);

  std::vector<double> range(double start, double end, int N, bool st_ex=false);
  std::vector<int> irange(int start, int end, int step=1, bool st_ex=false);
}
