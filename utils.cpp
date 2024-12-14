#include <iostream>
#include <vector>
#include "utils.h"

namespace utils {
  //lol I should figure out generics in c++ at some point
  void print_arr(std::vector<double> arr){
    auto prn = [](double el) { std::cout << el <<", " ; } ;
    std::cout << "[";
    std::for_each(arr.begin(), arr.end()-1, prn);
    std::cout << arr[arr.size() -1] << "]\n";
  }

  void print_i_arr(std::vector<int> arr){
    auto prn = [](int el) { std::cout << el <<", " ; } ;
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

  std::vector<int> irange(int start, int end, int step/*=1*/, bool st_ex/*=false*/) {
    int const out_len = (int)((end - (st_ex ? start+1 : start))/step);
    int const first_step = (st_ex  ? start+step : start);
    std::vector<int> out(out_len);
    for (int i = 0; i < out.size(); i++) out[i] = (first_step + i*step);
    return out;
  }

  void prn_c_arr(double* arr, int length) {
    std::cout << "[";
    for(int i = 0 ; i < length; i++) {
      std::cout << *(arr+i);
      if(i != length-1) std::cout << ",";
    }
    std::cout << "]\n";
  }
}
