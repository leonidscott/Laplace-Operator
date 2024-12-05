#include <iostream>
#include <vector>
#include "utils.h"

int main() {
  //* range(0,1,N) tests *//
  // Start Inclusive - End Exclusive
  std::vector<double> t0 = {0, 0.25, 0.5, 0.75};
  if(utils::range(0, 1, 4) != t0) std::cout<< "t0 failure\n";

  // Start Exclusive - End Exclusive
  std::vector<double> t1 = {0.2, 0.4, 0.6, 0.8};
  if(utils::range(0,1,4,/*st_ex=*/true) != t1) {
    std::cout << "t1 failure\n";
    std::cout << "expected: ";
    utils::print_arr(t1);
    std::cout << "got:      ";
    utils::print_arr(utils::range(0,1,4,/*st_ex=*/true));
    std::cout << std::endl;
      }

  //* range(1,5,N) tests  *//
  // Start Inclusive - End Exclusive
  std::vector<double> t2 = { 1, 7.0/3.0, 11.0/3.0};
  if(utils::range(1,5,3,/*st_ex=*/false) != t2) {
    std::cout << "t2 failure\n";
    std::cout << "expected: ";
    utils::print_arr(t2);
    std::cout << "got:      ";
    utils::print_arr(utils::range(1,5,3,/*st_ex=*/false));
    std::cout << std::endl;
  }

  // Start Exclusive - End Exclusive
  std::vector<double> t3 = { 2, 3, 4 };
  if(utils::range(1,5,3,/*st_ex=*/true) != t3) std::cout<< "t3 failure\n";
}
