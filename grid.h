#define _USE_MATH_DEFINES
#include <Eigen/SparseCore>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include "utils.h"
using namespace std;

#pragma once

namespace grid {
  /*CSV File Writing Functions*/
  extern const string fluid_header;
  void write_csv_header(string csv_name, string header);
  void write_fluid_vars(string csv_name,
                        tuple<Eigen::VectorXd, Eigen::VectorXd,
                              Eigen::VectorXd, Eigen::VectorXd> data,
                        bool overwrite);
  void write_grid(string csv_name,
                  tuple<Eigen::VectorXd, Eigen::VectorXd> data,
                  bool overwrite);

  /*Point Generators*/
  std::pair<std::vector<double>, std::vector<double>>
  polar_points(std::vector<double> r_locs, std::vector<double> t_locs);
  tuple<Eigen::VectorXd, Eigen::VectorXd> all_points(double R_min, double R_max,
                                                     int N);

  /*Grid Metrics*/
  double calc_dr(std::vector<double> r_pts);
  double calc_dr(Eigen::VectorXd r_pts);
  }
