#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <functional>
#include "utils.h"
#include "grid.h"

using namespace std;

namespace numerical_methods {

  // Public Matrix Functions
  const auto zero = [](int row){return 0.0;};
  const auto one  = [](int row){return 1.0;};
  const auto n_constructor = [](double n){ return [n](int row){ return n; }; };

  // Vector Factories
  Eigen::VectorXd def_vector(int dim, std::function<double(int)> fn);
  vector<double> def_std_vector(int dim, std::function<int(double)> fn);

  // Fluid Equations (Arguably Should go in a fluids.cpp)
  double exact_potential(double r, double t);

  // Derivative Stencils
  Eigen::SparseMatrix<double> D1(int const dim);
  Eigen::SparseMatrix<double> D1_periodic(int const dim);
  Eigen::SparseMatrix<double> L2_stencil(Eigen::VectorXd r_ups,
                                         Eigen::VectorXd t_ups);
  Eigen::VectorXd L2_polar(Eigen::VectorXd r_ups, Eigen::VectorXd t_ups,
                           Eigen::VectorXd f,     Eigen::VectorXd bcs);
}
