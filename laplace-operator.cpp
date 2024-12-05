#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <cmath>
#include <functional>
#include "utils.h"

// Public Matrix Functions
auto zero = [](int row){return 0;};
auto one  = [](int row){return 1;};

Eigen::MatrixXd tridiagonal(int dim,
                            std::function<double(int)> l_left,
                            std::function<double(int)> l_center,
                            std::function<double(int)> l_right) {
  Eigen::MatrixXd tdm(dim, dim); //tdm = Tri-Diagonal Matrix

  for (int i = 0; i < dim; i++) {
    //Compute Center
    tdm(i,i) = l_center(i);
    //Compute Left
    if(i > 0) tdm(i,i-1) = l_left(i);
    //Compute Right
    if(i < dim-1) tdm(i,i+1) = l_right(i);
  }

  return tdm;
}

Eigen::VectorXd vector(int dim, std::function<double(int)> fn) {
  Eigen::VectorXd v(dim);
  for (int i = 0; i< dim; i++){
    v(i) = fn(i);
  }
  return v;
}

Eigen::MatrixXd stencil(int nip, std::vector<double> r_pts, std::vector<double> t_pts) {
  double dr = r_pts[1] - r_pts[0];
  double dt = t_pts[1] - t_pts[0];
  auto center = [r_pts, t_pts, dr, dt](int row) {
    // -2 * (1/dr^2 + 1/(r^2 dt^2))
    return -2 * (1/pow(dr, 2) + 1/(pow(r_pts[row%t_pts.size()],2) * pow(dt,2)));
  };
  auto off_center = [r_pts, t_pts, dt](int row) {
    // 1 / (r^2 dt^2)
    int r = r_pts[row%t_pts.size()];
    return 1 / (pow(r, 2) * pow(dt, 2));
  };
  Eigen::MatrixXd T = tridiagonal(nip, off_center, center, off_center);
  return T;
}

Eigen::MatrixXd L2_polar(std::vector<double> r_points, std::vector<double> t_points) { //TODO: Fix return type
  /*
   * L2 = (I(ss) ⊗ T) + stuff
   */

  //Assumes r_points.size() == t_points.size()
  //Assumes r_points.size() >= 2
  //Assumes t_points.size() >= 2
  //Could derive system_size from r_points, and made sure that r_points
  int nip = r_points.size(); //Number of internal points

  //I(ss)
  Eigen::MatrixXd I_ss = Eigen::MatrixXd::Identity(nip, nip);

  //T
  //Need to create logic to filter out bc from stencil?
  Eigen::MatrixXd Stencil = stencil(nip , r_points, t_points);

  //Off Diagonal Terms
  // --> Lower Diagonal Term
  int ntp = t_points.size(); // Number of theta points
  double dr = r_points[1] - r_points[0];
  auto lower_off = [r_points, ntp, dr](int row) {
    // (1/dr^2 - 1/(r 2 dr))
    int r = r_points[row%ntp];
    return ( 1/pow(dr,2) - (1/(r * 2 * dr)) );
  };
  //Eigen::MatrixXd lower_diag = tridiagonal(pow(nip,2), lower_off, zero, zero);
  Eigen::MatrixXd lower_diag = Eigen::kroneckerProduct
    (tridiagonal(nip, one, zero, zero),
     tridiagonal(nip, zero, lower_off, zero)
     );

  // --> Upper Diagonal Term
  auto upper_off = [r_points, ntp, dr](int row) {
    // (1/dr^2 - 1/(r 2 dr))
    int r = r_points[row%ntp];
    return ( 1/pow(dr,2) + (1/(r * 2 * dr)) );
  };
  //Eigen::MatrixXd upper_diag = tridiagonal(pow(nip,2), zero, zero, upper_off);
  Eigen::MatrixXd upper_diag = Eigen::kroneckerProduct
    (tridiagonal(nip, zero , zero, one),
    tridiagonal(nip, zero, upper_off, zero)
     );

  //(I(ss) ⊗ T)
  return Eigen::kroneckerProduct(I_ss, Stencil) + lower_diag + upper_diag;
}

int main() {
  std::cout << "Hello World" << std::endl;

  /*
  int const N = 10;
  double const step = (double)1/(double)N;

  std::vector<double> x = utils::range(0, 1+step, step);
  std::cout << "x: " ;
  utils::print_arr(x);
  std::vector<double> y = utils::range(0, 1+step, step);
  std::cout << "y: " ;
  utils::print_arr(y);
  */

  double R_max = 5;
  double R_min = 1;
  double const N = 3;
  double const r_step = (R_max - R_min)/(N+1);
  double const t_step = (2*M_PI- 0)/N;
  //Excluding Boundary Conditions for now
  std::vector<double> r = utils::range(1+r_step, R_max, r_step);
  utils::print_arr(r);
  std::cout << "r.size() " << r.size() << std::endl;
  std::vector<double> t = utils::range(0, 2*M_PI, t_step);
  utils::print_arr(t);
  std::cout << "t.size() " << t.size() << std::endl;

  // Construct L2
  Eigen::MatrixXd L2 = L2_polar(r, t); //TODO: Correct for boundary conditions

  return 0;
}
