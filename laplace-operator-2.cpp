#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <cmath>
#include <functional>
#include "utils.h"

using namespace std;

void prn_c_arr(double* arr, int length) {
  std::cout << "[";
  for(int i = 0 ; i < length; i++) {
    std::cout << *(arr+i);
    if(i != length-1) std::cout << ",";
  }
  std::cout << "]\n";
}


// Public Matrix Functions
auto zero = [](int row){return 0;};
auto one  = [](int row){return 1;};

Eigen::VectorXd def_vector(int dim, std::function<double(int)> fn) {
  Eigen::VectorXd v(dim);
  for (int i = 0; i< dim; i++){
    v(i) = fn(i);
  }
  return v;
}

double calc_dr(vector<double> r_pts){
  double r_init = r_pts[0];
  for(int i = 0; i < r_pts.size() ; i++) {
    if (r_pts[i] != r_init) return (r_pts[i] - r_init);
  }
  return 0.0;
}

/**
   (For now) takes a list unknown points specified in radius and theta values
**/
Eigen::MatrixXd L2_matrix(std::vector<double> r_ups, std::vector<double> t_ups) {
  //TODO: Add a check that r_ups and t_ups are the same length
  int nup = r_ups.size(); //Number of unknown points
  int dim = std::sqrt(nup);
  Eigen::MatrixXd L2(nup, nup);

  //Start Lambda Definitions
  double dr = calc_dr(r_ups);
  double dt = t_ups[1] - t_ups[0];

  auto center = [r_ups, t_ups, dr, dt](int row) {
    // -2 * (1/dr^2 + 1/(r^2 dt^2))
    return -2 * (1/pow(dr, 2) + 1/(pow(r_ups[row%t_ups.size()],2) * pow(dt,2)));
  };

  auto l_and_r= [r_ups, dt](int row) {
    // 1 / (r^2 dt^2)
    int r = r_ups[row];
    return 1 / (pow(r, 2) * pow(dt, 2));
  };

  auto upper= [r_ups, dr](int row) {
    // (1/dr^2 - 1/(r 2 dr))
    int r = r_ups[row];
    return ( 1/pow(dr,2) + (1/(r * 2 * dr)) );
  };

  auto lower= [r_ups, dr](int row) {
    // (1/dr^2 - 1/(r 2 dr))
    int r = r_ups[row];
    return ( 1/pow(dr,2) - (1/(r * 2 * dr)) );
  };
  //End Lambda Definitions

  for (int i = 0; i < nup; i++){
    //Center Function
    L2(i,i) = center(i);
    //Left & Right
    switch((i+1) % dim) {
    case 1: //First row of T Matrix
      L2(i, i+1) = l_and_r(i);
      L2(i, i + (dim-1)) = l_and_r(i);
      break;
    case 0: //Last row of T matrix
      L2(i, i-1) = l_and_r(i);
      L2(i, i - (dim-1)) = l_and_r(i);
      break;
    default: //All ofther rows of T matrix
      L2(i, i-1) = l_and_r(i);
      L2(i, i+1) = l_and_r(i);
    }
    //Upper Function
    if(i - dim >= 0)  L2(i, i - dim) = upper(i);
    //Lower Function
    if(i + dim < nup) L2(i, i + dim) = lower(i);
  }

  return L2;
}

Eigen::VectorXd L2_polar(std::vector<double> r_ups, std::vector<double> t_ups, Eigen::VectorXd f, Eigen::VectorXd bcs) {
  Eigen::MatrixXd L2_mat = L2_matrix(r_ups, t_ups);
  Eigen::VectorXd RHS = f - bcs;

  //Solve L2_mat u = RHS
  Eigen::VectorXd unknowns = L2_mat.lu().solve(RHS);
  return unknowns;
}

/** Radius Major Points given radius and angle locations**/
pair<vector<double>, vector<double>> polar_points(vector<double> r_locs, vector<double> t_locs) {
  vector<double> r_pts(r_locs.size() * t_locs.size());
  vector<double> t_pts(r_locs.size() * t_locs.size());

  for (int r = 0; r < r_locs.size() ; r++) {
    for (int t = 0 ; t < t_locs.size() ; t++) {
      r_pts[r*t_locs.size() + t] = r_locs[r];
      t_pts[r*t_locs.size() + t] = t_locs[t];
    }
  }

  utils::print_arr(r_pts);
  utils::print_arr(t_pts);

  pair<vector<double>, vector<double>> polar_pts = { r_pts, t_pts };
  return polar_pts;
}

extern "C" {
  /*double**/ void polar_laplace(int R_min, int R_max, int N, double* out){
    // Create Internal Points
    double const r_step = (R_max - R_min)/(N+1);
    double const t_step = (2*M_PI- 0)/N;

    std::vector<double> rlocs = utils::range(1+r_step, R_max, r_step);
    std::vector<double> tlocs = utils::range(0, 2*M_PI, t_step);

    vector<double> r;
    vector<double> t;
    std::tie(r,t) = polar_points(rlocs, tlocs);

    // F Function
    Eigen::VectorXd f = def_vector(rlocs.size() * tlocs.size(), zero);

    // Boundary Conditions
    // 1. R = R_min
    auto rmin_bc = [rlocs, R_min](int i) {
      return (int(i/rlocs.size()) == 0 ? 0.0 : 0.0);
    };
    // 2. R = R_max
    auto rmax_bc = [rlocs, tlocs, t, R_max](int i) {
      int dim = rlocs.size() * tlocs.size();
      int t_val = t[i];
      return (int((dim-1 -i)/rlocs.size()) == 0 ? R_max*sin(t_val) : 0.0);
    };
    Eigen::VectorXd bcs = \
      def_vector(rlocs.size() * tlocs.size(), rmin_bc) +
      def_vector(rlocs.size() * tlocs.size(), rmax_bc);

    Eigen::VectorXd L2 = L2_polar(r, t, f, bcs);

    double* result;
    Eigen::Map<Eigen::VectorXd>(result, L2.size())= L2;
    std::cout << "hi\n";
    prn_c_arr(result, 9);
    out = result;
    //return out;
  }
}

int main() {
  std::cout << "Hello world\n";
  double const R_max = 5;
  double const R_min = 1;
  double const N = 3;
  double const r_step = (R_max - R_min)/(N+1);
  double const t_step = (2*M_PI- 0)/N;

  //Excluding Boundary Conditions for now
  std::vector<double> rlocs = utils::range(1+r_step, R_max, r_step);
  utils::print_arr(rlocs);
  std::cout << "rlocs.size() " << rlocs.size() << std::endl;
  std::vector<double> tlocs = utils::range(0, 2*M_PI, t_step);
  utils::print_arr(tlocs);
  std::cout << "tlocs.size() " << tlocs.size() << std::endl;

  vector<double> r;
  vector<double> t;
  std::tie(r,t) = polar_points(rlocs, tlocs);

  // F function
  Eigen::VectorXd f = def_vector(rlocs.size() * tlocs.size(), zero);

  // Boundary Conditions
  // 1. R = R_min
  auto rmin_bc = [rlocs, R_min](int i) {
    return (int(i/rlocs.size()) == 0 ? 0.0 : 0.0);
  };
  // 2. R = R_max
  auto rmax_bc = [rlocs, tlocs, t, R_max](int i) {
    int dim = rlocs.size() * tlocs.size();
    int t_val = t[i];
    return (int((dim-1 -i)/rlocs.size()) == 0 ? R_max*sin(t_val) : 0.0);
  };
  Eigen::VectorXd bcs = \
    def_vector(rlocs.size() * tlocs.size(), rmin_bc) +
    def_vector(rlocs.size() * tlocs.size(), rmax_bc);

  Eigen::MatrixXd L2 = L2_polar(r, t, f, bcs);
  std::cout << L2 << std::endl;

  double* out;
  polar_laplace(1, 5, 3, out);
  prn_c_arr(out, 9);
}
