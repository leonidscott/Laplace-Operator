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
auto zero = [](int row){return 0.0;};
auto one  = [](int row){return 1;};

Eigen::VectorXd def_vector(int dim, std::function<double(int)> fn) {
  Eigen::VectorXd v(dim);
  v.setZero();
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
Eigen::MatrixXd L2_stencil(std::vector<double> r_ups, std::vector<double> t_ups) {
  //TODO: Add a check that r_ups and t_ups are the same length
  int nup = r_ups.size(); //Number of unknown points
  int dim = std::sqrt(nup);
  Eigen::MatrixXd L2(nup, nup);
  L2.setZero();

  //Start Lambda Definitions
  double dr = calc_dr(r_ups);
  double dt = t_ups[1] - t_ups[0];

  auto center = [r_ups, dr, dt](int row) {
    // -2 * (1/dr^2 + 1/(r^2 dt^2))
    if(std::isnan(2 * (1/pow(dr, 2) + 1/(pow(r_ups[row],2) * pow(dt,2))))) std::cout << "center NAN\n";
    return -2.0 * (1/pow(dr, 2) + 1/(pow(r_ups[row],2) * pow(dt,2)));
  };

  auto l_and_r= [r_ups, dt](int row) {
    // 1 / (r^2 dt^2)
    double r = r_ups[row];
    if (std::isnan(1 / (pow(r, 2) * pow(dt, 2)))) std::cout << "l_and_r NAN\n";
    return 1 / (pow(r, 2) * pow(dt, 2));
  };

  auto upper= [r_ups, dr](int row) {
    // (1/dr^2 + 1/(r 2 dr))
    double r = r_ups[row];
    if (std::isnan(( 1/pow(dr,2) + (1/(r * 2 * dr)) ))) std::cout << "upper NAN\n";
    return ( 1/pow(dr,2) + (1/(2 * r * dr)) );
  };

  auto lower= [r_ups, dr](int row) {
    // (1/dr^2 - 1/(r 2 dr))
    double r = r_ups[row];
    if(std::isnan( 1/pow(dr,2) - (1/(r * 2 * dr)) )) std::cout << "lower NAN\n";
    return ( 1/pow(dr,2) - (1/(2 * r * dr)) );
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
  Eigen::MatrixXd stencil= L2_stencil(r_ups, t_ups);
  Eigen::VectorXd RHS = f - bcs;

  //Solve L2_stencil u = RHS
  Eigen::VectorXd unknowns = stencil.lu().solve(RHS);
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
  std::cout << "r_pts.size(): " << r_pts.size() << std::endl;
  std::cout << "t_pts.size(): " << t_pts.size() << std::endl;

  pair<vector<double>, vector<double>> polar_pts = { r_pts, t_pts };
  return polar_pts;
}

extern "C" {
  double* polar_laplace(double R_min, double R_max, int N) {
    std::cout << "R_min: " << R_min << std::endl;
    std::cout << "R_max: " << R_max<< std::endl;
    std::cout << "N: " << N << std::endl;

    // Setup Internal Points
    std::vector<double> rlocs = utils::range(R_min, R_max, N, true);
    std::vector<double> tlocs = utils::range(0, 2*M_PI, N);

    vector<double> r;
    vector<double> t;
    std::tie(r,t) = polar_points(rlocs, tlocs);

    // F Function
    Eigen::VectorXd f = def_vector(rlocs.size() * tlocs.size(), zero);

    // Boundary Conditions
    // 1. R = R_min
    auto rmin_bc = [rlocs, R_min](int i) {
      if(std::isnan((int(i/rlocs.size()) == 0 ? 0.0 : 0.0))) std::cout << "rmin_bc nan\n";
      return (int(i/rlocs.size()) == 0 ? 0.0 : 0.0);
    };
    // 2. R = R_max
    auto rmax_bc = [rlocs, tlocs, t, R_max](int i) {
      int dim = rlocs.size() * tlocs.size();
      double t_val = t[i];
      if(std::isnan(int((dim-1 -i)/rlocs.size()) == 0 ? R_max*sin(t_val) : 0.0)) std::cout << "rmax_bc nan\n";
      return (int((dim-1 -i)/rlocs.size()) == 0 ? R_max*sin(t_val) : 0.0);
    };
    Eigen::VectorXd bcs = \
      def_vector(rlocs.size() * tlocs.size(), rmin_bc) +
      def_vector(rlocs.size() * tlocs.size(), rmax_bc);

    // Compute Laplacian
    Eigen::VectorXd L2 = L2_polar(r, t, f, bcs);
    std::cout << "L2 Size: " << L2.size() << std::endl;

    double* tmp = new double;
    Eigen::Map<Eigen::VectorXd>(tmp, L2.size())= L2;
    return tmp;
  }

  void freeme(char *ptr)
  {
    printf("freeing address: %p\n", ptr);
    free(ptr);
  }
}



int main() {
  std::cout << "Hello world\n";
  /*
  double const R_max = 5;
  double const R_min = 1;
  double const N = 10;
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
  std::cout << "L2.size(): " << L2.size() << std::endl;
  */

  double* sub_out = polar_laplace(1.0, 5.0, 50);
}
