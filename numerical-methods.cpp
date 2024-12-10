#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
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
auto one  = [](int row){return 1.0;};

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
Eigen::SparseMatrix<double> L2_stencil(std::vector<double> r_ups, std::vector<double> t_ups) {
  //TODO: Add a check that r_ups and t_ups are the same length
  int nup = r_ups.size(); //Number of unknown points
  int dim = std::sqrt(nup);

  Eigen::SparseMatrix<double> Stencil(nup,nup);
  Stencil.reserve(Eigen::VectorXi::Constant(nup, 5));
  // ↑ Because we have a 5 point stencil, most columns will have 5 non zero values

  //Start Lambda Definitions
  double dr = calc_dr(r_ups);
  double dt = t_ups[1] - t_ups[0];

  auto center = [r_ups, dr, dt](int row) {
    // -2 * (1/dr^2 + 1/(r^2 dt^2))
    return (double)(-2.0 * ((1.0/pow(dr, 2)) + (1.0/(pow(r_ups[row],2) * pow(dt,2)))));
  };

  auto l_and_r= [r_ups, dt](int row) {
    // 1 / (r^2 dt^2)
    double r = r_ups[row];
    return (double)(1.0 / (pow(r, 2) * pow(dt, 2)));
  };

  auto upper= [r_ups, dr](int row) {
    // (1/dr^2 + 1/(r 2 dr))
    double r = r_ups[row];
    return (double)( (1.0/pow(dr,2)) + (1.0/(2.0 * r * dr)) );
  };

  auto lower= [r_ups, dr](int row) {
    // (1/dr^2 - 1/(r 2 dr))
    double r = r_ups[row];
    return (double)( (1.0/pow(dr,2)) - (1.0/(2 * r * dr)) );
  };
  //End Lambda Definitions

  for (int i = 0; i < nup; i++){
    //Center Function
    Stencil.coeffRef(i,i) = center(i);
    //Left & Right
    switch((i+1) % dim) {
    case 1: //First row of T Matrix
      Stencil.coeffRef(i, i+1) = l_and_r(i);
      Stencil.coeffRef(i, i + (dim-1)) = l_and_r(i);
      break;
    case 0: //Last row of T matrix
      Stencil.coeffRef(i, i-1) = l_and_r(i);
      Stencil.coeffRef(i, i - (dim-1)) = l_and_r(i);
      break;
    default: //All ofther rows of T matrix
      Stencil.coeffRef(i, i-1) = l_and_r(i);
      Stencil.coeffRef(i, i+1) = l_and_r(i);
    }
    //Upper Function
    if(i + dim < nup) Stencil.coeffRef(i, i + dim) = upper(i);
    //Lower Function
    if(i - dim >= 0)  Stencil.coeffRef(i, i - dim) = lower(i);
  }

  Stencil.makeCompressed();
  return Stencil;
}

Eigen::VectorXd L2_polar(std::vector<double> r_ups, std::vector<double> t_ups, Eigen::VectorXd f, Eigen::VectorXd bcs) {
  Eigen::SparseMatrix<double> stencil= L2_stencil(r_ups, t_ups);
  Eigen::VectorXd RHS = f - bcs;

  //Solve L2_stencil u = RHS
  Eigen::SparseLU<Eigen::SparseMatrix<double>> lu_decomp;
  lu_decomp.compute(stencil);
  Eigen::VectorXd unknowns = lu_decomp.solve(RHS);
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

double exact_potential(double r, double t) {
  double const R = 1.0;
  double const U = 1.0;
  return ( (r - (pow(R,2)/r)) * U * sin(t) );
}

extern "C" {
  /*
   * BC_type selects the BC at R_max:
   *   0 = R_max * sin(t)
   *   1 = Exact Solution
   */
  double* polar_laplace(double R_min, double R_max, int N, int BC_type) {
    std::cout << "Laplace Settings: \n";
    std::cout << "  R_min: " << R_min << std::endl;
    std::cout << "  R_max: " << R_max<< std::endl;
    std::cout << "  N: " << N << std::endl;
    std::cout << "  BC_type: " << (BC_type == 0 ? "sin" : "exact") << std::endl << std::endl;

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
      return (int(i/rlocs.size()) == 0 ? 0.0 : 0.0);
    };
    // 2. R = R_max
    double const dr = calc_dr(r);
    double scaling_value = ((1.0/pow(dr, 2) + ((1.0/(R_max - dr) * (1.0/(2.0*dr))))));

    std::function<double(int)> sin_bc = [rlocs, tlocs, t, R_max, scaling_value](int i) {
      int dim = rlocs.size() * tlocs.size();
      double t_val = t[i];
      if(int((dim-1 -i)/rlocs.size()) == 0) {
      }
      return (int((dim-1 -i)/rlocs.size()) == 0 ? scaling_value * R_max*sin(t_val) : 0.0);
      };

    std::function<double(int)> exact_bc = [rlocs, tlocs, t , scaling_value, R_max](int i) {
      int dim = rlocs.size() * tlocs.size();
      double t_val = t[i];
      return (int((dim-1 -i)/rlocs.size()) == 0 ? scaling_value * exact_potential(R_max,t_val) : 0.0);
    };
    std::function<double(int)> rmax_bc = (BC_type == 1 ? exact_bc : sin_bc);

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

  double* sub_out = polar_laplace(1.0, 5.0, 100,1);
}
