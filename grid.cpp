#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <Eigen/SparseCore>
#include "utils.h"
using namespace std;

namespace grid {
  /*CSV File Writing Functions*/
  extern const string fluid_header = "omega,phi,U_r,U_theta";

  void write_csv_header(string csv_name, string header) {
    std::ofstream csv(csv_name);
    csv << header << std::endl;
    csv.close();
  }

  /** CAUTION: Make sure to pass in absolute paths of csv file **/
  //TODO: Remove grid info, eventually add in iteration + time info
  void write_fluid_vars(string csv_name,
                        tuple< Eigen::VectorXd, Eigen::VectorXd,
                        Eigen::VectorXd, Eigen::VectorXd> data,
                        bool overwrite) {
    Eigen::VectorXd omega, phi, U_r, U_theta;
    std::tie(omega, phi, U_r, U_theta) = data;

    vector<double> vomega(omega.data(), omega.data() + omega.rows() * omega.cols());
    vector<double> vphi(phi.data(), phi.data() + phi.rows() * phi.cols());
    vector<double> vU_r(U_r.data(), U_r.data() + U_r.rows() * U_r.cols());
    vector<double> vU_theta(U_theta.data(), U_theta.data() + U_theta.rows() * U_theta.cols());

    std::ofstream csv(csv_name, (overwrite? std::ios::out : std::ios::app));
    // Write Header Row
    if(overwrite) csv << grid::fluid_header << std::endl;

    for(int i = 0; i < vomega.size() ; i++) {
      string line =
        std::to_string(vomega[i]) + ", " + std::to_string(vphi[i]) + ", " +
        std::to_string(vU_r[i])   + ", " + std::to_string(vU_theta[i]);
      csv << line << std::endl;
    }

    csv.close();
  }

  void write_grid(string csv_name,
                  tuple<Eigen::VectorXd, Eigen::VectorXd> data,
                  bool overwrite) {
    Eigen::VectorXd r, theta;
    std::tie(r, theta) = data;

    vector<double> vr(r.data(), r.data() + r.rows() * r.cols());
    vector<double> vtheta(theta.data(), theta.data() + theta.rows() * theta.cols());

    std::ofstream csv(csv_name, (overwrite? std::ios::out : std::ios::app));
    // Write Header Row
    const string grid_header="r,theta";
    if(overwrite) csv << grid_header << std::endl;

    for(int i = 0; i < r.rows() ; i++) {
      string line = std::to_string(vr[i]) + ", " + std::to_string(vtheta[i]);
      csv << line << std::endl;
    }
    csv.close();
  }

  /*CSV File Writing Functions*/

  /** Radius Major Points given radius and angle locations**/
  pair<vector<double>, vector<double>> polar_points(vector<double> r_locs,
                                                    vector<double> t_locs) {
    vector<double> r_pts(r_locs.size() * t_locs.size());
    vector<double> t_pts(r_locs.size() * t_locs.size());

    for (int r = 0; r < r_locs.size(); r++) {
      for (int t = 0; t < t_locs.size(); t++) {
        r_pts[r * t_locs.size() + t] = r_locs[r];
        t_pts[r * t_locs.size() + t] = t_locs[t];
      }
    }
    std::cout << "r_pts.size(): " << r_pts.size() << std::endl;
    std::cout << "t_pts.size(): " << t_pts.size() << std::endl;

    pair<vector<double>, vector<double>> polar_pts = {r_pts, t_pts};
    return polar_pts;
  }

  double calc_dr(vector<double> r_pts) {
    double r_init = r_pts[0];
    for (int i = 0; i < r_pts.size(); i++) {
      if (r_pts[i] != r_init)
        return (r_pts[i] - r_init);
    }
    return 0.0;
  }

  double calc_dr(Eigen::VectorXd r_pts) {
    double r_init = r_pts(0);
    for (int i = 0; i < r_pts.rows(); i++) {
      if (r_pts(i) != r_init) { return (r_pts(i) - r_init); }
    }
    return 0.0;
  }

  /*Internal Functions*/
  pair<vector<double>, vector<double>> lower_bc_points(double R_min, int N) {
    std::vector<double> rlocs(N);
    std::transform(rlocs.begin(), rlocs.end(), rlocs.begin(),
                   [R_min](double rval) { return R_min; });
    std::vector<double> tlocs = utils::range(0, 2 * M_PI, N);
    return {rlocs, tlocs};
  }
  std::tuple<vector<double>, vector<double>> interior_points(double R_min,
                                                        double R_max, int N) {
    std::vector<double> rlocs = utils::range(R_min, R_max, N, true);
    std::vector<double> tlocs = utils::range(0, 2 * M_PI, N);

    vector<double> r;
    vector<double> t;
    std::tie(r, t) = polar_points(rlocs, tlocs);

    return {r, t};
  }
  pair<vector<double>, vector<double>> upper_bc_points(double R_max, int N) {
    std::vector<double> rlocs(N);
    std::transform(rlocs.begin(), rlocs.end(), rlocs.begin(),
                   [R_max](double rval){ return R_max; });
    std::vector<double> tlocs = utils::range(0, 2*M_PI, N);
    return {rlocs, tlocs};
  }
  /*Internal Functions*/

  /** Returns a tuple of r and theta points in Eigen Vectors.
     These points include interior points, lower, and upper boundary conditions
     **/
  tuple<Eigen::VectorXd, Eigen::VectorXd> all_points(double R_min, double R_max,
                                                     int N) {
    // Get R, θ Points
    vector<double> lb_r, lb_theta;
    vector<double> in_r, in_theta;
    vector<double> ub_r, up_theta;
    std::tie(lb_r, lb_theta) = lower_bc_points(R_min, N);
    std::tie(in_r, in_theta) = interior_points(R_min, R_max, N);
    std::tie(ub_r, up_theta) = upper_bc_points(R_max, N);

    // Assemble full R, θ vectors
    std::vector<double> full_r;
    full_r.insert(full_r.end(), lb_r.begin(), lb_r.end());
    full_r.insert(full_r.end(), in_r.begin(), in_r.end());
    full_r.insert(full_r.end(), ub_r.begin(), ub_r.end());

    std::vector<double> full_theta;
    full_theta.insert(full_theta.end(), lb_theta.begin(), lb_theta.end());
    full_theta.insert(full_theta.end(), in_theta.begin(), in_theta.end());
    full_theta.insert(full_theta.end(), up_theta.begin(), up_theta.end());

    // Convert to eigen
    Eigen::VectorXd eigen_full_r =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(full_r.data(), full_r.size());
    Eigen::VectorXd eigen_full_theta =
        Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(full_theta.data(),
                                                      full_theta.size());
    return { eigen_full_r, eigen_full_theta };
  }
}

