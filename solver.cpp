#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <functional>
#include <string>

#include <Eigen/SparseCore>

#include "utils.h"
#include "grid.h"
#include "numerical-methods.h"

#pragma clang diagnostic ignored "-Wc++14-extensions"

using namespace std;
namespace nm = numerical_methods;

/** Initial Conditions for omega**/
double omega_bc(double const R_min, double const R_max,
                double const dr, int const N,
                Eigen::VectorXd full_r, Eigen::VectorXd full_phi,
                int row) {
  double const r = full_r[row];
  if (r == R_min) {
    // phi_1 = phi(R_min + dr, theta);
    double phi_1 = full_phi[row + N];
    return (double)((-2.0 * phi_1) / pow(dr, 2));
  } else if (r == R_max) {
    return 0.0;
  } else {
    return 0.0;
  }
}

/** Initial Conditions for Phi**/
double phi_bc(double const R_min, double r, double t) {
  //r = R_max -> exact solution
  //  r = R_min -> 0
  //  r = else  -> 0*/
  if (r == R_min){
    return 0.0;
  } else {
    return nm::exact_potential(r,t);
  }
}

double gaussian_omega_ic(double r, double theta) {
  return pow(M_E,
             (-1.0 * (pow((r - 2.5), 2) +
                      pow((theta - M_PI), 2))));
}

int closest_val_idx(Eigen::VectorXd vec, double target) {
    auto err = [](double val, double target) {return abs(target - val);};
    vector<double> err_vec(vec.rows());
    std::transform(vec.data(), vec.data() + vec.size(), err_vec.begin(),
                   [err, target](double val){ return err(val, target); });
    std::vector<double>::iterator err_min = std::min_element(err_vec.begin(), err_vec.end());
    int min_idx = std::distance(err_vec.begin(), err_min);
    return min_idx;
}


bool nan_trap(Eigen::VectorXd v, string vname) {
  vector<double> nans = {};
  std::copy_if(v.data(), v.data()+v.size(), back_inserter(nans),
               [vname](double el) {
                 return isnan(el);
               });
  if(nans.size() > 1) {std::cout << nans.size() << " nans in " << vname << std::endl;}
  return (nans.size() > 1);
}

class Solver {
  // Grid Variables
  Eigen::VectorXd full_r, full_theta;
  double const R_min;
  double const R_max;
  double dr;
  int const N;
  int const indim;
  int const fdim;

  // Fluid Variables
  Eigen::VectorXd ic_omega, ic_phi;
  Eigen::VectorXd spg_layer;

  // Derivative Stencils
  Eigen::SparseMatrix<double> D1, D1_periodic, L2_stencil;
  // ↓ LU Decomp of the L2 Stencil Matrix on interior points
  Eigen::SparseLU<Eigen::SparseMatrix<double>> L2_in_decomp;

  // Solver Variables
  int max_iters;
  int nplot;
  double dt;
  double Redh;
  string solf_name;

  pair<Eigen::VectorXd, Eigen::VectorXd> create_grid
  (double const R_min, double const R_max, int const N) {
    std::cout << "== Assembling Grid ==\n";
    return grid::all_points(R_min, R_max, N);
  }

  void write_grid(Eigen::VectorXd r, Eigen::VectorXd theta, string const csv_name) {
    std::cout << "Writing Grid\n\n";
    string abs_path =  std::__fs::filesystem::current_path();
    grid::write_grid(abs_path+"/"+csv_name, {r, theta}, true);
  }

  pair<Eigen::VectorXd, Eigen::VectorXd> apply_bcs
  (double const R_min, double const R_max, int const N,
   Eigen::VectorXd const full_r, Eigen::VectorXd const full_theta)
  {
    int const indim = N*N; //Dimension of interior points
    int const fdim = 2*N + indim; //Dimension of full domain
    // ↓ Populate phi
    Eigen::VectorXd full_phi(fdim);
    std::transform(full_r.data(), full_r.data() + full_r.size(),
                   full_theta.data(), full_phi.data(),
                   [R_min, R_max](double r, double theta) {
                     return phi_bc(R_min, r, theta);
                   });
    // ↓ Populate Omega
    //   NOTE: Omega BC's depend on the phi values,
    //         Applying omega BC's must happen AFTER appling phi BC's for correct results.

    /*
    auto marshalled_omega_bc = [R_min, R_max, dr=this->dr, N, full_r, full_phi](int row) {
      return omega_bc(R_min, R_max, dr, N, full_r, full_phi, row);
    };
    Eigen::VectorXd full_omega = nm::def_vector(fdim, marshalled_omega_bc);
    */

    Eigen::VectorXd full_omega(fdim);
    std::transform(full_r.data(), full_r.data() + full_r.size(),
                   full_theta.data(), full_omega.data(),
                   [](double r, double theta){ return gaussian_omega_ic(r, theta); });

    // Add purtabation at omega(R=2.5, θ=Pi)
    /*
    double const target_r = 2.5;
    double const target_theta = M_PI;
    int r_idx = closest_val_idx(full_r, target_r);
    int theta_idx = closest_val_idx(full_theta, target_theta);
    full_omega(r_idx + theta_idx) = 12;
    */

    return {full_omega, full_phi};
  }

  Eigen::VectorXd sponge_layer(Eigen::VectorXd full_r, double R_max) {
    auto sl = [R_max](int r) {
      return pow(fmax(r - (R_max - 2), 0)/2.0, 2);
    };
    Eigen::VectorXd s_layer(full_r.rows());
    std::transform(full_r.data(), full_r.data() + full_r.size(), s_layer.data(), sl);
    return s_layer;
  }

  // Leaving out BC functions until we know we need them
  // f is a vector in indim size
  Eigen::VectorXd L2_solve(Eigen::VectorXd f, Eigen::VectorXd bcs) {
    double scaling_value = ((1.0/pow(dr, 2) + ((1.0/(R_max - dr) * (1.0/(2.0*dr))))));
    return L2_in_decomp.solve(scaling_value * f - scaling_value * bcs);
  }

void solver_loop (int iters_2go, const Eigen::MatrixXd& full_omega, const Eigen::MatrixXd& full_phi) {
    std::cout << "-> Starting iteration: " << max_iters - iters_2go << std::endl;

    /** 1. Calculate radial and theta components of velocity
     * U_r = 1/r * dphi/dtheta, U_theta = -1 * dphi/dr
     **/
    Eigen::VectorXd r_inv(this->fdim);
    std::transform(full_r.data(), full_r.data() + full_r.size(), r_inv.data(),
                   [](double r){ return (double)(1.0/r); });
    Eigen::VectorXd Dphi_Dtheta= D1_periodic * full_phi;
    Eigen::VectorXd U_r = r_inv.cwiseProduct(Dphi_Dtheta);
    Eigen::VectorXd Dphi_Dr = D1 * full_phi;
    Eigen::VectorXd U_theta = -1.0 * Dphi_Dr;

    /** 2. Use U_r and U_theta to calculate omega@t=n+1
     * omega@t=n+1 = omega@t=n + dt * [1/Re_dh * L2(omega) - U_r * dw/dr - 1/r *
     * U_theta * dw/dtheta]
     **/
    // 2.a) domega terms
    Eigen::VectorXd Domega_Dtheta = D1_periodic * full_omega;
    Eigen::VectorXd Domega_Dr = D1 * full_omega;

    // 2.b) Laplacian of omega (L2*omega)
    Eigen::VectorXd grad2_omega = L2_stencil * full_omega;

    // 2.c) Calculate omega@n+1
    Eigen::VectorXd full_omega_np1 = full_omega + dt *
      (((1 / Redh) * grad2_omega) -
       U_r.cwiseProduct(Domega_Dr) -
       r_inv.cwiseProduct(U_theta).cwiseProduct(Domega_Dtheta) -
       (spg_layer.cwiseProduct(full_omega))
       );

    nan_trap(full_omega_np1, "full_omega_np1 - right after assignment");
    /** 3. Solve for in_phi@n+1
     *  grad2(in_phi@n+1) = -in_full_omega_np1
     **/
    Eigen::VectorXd phi_bcs(indim);
    phi_bcs.setZero();
    Eigen::VectorXd ub_r      = full_r(Eigen::seq(Eigen::last-N+1, Eigen::last));     //r, upper boundary points
    Eigen::VectorXd ub_theta  = full_theta(Eigen::seq(Eigen::last-N+1, Eigen::last)); //theta, upper boundary points
    std::transform(ub_r.data(), ub_r.data()+ ub_r.size(),
                   ub_theta.data(), phi_bcs.data()+phi_bcs.size() - N,
                   [R_min = this->R_min](double r, double theta) {
                     return phi_bc(R_min, r, theta);
                   });
    Eigen::VectorXd in_phi_np1 = L2_solve(-1.0*full_omega_np1(Eigen::seq(N,Eigen::last-N)), phi_bcs);
    nan_trap(in_phi_np1, "in_phi_np1");


    /** 4. Apply BC's **/
    // 4.a) Phi
    Eigen::VectorXd full_phi_np1(fdim);
    //    -> Lower BC
    std::transform(full_r.data(), full_r.data()+N,
                   full_theta.data(), full_phi_np1.data(),
                   [R_min = this->R_min](double r, double theta){
                     return phi_bc(R_min, r, theta);
                   });
    //    -> Interior Points
    full_phi_np1(Eigen::seq(N, Eigen::last-N)) = in_phi_np1;
    //    -> Upper BC
    //Eigen::VectorXd ub_r      = full_r(Eigen::seq(Eigen::last-N+1, Eigen::last));     //r, upper boundary points
    //Eigen::VectorXd ub_theta  = full_theta(Eigen::seq(Eigen::last-N+1, Eigen::last)); //theta, upper boundary points
    std::transform(ub_r.data(), ub_r.data() + ub_r.size(),
                   ub_theta.data(), full_phi_np1.data()+full_phi_np1.size()-N,
                   [R_min = this->R_min](double r, double theta){
                     return phi_bc(R_min, r, theta);
                   });

    // 4.b) Omega
    auto apply_omega_bc = [R_min = this->R_min, R_max=this->R_max, dr=this->dr,
                           N=this->N, full_r = this->full_r, full_phi] (int row) {
      return omega_bc(R_min, R_max, dr, N, full_r, full_phi, row);
    };
    //    -> Lower BC
    vector<int> rng = utils::irange(0, N);
    std::transform(rng.begin(), rng.end(), full_omega_np1.data(), apply_omega_bc);
    //    -> Upper BC
    rng = utils::irange(fdim-N, fdim);
    std::transform(rng.begin(), rng.end(),
                   full_omega_np1.data()+full_omega_np1.size()-N,
                   apply_omega_bc);

    /** 5. Write Solution
     * NOTE: Time/iteration needs to be introduced into the csv
     **/
    // Fix this so it points to iterations
    if ((max_iters-iters_2go)%nplot == 0 || iters_2go == 0) {
      // Path should probs be passed in as an absolute path but later
      std::cout << "Writing Solution @ iteration " << max_iters-iters_2go << std::endl;
      string absolute_path =  std::__fs::filesystem::current_path();
      grid::write_fluid_vars(absolute_path + "/" + solf_name,
                             {full_omega_np1, full_phi_np1, U_r, U_theta},
                             false);
    }

    if(nan_trap(full_omega_np1, "full_omega_np1")) return;

    /** Conditionally Recurse **/
    if (iters_2go == 0) {
      return;
    }
    else {
      solver_loop(iters_2go - 1, full_omega_np1, full_phi_np1);
    }
  }

public:
  Solver(double const R_min, double const R_max, int const N, string const gf_name) :
    R_min(R_min), R_max(R_max), N(N), indim(N*N), fdim(2*N + this->indim) //Setting Const Variables
  {
    // Create Grid
    std::tie(this->full_r, this->full_theta) = create_grid(R_min, R_max, N);
    write_grid(this->full_r, this->full_theta, gf_name);
    this->dr = grid::calc_dr(this->full_r);

    // Populate stencils once instead of once per iteration
    this->D1 = nm::D1(fdim);
    this->D1_periodic = nm::D1_periodic(this->fdim);
    this->L2_stencil = nm::L2_stencil(this->full_r, this->full_theta);
    Eigen::SparseMatrix<double> L2_in_stencil = nm::L2_stencil(this->full_r(Eigen::seq(N,Eigen::last-N)),
                                                               this->full_theta(Eigen::seq(N,Eigen::last-N)));
    this->L2_in_decomp.compute(L2_in_stencil);

    //Compute Sponge Layer
    this->spg_layer = sponge_layer(this->full_r, this->R_max);

    // Apply BC's
    std::tie(this->ic_omega, this-> ic_phi) = apply_bcs(R_min, R_max, N, this->full_r, this->full_theta);
  }

  ~Solver() {
    //IDK I think things are fine
  }

  void run(int const iterations, int const nplot,
           double const dt, double const Redh,
           string const solf_name) {
    this->max_iters = iterations;
    this->nplot     = nplot;
    this->dt        = dt;
    this->Redh      = Redh;
    this->solf_name = solf_name;

    string abs_path =  std::__fs::filesystem::current_path();
    grid::write_csv_header(abs_path + "/" + solf_name, grid::fluid_header);
    solver_loop(iterations-1, this->ic_omega, this->ic_phi);
  }

};

void cpp_version() {
  if (__cplusplus == 202101L) std::cout << "C++23";
  else if (__cplusplus == 202002L) std::cout << "C++20";
  else if (__cplusplus == 201703L) std::cout << "C++17";
  else if (__cplusplus == 201402L) std::cout << "C++14";
  else if (__cplusplus == 201103L) std::cout << "C++11";
  else if (__cplusplus == 199711L) std::cout << "C++98";
  else std::cout << "pre-standard C++." << __cplusplus;
  std::cout << "\n";
}

void solver_controller() {
  // Grid Variables
  int const N = 100;
  double const R_min = 1.0;
  double const R_max = 20.0;

  // Solver Variables
  double total_run_time = 60;
  double const dt = 0.01;
  int const iterations = total_run_time/dt;

  int const nplot = 100; //Just print at the end
  double const Redh = 23.0;

  Solver solver(R_min, R_max, N, "grid.csv");
  solver.run(iterations, nplot, dt, Redh, "solver-data.csv");
}

int main(){
  std::cout << "Hello world\n";
  cpp_version();
  solver_controller();
}
