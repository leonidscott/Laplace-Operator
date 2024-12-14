
void solver(){
  // Create Grid
  int const N = 3;
  double const R_min = 1.0;
  double const R_max = 15.0;

  // Construct Interior Points
  vector<double> irp; //Interior Radius Points
  vector<double> ithetap; //Interior Theta Points
  std::tie(irp,ithetap) = interior_points(R_min, R_max, N);
  int const idim = irp.size(); //Interior Problem Dimension

  // Construct Boundary Points
  vector<double> bl_rp; //Lower Boundary Radius Points
  vector<double> bl_thetap; //Lower Boundary Radius Points
  std::tie(bl_rp, bl_thetap) = lower_bc_points(R_min, N);

  vector<double> bu_rp;
  vector<double> bu_thetap;
  std::tie(bu_rp, bu_thetap) = upper_bc_points(R_max, N);

  // Set Initial Conditions
  // -> Set at interior points
  std::vector<double> i_omega(idim); //Interior Omega Initial Condition Points
  std::transform(irp.begin(), irp.end(),
                 ithetap.begin(), i_omega.begin(),
                 [R_min, R_max](double r_val, double t_val){
                   return omega_bc(R_min, R_max, r_val, t_val);
                 });
  std::vector<double> i_phi(idim); //Interior Phi Initial Condition Points
  std::transform(irp.begin(), irp.end(),
                 ithetap.begin(), i_phi.begin(),
                 [R_min, R_max](double r_val, double t_val){
                   return omega_bc(R_min, R_max, r_val, t_val);
                   });
  // -> Set at boundary points
  std::vector<double> bl_omega(N); //Lower Boundary Omega IC Points
  std::vector<double> bl_phi(N);   //Lower Boundary Phi IC Points
  std::transform(bl_rp.begin(), bl_rp.end(),
                 bl_thetap.begin(), bl_omega.begin(),
                 [R_min, R_max](double rv, double thetav){
                   return omega_bc(R_min, R_max, rv, thetav);
                 });
  std::transform(bl_rp.begin(), bl_rp.end(),
                 bl_thetap.begin(), bl_phi.begin(),
                 [R_min, R_max](double rv, double thetav){
                   return phi_ic(R_min, R_max, rv, thetav);
                 });
  std::vector<double> bu_omega(N); //Lower Boundary Omega IC Points
  std::vector<double> bu_phi(N);   //Lower Boundary Phi IC Points
  std::transform(bu_rp.begin(), bu_rp.end(),
                 bu_thetap.begin(), bu_omega.begin(),
                 [R_min, R_max](double rv, double thetav){
                   return omega_bc(R_min, R_max, rv, thetav);
                 });
  std::transform(bu_rp.begin(), bu_rp.end(),
                 bu_thetap.begin(), bu_phi.begin(),
                 [R_min, R_max](double rv, double thetav){
                   return phi_ic(R_min, R_max, rv, thetav);
                 });

  // Assemble Full Solution
  std::vector<double> all_r;
  all_r.insert(all_r.end(), bl_rp.begin(), bl_rp.end());
  all_r.insert(all_r.end(), irp.begin(), irp.end());
  all_r.insert(all_r.end(), bu_rp.begin(), bu_rp.end());

  std::vector<double> all_theta;
  all_theta.insert(all_theta.end(), bl_thetap.begin(), bl_thetap.end());
  all_theta.insert(all_theta.end(), ithetap.begin(), ithetap.end());
  all_theta.insert(all_theta.end(), bu_thetap.begin(), bu_thetap.end());

  std::vector<double> all_omega;
  all_omega.insert(all_omega.end(), bl_omega.begin(), bl_omega.end());
  all_omega.insert(all_omega.end(), i_omega.begin(), i_omega.end());
  all_omega.insert(all_omega.end(), bu_omega.begin(), bu_omega.end());

  std::vector<double> all_phi;
  all_phi.insert(all_phi.end(), bl_phi.begin(), bl_phi.end());
  all_phi.insert(all_phi.end(), i_phi.begin(), i_phi.end());
  all_phi.insert(all_phi.end(), bu_phi.begin(), bu_phi.end());

  int const fdim = all_r.size();

  // Begin Solver Loop
  solver_loop(R_min, R_max, idim,
              all_r, all_theta, all_omega, all_phi);

  // Write Solution
  //TODO: Delete U_r and U_theta?
  vector<double> U_r = nm::def_std_vector(fdim, nm::zero);
  vector<double> U_theta = nm::def_std_vector(fdim, nm::zero);

  string abs_path =  std::__fs::filesystem::current_path();
  write_to_csv(abs_path + "/" + "solver-data.csv",
               {all_r, all_theta, all_omega, all_phi, U_r, U_theta},
               true);
}
