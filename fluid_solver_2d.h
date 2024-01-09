#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

double square(double x) { return x * x; }

class fluid_solver_2d {
public:
  fluid_solver_2d(double Lx0, double Ly0, int Nx0, int Ny0, double gamma0,
                  double output_dt0)
      : gamma(gamma0), Lx(Lx0), Ly(Ly0), Nx(Nx0), Ny(Ny0) {
    output_dt = output_dt0;
    dx = Lx / (Nx - 2);
    dy = Ly / (Ny - 2);

    rho.resize(Nx * Ny);
    vx.resize(Nx * Ny);
    vy.resize(Nx * Ny);
    P.resize(Nx * Ny);

    mass.resize(Nx * Ny);
    mom_x.resize(Nx * Ny);
    mom_y.resize(Nx * Ny);
    energy.resize(Nx * Ny);

    rho_tmp.resize(Nx * Ny);
    vx_tmp.resize(Nx * Ny);
    vy_tmp.resize(Nx * Ny);
    P_tmp.resize(Nx * Ny);

    rho_Lx.resize(Nx * Ny);
    rho_Rx.resize(Nx * Ny);
    rho_Ly.resize(Nx * Ny);
    rho_Ry.resize(Nx * Ny);

    vx_Lx.resize(Nx * Ny);
    vx_Rx.resize(Nx * Ny);
    vx_Ly.resize(Nx * Ny);
    vx_Ry.resize(Nx * Ny);
    vy_Lx.resize(Nx * Ny);
    vy_Rx.resize(Nx * Ny);
    vy_Ly.resize(Nx * Ny);
    vy_Ry.resize(Nx * Ny);
    P_Lx.resize(Nx * Ny);
    P_Rx.resize(Nx * Ny);
    P_Ly.resize(Nx * Ny);
    P_Ry.resize(Nx * Ny);

    mass_flux_x.resize(Nx * Ny);
    mass_flux_y.resize(Nx * Ny);
    momx_flux_x.resize(Nx * Ny);
    momx_flux_y.resize(Nx * Ny);
    momy_flux_x.resize(Nx * Ny);
    momy_flux_y.resize(Nx * Ny);
    energy_flux_x.resize(Nx * Ny);
    energy_flux_y.resize(Nx * Ny);
  }

  ~fluid_solver_2d() {}

  
  //They are related to mass
  void primitive_to_conserved() {
    for(int i=1; i<Nx-1; i++){
      for(int j =1; j<Ny-1;j++){
        int idx = i+Nx*j;
        mass[idx] = rho[idx] * dx * dy;
        mom_x[idx]= rho[idx]*vx[idx]*dx*dy;
        mom_y[idx]=rho[idx]*vy[idx]*dx*dy;
        double u = P[idx]/((gamma-1)*rho[idx]); 
        energy[idx] = rho[idx]*dx*dy*(0.5*(vx[idx]*vx[idx]+vy[idx]*vy[idx])+u);
      }
    } 
    periodic_boundary(mass);
    periodic_boundary(mom_x);
    periodic_boundary(mom_y);
    periodic_boundary(energy);
  }

  //primative values are related to density
  void conserved_to_primitive() {
    for(int i=1; i<Nx-1; i++){
      for(int j =1; j<Ny-1;j++){
        int idx = i+Nx*j;
        rho[idx] = mass[idx] / (dx * dy);
        vx[idx] = mom_x[idx]/(rho[idx]*dx*dy);
        vy[idx] = mom_y[idx]/(rho[idx]*dx*dy);
        double u = (energy[idx]/(dy*dx))/rho[idx] - 0.5*(vx[idx]*vx[idx]+vy[idx]*vy[idx]);
        P[idx] = (gamma-1)*rho[idx]*u;
      }
    }
    periodic_boundary(rho);
    periodic_boundary(vx);
    periodic_boundary(vy);
    periodic_boundary(P);
  }

  void init(const std::vector<double> &rho0, const std::vector<double> &vx0,
            const std::vector<double> &vy0, const std::vector<double> &P0) {
    // TODO: Initialize the primitive variables using the given initial
    // condition
    rho = rho0;
    vx=vx0;
    vy=vy0;
    P=P0;
    primitive_to_conserved();

  }

  double find_dt() {
    // TODO: Find the optimal dt that satisfies the CFL condition, and return
    // its value
    double Ccfl = 0.4;
    double numerator = Ccfl * std::min(dx, dy);
    double denominator = std::sqrt(gamma * P[0] / rho[0]) + std::sqrt(vx[0] * vx[0] + vy[0] * vy[0]);
    for(int j = 1; j < Ny-1; j++) {
      for(int i = 1; i < Nx-1; i++) {
        int idx = i + j * Nx;
        double new_denom = std::sqrt(gamma * P[idx] / rho[idx]) + std::sqrt(vx[idx] * vx[idx] + vy[idx] * vy[idx]);
        denominator = std::max(denominator,new_denom);
      }
    }
    return numerator/denominator;
  }

  void solve(double t0, double t_end) {
    // Solve the fluid equations, starting from t0 and stoping at t_end
    double t = t0;
    int n = 0; // n labels the output file
    while (t < t_end) {
      //std::cout<<n<<std::endl;
      if (t >= output_dt * n) {
        //std::cout<<"hey"<<std::endl;
        output(n);
        n += 1;
      }
      double dt = find_dt();
      //std::cout<<dt<<std::endl;
      step(dt);
      t += dt;
    }
  }

  void step(double dt) {
    // extrapolate a half step in time using primitive equations

    primitive_update(0.5 * dt);

    // compute fluxes
    compute_fluxes();

    // update solultion from fluxes
    update_conserved(dt);

    // update primitive variables
    conserved_to_primitive();
  }

  void periodic_boundary(std::vector<double> &f) {
    // TODO: apply periodic boundary conditions to an array f
    //f(0, j) = f (Nx − 2, j), f (Nx − 1, j) = f (1, j)
    //f(i, 0) = f (i, Ny − 2), f (i, Ny − 1) = f (i, 1)
    for(int j = 0; j < Ny; j++) {
      f[0 + j * Nx] = f[(Nx - 2) + j * Nx];
      f[(Nx - 1) + j * Nx] = f[1 + j * Nx];
    }

    for(int i = 0; i < Nx; i++) {
      f[i + 0 * Nx] = f[i + (Ny - 2) * Nx];
      f[i + (Ny - 1) * Nx] = f[i + 1 * Nx];
    }
  }

  void primitive_update(double dt) {
    // TODO: update the primitive variables using Euler equations in primitive
    // form using an FTCS scheme

    //first half time step
    for(int  i = 1; i < Nx - 1; i++) {
      for(int j = 1; j < Ny - 1; j++) {
        int idx = i + j * Nx;
        double drho_dx = (rho[idx+1] - rho[idx-1])/(2.0*dx);
        double drho_dy = (rho[idx + Nx] - rho[idx - Nx])/(2.0*dy);
        double dvx_dx = (vx[idx+1] - vx[idx-1])/(2.0*dx);
        double dvy_dy = (vy[idx + Nx] - vy[idx - Nx])/(2.0*dy); 
        rho_tmp[idx] = rho[idx] - dt *(vx[idx]*drho_dx + vy[idx]*drho_dy + rho[idx]*(dvx_dx + dvy_dy));
        double dP_dx = (P[idx+1] - P[idx-1])/(2.0*dx);
        double dvx_dy = (vx[idx+Nx] - vx[idx-Nx])/(2.0*dy);
        vx_tmp[idx] = vx[idx] - dt/(rho[idx])*(dP_dx + vx[idx]*dvx_dx + vy[idx]*dvx_dx); 
        double dvy_dx = (vy[idx+1] - vy[idx-1])/(2.0*dx);
        double dP_dy = (P[idx+Nx] - P[idx-Nx])/(2.0*dy);
        vy_tmp[idx] = vy[idx]- dt/(rho[idx])*(dP_dy+ vy[idx]*dvy_dy + vx[idx]*dvy_dx);
        P_tmp[idx] = P[idx] - dt*(gamma*P[idx]*(dvx_dx + dvy_dy)+ vy[idx]*dP_dy + vx[idx]*dP_dx); 
      }
    }
    periodic_boundary(rho_tmp);
    periodic_boundary(vx_tmp);
    periodic_boundary(vy_tmp);
    periodic_boundary(P_tmp);
    rho=rho_tmp;
    vx=vx_tmp;
    vy=vy_tmp;
    P=P_tmp;

  }

  void extrapolate_to_interface() {
    // TODO: compute rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, and P_R here
    
    for(int i = 1; i < Nx - 1; i++) {
      for(int j = 1; j < Ny - 1; j++) {
        int index = i + j * Nx;
        //compute rho_Lx, rho_Rx, rho_Ly, rho_Ry
        rho_Lx[index] = rho[index] - 0.25 * (rho[index+1] - rho[index-1]);
        rho_Rx[index] = rho[index] + 0.25 * (rho[index+1] - rho[index-1]);
        rho_Ly[index] = rho[index] - 0.25 * (rho[index+Nx] - rho[index-Nx]);
        rho_Ry[index] = rho[index] + 0.25 * (rho[index+Nx] - rho[index-Nx]);

        //compute vx_Lx, vx_Rx, vx_Ly, vx_Ry
        vx_Lx[index] = vx[index] - 0.25 * (vx[index+1] - vx[index-1]);
        vx_Rx[index] = vx[index] + 0.25 * (vx[index+1] - vx[index-1]);
        vx_Ly[index] = vx[index] - 0.25 * (vx[index+Nx] - vx[index-Nx]);
        vx_Ry[index] = vx[index] + 0.25 * (vx[index+Nx] - vx[index-Nx]);

        //compute vy_Lx, vy_Rx, vy_Ly, vy_Ry
        vy_Lx[index] = vy[index] - 0.25 * (vy[index+1] - vy[index-1]);
        vy_Rx[index] = vy[index] + 0.25 * (vy[index+1] - vy[index-1]);
        vy_Ly[index] = vy[index] - 0.25 * (vy[index+Nx] - vy[index-Nx]);
        vy_Ry[index] = vy[index] + 0.25 * (vy[index+Nx] - vy[index-Nx]);

        //compute P_Lx, P_Rx, P_Ly, P_Ry
        P_Lx[index] = P[index] - 0.25 * (P[index+1] - P[index-1]);
        P_Rx[index] = P[index] + 0.25 * (P[index+1] - P[index-1]);
        P_Ly[index] = P[index] - 0.25 * (P[index+Nx] - P[index-Nx]);
        P_Ry[index] = P[index] + 0.25 * (P[index+Nx] - P[index-Nx]);
      }
    }
    periodic_boundary(rho_Lx);
    periodic_boundary(rho_Rx);
    periodic_boundary(rho_Ly);
    periodic_boundary(rho_Ry);
    periodic_boundary(vx_Lx);
    periodic_boundary(vx_Rx);
    periodic_boundary(vx_Ly);
    periodic_boundary(vx_Ry);
    periodic_boundary(vy_Lx);
    periodic_boundary(vy_Rx);
    periodic_boundary(vy_Ly);
    periodic_boundary(vy_Ry);
    periodic_boundary(P_Lx);
    periodic_boundary(P_Ly);
    periodic_boundary(P_Rx);
    periodic_boundary(P_Ry);

    
  }

  void compute_fluxes() {
    // TODO: compute the fluxes
    extrapolate_to_interface();
    for(int i = 1; i < Nx-1; i++){
      for(int j = 1; j < Ny-1; j++){
      int idx= i + j*Nx;
      //finding vmax
      // double vmax_Rx = std::sqrt(gamma*P_Lx[idx+1]/rho_Lx[idx+1])+std::sqrt(vx_Lx[idx+1]*vx_Lx[idx+1]+vy_Lx[idx+1]*vy_Lx[idx+1]);
      // double vmax_Lx = std::sqrt(gamma*P_Rx[idx]/rho_Rx[idx])+std::sqrt(vx_Rx[idx]*vx_Rx[idx]+vy_Rx[idx]*vy_Rx[idx]);
      double vmax_Rx = std::sqrt(gamma*P_Lx[idx+1]/rho_Lx[idx+1]) + std::abs(vx_Lx[idx+1]);
      double vmax_Lx = std::sqrt(gamma*P_Rx[idx]/rho_Rx[idx]) + std::abs(vx_Rx[idx]);
      double vmax_x= std::max(vmax_Lx,vmax_Rx);
      // double vmax_Ry = std::sqrt(gamma*P_Ly[idx+Nx]/rho_Ly[idx+Nx])+std::sqrt(vx_Ly[idx+Nx]*vx_Ly[idx+Nx]+vy_Ly[idx+Nx]*vy_Ly[idx+Nx]);
      // double vmax_Ly = std::sqrt(gamma*P_Ry[idx]/rho_Ry[idx])+std::sqrt(vx_Ry[idx]*vx_Ry[idx+Nx]+vy_Ry[idx]*vy_Ry[idx]);
      double vmax_Ry = std::sqrt(gamma*P_Ly[idx+Nx]/rho_Ly[idx+Nx]) + std::abs(vy_Ly[idx+Nx]);
      double vmax_Ly = std::sqrt(gamma*P_Ry[idx]/rho_Ry[idx]) + std::abs(vy_Ry[idx]);
      double vmax_y= std::max(vmax_Ly,vmax_Ry);
      // mass flux in x
      double F_Lx = rho_Rx[idx]*vx_Rx[idx];
      double F_Rx = rho_Lx[idx+1]*vx_Lx[idx+1];
      double Q_Lx = rho_Rx[idx];
      double Q_Rx = rho_Lx[idx+1];
      mass_flux_x[idx]= (0.5)*(F_Lx+F_Rx) - (vmax_x/2.0)*(Q_Rx-Q_Lx);
      //std::cout<<mass_flux_x[idx]<<std::endl;
      // mass flux in y
      double F_Ly = rho_Ry[idx]*vy_Ry[idx];
      double F_Ry = rho_Ly[idx+Nx]*vy_Ly[idx+Nx];
      double Q_Ly = rho_Ry[idx];
      double Q_Ry = rho_Ly[idx+Nx];
      mass_flux_y[idx]= (0.5)*(F_Ly+F_Ry) - (vmax_y/2.0)*(Q_Ry-Q_Ly);
      // momentum_x flux in x
      F_Lx = rho_Rx[idx]*vx_Rx[idx]*vx_Rx[idx] + P_Rx[idx];
      F_Rx = rho_Lx[idx+1]*vx_Lx[idx+1]*vx_Lx[idx+1] + P_Lx[idx+1];
      Q_Lx = rho_Rx[idx]*vx_Rx[idx]; //Q_Lx = mom_x to the right
      Q_Rx= rho_Lx[idx+1]*vx_Lx[idx+1]; //Q_Rx = mom_x to the left
      momx_flux_x[idx]= (0.5)*(F_Lx+F_Rx) - (vmax_x/2.0)*(Q_Rx-Q_Lx);
      //momentum_x flux in y
      F_Ly = rho_Ry[idx]*vx_Ry[idx]*vy_Ry[idx];
      F_Ry = rho_Ly[idx+Nx]*vx_Ly[idx+Nx] *vy_Ly[idx+Nx];
      Q_Ly = rho_Ry[idx]*vx_Ry[idx];
      Q_Ry= rho_Ly[idx+Nx]*vx_Ly[idx+Nx];
      momx_flux_y[idx]= (0.5)*(F_Ly+F_Ry) - (vmax_y/2.0)*(Q_Ry-Q_Ly);
      //momentum_y flux in x
      F_Lx = rho_Rx[idx]*vx_Rx[idx]*vy_Rx[idx]; 
      F_Rx = rho_Lx[idx+1]*vx_Lx[idx+1]*vy_Lx[idx+1];
      Q_Lx = rho_Rx[idx]*vy_Rx[idx]; //Q_Lx = mom_y to the right
      Q_Rx = rho_Lx[idx+1]*vy_Lx[idx+1]; // Q_Rx = mom_y to the left 
      momy_flux_x[idx] = (0.5)*(F_Lx +F_Rx) -(vmax_x/2.0)*(Q_Rx - Q_Lx); 
      // momentum_y flux in y
      F_Ly = (rho_Ry[idx]*vy_Ry[idx]*vy_Ry[idx]) + P_Ry[idx]; 
      F_Ry = (rho_Ly[idx+Nx]*vy_Ly[idx+Nx] * vy_Ly[idx+Nx])+ P_Ly[idx+Nx]; 
      Q_Ly = rho_Ry[idx]*vy_Ry[idx]; //Q_Lx = mom_y to the right
      Q_Ry = rho_Ly[idx+Nx]*vy_Ly[idx+Nx]; // Q_Rx = mom_y to the left
      momy_flux_y[idx] = (0.5)*(F_Ly+F_Ry) - (vmax_y/2.0)*(Q_Ry - Q_Ly);
      //energy flux in x --> stuck-ish here 
      double u_Lx = P_Rx[idx]/((gamma-1.0)*rho_Rx[idx]);
      double u_Rx = P_Lx[idx+1]/((gamma-1.0)*rho_Lx[idx+1]);
      double U_Lx = rho_Rx[idx]*(0.5*(vx_Rx[idx]*vx_Rx[idx]+vy_Rx[idx]*vy_Rx[idx])+ u_Lx); 
      double U_Rx = rho_Lx[idx+1]*(0.5*(vx_Lx[idx+1]*vx_Lx[idx+1]+vy_Lx[idx+1]*vy_Lx[idx+1])+u_Rx); 
      F_Lx = (U_Lx + P_Rx[idx])*vx_Rx[idx];
      F_Rx = (U_Rx + P_Lx[idx+1])*vx_Lx[idx+1];
      Q_Lx = U_Lx; 
      Q_Rx = U_Rx;
      energy_flux_x[idx]= (0.5)*(F_Lx +F_Rx) -(vmax_x/2.0)*(Q_Rx - Q_Lx); 
      //energy flux in y 
      double u_Ly = P_Ry[idx]/((gamma-1.0)*rho_Ry[idx]);
      double u_Ry = P_Ly[idx+Nx]/((gamma-1.0)*rho_Ly[idx+Nx]);
      double U_Ly = rho_Ry[idx]*(0.5*(vx_Ry[idx]*vx_Ry[idx]+vy_Ry[idx]*vy_Ry[idx])+u_Ly);
      double U_Ry = rho_Ly[idx+Nx]*(0.5*(vx_Ly[idx+Nx]*vx_Ly[idx+Nx]+vy_Ly[idx+Nx]*vy_Ly[idx+Nx])+u_Ry);
      F_Ly = (U_Ly + P_Ry[idx])*vy_Ry[idx]; 
      F_Ry = (U_Ry + P_Ly[idx+Nx])*vy_Ly[idx+Nx];
      Q_Ly = U_Ly;
      Q_Ry = U_Ry;
      energy_flux_y[idx] = (0.5)*(F_Ly+F_Ry) - (vmax_y/2.0)*(Q_Ry - Q_Ly);
      }
    }
    periodic_boundary(mass_flux_x);
    periodic_boundary(mass_flux_y);
    periodic_boundary(momx_flux_x);
    periodic_boundary(momx_flux_y);
    periodic_boundary(momy_flux_x);
    periodic_boundary(momy_flux_y);
    periodic_boundary(energy_flux_x);
    periodic_boundary(energy_flux_y);
  }

  void update_conserved(double dt) {
    for(int i = 1; i < Nx-1; i++){
      for(int j = 1; j < Ny-1; j++){
      int idx= i + j*Nx;
      //updating mass, keep factors of dx and dy in because we didn't in flux calculations
      //std::cout<<mass[idx]<<std::endl;
      mass[idx]= mass[idx] - (mass_flux_x[idx]-mass_flux_x[idx-1])*dy*dt - (mass_flux_y[idx]-mass_flux_y[idx-Nx])*dx*dt;
      //std::cout<<mass[idx]/(dy*dx)<<std::endl;
      //updating momx
      mom_x[idx]= mom_x[idx]-(momx_flux_x[idx]-momx_flux_x[idx-1])*dt*dy - (momx_flux_y[idx]-momx_flux_y[idx-Nx])*dt*dx;
      //updating momy
      mom_y[idx]= mom_y[idx]-(momy_flux_x[idx]-momy_flux_x[idx-1])*dt*dy - (momy_flux_y[idx]-momy_flux_y[idx-Nx])*dt*dx;
      //updating energy
      energy[idx]=energy[idx]-(energy_flux_x[idx]-energy_flux_x[idx-1])*dt*dy - (energy_flux_y[idx]-energy_flux_y[idx-Nx])*dx*dt;
      }
    }
    periodic_boundary(mass);
    periodic_boundary(mom_x);
    periodic_boundary(mom_y);
    periodic_boundary(energy);
  }

  void output(int n) {
    std::ofstream outfile("output_rho_" + std::to_string(n) + ".csv");
    for (int j = 1; j < Ny - 1; j++) {
      for (int i = 1; i < Nx - 1; i++) {
        int idx = i + j * Nx;
        outfile << rho[idx];
        if (i != Nx - 2)
          outfile << ", ";
        else
          outfile << std::endl;
      }
    }
    outfile.close();
  }

  int Nx, Ny;
  double Lx, Ly;
  double dx, dy;
  double gamma, output_dt;
  std::vector<double> rho, vx, vy, P;                 // primitive variables
  std::vector<double> mass, mom_x, mom_y, energy;     // conserved variables
  // arrays to hold the results during primitive_update
  std::vector<double> rho_tmp, vx_tmp, vy_tmp, P_tmp;
  // arrays of fluxes for each conserved variable:
  std::vector<double> mass_flux_x, mass_flux_y;
  std::vector<double> momx_flux_x, momx_flux_y;
  std::vector<double> momy_flux_x, momy_flux_y;
  std::vector<double> energy_flux_x, energy_flux_y;
  // arrays for extrapolating to cell interfaces:
  std::vector<double> rho_Lx, rho_Ly, rho_Rx, rho_Ry;
  std::vector<double> vx_Lx, vx_Ly, vx_Rx, vx_Ry;
  std::vector<double> vy_Lx, vy_Ly, vy_Rx, vy_Ry;
  std::vector<double> P_Lx, P_Ly, P_Rx, P_Ry;
};
