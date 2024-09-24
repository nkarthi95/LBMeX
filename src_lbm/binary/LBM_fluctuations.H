#ifndef LBM_FLUCTUATIONS_H_
#define LBM_FLUCTUATIONS_H_

#ifndef AMREX_USE_CUDA
  #include <lapacke.h>
#endif

#include <AMReX_GpuComplex.H>
#include <math.h>
#include <vector>
#include "LBM_d3q19.H"
#include "LBM_binary.H"
#include "LBM_FFT.H"

const int ncons = 2 + AMREX_SPACEDIM;
const int ndof = 2*nvel;

// Cholesky decomposition of matrix A
// result is stored in lower triangle of A
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cholesky_decomp(GpuArray<Real,ndof*ndof>& A, const int n, const int bstart) {
  // Cholesky-Banachiewicz algorithm
  Real sum;
  for (int i=bstart; i<n; ++i) {
    for (int j=bstart; j<=i; ++j) {
      sum = A[i*n+j];
      for (int k=j-1; k>=bstart; --k) {
	      sum -= A[i*n+k]*A[j*n+k];
      }
      if (i==j) {
	      if (sum>=0) {
	        A[i*n+j] = std::sqrt(sum);
	      } else {
	        A[i*n+j] = 0.0;
          Print() << "Row " << i << " matrix not positive definite! " << sum << std::endl;
          exit(-1);
	      }
      } else {
	      if (A[j*n+j]>0) {
	        A[i*n+j] = sum/A[j*n+j];
	      } else {
          Print() << "Cholesky decomposition should not reach " << __FILE__ <<":"<< __LINE__ << std::endl;
	        exit(-1);
	      }
      }
    }
  }
  for (int i=0; i<n; ++i) {
    for (int j=i+1; j<n; ++j) {
      A[i*n+j] = 0.0;
    }
  }
}

#if 0
GpuArray<Real,ndof*ndof> noise_covariance(Real rho0, Real phi0, Real k2) {
  BL_PROFILE_VAR("noise_covariance()",noise_covariance);
  const int Q = nvel;
  const Real kT = temperature;
  // const Real k2 = fourier_laplace_operator(kx, ky, kz, domain);

  // TODO: generalize to MRT
  const Real lambda_r = 1./tau_r;
  const Real lambda_p = 1./tau_p;

  const Real cs2k = T + kappa*k2*rho0;// + kappa*k2*phi0; //checked after conversion
  const Real mu_rho = (-2.*T*pow(rho0, 4) - chi*pow(phi0, 4) + chi*pow(phi0*rho0, 2))/(2.*pow(rho0, 3)*(pow(phi0, 2) - pow(rho0, 2))) + k2*kappa;
  const Real mu_phi = (-2.*T*pow(rho0, 2) - chi*pow(phi0, 2) + chi*pow(rho0, 2))/(2.*rho0*(pow(phi0, 2) - pow(rho0, 2))) + k2*kappa;
  const Real p_C = k2*kappa*phi0; //checked after conversion
  
  GpuArray<Real,ndof*ndof> Xi = {};
  Xi.fill(0.);

  // diagonal part sector
  Xi[195] = 2.*Gamma*kT*lambda_p/rho0;
  Xi[234] = 2.*Gamma*kT*lambda_p/rho0;
  Xi[273] = 2.*Gamma*kT*lambda_p/rho0;
  Xi[312] = 2.*kT*lambda_r*rho0*(5. - 9.*cs2k);
  Xi[351] = 8.*kT*lambda_r*rho0;
  Xi[390] = (8.0/3.0)*kT*lambda_r*rho0;
  Xi[429] = (2.0/3.0)*kT*lambda_r*rho0;
  Xi[468] = (2.0/3.0)*kT*lambda_r*rho0;
  Xi[507] = (2.0/3.0)*kT*lambda_r*rho0;
  Xi[546] = 4.*kT*lambda_r*rho0;
  Xi[585] = 4.*kT*lambda_r*rho0;
  Xi[624] = 4.*kT*lambda_r*rho0;
  Xi[663] = (4.0/3.0)*kT*lambda_r*rho0;
  Xi[702] = (4.0/3.0)*kT*lambda_r*rho0;
  Xi[741] = (4.0/3.0)*kT*lambda_r*rho0;
  Xi[780] = 18.*kT*lambda_r*rho0*(1. - cs2k);
  Xi[819] = 8.*kT*lambda_r*rho0;
  Xi[858] = (8.0/3.0)*kT*lambda_r*rho0;
  Xi[897] = 2.*Gamma*kT*lambda_p*(-9.*Gamma*mu_phi + 5.)/rho0;
  Xi[936] = 8.*Gamma*kT*lambda_p/rho0;
  Xi[975] = (8.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1014] = (2.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1053] = (2.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1092] = (2.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1131] = 4.*Gamma*kT*lambda_p/rho0;
  Xi[1170] = 4.*Gamma*kT*lambda_p/rho0;
  Xi[1209] = 4.*Gamma*kT*lambda_p/rho0;
  Xi[1248] = (4.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1287] = (4.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1326] = (4.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1365] = 18.*Gamma*kT*lambda_p*(-Gamma*mu_phi + 1.)/rho0;
  Xi[1404] = 8.*Gamma*kT*lambda_p/rho0;
  Xi[1443] = (8.0/3.0)*Gamma*kT*lambda_p/rho0;

  // rho-phi or phi-rho sector
  Xi[23] = -3.*Gamma*kT*lambda_p*mu_rho*rho0/cs2k;
  Xi[35] = 3.*Gamma*kT*lambda_p*mu_rho*rho0/cs2k;
  Xi[58] = 3.*kT*lambda_r*(phi0*k2*kappa + p_C)/(mu_phi*rho0);
  Xi[324] = 6.*kT*lambda_r*rho0*(3*cs2k - 1);
  Xi[327] = -3.*kT*(Gamma*mu_phi*lambda_p*mu_rho*pow(rho0, 2)*(3.*cs2k - 1.) + cs2k*lambda_r*(3.*Gamma*mu_phi - 1.)*(phi0*k2*kappa + p_C))/(mu_phi*cs2k*rho0);
  Xi[339] = 3.*kT*(Gamma*mu_phi*lambda_p*mu_rho*pow(rho0, 2)*(3.*cs2k - 1.) + cs2k*lambda_r*(3.*Gamma*mu_phi - 1.)*(phi0*k2*kappa + p_C))/(mu_phi*cs2k*rho0);
  Xi[761] = 3.*kT*lambda_r*(phi0*k2*kappa + p_C)/(mu_phi*rho0);
  Xi[768] = 6.*kT*lambda_r*rho0*(3.*cs2k - 1.);
  Xi[874] = -3.*Gamma*kT*lambda_p*mu_rho*rho0/cs2k;
  Xi[882] = -3.*kT*(Gamma*mu_phi*lambda_p*mu_rho*pow(rho0, 2)*(3.*cs2k - 1.) + cs2k*lambda_r*(3*Gamma*mu_phi - 1.)*(phi0*k2*kappa + p_C))/(mu_phi*cs2k*rho0);
  Xi[1330] = 3.*Gamma*kT*lambda_p*mu_rho*rho0/cs2k;
  Xi[1338] = 3.*kT*(Gamma*mu_phi*lambda_p*mu_rho*pow(rho0, 2)*(3.*cs2k - 1.) + cs2k*lambda_r*(3*Gamma*mu_phi - 1.)*(phi0*k2*kappa + p_C))/(mu_phi*cs2k*rho0);

  // rho-rho and phi-phi off diagonal
  Xi[46] = -3.*kT*lambda_r*(phi0*k2*kappa + p_C)/(mu_phi*rho0);
  Xi[81] = -phi0*kT*lambda_p;
  Xi[120] = -phi0*kT*lambda_p;
  Xi[159] = -phi0*kT*lambda_p;
  Xi[192] = -phi0*kT*lambda_p;
  Xi[231] = -phi0*kT*lambda_p;
  Xi[270] = -phi0*kT*lambda_p;
  Xi[305] = -3.*kT*lambda_r*(phi0*k2*kappa + p_C)/(mu_phi*rho0);
  Xi[783] = 3.*kT*(Gamma*mu_phi*lambda_p*mu_rho*pow(rho0, 2)*(3*cs2k - 1.) + cs2k*lambda_r*(3.*Gamma*mu_phi - 1.)*(phi0*k2*kappa + p_C))/(mu_phi*cs2k*rho0);
  Xi[795] = -3.*kT*(Gamma*mu_phi*lambda_p*mu_rho*pow(rho0, 2)*(3*cs2k - 1.) + cs2k*lambda_r*(3.*Gamma*mu_phi - 1.)*(phi0*k2*kappa + p_C))/(mu_phi*cs2k*rho0);
  Xi[894] = 3.*kT*(Gamma*mu_phi*lambda_p*mu_rho*pow(rho0, 2)*(3*cs2k - 1.) + cs2k*lambda_r*(3.*Gamma*mu_phi - 1.)*(phi0*k2*kappa + p_C))/(mu_phi*cs2k*rho0);
  Xi[909] = 6.*Gamma*kT*lambda_p*(3.*Gamma*mu_phi - 1.)/rho0;
  Xi[1350] = -3.*kT*(Gamma*mu_phi*lambda_p*mu_rho*pow(rho0, 2)*(3*cs2k - 1.) + cs2k*lambda_r*(3.*Gamma*mu_phi - 1.)*(phi0*k2*kappa + p_C))/(mu_phi*cs2k*rho0);
  Xi[1353] = 6.*Gamma*kT*lambda_p*(3.*Gamma*mu_phi - 1.)/rho0;

  return Xi;
}
#endif

#if 0
GpuArray<Real,ndof*ndof> noise_covariance(Real rho0, Real phi0, Real k2) {
  BL_PROFILE_VAR("noise_covariance()",noise_covariance);
  const int Q = nvel;
  const Real kT = temperature;
  // const Real k2 = fourier_laplace_operator(kx, ky, kz, domain);

  // TODO: generalize to MRT
  const Real lambdaLB_r = -1./tau_r;
  const Real lambdaLB_p = -1./tau_p;
  const Real lambda_r  = -lambdaLB_r*(2.+lambdaLB_r)/2.;
  const Real lambda_p  = -lambdaLB_p*(2.+lambdaLB_p)/2.;
  const Real lambda_rp = -lambdaLB_r*(2.+lambdaLB_p)/2.;
  const Real lambda_pr = -lambdaLB_p*(2.+lambdaLB_r)/2.;

  const Real cs2k = T + kappa*k2*rho0; //\frac{d pr}{d \rho } Eq 40 of swift et al
  const Real p_phi = kappa*k2*phi0;//\frac{d pr}{d \phi } Eq 40 of swift et al
  const Real mu_rho = -T*phi0/(rho0*rho0-phi0*phi0) + chi/2*phi0/(rho0*rho0); //\frac{d^2 f}{d \phi d\rho} Eq 37 of Swift et al
  const Real mu_phi = T*rho0/(rho0*rho0-phi0*phi0) - chi/2/rho0 + kappa*k2; //\frac{d^2 f}{d \phi d\phi} Eq 37 of Swift et al

  // const Real cs2k = T + kappa*k2*rho0;// + kappa*k2*phi0; //checked after conversion
  // const Real mu_rho = (-2.*T*pow(rho0, 4) - chi*pow(phi0, 4) + chi*pow(phi0*rho0, 2))/(2.*pow(rho0, 3)*(pow(phi0, 2) - pow(rho0, 2))) + k2*kappa;
  // const Real mu_phi = (-2.*T*pow(rho0, 2) - chi*pow(phi0, 2) + chi*pow(rho0, 2))/(2.*rho0*(pow(phi0, 2) - pow(rho0, 2))) + k2*kappa;
  // const Real p_phi = k2*kappa*phi0; //checked after conversion
  
  GpuArray<Real,ndof*ndof> Xi = {};
  Xi.fill(0.);

  // diagonal part sector
  Xi[195] = 2*Gamma*kT*lambda_p/rho0;
  Xi[234] = 2*Gamma*kT*lambda_p/rho0;
  Xi[273] = 2*Gamma*kT*lambda_p/rho0;
  Xi[312] = 2*kT*lambda_r*rho0*(5 - 9*cs2k);
  Xi[351] = 8*kT*lambda_r*rho0;
  Xi[390] = (8.0/3.0)*kT*lambda_r*rho0;
  Xi[429] = (2.0/3.0)*kT*lambda_r*rho0;
  Xi[468] = (2.0/3.0)*kT*lambda_r*rho0;
  Xi[507] = (2.0/3.0)*kT*lambda_r*rho0;
  Xi[546] = 4*kT*lambda_r*rho0;
  Xi[585] = 4*kT*lambda_r*rho0;
  Xi[624] = 4*kT*lambda_r*rho0;
  Xi[663] = (4.0/3.0)*kT*lambda_r*rho0;
  Xi[702] = (4.0/3.0)*kT*lambda_r*rho0;
  Xi[741] = (4.0/3.0)*kT*lambda_r*rho0;
  Xi[780] = 18*kT*lambda_r*rho0*(1 - cs2k);
  Xi[819] = 8*kT*lambda_r*rho0;
  Xi[858] = (8.0/3.0)*kT*lambda_r*rho0;
  Xi[897] = 2*Gamma*kT*lambda_p*(-9*Gamma*mu_phi + 5)/rho0;
  Xi[936] = 8*Gamma*kT*lambda_p/rho0;
  Xi[975] = (8.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1014] = (2.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1053] = (2.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1092] = (2.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1131] = 4*Gamma*kT*lambda_p/rho0;
  Xi[1170] = 4*Gamma*kT*lambda_p/rho0;
  Xi[1209] = 4*Gamma*kT*lambda_p/rho0;
  Xi[1248] = (4.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1287] = (4.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1326] = (4.0/3.0)*Gamma*kT*lambda_p/rho0;
  Xi[1365] = 18*Gamma*kT*lambda_p*(-Gamma*mu_phi + 1)/rho0;
  Xi[1404] = 8*Gamma*kT*lambda_p/rho0;
  Xi[1443] = (8.0/3.0)*Gamma*kT*lambda_p/rho0;

  // rho-phi or phi-rho sector
  Xi[23] = -3*Gamma*kT*lambda_pr*mu_rho*rho0/cs2k;
  Xi[35] = 3*Gamma*kT*lambda_pr*mu_rho*rho0/cs2k;
  Xi[58] = 3*kT*lambda_rp*p_phi/(mu_phi*rho0);
  Xi[324] = 6*kT*lambda_rp*rho0*(3*cs2k - 1);
  Xi[327] = -9*Gamma*kT*lambda_pr*mu_rho*rho0 - 9*Gamma*kT*lambda_rp*p_phi/rho0 + 3*Gamma*kT*lambda_pr*mu_rho*rho0/cs2k + 3*kT*lambda_rp*p_phi/(mu_phi*rho0);
  Xi[339] = 3*kT*(Gamma*lambda_pr*mu_phi*mu_rho*pow(rho0, 2)*(3*cs2k - 1) + cs2k*lambda_rp*p_phi*(3*Gamma*mu_phi - 1))/(cs2k*mu_phi*rho0);
  Xi[761] = 3*kT*lambda_rp*p_phi/(mu_phi*rho0);
  Xi[768] = 6*kT*lambda_rp*rho0*(3*cs2k - 1);
  Xi[874] = -3*Gamma*kT*lambda_pr*mu_rho*rho0/cs2k;
  Xi[882] = -9*Gamma*kT*lambda_pr*mu_rho*rho0 - 9*Gamma*kT*lambda_rp*p_phi/rho0 + 3*Gamma*kT*lambda_pr*mu_rho*rho0/cs2k + 3*kT*lambda_rp*p_phi/(mu_phi*rho0);
  Xi[1330] = 3*Gamma*kT*lambda_pr*mu_rho*rho0/cs2k;
  Xi[1338] = 3*kT*(Gamma*lambda_pr*mu_phi*mu_rho*pow(rho0, 2)*(3*cs2k - 1) + cs2k*lambda_rp*p_phi*(3*Gamma*mu_phi - 1))/(cs2k*mu_phi*rho0);

  // rho-rho and phi-phi off diagonal
  Xi[46] = -3*kT*lambda_rp*p_phi/(mu_phi*rho0);
  Xi[81] = -kT*lambda_pr*phi0;
  Xi[120] = -kT*lambda_pr*phi0;
  Xi[159] = -kT*lambda_pr*phi0;
  Xi[192] = -kT*lambda_pr*phi0;
  Xi[231] = -kT*lambda_pr*phi0;
  Xi[270] = -kT*lambda_pr*phi0;
  Xi[305] = -3*kT*lambda_rp*p_phi/(mu_phi*rho0);
  Xi[783] = 3*kT*(Gamma*lambda_pr*mu_phi*mu_rho*pow(rho0, 2)*(3*cs2k - 1) + cs2k*lambda_rp*p_phi*(3*Gamma*mu_phi - 1))/(cs2k*mu_phi*rho0);
  Xi[795] = -9*Gamma*kT*lambda_pr*mu_rho*rho0 - 9*Gamma*kT*lambda_rp*p_phi/rho0 + 3*Gamma*kT*lambda_pr*mu_rho*rho0/cs2k + 3*kT*lambda_rp*p_phi/(mu_phi*rho0);
  Xi[894] = 3*kT*(Gamma*lambda_pr*mu_phi*mu_rho*pow(rho0, 2)*(3*cs2k - 1) + cs2k*lambda_rp*p_phi*(3*Gamma*mu_phi - 1))/(cs2k*mu_phi*rho0);
  Xi[909] = 6*Gamma*kT*lambda_pr*(3*Gamma*mu_phi - 1)/rho0;
  Xi[1350] = -9*Gamma*kT*lambda_pr*mu_rho*rho0 - 9*Gamma*kT*lambda_rp*p_phi/rho0 + 3*Gamma*kT*lambda_pr*mu_rho*rho0/cs2k + 3*kT*lambda_rp*p_phi/(mu_phi*rho0);
  Xi[1353] = 6*Gamma*kT*lambda_pr*(3*Gamma*mu_phi - 1)/rho0;

  return Xi;
}
#endif

// noise covariance matrix of the LB modes
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<Real,ndof*ndof> noise_covariance(const Real rho0, const Real phi0, const Real k2) {
  const uint Q = nvel;
  const Real kT = temperature;

  const Real cs2 = T + kappa*k2*rho0; // p_rho => modified speed of sound
  const Real p_phi = kappa*k2*phi0;
  const Real mu_rho = -T*phi0/(rho0*rho0-phi0*phi0) + chi/2*phi0/(rho0*rho0);
  const Real mu_phi = T*rho0/(rho0*rho0-phi0*phi0) - chi/2/rho0 + kappa*k2;

  // Noise covariance of the modes in k-space
  GpuArray<Real,ndof*ndof> Xi = {};
  Xi.fill(0);

  // TODO: generalize to MRT
  const Real lambdaLB_r = -1./tau_r;
  const Real lambdaLB_p = -1./tau_p;

  // TODO: discretized needs check!
  Real lambda_r = -lambdaLB_r*(2+lambdaLB_r)/2;
  Real lambda_p = -lambdaLB_p*(2+lambdaLB_p)/2;
  Real lambda_rp = -lambdaLB_r*(2+lambdaLB_p)/2;
  Real lambda_pr = -lambdaLB_p*(2+lambdaLB_r)/2;

  // diagonal part
  Xi[(   5)*ndof+(   5)] = 2.*Gamma*kT/rho0*lambda_p;
  Xi[(   6)*ndof+(   6)] = 2.*Gamma*kT/rho0*lambda_p;
  Xi[(   7)*ndof+(   7)] = 2.*Gamma*kT/rho0*lambda_p;
  Xi[(   8)*ndof+(   8)] = 2.*kT*rho0*(5 - 9*cs2)*lambda_r;
  Xi[(   9)*ndof+(   9)] = 8.*kT*rho0*lambda_r;
  Xi[(  10)*ndof+(  10)] = (8.0/3.0)*kT*rho0*lambda_r;
  Xi[(  11)*ndof+(  11)] = (2.0/3.0)*kT*rho0*lambda_r;
  Xi[(  12)*ndof+(  12)] = (2.0/3.0)*kT*rho0*lambda_r;
  Xi[(  13)*ndof+(  13)] = (2.0/3.0)*kT*rho0*lambda_r;
  Xi[(  14)*ndof+(  14)] = 4.*kT*rho0*lambda_r;
  Xi[(  15)*ndof+(  15)] = 4.*kT*rho0*lambda_r;
  Xi[(  16)*ndof+(  16)] = 4.*kT*rho0*lambda_r;
  Xi[(  17)*ndof+(  17)] = (4.0/3.0)*kT*rho0*lambda_r;
  Xi[(  18)*ndof+(  18)] = (4.0/3.0)*kT*rho0*lambda_r;
  Xi[(Q+ 0)*ndof+(Q+ 0)] = (4.0/3.0)*kT*rho0*lambda_r;
  Xi[(Q+ 1)*ndof+(Q+ 1)] = 18.*kT*rho0*(1 - cs2)*lambda_r;
  Xi[(Q+ 2)*ndof+(Q+ 2)] = 8.*kT*rho0*lambda_r;
  Xi[(Q+ 3)*ndof+(Q+ 3)] = (8.0/3.0)*kT*rho0*lambda_r;
  Xi[(Q+ 4)*ndof+(Q+ 4)] = 2.*Gamma*kT/rho0*(-9*Gamma*mu_phi + 5)*lambda_p;
  Xi[(Q+ 5)*ndof+(Q+ 5)] = 8.*Gamma*kT/rho0*lambda_p;
  Xi[(Q+ 6)*ndof+(Q+ 6)] = (8.0/3.0)*Gamma*kT/rho0*lambda_p;
  Xi[(Q+ 7)*ndof+(Q+ 7)] = (2.0/3.0)*Gamma*kT/rho0*lambda_p;
  Xi[(Q+ 8)*ndof+(Q+ 8)] = (2.0/3.0)*Gamma*kT/rho0*lambda_p;
  Xi[(Q+ 9)*ndof+(Q+ 9)] = (2.0/3.0)*Gamma*kT/rho0*lambda_p;
  Xi[(Q+10)*ndof+(Q+10)] = 4.*Gamma*kT/rho0*lambda_p;
  Xi[(Q+11)*ndof+(Q+11)] = 4.*Gamma*kT/rho0*lambda_p;
  Xi[(Q+12)*ndof+(Q+12)] = 4.*Gamma*kT/rho0*lambda_p;
  Xi[(Q+13)*ndof+(Q+13)] = (4.0/3.0)*Gamma*kT/rho0*lambda_p;
  Xi[(Q+14)*ndof+(Q+14)] = (4.0/3.0)*Gamma*kT/rho0*lambda_p;
  Xi[(Q+15)*ndof+(Q+15)] = (4.0/3.0)*Gamma*kT/rho0*lambda_p;
  Xi[(Q+16)*ndof+(Q+16)] = 18.*Gamma*kT/rho0*(-Gamma*mu_phi + 1)*lambda_p;
  Xi[(Q+17)*ndof+(Q+17)] = 8.*Gamma*kT/rho0*lambda_p;
  Xi[(Q+18)*ndof+(Q+18)] = (8.0/3.0)*Gamma*kT/rho0*lambda_p;

  // rho-rho sector [0, 2..4, 8..(Q+3)]
  Xi[(   8)*ndof+(Q+ 1)] = 6.*kT*rho0*(3*cs2 - 1)*lambda_r;
  Xi[(Q+ 1)*ndof+(   8)] = 6.*kT*rho0*(3*cs2 - 1)*lambda_r;

  // phi-phi sector [1, 5..7, (Q+4)..(Q-1)]
  Xi[(Q+ 4)*ndof+(Q+16)] = 6.*Gamma*kT/rho0*(3*Gamma*mu_phi - 1)*lambda_p;
  Xi[(Q+16)*ndof+(Q+ 4)] = 6.*Gamma*kT/rho0*(3*Gamma*mu_phi - 1)*lambda_p;

  // rho-phi sector
  Xi[(   8)*ndof+(Q+ 4)] = -3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0);
  Xi[(Q+ 4)*ndof+(   8)] = -3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0);
  Xi[(   8)*ndof+(Q+16)] = 3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0);
  Xi[(Q+16)*ndof+(   8)] = 3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0);

  Xi[(Q+ 1)*ndof+(Q+ 4)] = 3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0);
  Xi[(Q+ 4)*ndof+(Q+ 1)] = 3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0);
  Xi[(Q+ 1)*ndof+(Q+16)] = -3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0);
  Xi[(Q+16)*ndof+(Q+ 1)] = -3.*kT*(Gamma*mu_phi*pow(rho0, 2)*(3*cs2 - 1)*mu_rho*lambda_pr + cs2*(3*Gamma*mu_phi - 1)*p_phi*lambda_rp)/(cs2*mu_phi*rho0);

  Xi[(   0)*ndof+(Q+ 4)] = -3.*Gamma*kT*rho0*mu_rho/cs2*lambda_pr;
  Xi[(Q+ 4)*ndof+(   0)] = -3.*Gamma*kT*rho0*mu_rho/cs2*lambda_pr;
  Xi[(   0)*ndof+(Q+16)] = 3.*Gamma*kT*rho0*mu_rho/cs2*lambda_pr;
  Xi[(Q+16)*ndof+(   0)] = 3.*Gamma*kT*rho0*mu_rho/cs2*lambda_pr;
  Xi[(   1)*ndof+(Q+ 1)] = 3.*kT*p_phi/(mu_phi*rho0)*lambda_rp;
  Xi[(Q+ 1)*ndof+(   1)] = 3.*kT*p_phi/(mu_phi*rho0)*lambda_rp;
  Xi[(   1)*ndof+(   8)] = -3.*kT*p_phi/(mu_phi*rho0)*lambda_rp;
  Xi[(   8)*ndof+(   1)] = -3.*kT*p_phi/(mu_phi*rho0)*lambda_rp;
  Xi[(   2)*ndof+(   5)] = -phi0*kT*lambda_pr;
  Xi[(   3)*ndof+(   6)] = -phi0*kT*lambda_pr;
  Xi[(   4)*ndof+(   7)] = -phi0*kT*lambda_pr;
  Xi[(   5)*ndof+(   2)] = -phi0*kT*lambda_pr;
  Xi[(   6)*ndof+(   3)] = -phi0*kT*lambda_pr;
  Xi[(   7)*ndof+(   4)] = -phi0*kT*lambda_pr;

  return Xi;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real fourier_laplace_operator(int ikx, int iky, int ikz, const Box& domain) {
  BL_PROFILE_VAR("fourier_laplace_operator()",fourier_laplace_operator);
  Real k2;

  IntVect n = domain.length();

  // FFTW convention for ordering of wave vectors
  Real kx = (ikx < (n[0]+1)/2) ? 2.*M_PI/n[0]*ikx : 2.*M_PI/n[0]*(ikx-n[0]);
  Real ky = (iky < (n[1]+1)/2) ? 2.*M_PI/n[1]*iky : 2.*M_PI/n[1]*(iky-n[1]);
  Real kz = (ikz < (n[2]+1)/2) ? 2.*M_PI/n[2]*ikz : 2.*M_PI/n[2]*(ikz-n[2]);

  Real cosx = cos(kx);
  Real cosy = cos(ky);
  Real cosz = cos(kz);

  Real expr1 = cosx + cosy + cosz;
  Real expr2 = cosx*cosy + cosy*cosz + cosx*cosz;
  k2 = -2./cs2*(1./9.*expr1 + 1./9.*expr2 - 2./3.);

  return k2;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<GpuComplex<Real>,ndof> kspace_white_noise(int kx, int ky, int kz, const Box& domain, RandomEngine const& engine) {
  BL_PROFILE_VAR("kspace_white_noise()",kspace_white_noise);
  GpuArray<GpuComplex<Real>,ndof> r = {};
  for (int i=ncons; i<ndof; ++i) {
    // symmetry points are purely real
    if (     ((kx == 0) || (kx == domain.length(0) - kx))
          && ((ky == 0) || (ky == domain.length(1) - ky))
          && ((kz == 0) || (kz == domain.length(2) - kz)) ) {
      // real Gaussian random variables with zero mean and variance 1
      r[i] = { RandomNormal(0., 1., engine), RandomNormal(0., 0., engine) };
    } else {
      // complex Gaussian random variables with zero mean and variance 0.5
      r[i] = { RandomNormal(0., std::sqrt(0.5), engine), RandomNormal(0., std::sqrt(0.5), engine) };
      // r[i] = { RandomNormal(0., 1., engine), RandomNormal(0., 0., engine) };
    }
  }
  return r;
}

// compute correlated noise vector from Gaussian random variables
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<GpuComplex<Real>,ndof> correlated_noise_vars(int kx, int ky, int kz, const Box& domain, const RandomEngine& engine, const Array4<Real>& ref) {
  
  BL_PROFILE_VAR("correlated_noise()",correlated_noise);
  
  const Real k2 = fourier_laplace_operator(kx, ky, kz, domain);
  const Real rho0 = ref(kx, ky, kz, 0);
  const Real phi0 = ref(kx, ky, kz, 1);
  // if(ky == 0 and kz == 0){Print() << "rho:" << rho0 << " phi0:" << phi0 << "\n";}
  GpuArray<GpuComplex<Real>,ndof> r, xi;
  GpuArray<Real,ndof*ndof> C;

  C = noise_covariance(rho0, phi0, k2);

  #ifdef AMREX_USE_CUDA
    cholesky_decomp(C,ndof,ncons);
  #else
    // // LAPACK CHOLESKY DECOMP. COMMENT OUT IF THIS DOESN'T WORK
    const char uplo = 'L';
    const int lda = ndof-ncons;
    const int order = lda;
    int info = 0;
    GpuArray<Real,lda*lda> SLC_C; SLC_C.fill(0.);
    for (int i = ncons; i < ndof; i++){
      for (int j = ncons; j < ndof; j++){
        SLC_C[(i-ncons)*lda+(j-ncons)] = C[i*ndof+j];
      }
    }
    size_t arraySize = SLC_C.size();
    dpotrf_(&uplo, &order, SLC_C.data(), &lda, &info, arraySize);
    C.fill(0.);
    for (int i = ncons; i < ndof; i++){
      for (int j = ncons; j < i+1; j++){
        C[i*ndof+j] = SLC_C[(j-ncons)*lda+(i-ncons)];
      }
    }
    // // LAPACK CHOLESKY DECOMP. COMMENT OUT IF THIS DOESN'T WORK
  #endif

  // need to generate the correct symmetries here? [uschill 07/25/2022]
  // possibly the case if C(k) != C(-k) [uschill 07/27/2022]
  r = kspace_white_noise(kx,ky,kz,domain,engine);

  // compute correlated noise vector from Gaussian random variables
  for (int i=0; i<ndof; ++i) {
    xi[i] = { 0, 0 };
    for (int j=0; j<=i; ++j) {
      xi[i] += C[i*ndof+j]*r[j];
    }
  }
  return xi;
}

// generate k-space noise for all non-conserved moments
// the required symmetries are not included here because of grid decomposition
// (this is to allow parallel generation of noise)
// the symmetries are handled when copying to one whole grid
inline void correlated_noise(const Geometry& geom,
				  MultiFab& kspace_noise_real,
				  MultiFab& kspace_noise_imag,
          MultiFab& ref_params) {
  BL_PROFILE_VAR("generate_kspace_noise()",generate_kspace_noise);
  
  const Box domain = geom.Domain(); BoxArray ba_onegrid(domain); DistributionMapping dm_onegrid(ba_onegrid); BoxArray ba(domain);

  // include density? [uschill 07/26/2022]
  
  // generate noise in whole box because of grid decomposition
  // (generating noise without grid decomposition may be faster)
  for (MFIter mfi(kspace_noise_real); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& xi_real = kspace_noise_real.array(mfi);
    const Array4<Real>& xi_imag = kspace_noise_imag.array(mfi);
    const Array4<Real>& h = ref_params.array(mfi);

    // construct noise in k-space
    ParallelForRNG(box, [=] AMREX_GPU_DEVICE(int kx, int ky, int kz, RandomEngine const& engine) {
      if (kx <= domain.length(0)/2) { // need only half of k-space for c2r FFT
        GpuArray<GpuComplex<Real>,ndof> xi = {};
	      // compute correlated noise in k-space
	      xi = correlated_noise_vars(kx,ky,kz,domain,engine,h);

        for (int i=0; i<ndof; ++i) {
          xi_real(kx,ky,kz,i) = xi[i].real();
          xi_imag(kx,ky,kz,i) = xi[i].imag();
        }
      }
    });
  }
}

// compute spatially uncorrelated noise vector from Gaussian random variables
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<Real,ndof> uncorrelated_noise(const Real rho0, const Real phi0, const RandomEngine& engine) {
  GpuArray<Real,ndof> r, xi;
  GpuArray<Real,ndof*ndof> C;

  // Cholesky decomposition of noise covariance matrix
  C = noise_covariance(rho0,phi0,0);
  #ifdef AMREX_USE_CUDA
    cholesky_decomp(C,ndof,ncons);
  #else
    // // LAPACK CHOLESKY DECOMP. COMMENT OUT IF THIS DOESN'T WORK
    const char uplo = 'L';
    const int lda = ndof-ncons;
    const int order = lda;
    int info = 0;
    GpuArray<Real,lda*lda> SLC_C; SLC_C.fill(0.);
    for (int i = ncons; i < ndof; i++){
      for (int j = ncons; j < ndof; j++){
        SLC_C[(i-ncons)*lda+(j-ncons)] = C[i*ndof+j];
      }
    }
    size_t arraySize = SLC_C.size();
    dpotrf_(&uplo, &order, SLC_C.data(), &lda, &info, arraySize);
    C.fill(0.);
    for (int i = ncons; i < ndof; i++){
      for (int j = ncons; j < i+1; j++){
        C[i*ndof+j] = SLC_C[(j-ncons)*lda+(i-ncons)];
      }
    }
    // // LAPACK CHOLESKY DECOMP. COMMENT OUT IF THIS DOESN'T WORK
  #endif

  // random white noise
  for (int i=ncons; i<ndof; ++i) {
    r[i] = RandomNormal(0., 1., engine);
  }

  // compute noise vector from Gaussian random variables
  for (int i=0; i<ndof; ++i) {
    xi[i] = 0;
    for (int j=0; j<=i; ++j) {
      xi[i] += C[i*ndof+j]*r[j];
    }
  }

  return xi;
}

inline void uncorrelated_noise(const Geometry& geom, MultiFab& hydrovs, MultiFab& noise) {
  for (MFIter mfi(noise); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& xi = noise.array(mfi);
    const Array4<Real>& h = hydrovs.array(mfi);
    ParallelForRNG(box, [=] AMREX_GPU_DEVICE(int x, int y, int z, RandomEngine const& engine) {
      const Real rho0 = h(x,y,z,0);
      const Real phi0 = h(x,y,z,1);

      // if(y == 0 and z == 0){Print() << "rho:" << rho0 << " phi0:" << phi0 << "\n";}
      GpuArray<Real,ndof> r = uncorrelated_noise(rho0,phi0,engine);
      for (int i=0; i<ndof; ++i) {
        xi(x,y,z,i) = r[i];
      }
    });
  }
}

// LB thermalization procedure for spatially correlated, non-diagonal noise
inline void generate_fluctuations(const Geometry& geom,
				  MultiFab& hydrovs,
				  MultiFab& noise,
          MultiFab& ref_params) {
  BL_PROFILE_VAR("generate_fluctuations()",generate_fluctuations);
  
  BoxArray ba = noise.boxArray();
  DistributionMapping dm = noise.DistributionMap();

  // generate noise in k-space
  if (use_correlated_noise){
  // Print() << "Correlated noise\n";
  MultiFab kspace_noise_real(ba, dm, ndof, 0);
  MultiFab kspace_noise_imag(ba, dm, ndof, 0);

  kspace_noise_real.setVal(0.);
  kspace_noise_imag.setVal(0.);

  correlated_noise(geom, kspace_noise_real, kspace_noise_imag, ref_params);
  compute_ifft(geom, noise, kspace_noise_real, kspace_noise_imag, ndof);
  }
  else{
    // Print() << "Uncorrelated noise\n";
    uncorrelated_noise(geom, ref_params, noise);
  }
}

#endif