#ifndef LBM_FLUCTUATIONS_H_
#define LBM_FLUCTUATIONS_H_

#ifdef AMREX_USE_CUDA
#include <cufft.h>
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
#endif

#include <lapacke.h>
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

GpuArray<Real,ndof*ndof> noise_covariance(const Real rho0, const Real phi0, const Real k2) { // TODO: add function arguments: rho, phi, k, others?
  BL_PROFILE_VAR("noise_covariance()",noise_covariance);
  const int Q = nvel;
  const Real kT = temperature;
  bulk_free_energy curr_state(rho0, phi0);

  // TODO: generalize to MRT
  const Real lambda_r = 1./tau_r;
  const Real lambda_p = 1./tau_p;

  const Real cs2k = curr_state.sound_speed_square() + kappa*k2*rho0;
  // const Real mu_rho = curr_state.calculate_dmur_drho() + k2*kappa;
  const Real mu_rho = curr_state.calculate_dmup_drho();
  const Real mu_phi = curr_state.calculate_dmup_dphi() + k2*kappa;
  const Real p_C = k2*kappa*phi0;
  
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

// compute correlated noise vector from Gaussian random variables
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<GpuComplex<Real>,ndof> correlated_noise(int kx, int ky, int kz,
						   const Box& domain,
						   const RandomEngine& engine,
               const Array4<Real>& hydrovs_ref) {
  BL_PROFILE_VAR("correlated_noise()",correlated_noise);
  GpuArray<GpuComplex<Real>,ndof> r, xi;
  GpuArray<Real,ndof*ndof> C;
  // GpuArray<GpuComplex<Real>,ndof*ndof> C;

  // // Cholesky decomposition of noise covariance matrix
  const Real rho0 = hydrovs_ref(kx, ky, kz, 0);
  const Real phi0 = hydrovs_ref(kx, ky, kz, 1);
  const Real k2 = fourier_laplace_operator(kx, ky, kz, domain);
  #if AMREX_DEBUG
    if (ky == 0 and kz == 0){Print() << "kx:" << kx << " rho:" << rho0 << " phi:" << phi0 << " k2:" << k2 << "\n"}
  #endif
    // if(ky == 0 and kz == 0){Print() << "kx:" << kx << " ky:" << ky << " kz:" << kz << " phi:" << phi0 << " convolution:" << phi_conv << "\n";
                          // Print() << "kx:" << kx << " ky:" << ky << " kz:" << kz << " rho:" << rho0 << " convolution:" << rho_conv << "\n";}
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
inline void generate_kspace_noise(const Geometry& geom,
				  MultiFab& kspace_noise_real,
				  MultiFab& kspace_noise_imag,
          MultiFab& ref_params) {
  BL_PROFILE_VAR("generate_kspace_noise()",generate_kspace_noise);
  
  const Box domain = geom.Domain(); 
  
  // generate noise in whole box because of grid decomposition
  // (generating noise without grid decomposition may be faster)
  for (MFIter mfi(kspace_noise_real); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& xi_real = kspace_noise_real.array(mfi);
    const Array4<Real>& xi_imag = kspace_noise_imag.array(mfi);
    const Array4<Real>& r = ref_params.array(mfi);

    // construct noise in k-space
    ParallelForRNG(box, [=] AMREX_GPU_DEVICE(int kx, int ky, int kz, RandomEngine const& engine) {
      if (kx <= domain.length(0)/2) { // need only half of k-space for c2r FFT
        GpuArray<GpuComplex<Real>,ndof> xi = {};
	      // compute correlated noise in k-space
	      xi = correlated_noise(kx,ky,kz,domain,engine,r);

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
  for (int i=ncons+1; i<ndof; ++i) {
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

// generate real-space noise for all non-conserved moments
inline void generate_realspace_noise(const Geometry& geom, MultiFab& hydrovs, MultiFab& noise) {
  for (MFIter mfi(noise); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& xi = noise.array(mfi);
    const Array4<Real>& h = hydrovs.array(mfi);
    ParallelForRNG(box, [=] AMREX_GPU_DEVICE(int x, int y, int z, RandomEngine const& engine) {
      const Real rho0 = h(x,y,z,0);
      const Real phi0 = h(x,y,z,1);
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
  MultiFab kspace_noise_real(ba, dm, ndof, 0);
  MultiFab kspace_noise_imag(ba, dm, ndof, 0);

  kspace_noise_real.setVal(0.);
  kspace_noise_imag.setVal(0.);

  // TODO: check if this takes care of ndof entries!?

  if (use_correlated_noise){
    // generate noise in k-space
    Print() << "Correlated noise\n";
    generate_kspace_noise(geom, kspace_noise_real, kspace_noise_imag, ref_params);
    // note that the k-space noise is generated without the required symmetries
    // (this is to allow for parallel generation of noise)
    // the k-space symmetries are handled when copying in compute_ifft

    // inverse Fourier transform noise vector to real space
    compute_ifft(geom, noise, kspace_noise_real, kspace_noise_imag, ndof);
  }
  else{
    Print() << "Uncorrelated noise\n";
    generate_realspace_noise(geom, ref_params, noise);
  }

}
#endif