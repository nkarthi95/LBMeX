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
#include <AMReX_Array.H>
#include "LBM_d3q19.H"
#include "LBM_FFT.H"
#include "LBM_thermodynamics.H"

const int ncons = 2 + AMREX_SPACEDIM;

#if 0
void conv_3d_full(const MultiFab& in_real_mf, const MultiFab& in_imag_mf, 
                  const MultiFab& kernel_mf, MultiFab& out_real_mf, MultiFab& out_imag_mf, int NVAR){
  int row, col, dep, krow, kcol, kdep, orow, ocol, odep, kcenterx, kcentery, kcenterz;
  int ri, rj, rk, system_volume;

  for (MFIter mfi(in_real_mf); mfi.isValid(); ++mfi) {
    Array4<Real const> const in_real = in_real_mf.const_array(mfi);
    Array4<Real const> const in_imag = in_imag_mf.const_array(mfi);
    Array4<Real const> const kernel =  kernel_mf.const_array(mfi);
    Array4<Real> const& out_real =  out_real_mf.array(mfi);
    Array4<Real> const& out_imag =  out_imag_mf.array(mfi);

    Dim3 ihi = ubound(in_imag);
    Dim3 khi = ubound(kernel);
    Dim3 ohi = ubound(out_imag);

    Print() << ihi.x << ihi.y << ihi.z << "\n";
    Print() << khi.x << khi.y << khi.z << "\n";
    Print() << ohi.x << ohi.y << ohi.z << "\n";

    orow = ohi.x; ocol = ohi.y; odep = ohi.z;
    row = ihi.x; col = ihi.y; dep = ihi.z;
    krow = khi.x; kcol = khi.y; kdep = khi.z;

    system_volume = row*col*dep;

    for(int idx = 0; idx < NVAR; idx++){
      for(int i = 0; i < orow; i++){
        for(int j = 0; j < ocol; j++){
          for(int k = 0; k < odep; k++){
            for(int m = 0; m < krow; m++){
              for(int n = 0; n < kcol; n++){
                for(int l = 0; l < kdep; l++){
                  ri = i - m;
                  rj = j - n;
                  rk = k - l;
                  if (ri >= 0 && ri < row && rj >= 0 && rj < col && rk >= 0 and rk < dep){
                    out_real(i,j,k,idx) += in_real(rj,rj,rk,idx)*kernel(m,n,l,idx)/system_volume;
                    out_imag(i,j,k,idx) += in_imag(rj,rj,rk,idx)*kernel(m,n,l,idx)/system_volume;
                  }
                }
              }
            }
          }
        }
      }
    }
    //MFIter loop end
  }
}

void conv_3d_same(const MultiFab& in_real_mf, const MultiFab& in_imag_mf, 
                  const MultiFab& kernel_mf, MultiFab& out_real_mf, MultiFab& out_imag_mf, int NVAR){
                    // to do
                    // adapt function to use multifabs to make steps later on easier
  int row, col, dep, krow, kcol, kdep, orow, ocol, odep, kcenterx, kcentery, kcenterz;
  int ri, rj, rk;
  Real system_volume;

  for (MFIter mfi(in_real_mf); mfi.isValid(); ++mfi) {
    Array4<Real const> const in_real = in_real_mf.const_array(mfi);
    Array4<Real const> const in_imag = in_imag_mf.const_array(mfi);
    Array4<Real const> const kernel =  kernel_mf.const_array(mfi);
    Array4<Real> const& out_real =  out_real_mf.array(mfi);
    Array4<Real> const& out_imag =  out_imag_mf.array(mfi);

    Dim3 ihi = ubound(in_imag);
    Dim3 khi = ubound(kernel);
    Dim3 ohi = ubound(out_imag);

    // Print() << ihi.x << ihi.y << ihi.z << "\n";
    // Print() << khi.x << khi.y << khi.z << "\n";
    // Print() << ohi.x << ohi.y << ohi.z << "\n";
    
    orow = ohi.x; ocol = ohi.y; odep = ohi.z;
    row = ihi.x; col = ihi.y; dep = ihi.z;
    krow = khi.x; kcol = khi.y; kdep = khi.z;

    system_volume = row*col*dep;
    kcenterx = krow/2; kcentery = kcol/2; kcenterz = kdep/2;

    for(int idx = 0; idx < NVAR; idx++){
      for(int i = 0; i <= orow; i++){
        for(int j = 0; j <= ocol; j++){
          for(int k = 0; k <= odep; k++){
            for(int m = 0; m <= krow; m++){
              for(int n = 0; n <= kcol; n++){
                for(int l = 0; l <= kdep; l++){
                  ri = i + (m - kcenterx);
                  rj = j + (n - kcentery);
                  rk = k + (l - kcenterz);
                  if (ri >= 0 && ri < row && rj >= 0 && rj < col && rk >= 0 and rk < dep){
                    out_real(i,j,k,idx) += in_real(rj,rj,rk,idx)*kernel(m,n,l,idx)/system_volume;
                    out_imag(i,j,k,idx) += in_imag(rj,rj,rk,idx)*kernel(m,n,l,idx)/system_volume;
                    
                  }
                }
              }
            }
          }
        }
      }
    }
    //MFIter loop end
  }
}

void convolve_3D(const MultiFab& in_real_mf, const MultiFab& in_imag_mf, 
                 const MultiFab& kernel_mf, MultiFab& out_real_mf, MultiFab& out_imag_mf, 
                 std::string mode, int NVAR){
    

    if(mode == "same"){
      conv_3d_same(in_real_mf,in_imag_mf, kernel_mf, 
                                    out_real_mf, out_imag_mf, NVAR);}
    else if(mode == "full"){conv_3d_full(in_real_mf,in_imag_mf, kernel_mf, 
                                    out_real_mf, out_imag_mf, NVAR);}
    else{throw std::invalid_argument( "convolution mode is invalid" );exit(-1);}
}
#endif

void conv_3d_full(const Array3D<Real>& in_real, const Array3D<Real>& in_imag, 
                  const Array3D<Real>& kernel, Array3D<Real>& out_real, Array3D<Real>& out_imag){
  
  const int orow = out_real.xlen(); const int ocol = out_real.ylen(); const int odep = out_real.zlen();
  const int krow = kernel.xlen(); const int kcol = kernel.ylen(); const int kdep = kernel.zlen();
  const int row = in_real.xlen(); const int col = in_real.ylen(); const int dep = in_real.zlen();
  
  const int kcenterx = krow/2; const int kcentery = kcol/2; const int kcenterz = kdep/2;
  const int system_volume = row*col*dep;

  int ri, rj, rk;

  for(int i = 0; i <= orow; i++){
    for(int j = 0; j <= ocol; j++){
      for(int k = 0; k <= odep; k++){
        for(int m = 0; m <= krow; m++){
          for(int n = 0; n <= kcol; n++){
            for(int l = 0; l <= kdep; l++){
              ri = i - m;
              rj = j - n;
              rk = k - l;
              if (ri >= 0 && ri < row && rj >= 0 && rj < col && rk >= 0 and rk < dep){
                out_real(i,j,k) += in_real(rj,rj,rk)*kernel(m,n,l)/system_volume;
                out_imag(i,j,k) += in_imag(rj,rj,rk)*kernel(m,n,l)/system_volume;
              }
            }
          }
        }
      }
    }
  }
}

void conv_3d_same(const Array3D<Real>& in_real, const Array3D<Real>& in_imag, 
                  const Array3D<Real>& kernel, Array3D<Real>& out_real, Array3D<Real>& out_imag){

  const int orow = out_real.xlen(); const int ocol = out_real.ylen(); const int odep = out_real.zlen();
  const int krow = kernel.xlen(); const int kcol = kernel.ylen(); const int kdep = kernel.zlen();
  const int row = in_real.xlen(); const int col = in_real.ylen(); const int dep = in_real.zlen();
  
  const int kcenterx = krow/2; const int kcentery = kcol/2; const int kcenterz = kdep/2;
  const int system_volume = row*col*dep;

  int ri, rj, rk;

  for(int i = 0; i <= orow; i++){
    for(int j = 0; j <= ocol; j++){
      for(int k = 0; k <= odep; k++){
        for(int m = 0; m <= krow; m++){
          for(int n = 0; n <= kcol; n++){
            for(int l = 0; l <= kdep; l++){
              ri = i + (m - kcenterx);
              rj = j + (n - kcentery);
              rk = k + (l - kcenterz);
              if (ri >= 0 && ri < row && rj >= 0 && rj < col && rk >= 0 and rk < dep){
                out_real(i,j,k) += in_real(rj,rj,rk)*kernel(m,n,l)/system_volume;
                out_imag(i,j,k) += in_imag(rj,rj,rk)*kernel(m,n,l)/system_volume;
              }
            }
          }
        }
      }
    }
  }
}

void convolve_3D(const MultiFab& in_real_mf, const MultiFab& in_imag_mf, 
                 const MultiFab& kernel_mf, MultiFab& out_real_mf, MultiFab& out_imag_mf, 
                 std::string mode, int comp, const Box& domain_in, const Box& domain_out){
    
    IntVect iBoxdims = domain_in.length();
    IntVect oBoxdims = domain_out.length();
    const in row = iBoxdims(0); const in col = iBoxdims(1); const in dep = iBoxdims(2);
    const in nrow = oBoxdims(0); const in ncol = oBoxdims(1); const in ndep = 2*oBoxdims(2);

    Array3D<Real, 0, row, 0, col, 0, dep, Order::C> in_real_arr;
    Array3D<Real, 0, row, 0, col, 0, dep, Order::C> in_imag_ar;
    Array3D<Real, 0, row, 0, col, 0, dep, Order::C> kernel_ar;
    Array3D<Real, 0, nrow, 0, ncol, 0, ndep, Order::C> out_real_arr;
    Array3D<Real, 0, nrow, 0, ncol, 0, ndep, Order::C> out_imag_arr;

    for (MFIter mfi(in_real_mf); mfi.isValid(); ++mfi) {
      Array4<Real const> const in_real = in_real_mf.const_array(mfi);
      Array4<Real const> const in_imag = in_imag_mf.const_array(mfi);
      Array4<Real const> const kernel =  kernel_mf.const_array(mfi);

      ParallelForRNG(box, [=] AMREX_GPU_DEVICE(int x, int y, int z, RandomEngine const& engine) {
        in_real_arr(x,y,z) = in_real(x,y,z,comp);
        in_imag_arr(x,y,z) = in_imag(x,y,z,comp);
        kernel_arr(x,y,z) = kernel(x,y,z,comp);
    });
    }

    if(mode == "same"){
      conv_3d_full(in_real_arr,in_imag_arr, kernel_arr, out_real_arr, out_imag_arr);
    }else if(mode == "full"){
      conv_3d_full(in_real_arr,in_imag_arr, kernel_arr, out_real_arr, out_imag_arr);
    } else{throw std::invalid_argument( "convolution mode is invalid" );exit(-1);}

    for (MFIter mfi(out_real_mf); mfi.isValid(); ++mfi) {
      Array4<Real const> const out_real = out_real_mf.const_array(mfi);
      Array4<Real const> const out_imag = out_imag_mf.const_array(mfi);

      ParallelForRNG(box, [=] AMREX_GPU_DEVICE(int x, int y, int z, RandomEngine const& engine) {
        out_real(x,y,z,comp) = out_real_arr(x,y,z);
        out_imag(x,y,z,comp) = out_imag_arr(x,y,z);
    });
    }

}


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
GpuArray<GpuComplex<Real>,ndof> correlated_noise(int kx, int ky, int kz,
						   const Box& domain,
						   const RandomEngine& engine,
               const Array4<Real>& hydrovs_ref) {
  BL_PROFILE_VAR("correlated_noise()",correlated_noise);
  GpuArray<GpuComplex<Real>,ndof> r, xi;
  GpuArray<Real,ndof*ndof> C;
  Real k2;
  // GpuArray<GpuComplex<Real>,ndof*ndof> C;
  if (use_correlated_noise){k2 = fourier_laplace_operator(kx, ky, kz, domain);}
  else{k2 = 0.;}

  // // Cholesky decomposition of noise covariance matrix
  const Real rho0 = hydrovs_ref(kx, ky, kz, 0);
  const Real phi0 = hydrovs_ref(kx, ky, kz, 1);
  // if (ky == 0 and kz == 0){Print() << "kx:" << kx << " rho:" << rho0 << " phi:" << phi0 << " k2:" << k2 << "\n";}
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
  const int nvar = 2; const int nghost = 0; const Real system_volume = domain.length(0)*domain.length(1)*domain.length(2);
  BoxArray ba2 = kspace_noise_real.boxArray();
  DistributionMapping dm2 = kspace_noise_real.DistributionMap();
  Box domain2(IntVect{0,0,0}, IntVect{domain.length(0), domain.length(1), domain.length(2)});
  // Print() << "x:" << domain.length(0) << " y:" << domain.length(1) << " z:" << domain.length(2) << "\n";

  // FFT //
  MultiFab fft_real(ba2, dm2, nvar, nghost);
  MultiFab fft_imag(ba2, dm2, nvar, nghost);  
  compute_fft(geom, ref_params, fft_real, fft_imag, nvar, false); //output is correct. Normalization is 1/boxVol
  // FFT //
  
  // convolution //
  std::string mode = "same";
  MultiFab kernel(ba2, dm2, nvar, nghost);
  kernel.setVal(1.);

  MultiFab kref(ba2, dm2, nvar, nghost); kref.setVal(0.);
  MultiFab kref_imag(ba2, dm2, nvar, nghost); kref_imag.setVal(0.);
  MultiFab conv_real(ba2, dm2, nvar, nghost); conv_real.setVal(0.);
  MultiFab conv_imag(ba2, dm2, nvar, nghost); conv_imag.setVal(0.);
  // if(mode == "full"){
  //   Print() << "Full mode convolution\n";
  //   IntVect n = domain.length();
  //   Box domain2(IntVect{0,0,0}, IntVect{2*domain.length(0)-1, 2*domain.length(1)-1, 2*domain.length(2)-1});
  //   BoxArray ba2(domain2);
  //   DistributionMapping dm2(ba2);
  //   MultiFab kref(ba2, dm2, nvar, nghost);
  //   MultiFab kref_imag(ba2, dm2, nvar, nghost);
  //   MultiFab conv_real(ba2, dm2, nvar, nghost);
  //   MultiFab conv_imag(ba2, dm2, nvar, nghost);
  // }
  // Print() << "x:" << domain2.length()[0] << " y:" << domain2.length()[1] << " z:" << domain2.length()[2] << "\n";
  // Print() << "x:" << domain2.length(0) << " y:" << domain2.length(1) << " z:" << domain2.length(2) << "\n";
  // convolve_3D(fft_real, fft_imag, kernel, conv_real, conv_imag, mode, nvar);
  convolve_3D(fft_real, fft_imag, kernel, conv_real, conv_imag, mode, 0, domain, domain2);

  // conv_real.mult(1./system_volume, 0, nvar); conv_imag.mult(1./system_volume, 0, nvar);

  // for (MFIter mfi(conv_real); mfi.isValid(); ++mfi){
  //   const Box& box = mfi.validbox();
  //   const Array4<Real>& xi_real = conv_real.array(mfi);
  //   const Array4<Real>& xi_imag = conv_imag.array(mfi);
  //   const Array4<Real>& r = ref_params.array(mfi);

  //   ParallelForRNG(box, [=] AMREX_GPU_DEVICE(int kx, int ky, int kz, RandomEngine const& engine) {
  //     if (kx <= domain.length(0)/2) { // need only half of k-space for c2r FFT
	//       // compute correlated noise in k-space
  //       if (ky == 0 && kz == 0){
  //         Print() << "kx:" << kx << " rho:" << r(kx, ky, kz, 0) << " rhok_r:" << xi_real(kx, ky, kz, 0) << " rhok_i:" << xi_imag(kx, ky, kz, 0) << "\n";
  //         Print() << "kx:" << kx << " phi:" << r(kx, ky, kz, 1) << " phik_r:" << xi_real(kx, ky, kz, 1) << " phik_i:" << xi_imag(kx, ky, kz, 1) << "\n";
  //       }
  //     }
  //   });
  // }

  MultiFab::Copy(kref, conv_real, 0, 0, nvar, nghost); //copy works fine
  MultiFab::Copy(kref_imag, conv_imag, 0, 0, nvar, nghost); //copy works fine
  MultiFab::Multiply(kref, conv_real, 0, 0, nvar, nghost); //squaring. Checked to work
  MultiFab::Multiply(kref_imag, conv_imag, 0, 0, nvar, nghost); // squaring. Checked to work
  MultiFab::Subtract(kref, kref_imag, 0, 0, nvar, nghost); // multiplying convolution by its conjugate for reference state which acts as (a + b)*(a - b). Works
  // convolution //
  
  // generate noise in whole box because of grid decomposition
  // (generating noise without grid decomposition may be faster)
  for (MFIter mfi(kspace_noise_real); mfi.isValid(); ++mfi) {
    const Box& box = mfi.validbox();
    const Array4<Real>& xi_real = kspace_noise_real.array(mfi);
    const Array4<Real>& xi_imag = kspace_noise_imag.array(mfi);
    const Array4<Real>& r = ref_params.array(mfi);
    const Array4<Real>& kr = kref.array(mfi);

    // construct noise in k-space
    ParallelForRNG(box, [=] AMREX_GPU_DEVICE(int kx, int ky, int kz, RandomEngine const& engine) {
      if (kx <= domain.length(0)/2) { // need only half of k-space for c2r FFT
        GpuArray<GpuComplex<Real>,ndof> xi = {};
	      // compute correlated noise in k-space
	      xi = correlated_noise(kx,ky,kz,domain,engine,r);
        // if (ky == 0 && kz == 0){Print() << "kx:" << kx << " rhok:" << kr(kx, ky, kz, 0) << " phik:" << kr(kx, ky, kz, 1) << "\n";}

        for (int i=0; i<ndof; ++i) {
          xi_real(kx,ky,kz,i) = xi[i].real();
          xi_imag(kx,ky,kz,i) = xi[i].imag();
        }
      }
    });
  }

}

#if 0
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
#endif

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
  generate_kspace_noise(geom, kspace_noise_real, kspace_noise_imag, ref_params);
  compute_ifft(geom, noise, kspace_noise_real, kspace_noise_imag, ndof);
  // if (use_correlated_noise){
  //   // generate noise in k-space
  //   // Print() << "Correlated noise\n";
  //   generate_kspace_noise(geom, kspace_noise_real, kspace_noise_imag, ref_params);
  //   // note that the k-space noise is generated without the required symmetries
  //   // (this is to allow for parallel generation of noise)
  //   // the k-space symmetries are handled when copying in compute_ifft

  //   // inverse Fourier transform noise vector to real space
  //   compute_ifft(geom, noise, kspace_noise_real, kspace_noise_imag, ndof);
  // }
  // else{
  //   // Print() << "Uncorrelated noise\n";
  //   generate_realspace_noise(geom, ref_params, noise);
  // }

}
#endif