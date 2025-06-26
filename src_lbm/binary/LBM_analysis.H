#ifndef LBM_analysis_H_
#define LBM_analysis_H_

#ifdef AMREX_USE_CUDA
#include <culapacke.h> // this does not work yet. GPU support not implemented
#include <cufft.h>
#else
#include <lapack.h>
#include <fftw3.h>
#endif

#include <math.h>

#include "LBM_binary.H"

MultiFab binarize_droplet(const MultiFab& hydrovars, int comp = 1, Real cutoff = 0.){
    BoxArray ba = hydrovars.boxArray();
    DistributionMapping dm = hydrovars.DistributionMap();
    MultiFab droplet(ba, dm, 1, 0);
    
    droplet.ParallelCopy(hydrovars, comp, 0, 1);
    auto const & d = droplet.arrays();
    ParallelFor(droplet, IntVect(0), [=] AMREX_GPU_DEVICE(int nbx, int x, int y, int z) {
        if (d[nbx](x, y, z, 0) <= cutoff){
            d[nbx](x, y, z, 0) = 0;}
    });
    return droplet;
}

Real droplet_radius(MultiFab& droplet){
    BoxArray ba = droplet.boxArray();
    DistributionMapping dm = droplet.DistributionMap();

    Box domain_box = ba.minimalBox();
    const IntVect domain_size2 = domain_box.size();
    const int nx = domain_size2[0]; const int ny = domain_size2[1]; const int nz = domain_size2[2];

    IntVect slc_lo, slc_hi;
    Real R = 0.; Real mass = 0.;
    Real rhod = 0.; Real rhom = 0.;

    slc_lo = {nx/2-1, ny/2-1, nz/2-1}; slc_hi = {nx/2+1, ny/2+1, nz/2+1};
    Box slice1(slc_lo, slc_hi);
    rhod = droplet.sum(slice1, 0);
    rhod /= 27.;

    slc_lo = {-1,-1,-1}; slc_hi = {1,1,1};
    Box slice2(slc_lo, slc_hi);
    rhom = droplet.sum(slice2, 0);
    rhom /= 8.;

    // Print() << "rhod: " << rhod << "\n";
    // Print() << "rhom: " << rhom << "\n";

    MultiFab mass_onegrid(ba, dm, 1, 0); // ghost cells is 0
    mass_onegrid.ParallelCopy(droplet, 0, 0, 1);

    mass_onegrid.plus(-rhom, 0);
    mass = mass_onegrid.sum(0);

    R = pow((3.*mass)/(4.*M_PI*(Real(rhod) - Real(rhom))), 1./3.);
    return R;
}

MultiFab get_indices(MultiFab& droplet){
    BoxArray ba = droplet.boxArray();
    DistributionMapping dm = droplet.DistributionMap();
    MultiFab indices(ba, dm, 3, 0);

    auto const & f = indices.arrays();
    ParallelFor(indices, IntVect(0), [=] AMREX_GPU_DEVICE(int nbx, int x, int y, int z) {
        f[nbx](x, y, z, 0) = x;
        f[nbx](x, y, z, 1) = y;
        f[nbx](x, y, z, 2) = z;
    });
    return indices;
}

GpuArray<Real,3> center_of_mass(MultiFab& droplet){
    GpuArray<Real,3> com;
    Real total_mass = droplet.sum(0);

    MultiFab indices;
    indices = get_indices(droplet);

    BoxArray ba = droplet.boxArray();
    DistributionMapping dm = droplet.DistributionMap();
    MultiFab com_onegrid(ba, dm, 1, 0); // ghost cells of hydrovars is 2

    for (int i = 0; i < 3; i++){
        com_onegrid.ParallelCopy(indices, i, 0, 1);
        com_onegrid.Multiply(com_onegrid, droplet, 0, 0, 1, 0);
        com[i] = com_onegrid.sum(0)/total_mass;
    }
    return com;
}

// Fitting expression to the density profile
#ifndef AMREX_USE_CUDA
void sort_by_r(const std::vector<Real>& r,
               const std::vector<Real>& rho,
               std::vector<Real>& r_sorted,
               std::vector<Real>& rho_sorted){
    int N = r.size();
    // std::vector<int> idxs(N);
    std::vector<int> idxs(N);
    std::iota(idxs.begin(), idxs.end(), 0); // idxs = [0, 1, ..., N-1]

    // Sort indices based on values in r
    std::sort(idxs.begin(), idxs.end(),
            [&r](int i1, int i2) { return r[i1] < r[i2]; });

    // Reorder r and rho based on sorted indices
    // r_sorted.resize(N);
    // rho_sorted.resize(N);
    for (int i = 0; i < N; ++i) {
        r_sorted[i] = r[idxs[i]];
        rho_sorted[i] = rho[idxs[i]];
    }
}

inline Real model(Real r, Real a, Real R, Real zeta, Real b) {
    return a*std::tanh((r - R)/(std::sqrt(2.0)*zeta)) + b;
}

Real compute_cost(const std::vector<Real>& r,
                    const std::vector<Real>& rho,
                    Real a, Real R, Real zeta, Real b) {
    Real cost = 0.0;
    for (size_t i = 0; i < r.size(); ++i) {
        Real res = model(r[i], a, R, zeta, b) - rho[i];
        cost += res * res;
    }
    return cost;
}

void compute_gradient(const std::vector<Real>& r,
                      const std::vector<Real>& rho,
                      Real a, Real R, Real zeta, Real b,
                      std::vector<Real>& grad,
                      Real epsilon = 1e-6) {
    Real base_cost = compute_cost(r, rho, a, R, zeta, b);

    // grad.resize(4);

    grad[0] = (compute_cost(r, rho, a + epsilon, R, zeta, b) - base_cost) / epsilon; // ∂/∂a
    grad[1] = (compute_cost(r, rho, a, R + epsilon, zeta, b) - base_cost) / epsilon; // ∂/∂R
    grad[2] = (compute_cost(r, rho, a, R, zeta + epsilon, b) - base_cost) / epsilon; // ∂/∂zeta
    grad[3] = (compute_cost(r, rho, a, R, zeta, b + epsilon) - base_cost) / epsilon; // ∂/∂b
}

std::vector<Real> fit_tanh_model(const std::vector<Real> r,
                    const std::vector<Real> rho,
                    std::vector<Real> init_params,
                    int max_iters = 1e4, Real lr = 1e-6) {
    
    Real a = init_params[0]; Real R = init_params[1]; Real zeta = init_params[2]; Real b = init_params[3];
    std::vector<Real> fit_params(4);

    Real previous_cost = 0.0;
    Real threshold = 1e-5;

    for (int iter = 0; iter < max_iters; ++iter) {
        std::vector<Real> grad(4);
        compute_gradient(r, rho, a, R, zeta, b, grad);

        // Update parameters
        a    -= lr * grad[0];
        R    -= lr * grad[1];
        zeta -= lr * grad[2];
        b    -= lr * grad[3];

        if (iter % 100 == 0) {
            Real cost = compute_cost(r, rho, a, R, zeta, b);
            if (Math::abs(cost - previous_cost) < threshold){break;}
            else{previous_cost = cost;}
        }
    }
    fit_params[0] = a; fit_params[1] = R; fit_params[2] = zeta; fit_params[3] = b;
    return fit_params;
}

std::vector<Real> copy_mf_to_vector(MultiFab& mf, std::vector<Real> vec){
    int offset = 0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();  // or tilebox()
        const auto& fab = mf[mfi];

        const Real* fab_ptr = fab.dataPtr(0);
        int npts = box.numPts();

        std::copy(fab_ptr, fab_ptr + npts, vec.begin() + offset);
        offset += npts;
    }

    return vec;
}

Real droplet_radius_profile_fit(MultiFab& droplet){
    BL_PROFILE_VAR("droplet_radius_profile_fit()",droplet_radius_profile_fit);
    Real fitted_radius;

    GpuArray<Real,3> com;
    com = center_of_mass(droplet);
        
    MultiFab indices;
    indices = get_indices(droplet);
    
    BoxArray ba = droplet.boxArray();
    DistributionMapping dm = droplet.DistributionMap();

    Box domain_box = ba.minimalBox();
    const IntVect domain_size = domain_box.size();
    const int nx = domain_size[0]; const int ny = domain_size[1]; const int nz = domain_size[2];
    int N = nx*ny*nz;
    // Print() << "nx: " << nx << ", ny: " << ny << ", nz: " << nz << ", N: "<< N << '\n';
    MultiFab R(ba, dm, 1, 0);

    // adding square of indices to radius before square root

    auto const & radius_arr = R.arrays();
    auto const & ind_arr = indices.arrays();
    ParallelFor(R, IntVect(0), [=] AMREX_GPU_DEVICE(int nbx, int x, int y, int z) {
        Real xdiff = Real(ind_arr[nbx](x,y,z,0)) - Real(com[0]);
        Real ydiff = Real(ind_arr[nbx](x,y,z,1)) - Real(com[1]);
        Real zdiff = Real(ind_arr[nbx](x,y,z,2)) - Real(com[2]);

        radius_arr[nbx](x,y,z,0) = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);
    });
    
    // need to convert from multifab to another data object for sorting using std::sort. NOT GPU COMPATIBLE FOR NOW
    std::vector<Real> R_flat(N), rho_flat(N), R_flat_sorted(N), rho_flat_sorted(N);

    R_flat = copy_mf_to_vector(R, R_flat);
    rho_flat = copy_mf_to_vector(droplet, rho_flat);

    // Print() << "Created object\n";
    // for (MFIter mfi(droplet); mfi.isValid(); ++mfi){
    //   const Array4<Real>& density = droplet.array(mfi);
    //   const Array4<Real>& f = R.array(mfi);

    //   for(int i = 0; i < nx; i++){
    //     for(int j = 0; j < ny; j++){
    //       for(int k = 0; k < nz; k++){
    //         int idx = i + nx*(j + ny*k);
    //         // Print() << density(i,j,k,0) << "\n";
    //         R_flat[idx] = f(i,j,k,0);
    //         rho_flat[idx] = density(i,j,k,0);
    //         Print() << "idx:" << idx << ", r_flat:" << R_flat[idx] << ", rho_flat:" << rho_flat[idx] << "\n";
    //                     }
    //                 }
    //             }
    //         }

    Print() << "Copy complete" << "\n";
    // Print() << "R at index 0: " << R_flat[0] << ", droplet at index 0: " << rho_flat[0] << "\n";
    // Print() << "R at index N/2: " << R_flat[N/2] << ", droplet at index N/2: " << rho_flat[N/2] << "\n";
    // Print() << "R at index -1: " << R_flat[N - 1] << ", droplet at index N/2: " << rho_flat[N - 1] << "\n";

    sort_by_r(R_flat, rho_flat, R_flat_sorted, rho_flat_sorted);
    Print() << "R at index 0: " << R_flat_sorted[0] << ", droplet at index 0: " << rho_flat_sorted[0] << "\n";
    Print() << "R at index N/2: " << R_flat_sorted[N/2] << ", droplet at index N/2: " << rho_flat_sorted[N/2] << "\n";
    Print() << "R at index -1: " << R_flat_sorted[N - 1] << ", droplet at index N/2: " << rho_flat_sorted[N - 1] << "\n";

    std::vector<Real> init_params = {(droplet.min(0) - droplet.max(0))/2., droplet_radius_prop*nx, sqrt(kappa/(chi/2. - T)),(droplet.min(0) + droplet.max(0))/2.}; // a, R, zeta, b
    std::vector<Real> fitted_params(4);

    fitted_params = fit_tanh_model(R_flat_sorted, rho_flat_sorted,init_params);
    fitted_radius = fitted_params[1];
    // fitted_radius = 0.3*nx;

    return fitted_radius;
}
#endif
// Fitting expression to the density profile

GpuArray<Real,9> gyration_tensor(GpuArray<Real,3> com, MultiFab& droplet){
    // Array1D<Real,0,9> S;
    GpuArray<Real, 9> S;
    Real total_mass = droplet.sum(0);

    MultiFab indices;
    indices = get_indices(droplet);

    BoxArray ba = droplet.boxArray();
    DistributionMapping dm = droplet.DistributionMap();
    MultiFab com_onegrid_i(ba, dm, 1, 0); // ghost cells of hydrovars is 2
    MultiFab com_onegrid_j(ba, dm, 1, 0); // ghost cells of hydrovars is 2
    MultiFab product(ba, dm, 1, 0); // ghost cells of hydrovars is 2

    Real com_i, com_j;

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            com_onegrid_i.ParallelCopy(indices, i, 0, 1);
            com_onegrid_j.ParallelCopy(indices, j, 0, 1);

            // curr_i -= cm[i]
            // curr_j -= cm[j]
            // com_i = com(i); com_j = com(j);
            com_i = com[i]; com_j = com[j];
            com_onegrid_i.plus(-com_i, 0);
            com_onegrid_j.plus(-com_j, 0);

            // curr_gyr = OutArray*curr_i*curr_j/total_mass
            com_onegrid_i.Multiply(com_onegrid_i, com_onegrid_j, 0, 0, 1, 0);
            com_onegrid_i.Multiply(com_onegrid_i, droplet, 0, 0, 1, 0);
            S[i*3+j] = com_onegrid_i.sum(0)/total_mass;
        }
    }
    return S;
}

GpuArray<Real,3> axial_radii(MultiFab& droplet){
    BL_PROFILE_VAR("axial_radii()",axial_radii);
    GpuArray<Real,3> com;
    GpuArray<Real,9> S;
    GpuArray<Real,3> radii;

    com = center_of_mass(droplet);
    S = gyration_tensor(com, droplet);

    const char JOBVL = 'N'; // const char* as input to eigenval calc
    const char JOBVR = 'N'; // const char* as input to eigenval calc
    const int orderS = 3; // Order of matrix S. const int* input to eigenval calc
    const int lda = 3; // leading dimension of matrix S. const int* input to eigenval calc
    GpuArray<Real,3> wr, wi; // Output data structures for eigenvalues. r is real, i is imaginary
    GpuArray<Real,9> vl, vr, WORK; // Output data structures for eigenvectors. r is real, i is imaginary
    const int eigvec_sz = 3; // Size of eigenvector
    const int LWORK = 3*orderS; // Span*size of eigenvector matrix
    int info = 0; // output condition of function

    //dgeev_ is from LAPACK and is used to calculate the eigenvalues of the gyration tensor
    LAPACK_dgeev(&JOBVL, &JOBVR, &orderS, S.data(), &lda, wr.data(), wi.data(), vl.data(), &eigvec_sz, vr.data(), &eigvec_sz, WORK.data(), &LWORK, &info);

    radii[0] = pow(wr[0], 1./3.)/pow(wr[1]*wr[2], 1./6.);
    radii[1] = pow(wr[1], 1./3.)/pow(wr[0]*wr[2], 1./6.);
    radii[2] = pow(wr[2], 1./3.)/pow(wr[0]*wr[1], 1./6.);

    return radii;
} 

#if 0
// implementing 1D FFT to directly calculate height fluctuation in AMReX
void forward_fft_1D(const GpuArray<Real>& height, GpuArray<Real>& real, GpuArray<Real>& imag, GpuArray<Real>& freqs, GpuArray<Real>& mags){
    // inspired by https://www.fftw.org/doc/Complex-One_002dDimensional-DFTs.html and https://www.fftw.org/doc/One_002dDimensional-DFTs-of-Real-Data.html
    
    const int N = height.size()
    const int FFTW_ESTIMATE;

    // create data structures for doing the FFT 
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = fftw_malloc ( sizeof ( fftw_complex ) * n);
    out = fftw_malloc ( sizeof ( fftw_complex ) * m * n );
    
    // Copy input data to FFTW input array
    for (int i = 0; i < n; ++i) {
        in[i] = height[i];
    }

    // Create FFTW plan
    fftw_plan plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    // Execute FFT plan
    fftw_execute(plan);

    // moving fftw output data to AMReX data structures for further processing
    for (int i = 0; i < n; i++){
        Real k = (i < (n+1)/2) ? 2.*M_PI/n*i : 2.*M_PI/n*(i-n);
        real[i] = out[i][0];
        imag[i] = out[i][1];
        freqs[i] = k;
        mags[i] = real[i]*real[i] + imag[i]*imag[i];
    }

    // Clean up
    fftw_destroy_plan(plan); fftw_free(in); fftw_free(out);
}

void conjugate_complex(GpuArray<Real>& real, GpuArray<Real>& imag, GpuArray<Real>& mag){
    const int N = real.size()
    for (int i = 0; i < N; i++){
        mag[i] = real[i]*real[i] + imag[i]*imag[i];
    }
}
#endif

#endif