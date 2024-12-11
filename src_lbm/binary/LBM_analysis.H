#ifndef LBM_analysis_H_
#define LBM_analysis_H_

#ifdef AMREX_USE_CUDA
#include <culapacke.h> // this does not work yet. GPU support not implemented
#include <cufft.h>
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <lapack.h>
#endif

#include <AMReX_Array.H>

#include "LBM_IO.H"
#include <math.h>


// def droplet_radius_mass(density, Vp = 0, np_sphere = 0, rho_sphere = 1):
//     nx, ny, nz = density.shape
//     center_slc = np.s_[nx//2-1:nx//2+2, ny//2-1:ny//2+2, nz//2-1:nz//2+2]
//     edge_slc = np.s_[0:nx:nx-1, 0:ny:ny-1, 0:nz:nz-1]
//     if isinstance(density, int):
//         return np.nan
//     else:
//         # center = tuple([ l//2 for l in density.shape ])
//         rho_d = density[center_slc].mean()
//         rho_m = density[edge_slc].mean()
//         # mass = np.sum(density - rho_m) + 0.5*Vp*np_sphere*rho_sphere
//         mass = droplet_mass(density - rho_m) + 0.5*Vp*np_sphere*rho_sphere
//         R = (3./4./np.pi*mass/(rho_d-rho_m))**(1./3.)
//         return R
Real droplet_radius(MultiFab& hydrovars, int comp = 1){
    Real R = 0; Real rhod = 0; Real rhom = 0; Real mass = 0;

  for (MFIter mfi(hydrovars); mfi.isValid(); ++mfi) {
    const Array4<Real>& h = hydrovars.array(mfi);
    // averaging a cubic box of size 3, centered around the center of the simulation domain for rhod
    for (int i = -1; i <= 1; i++){
        for (int j = -1; j <= 1; j++){
            for (int k = -1; k <= 1; k++){
                rhod += h(nx/2 + i, ny/2 + j, nz/2 + k, comp);
            }
        }
    }
    // rhod += h(nx/2, ny/2, nz/2, comp);
    // averaging a corner points, for rhom
    rhom += h(0, 0, 0, comp); 
    rhom += h(nx - 1, 0, 0, comp); 
    rhom += h(0, ny - 1, 0, comp); 
    rhom += h(0, 0, nz - 1, comp);
    rhom += h(nx - 1, ny - 1, 0, comp); 
    rhom += h(0, ny - 1, nz - 1, comp); 
    rhom += h(nx - 1, 0, nz - 1, comp);
    rhom += h(nx - 1, ny - 1, nz - 1, comp);
  }
    rhod /= 27;
    rhom /= 8;
    
    // Print() << "rhod:" << rhod << "\n";
    // Print() << "rhom:" << rhom << "\n";

    BoxArray ba = hydrovars.boxArray();
    DistributionMapping dm = hydrovars.DistributionMap();
    MultiFab mass_onegrid(ba, dm, 1, 0); // ghost cells of hydrovars is 2
    mass_onegrid.ParallelCopy(hydrovars, comp, 0, 1);

    mass_onegrid.plus(-rhom, 0);
    mass = mass_onegrid.sum(0);

    R = pow((3.*mass)/(4.*M_PI*(rhod - rhom)), 1./3.);
    return R;
}


// def get_indices(nx, ny, nz):
//     idxs_test = np.zeros((3, nx, ny, nz))
//     for x in range(nx):
//         for y in range(ny):
//             for z in range(nz):
//                 idxs_test[0, x, y, z] = x
//                 idxs_test[1, x, y, z] = y
//                 idxs_test[2, x, y, z] = z
    
//     return idxs_test
MultiFab get_indices(MultiFab& hydrovars){
    BoxArray ba = hydrovars.boxArray();
    DistributionMapping dm = hydrovars.DistributionMap();
    MultiFab indices(ba, dm, 3, 0);

    auto const & f = indices.arrays();
    ParallelFor(indices, IntVect(0), [=] AMREX_GPU_DEVICE(int nbx, int x, int y, int z) {
        f[nbx](x, y, z, 0) = x;
        f[nbx](x, y, z, 1) = y;
        f[nbx](x, y, z, 2) = z;
    });
    return indices;
}

// def com(OutArray):
//     total_mass = np.sum(OutArray)
//     com = np.zeros(3)
//     idxs = get_indices(nx, ny, nz)
    
//     for i in range(3):
//         field1 = idxs[0]
//         com[i] = np.sum(field1*OutArray)/total_mass
//     return com
GpuArray<Real,3> center_of_mass(MultiFab& hydrovars, int comp = 1){
    GpuArray<Real,3> com;
    Real total_mass = hydrovars.sum(comp);

    MultiFab indices;
    indices = get_indices(hydrovars);

    BoxArray ba = hydrovars.boxArray();
    DistributionMapping dm = hydrovars.DistributionMap();
    MultiFab com_onegrid(ba, dm, 1, 0); // ghost cells of hydrovars is 2

    for (int i = 0; i < 3; i++){
        com_onegrid.ParallelCopy(indices, i, 0, 1);
        com_onegrid.Multiply(com_onegrid, hydrovars, comp, 0, 1, 0);
        com[i] = com_onegrid.sum(0)/total_mass;
    }
    return com;
}

// def gyration_tensor_new(cm, OutArray):
//     nx, ny, nz = OutArray.shape
//     idxs = get_indices(nx, ny, nz)
//     total_mass = np.sum(OutArray)
//     S = np.zeros(9)

//     for i in range(3):
//         for j in range(3):
//             curr_i = idxs[i].copy()
//             curr_j = idxs[j].copy()

//             curr_i -= cm[i]
//             curr_j -= cm[j]

//             curr_gyr = OutArray*curr_i*curr_j/total_mass

//             S[i*3+j] = curr_gyr.sum()
//     return S
GpuArray<Real,9> gyration_tensor(GpuArray<Real,3> com, MultiFab& hydrovars, int comp = 1){
    // Array1D<Real,0,9> S;
    GpuArray<Real, 9> S;
    Real total_mass = hydrovars.sum(comp);

    MultiFab indices;
    indices = get_indices(hydrovars);

    BoxArray ba = hydrovars.boxArray();
    DistributionMapping dm = hydrovars.DistributionMap();
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
            com_onegrid_i.Multiply(com_onegrid_i, hydrovars, comp, 0, 1, 0);
            S[i*3+j] = com_onegrid_i.sum(0)/total_mass;
        }
    }
    return S;
}

GpuArray<Real,3> axial_radii(MultiFab& hydrovars, int comp = 1){
    BL_PROFILE_VAR("axial_radii()",axial_radii);
    GpuArray<Real,3> com;
    GpuArray<Real,9> S;
    GpuArray<Real,3> radii;

    com = center_of_mass(hydrovars, comp);
    S = gyration_tensor(com, hydrovars, comp);

    const char JOBVL = 'N'; // const char* as input to eigenval calc
    const char JOBVR = 'N'; // const char* as input to eigenval calc
    GpuArray<Real,3> wr, wi; // Output data structures for eigenvalues. r is real, i is imaginary
    GpuArray<Real,9> vl, vr, WORK; // Output data structures for eigenvectors. r is real, i is imaginary
    const int orderS = 3; // Order of matrix S. const int* input to eigenval calc
    const int lda = 3; // leading dimension of matrix S. const int* input to eigenval calc
    const int eigvec_sz = 3; // Size of eigenvector
    const int LWORK = 3*orderS; // Span*size of eigenvector matrix
    int info = 0; // output condition of function
    size_t SSIZE = S.size(); // memory allocation size of array for eigvals to be found
    size_t WSIZE = WORK.size(); // Size of array where eigen vectors are stored.

    //dgeev_ is from LAPACK and is used to calculate the eigenvalues of the gyration tensor
    dgeev_(&JOBVL, &JOBVR, &orderS, S.data(), &lda, wr.data(), wi.data(), vl.data(), &eigvec_sz, vr.data(), &eigvec_sz, WORK.data(), &LWORK, &info, SSIZE, WSIZE);
    // const char*  = "N"
    // const char*  = "N"
    // const int*   = &orderS
    // double*      = S.data()
    // const int*   = &lda
    // double*      = wr.data()
    // double*      = wi.data()
    // double*      = vl.data()
    // const int*   = &eigvec_sz
    // double*      = vr.data()
    // const int*   = &eigvec_sz
    // double*      = WORK.data()
    // const int*   = &LWORK
    // int*         = info
    // size_t       = S.size()
    // size_t       = WORK.size()

    radii[0] = pow(wr[0], 1./3.)/pow(wr[1]*wr[2], 1./6.);
    radii[1] = pow(wr[1], 1./3.)/pow(wr[0]*wr[2], 1./6.);
    radii[2] = pow(wr[2], 1./3.)/pow(wr[0]*wr[1], 1./6.);

    return radii;
} 

// void write_data()

#endif