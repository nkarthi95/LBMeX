#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include "StructFact.H"
using namespace amrex;
#include "LBM_binary.H"
#include "tests.H"
#include "LBM_IO.H"
#include "LBM_thermodynamics.H"
#include "LBM_analysis.H"

#include <iostream>
#include <fstream>

void main_driver(const char* argv) {
  // test_case_fft();
  // if (!cholesky_test(100)) exit(-1);

  // store the current time so we can later compute total run time.
  Real strt_time = ParallelDescriptor::second();
  
  const std::string dist_plt = "xi_plt_";
  const std::string hydro_plt = "hydro_plt_";
  const std::string SF_plt = "SF_plt";
  const std::string hydro_chk = "chk_hydro_";
  const std::string ref_plt = "ref_plt_";

  // default grid parameters
  // int nx = 16; int ny = 16; int nz = 16;
  IntVect max_grid_size = {nx, ny, nz};
  // int ic = 0;
  // Real R = 0.3;

  // default time stepping parameters
  int n_sci_start = 0;
  int nsteps = 100;
  int dump_hydro = 1;
  int n_hydro = 10;
  int dump_SF = 1;
  int n_SF = 10;
  int start_step = 0;
  int output_hdf = 0;

  // run output
  int record_start = 550000;
  int record_interval = 10;
  int binarize = 1;

  // input parameters
  ParmParse pp;
  // box parameters
  pp.query("nx", nx);
  pp.query("ny", ny);
  pp.query("nz", nz);
  pp.query("max_grid_size_x", max_grid_size[0]);
  pp.query("max_grid_size_y", max_grid_size[1]);
  pp.query("max_grid_size_z", max_grid_size[2]);
  pp.query("init_cond", ic);
  pp.query("nz", nz);
  pp.query("R", R);

  // plot parameters
  pp.query("nsteps", nsteps);
  int n_checkpoint = nsteps;
  pp.query("n_sci_start", n_sci_start);
  pp.query("dump_hydro", dump_hydro);
  pp.query("n_hydro", n_hydro);
  pp.query("dump_SF", dump_SF);
  pp.query("n_SF", n_SF);
  pp.query("n_checkpoint", n_checkpoint);
  pp.query("start_step", start_step);
  pp.query("output_hdf5", output_hdf);

  // model parameters
  pp.query("kappa", kappa);
  pp.query("lambda", chi);
  pp.query("T", T);
  pp.query("gamma", Gamma);

  //fluctuation parameters
  pp.query("temperature", temperature);
  pp.query("correlated_noise", use_correlated_noise);

  // set up Box and Geomtry
  IntVect dom_lo(0, 0, 0);
  IntVect dom_hi(nx-1, ny-1, nz-1);
  Array<int,3> periodicity({1,1,1});
  Box domain(dom_lo, dom_hi);
  RealBox real_box({0.,0.,0.},{1.,1.,1.});
  Geometry geom(domain, real_box, CoordSys::cartesian, periodicity);
  BoxArray ba(domain);
  // split BoxArray into chunks no larger than "max_grid_size" along a direction
  ba.maxSize(max_grid_size);
  DistributionMapping dm(ba);
  // need two halo layers for gradients
  Real time = start_step;
  const int nghost = 2;
  
  // set up MultiFabs
  MultiFab fold(ba, dm, nvel, nghost);
  MultiFab fnew(ba, dm, nvel, nghost);
  MultiFab gold(ba, dm, nvel, nghost);
  MultiFab gnew(ba, dm, nvel, nghost);
  MultiFab hydrovs(ba, dm, 2*nvel, nghost);
  MultiFab noise(ba, dm, 2*nvel, nghost);
  MultiFab ref_params(ba, dm, 2, nghost); //reference rho and C for each point of the lattice
  MultiFab mf_checkpoint(ba, dm, 2*nvel+2, nghost);

  // structure factor stuff
  int nStructVars = 5;
  const Vector<std::string> var_names = VariableNames(nStructVars);
  Vector<Real> var_scaling(nStructVars*(nStructVars+1)/2); var_scaling.assign(var_scaling.size(), 1.);
  StructFact structFact(ba, dm, var_names, var_scaling);

  // INITIALIZE
  switch(ic){
    case 0:
      LBM_init_mixture(fold, gold, hydrovs, ref_params);
      break;
    case 1:
      LBM_init_flat_interface(geom, fold, gold, hydrovs, ref_params);
      break;
    case 2:
      LBM_init_droplet(R, geom, fold, gold, hydrovs, ref_params);
      break;
    case 30:
      LBM_init_cylinder_x(R, geom, fold, gold, hydrovs, ref_params);
      break;
    case 31:
      LBM_init_cylinder_y(R, geom, fold, gold, hydrovs, ref_params);
      break;
    case 32:
      LBM_init_cylinder_z(R, geom, fold, gold, hydrovs, ref_params);
      break;
    case 10:
      checkpointRestart(start_step, mf_checkpoint, hydro_chk, fold, gold, ba, dm);--start_step;
      hydrovs.ParallelCopy(mf_checkpoint, 0, 0, 2*nvel);
      ref_params.ParallelCopy(mf_checkpoint, 2*nvel, 0, 2);
      break;
  }

  if (n_checkpoint > 0 && ic != 10){
    // WriteCheckPoint(start_step, hydrovs, hydro_chk);start_step = 0;
    mf_checkpoint.ParallelCopy(hydrovs,0,0,2*nvel);
    mf_checkpoint.ParallelCopy(ref_params,0,2*nvel,2);
    WriteCheckPoint(start_step, mf_checkpoint, hydro_chk);start_step = 0;
    WriteOutput(start_step, ref_params, geom, ref_plt, 2, output_hdf);
    }
  // checkpoint read of hydrovs to generate fold and gold to be used for further simulations

  // Write a plotfile of the initial data if plot_int > 0
  if (dump_hydro == 1){WriteOutput(start_step, hydrovs, geom, hydro_plt, 2*nvel, output_hdf);}
  Print() << "LB initialized\n";
  start_step++;

  MultiFab droplet(ba, dm, 1, 0);
  std::ofstream outfile("example.csv");
  outfile << "Timestep,Radius,cx,cy,cz,dx,dy,dz" << std::endl;
  if (binarize == 1){droplet = binarize_droplet(hydrovs, 1, 0.);} else{droplet.ParallelCopy(hydrovs, 1, 0, 1);}
  Real R = droplet_radius(droplet); // AMReX implementation assumes binarized droplet as input.
  GpuArray<Real, 3> com = center_of_mass(droplet);
  GpuArray<Real, 3> dr = axial_radii(droplet);
  outfile << start_step-1 << "," << R << "," << com[0] << "," << com[1] << "," << com[2] << "," << dr[0] << "," << dr[1] << "," << dr[2] << std::endl; //"Timestep, Radius, dx, dy, dz\n" 
  // myfile.open("example.txt");

  // TIMESTEP
  for (int step=start_step; step <= nsteps; ++step) {
    LBM_timestep(geom, fold, gold, fnew, gnew, hydrovs, noise, ref_params);
    
    // droplet.ParallelCopy(hydrovars, comp, 0, 1);
    // droplet = binarize_droplet(hydrovs, 1, 0.);
    // Real R = droplet_radius(droplet); // AMReX implementation assumes binarized droplet as input.
    // GpuArray<Real, 3> dr = axial_radii(droplet);
    // // testing //
    // Print() << "Droplet Radius: " << R << "\n";
    // GpuArray<Real, 3> com = center_of_mass(droplet);
    // Print() << "COM [x, y, z]: " << "[" << com[0] << ", " << com[1] << ", " << com[2] << "]\n";
    // GpuArray<Real, 9> S = gyration_tensor(com, droplet);
    // for (int i = 0; i < 3; i++){
    //   for (int j = 0; j < 3; j++){
    //     Print() <<"idx: " << i*3 + j << " " << S[i*3 + j] << "\n";
    //   }
    //   // Print() << "]\n";
    // }
    // Print() << "dr [dx, dy, dz]: " << "[" << dr[0] << ", " << dr[1] << ", " << dr[2] << "]\n";
    // // testing //

    if (step >= record_start and step%record_interval == 0){
      // droplet.ParallelCopy(hydrovs, 1, 0, 1);
      // droplet = binarize_droplet(hydrovs, 1, 0.);
      if (binarize == 1){droplet = binarize_droplet(hydrovs, 1, 0.);} else{droplet.ParallelCopy(hydrovs, 1, 0, 1);}
      Real R = droplet_radius(droplet); // AMReX implementation assumes binarized droplet as input.
      GpuArray<Real, 3> com = center_of_mass(droplet);
      GpuArray<Real, 3> dr = axial_radii(droplet);
      outfile << step << "," << R << "," << com[0] << "," << com[1] << "," << com[2] << "," << dr[0] << "," << dr[1] << "," << dr[2] << std::endl; //"Timestep, Radius, dx, dy, dz\n" 
    }

    if (step >= n_sci_start){

    if (dump_SF == 1 && temperature > 0){structFact.FortStructure(hydrovs, geom);}
    
    if (dump_hydro == 1 && step%n_hydro == 0){
      WriteOutput(step, hydrovs, geom, hydro_plt, 2*nvel, output_hdf);
      // WriteOutput(step, noise, geom, dist_plt, 2*nvel, output_hdf);
      }

    if(dump_SF == 1 && step%n_SF == 0 && temperature > 0){
      // structFact.WritePlotFile(step, static_cast<Real>(step), geom, SF_plt, 0);
      structFact.WritePlotFile(step, static_cast<Real>(step), geom, SF_plt);
      StructFact structFact(ba, dm, var_names, var_scaling);
      }
    }

    if (n_checkpoint > 0 && step%n_checkpoint == 0){
      mf_checkpoint.ParallelCopy(hydrovs,0,0,2*nvel);
      mf_checkpoint.ParallelCopy(ref_params,0,2*nvel,2);
      WriteCheckPoint(step, mf_checkpoint, hydro_chk);
    }
    Print() << "LB step " << step << " completed\n";
  }
  // Call the timer again and compute the maximum difference between the start time 
  // and stop time over all processors
  outfile.close();
  Real stop_time = ParallelDescriptor::second() - strt_time;
  ParallelDescriptor::ReduceRealMax(stop_time);
  amrex::Print() << "Run time = " << stop_time << std::endl;
}