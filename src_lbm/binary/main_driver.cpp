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

void main_driver(const char* argv) {

  // if (!cholesky_test(100)) exit(-1);

  // store the current time so we can later compute total run time.
  Real strt_time = ParallelDescriptor::second();
    
  // default grid parameters
  int nx = 16;
  int max_grid_size = 8;
  int ic = 0;

  // default time stepping parameters
  int nsteps = 100;
  int plot_int = 10;
  int n_checkpoint = nsteps;
  int start_step = 0;

  // input parameters
  ParmParse pp;
  // box parameters
  pp.query("nx", nx);
  pp.query("max_grid_size", max_grid_size);
  pp.query("init_cond", ic);

  // plot parameters
  pp.query("nsteps", nsteps);
  pp.query("plot_int", plot_int);
  pp.query("n_checkpoint", n_checkpoint);
  pp.query("start_step", start_step);

  // model parameters
  pp.query("kappa", kappa);
  pp.query("lambda", chi);
  pp.query("T", T);
  pp.query("temperature", temperature);
  pp.query("gamma", Gamma);

  // set up Box and Geomtry
  IntVect dom_lo(0, 0, 0);
  IntVect dom_hi(nx-1, nx-1, nx-1);
  Array<int,3> periodicity({1,1,1});
  Box domain(dom_lo, dom_hi);
  RealBox real_box({0.,0.,0.},{1.,1.,1.});
  Geometry geom(domain, real_box, CoordSys::cartesian, periodicity);
  BoxArray ba(domain);
  // split BoxArray into chunks no larger than "max_grid_size" along a direction
  ba.maxSize(max_grid_size);
  DistributionMapping dm(ba);
  // need two halo layers for gradients
  int nghost = 2;
  Real time = start_step;

  // set up MultiFabs
  MultiFab fold(ba, dm, nvel, nghost);
  MultiFab fnew(ba, dm, nvel, nghost);
  MultiFab gold(ba, dm, nvel, nghost);
  MultiFab gnew(ba, dm, nvel, nghost);
  MultiFab hydrovs(ba, dm, 2*nvel, nghost);
  MultiFab noise(ba, dm, 2*nvel, nghost);
  MultiFab ref_params(ba, dm, 2, nghost); //reference rho and C for each point of the lattice

  StructFact structFact;
  int nStructVars = 5;
  const Vector<std::string> var_names = VariableNames(nStructVars);
  Vector<Real> var_scaling(nStructVars*(nStructVars+1)/2);
  for (int i=0; i<var_scaling.size(); ++i) {
    if (temperature>0) var_scaling[i] = temperature; else var_scaling[i] = 1.;
  }

  if (ic != 10){StructFact structFact(ba, dm, var_names, var_scaling);}

  Print() << "Struct fact object initialized\n";
  // INITIALIZE
  switch(ic){
    case 0:
      LBM_init_mixture(fold, gold, hydrovs);
      break;
    case 1:
      LBM_init_flat_interface(geom, fold, gold, hydrovs);
      break;
    case 2:
      LBM_init_droplet(0.3, geom, fold, gold, hydrovs);
      break;
    case 10:
      checkpointRestart(start_step, hydrovs, fold, gold, ba, dm);
      --start_step;
      // structFact.ReadCheckPoint(step_start,time,"chk_SF",ba,dmap);
      if (temperature != 0){structFact.ReadCheckPoint(start_step, time, "chk_SF_", ba, dm);--start_step;}
      break;
  }
  if (ic != 10){
    WriteCheckPoint(start_step, hydrovs);
    if (temperature > 0){structFact.WriteCheckPoint(start_step,time,"chk_SF_");}
    }

  hydrovs.Copy(ref_params, hydrovs, 0, 0, 2, nghost);
  // Write a plotfile of the initial data if plot_int > 0
  if (plot_int > 0) WriteOutput(start_step, hydrovs, geom, "hydro_plt");
  Print() << "LB initialized\n";
  start_step++;

  // TIMESTEP
  for (int step=start_step; step <= nsteps; ++step) {
    LBM_timestep(geom, fold, gold, fnew, gnew, hydrovs, noise, ref_params);
    Print() << "LBM time integration complete" << "\n";
    if (plot_int > 0 && step%plot_int ==0) {
      WriteOutput(step, hydrovs, geom, "hydro_plt");
      Print() << "hydrovars plotted" << "\n";
      if (step%n_checkpoint == 0){WriteCheckPoint(step, hydrovs);}
      if (temperature > 0){
        // WriteOutput(step, noise, geom, "xi_plt"); 
        structFact.WritePlotFile(step, time, geom, "SF_plt", 0); Print() << "SF plotted" << "\n"; // remove 0 if k = 0 point is to be zeroed in output
        StructFact structFact(ba, dm, var_names, var_scaling); Print() << "SF remade" << "\n";;
        if (step%n_checkpoint == 0){structFact.WriteCheckPoint(step,time,"chk_SF_");}
        }
    }
    if (temperature > 0){structFact.FortStructure(hydrovs, geom, 0);Print() << "SF added" << "\n";}
    Print() << "LB step " << step << "\n";
    time++;
    }
  
  // Call the timer again and compute the maximum difference between the start time 
  // and stop time over all processors
  Real stop_time = ParallelDescriptor::second() - strt_time;
  ParallelDescriptor::ReduceRealMax(stop_time);
  amrex::Print() << "Run time = " << stop_time << std::endl;
  
}