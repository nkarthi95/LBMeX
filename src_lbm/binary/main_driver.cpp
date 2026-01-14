#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <StructFact.H>

using namespace amrex;

#include "LBM_binary.H"
#include "LBM_fluctuations.H"
#include "LBM_IO.H"

#ifndef AMREX_USE_CUDA
#include "LBM_tests.H"
#endif

/// Default grid parameters.
// domain_size : Number of cells in each spatial direction.
/// max_box_size : Maximum box size for AMReX grid decomposition.
IntVect domain_size(16);
IntVect max_box_size(32);

/// Default time stepping parameters.
/// nsteps          : Total number of timesteps to run.
/// checkpoint_int  : Interval (in timesteps) between checkpoint outputs.
/// start_time      : Restart time index (used when restoring from checkpoint).
int nsteps = 10;
int checkpoint_int = nsteps;
int start_time = 0;

/// Default output parameters.
/// nvars_dump        : Number of variables written to plotfiles.
/// dump_start        : Timestep at which plotfile output begins.
/// hydrovars_int     : Interval for writing hydrodynamic variables.
/// SF_int            : Interval for writing structure factor data (0 disables).
/// distribution_int  : Interval for writing distribution functions (0 disables).
int nvars_dump = 2;
int dump_start = 0;
int hydrovars_int = nsteps;
int SF_int = 0;
int distribution_int = 0;



// fluctuation_start in LBM_binary.H
/* binary fluid parameter defaults in LBM_binary.H */
/* Default initial conditions all in LBM_binary.H */
/* fluctuations parameter defaults in LBM_fluctuations.H */

/**
 * \brief Read runtime parameters from the input file.
 *
 * This function uses AMReX ParmParse to read simulation parameters
 * from the input file and override default values. Parameters are
 * grouped into grid setup, initial conditions, time stepping,
 * output control, binary fluid physics, and fluctuation settings.
 *
 * All parameters are optional; if a parameter is not found in the
 * input file, its default value remains unchanged.
 */
inline void ReadInput()
{
  ParmParse pp;

  /** \name Grid parameters */
  ///@{

  /// nx : Number of cells in the x-direction.
  /// ny : Number of cells in the y-direction (defaults to nx).
  /// nz : Number of cells in the z-direction (defaults to nx).
  pp.query("nx", domain_size[0]);
  domain_size[2] = domain_size[1] = domain_size[0]; // default cubic domain
  pp.query("ny", domain_size[1]);
  pp.query("nz", domain_size[2]);

  /// max_grid_size_x : Maximum box size in the x-direction.
  /// max_grid_size_y : Maximum box size in the y-direction (defaults to x).
  /// max_grid_size_z : Maximum box size in the z-direction (defaults to x).
  pp.query("max_grid_size_x", max_box_size[0]);
  max_box_size[2] = max_box_size[1] = max_box_size[0]; // uniform max box size
  pp.query("max_grid_size_y", max_box_size[1]);
  pp.query("max_grid_size_z", max_box_size[2]);

  ///@}

  /** \name Initial condition parameters */
  ///@{

  /// init_cond : Initial condition selector.
  ///   0 = Homogeneous mixture
  ///   1 = Flat interface in the yz-plane
  ///   2 = Spherical droplet of C1 at domain center
  ///   3 = Cylindrical C1 region aligned with z-axis
  pp.query("init_cond", init_cond);

  /// droplet_radius_prop : Droplet radius as a fraction of domain size.
  pp.query("droplet_radius_prop", droplet_radius_prop);

  /// C1_in  : Concentration of component C1 inside droplet/region.
  /// C1_out : Concentration of component C1 outside droplet/region.
  /// C2_in  : Concentration of component C2 inside droplet/region.
  /// C2_out : Concentration of component C2 outside droplet/region.
  pp.query("C1_in", C1_in);
  pp.query("C1_out", C1_out);
  pp.query("C2_in", C2_in);
  pp.query("C2_out", C2_out);

  /// tau_r : Relaxation time for order parameter dynamics.
  /// tau_p : Relaxation time for momentum dynamics.
  pp.query("tau_r", tau_r);
  pp.query("tau_p", tau_p);

  ///@}

  /** \name Time stepping parameters */
  ///@{

  /// nsteps            : Total number of timesteps.
  /// n_checkpoint      : Interval between checkpoint writes.
  /// restore_string    : Restart timestep index.
  /// fluctuation_start : Timestep at which thermal fluctuations are enabled.
  pp.query("nsteps", nsteps);
  pp.query("n_checkpoint", checkpoint_int);
  pp.query("restore_string", start_time);
  pp.query("fluctuation_start", fluctuation_start);

  ///@}

  /** \name Output control parameters */
  ///@{

  /// nvars_dump        : Number of fields written to plotfiles.
  /// dump_start        : Timestep to begin plotfile output.
  /// hydrovars_int     : Interval for hydrodynamic variable output.
  /// SF_int            : Interval for structure factor output.
  /// distribution_int  : Interval for distribution function output.
  pp.query("nvars_dump", nvars_dump);
  pp.query("dump_start", dump_start);
  pp.query("hydrovars_int", hydrovars_int);
  pp.query("SF_int", SF_int);
  pp.query("distribution_int", distribution_int);

  ///@}

  /** \name Binary fluid parameters */
  ///@{

  /// chi   : Flory–Huggins interaction parameter.
  /// T     : Thermodynamic temperature.
  /// kappa : Gradient energy coefficient.
  /// gamma : Mobility or coupling parameter.
  pp.query("chi", chi);
  pp.query("T", T);
  pp.query("kappa", kappa);
  pp.query("gamma", Gamma);

  ///@}

  /** \name Fluctuation parameters */
  ///@{

  /// temperature : Noise temperature used in fluctuating hydrodynamics.
  pp.query("temperature", temperature);

  ///@}
}

void main_driver(const char* argv) {

  // store the current time so we can later compute total run time.
  Real strt_time = ParallelDescriptor::second();

  // read input parameters
  ReadInput();

  // set up Box and Geomtry
  RealBox real_box({0.,0.,0.},{1.,1.,1.});
  IntVect dom_lo(0, 0, 0);
  IntVect dom_hi(domain_size-1);
  Array<int,3> periodicity({1,1,1});

  Box domain(dom_lo, dom_hi);
  Geometry geom(domain, real_box, CoordSys::cartesian, periodicity);
  BoxArray ba(domain);
  ba.maxSize(max_box_size); // chop domain into boxes
  DistributionMapping dm(ba);

  // set up MultiFabs
  MultiFab fold(ba, dm, nvel, nghost);
  MultiFab fnew(ba, dm, nvel, nghost);
  MultiFab gold(ba, dm, nvel, nghost);
  MultiFab gnew(ba, dm, nvel, nghost);
  MultiFab hydrovs(ba, dm, 2*nvel, nghost);
  MultiFab noise(ba, dm, 2*nvel, nghost);

  // set up StructFactF
  const Vector<std::string> var_names = hydrovars_names(ndof);
  Vector<int> pairA(ndof); std::iota(pairA.begin(), pairA.end(), 0); pairA.emplace_back(1); // idxs = [0, 1, ..., N-1, 1]
  Vector<int> pairB(ndof); std::iota(pairB.begin(), pairB.end(), 0); pairB.emplace_back(0); // idxs = [0, 1, ..., N-1, 0]
  const Vector<Real> var_scaling(pairA.size(), 1.0);
  StructFact structFact(ba, dm, var_names, var_scaling, pairA, pairB);

  // INITIALIZE
  switch(init_cond){
    case 0:
      LBM_init_mixture(fold, gold, hydrovs);
      start_time = 0;
      break;
    case 1:
      LBM_init_flat_interface(geom, fold, gold, hydrovs);
      start_time = 0;
      break;
    case 2:
      LBM_init_droplet(droplet_radius_prop, geom, fold, gold, hydrovs);
      start_time = 0;
      break;
    case 3:
      LBM_init_cylinder(droplet_radius_prop, geom, fold, gold, hydrovs);
      start_time = 0;
      break;
    case 7:
      checkpointRestart(start_time, hydrovs, fold, gold, ba, dm); start_time--; //start_time is increased by 1 when checkpoint restart is done.
      break;
    default:
      Print() << "Initial condition specified does not exist. Please enter a difference choice" << std::endl;
  }
  // if (hydrovars_int > 0) WriteOutput(start_time, geom, hydrovs, structFact);
  if (hydrovars_int > 0) WriteHydrovars(start_time, geom, hydrovs, nvars_dump);
  if (checkpoint_int > 0) WriteCheckPoint(start_time, hydrovs); start_time++;
  Print() << "LB initialized lattice " << domain <<"\n" << ba << dm << std::endl;

  #ifndef AMREX_USE_CUDA
    #if AMREX_DEBUG
      unit_tests(geom, hydrovs);
    #endif
  #endif

  // TIMESTEP
  for (int step=start_time; step <= nsteps; ++step) {
      LBM_timestep(geom, fold, gold, fnew, gnew, hydrovs, noise, step);
    structFact.FortStructure(hydrovs);

    if (hydrovars_int > 0 && step%hydrovars_int == 0 && step >= dump_start){
      WriteHydrovars(step, geom, hydrovs, nvars_dump);
      Print() << "LB step " << step << std::endl;
    }
    
    if (SF_int > 0 && step%SF_int == 0 && step >= dump_start){
      WriteSF(step, structFact);
    }

    if(distribution_int > 0 && step%distribution_int == 0 && step >= dump_start){
        WriteDists(step, geom, fold, gold);
    }

    if (checkpoint_int > 0 && step%checkpoint_int ==0){
      WriteCheckPoint(step, hydrovs);
    }
  }

  Print() << "LB completed " << nsteps << " time steps" << std::endl;

  // Call the timer again and compute the maximum difference between the start time 
  // and stop time over all processors
  Real stop_time = ParallelDescriptor::second() - strt_time;
  ParallelDescriptor::ReduceRealMax(stop_time);
  amrex::Print() << "Run time = " << stop_time << " s (" << domain.numPts()*nsteps/stop_time << " LUP/s)" << std::endl;
}
