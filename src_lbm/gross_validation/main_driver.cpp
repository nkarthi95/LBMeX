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

inline void WriteDist(int step, 
      const MultiFab& fold,
      const MultiFab& gold, 
      const Vector<std::string>& var_names,
      const Geometry& geom){
      
      const Real time = step;
      std::string pltfile = amrex::Concatenate("f_plt_",step,5);
      WriteSingleLevelPlotfile(pltfile, fold, var_names, geom, time, step);

      pltfile = amrex::Concatenate("g_plt_",step,5);
      WriteSingleLevelPlotfile(pltfile, gold, var_names, geom, time, step);
      }

inline Vector<std::string> VariableNames(const int numVars) {
  // set variable names for output
  Vector<std::string> var_names(numVars);
  std::string name;
  int cnt = 0;
  var_names[cnt++] = "rho";
  // Print() << "Density specified\n";
  // velx, vely, velz
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    name = "u";
    name += (120+d);
    var_names[cnt++] = name;
  }
  // Print() << "velocity specified\n";
  if(cnt >= numVars){return var_names;}
  for (int i=0; i<AMREX_SPACEDIM; ++i) {
    for (int j=i; j<AMREX_SPACEDIM; ++j) {
      name = "p";
      name += (120+i);
      name += (120+j);
      var_names[cnt++] = name;
    }
  }
  // Print() << "pressure specified\n";
  if(cnt >= numVars){return var_names;}
  for (; cnt < numVars;) {
    name = "m";
    name += std::to_string(cnt);
    var_names[cnt++] = name;
  }
  return var_names;
}

inline void WriteOutput(int step,
			const MultiFab& hydrovs,
			const Geometry& geom) {
  // set up variable names for output
  const Vector<std::string> var_names = VariableNames(nvel);
  const std::string& pltfile = amrex::Concatenate("hydro_plt",step,7);
  WriteSingleLevelPlotfile(pltfile, hydrovs, var_names, geom, Real(step), step);
}

void main_driver(const char* argv) {

  // if (!cholesky_test(100)) exit(-1);

  // store the current time so we can later compute total run time.
  Real strt_time = ParallelDescriptor::second();
    
  // default grid parameters
  int nx = 16;
  int max_grid_size = 8;
  int init_cond = 0;
  Real radius = 0.3;

  // default time stepping parameters
  int nsteps = 100;
  int plot_int = 10;

  // fft test input
  int reps = 1;

  // input parameters
  ParmParse pp;
  pp.query("nx", nx);
  pp.query("max_grid_size", max_grid_size);
  pp.query("nsteps", nsteps);
  pp.query("plot_int", plot_int);
  pp.query("kappa", kappa);
  pp.query("rhov", rhov);
  pp.query("rhol", rhol);
  pp.query("B", Beta);
  pp.query("temperature", temperature);
  pp.query("init_cond", init_cond);
  pp.query("droplet_radius", radius);
  // pp.query("reps", reps);
  Print() << "parameters parsed\n";
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

  // Print() << "Geometry generated\n";

  // need two halo layers for gradients
  int nghost = 2;

  // set up MultiFabs
  MultiFab fold(ba, dm, nvel, nghost);
  MultiFab fnew(ba, dm, nvel, nghost);
  MultiFab hydrovs(ba, dm, nvel, nghost);
  MultiFab noise(ba, dm, nvel, nghost);
  // MultiFab ic_setup(ba, dm, nvel, nghost);
  // MultiFab test_noise(ba, dm, 2*nvel, nghost);
  Real rho0;
  if (init_cond == 0) {
    LBM_init_liquid(fold, hydrovs);
    rho0 = rhol;
    }
  else if (init_cond == 1){
    LBM_init_mixture(fold, hydrovs);
    rho0 = 0.5*rhol + 0.5*rhov;
    }
  else if (init_cond == 2){
    LBM_init_flat_interface(geom, fold, hydrovs);
    rho0 = 0.5*rhol + 0.5*rhov;
    }
  else if (init_cond == 3){
    LBM_init_droplet(radius, geom, fold, hydrovs);
    rho0 = 0.5*rhol + 0.5*rhov;
    }
  // Print() << "Data structures generated\n";

  int nStructVars = 10;
  const Vector<std::string> var_names = VariableNames(nStructVars);
  // Print() << "SF names specified\n";
  Vector<Real> var_scaling(nStructVars*(nStructVars+1)/2);
  for (int i=0; i<var_scaling.size(); ++i) {
    if (temperature>0) var_scaling[i] = rho0*temperature; else var_scaling[i] = 1.;
  }
  StructFact structFact(ba, dm, var_names, var_scaling);
  // Print() << "StructFact object generated\n";
  // INITIALIZE
  // Print() << "Initial condition created\n";
  // Print() << "thermodynamic cs2: " << system_cs2 << "\n";

  // Write a plotfile of the initial data if plot_int > 0
  if (plot_int > 0) {WriteOutput(0, hydrovs, geom);}
  // Print() << "LB initialized\n";

  // TIMESTEP
  for (int step=1; step <= nsteps; ++step) {
    LBM_timestep(geom, fold, fnew, hydrovs, noise);
    if (plot_int > 0 && step%plot_int ==0) {WriteOutput(step, hydrovs, geom);}
    if (step > 0.9*nsteps){structFact.FortStructure(hydrovs, geom);}
    Print() << "LB step " << step << "\n";
  }

  // structFact.WritePlotFile(nsteps, nsteps, geom, "SF_plt");
  structFact.WritePlotFile(nsteps, static_cast<Real>(nsteps), geom, "SF_plt");
  // Call the timer again and compute the maximum difference between the start time 
  // and stop time over all processors
  Real stop_time = ParallelDescriptor::second() - strt_time;
  ParallelDescriptor::ReduceRealMax(stop_time);
  amrex::Print() << "Run time = " << stop_time << std::endl;
  
}