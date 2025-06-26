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
#include "LBM_analysis.H"
#include "LBM_tests.H"
#endif

// default grid parameters
IntVect domain_size(16);
IntVect max_box_size(32);

// default time stepping parameters
int nsteps = 10;
int checkpoint_int = nsteps;
int start_time = 0;
int dump_SF = 0;
int dump_hydro = 1;
int dump_start = 0;
std::string analysis_filePath = "droplet_analysis.csv";
int analysis_int = 10;
std::vector<std::string> col_headers;

inline void ReadInput() {
  ParmParse pp;

  /* grid parameters */
  pp.query("nx", domain_size[0]);
  domain_size[2] = domain_size[1] = domain_size[0]; // default to cubic box
  pp.query("ny", domain_size[1]);
  pp.query("nz", domain_size[2]);

  pp.query("max_grid_size_x", max_box_size[0]);
  max_box_size[2] = max_box_size[1] = max_box_size[0]; // default to same maxSize in all directions
  pp.query("max_grid_size_y", max_box_size[1]);
  pp.query("max_grid_size_z", max_box_size[2]);
  pp.query("init_cond", init_cond);
  pp.query("droplet_radius_prop", droplet_radius_prop);
  pp.query("C1", C1);

  /* time stepping and output parameters */
  pp.query("nsteps", nsteps);
  pp.query("plot_int", plot_int);
  pp.query("n_checkpoint", checkpoint_int);
  pp.query("restore_string", start_time);
  pp.query("fluctuation_start", fluctuation_start);
  pp.query("dump_start", dump_start);
  pp.query("dump_SF", dump_SF);
  pp.query("dump_hydro", dump_hydro);
  pp.query("analysis_int", analysis_int);
  pp.query("analysis_file", analysis_filePath);

  /* binary fluid parameters */
  pp.query("chi", chi);
  pp.query("T", T);
  pp.query("kappa", kappa);
  pp.query("gamma", Gamma);

  /* noise parameters */
  pp.query("temperature", temperature);
  // pp.dumpTable()
  pp.dumpTable(amrex::OutStream(), true);
}

inline void WriteOutput(int step,
      const Geometry& geom,
			const MultiFab& hydrovs,
      StructFact& structFact) {
  // set up variable names for output
  const int zero_avg = 1;
  const int nvars = 14;
  const Vector<std::string> var_names = hydrovars_names(nvars);
  const std::string& pltfile = amrex::Concatenate("hydro_plt",step,9);
  if (dump_hydro) {WriteSingleLevelPlotfile(pltfile, hydrovs, var_names, geom, Real(step), step);}
  if (dump_SF) {structFact.WritePlotFile(step, static_cast<Real>(step), "SF_plt", zero_avg);}
}

#ifndef AMREX_USE_CUDA
inline void write_csv(const std::string analysis_filePath, Array1D<Real, 0, 8> data_to_append){
  for (int i = 0; i < 8; i++){
    PrintToFile(analysis_filePath, 0) << data_to_append(i) << ",";
  }
  PrintToFile(analysis_filePath, 0) << "\n";
}

inline void write_csv(const std::string analysis_filePath, std::vector<std::string> data_to_append){
  for (int i = 0; i < 8; i++){
    PrintToFile(analysis_filePath, 0) << data_to_append[i] << ",";
  }
  PrintToFile(analysis_filePath, 0) << "\n";
}

inline void droplet_analysis(const std::string analysis_filePath, const int step, const MultiFab& hydrovs, MultiFab& droplet){
    Array1D<Real, 0, 8> droplet_data; //"Timestep, Radius, com_x, com_y, com_z, dx, dy, dz\n" 
    MultiFab::Copy(droplet, hydrovs, 1, 0, 1, 0);
    // droplet.ParallelCopy(hydrovs, 1, 0, 1);
    MultiFab::Add(droplet, hydrovs, 0, 0, 1, 0);
    droplet.mult(0.5, 0);

    droplet = binarize_droplet(droplet, 0, 0.5);
    
    Real R = droplet_radius(droplet);
    // Real R = droplet_radius_profile_fit(droplet);
    // Print() << R << "\n";
    GpuArray<Real, 3> com = center_of_mass(droplet);
    GpuArray<Real, 3> dr = axial_radii(droplet);
    droplet_data(0) = step;
    droplet_data(1) = R;
    droplet_data(2) = com[0]; droplet_data(3) = com[1]; droplet_data(4) = com[2];
    droplet_data(5) = dr[0]; droplet_data(6) = dr[1]; droplet_data(7) = dr[2];
    write_csv(analysis_filePath, droplet_data);
}
#endif

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
  int nghost = 2; // need two halo layers for gradients

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
  MultiFab test_noise(ba, dm, 2*nvel, nghost);
  MultiFab reference(ba, dm, 2, nghost);
  reference.setVal(1.0, 0, 1, nghost); // set reference density to a constant value rho = 1.0
  reference.setVal(0., 1, 1, nghost); // set reference order parameter to a constant value phi = 0.0

  // droplet analysis
  MultiFab droplet(ba, dm, 1, 0);

  // set up StructFact
  int nStructVars = 14;
  const Vector<std::string> var_names = hydrovars_names(nStructVars);
  // const Vector<int> pairA = { 0, 1, 2, 3, 4 };
  // const Vector<int> pairB = { 0, 1, 2, 3, 4 };
  // const Vector<Real> var_scaling = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  Vector<int> pairA(nStructVars); std::iota(pairA.begin(), pairA.end(), 0); // idxs = [0, 1, ..., N-1]
  Vector<int> pairB(nStructVars); std::iota(pairB.begin(), pairB.end(), 0); // idxs = [0, 1, ..., N-1]
  const Vector<Real> var_scaling(nStructVars, 1.0);
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
      #ifndef AMREX_USE_CUDA
      col_headers = {"Timestep", "Radius", "cx", "cy", "cz", "dx", "dy", "dz"};
      write_csv(analysis_filePath, col_headers);
      droplet_analysis(analysis_filePath, start_time, hydrovs, droplet);
      #endif
      break;
    case 7:
      checkpointRestart(start_time, hydrovs, fold, gold, ba, dm); start_time--; //start_time is increased by 1 when checkpoint restart is done.
      reference.Copy(reference, hydrovs, 0, 0, 2, 2); // need 2 ghost cells
      break;
    default:
      Print() << "Initial condition specified does not exist. Please enter a difference choice" << std::endl;
  }
  if (plot_int > 0) WriteOutput(start_time, geom, hydrovs, structFact);
  if (checkpoint_int > 0) WriteCheckPoint(start_time, hydrovs); start_time++;
  Print() << "LB initialized lattice " << domain <<"\n" << ba << dm << std::endl;

  #ifndef AMREX_USE_CUDA
    #if AMREX_DEBUG
      unit_tests(geom, hydrovs);
    #endif
  #endif

  // TIMESTEP
  for (int step=start_time; step <= nsteps; ++step) {
    LBM_timestep(geom, fold, gold, fnew, gnew, hydrovs, noise, reference, step);
    structFact.FortStructure(hydrovs);
    if (plot_int > 0 && step%plot_int == 0 && step >= dump_start) {
      WriteOutput(step, geom, hydrovs, structFact);
      Print() << "LB step " << step << std::endl;
    }
    if (checkpoint_int > 0 && step%checkpoint_int ==0){
      WriteCheckPoint(step, hydrovs);
    }
    #ifndef AMREX_USE_CUDA
    if (analysis_int > 0 && step%analysis_int == 0){droplet_analysis(analysis_filePath, step, hydrovs, droplet);}
    #endif
  }

  Print() << "LB completed " << nsteps << " time steps" << std::endl;

  // Call the timer again and compute the maximum difference between the start time 
  // and stop time over all processors
  Real stop_time = ParallelDescriptor::second() - strt_time;
  ParallelDescriptor::ReduceRealMax(stop_time);
  amrex::Print() << "Run time = " << stop_time << " s (" << domain.numPts()*nsteps/stop_time << " LUP/s)" << std::endl;
}