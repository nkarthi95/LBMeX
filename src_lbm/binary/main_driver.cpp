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

/* Default initial conditions all in LBM_binary.H */

// time stepping
int nsteps = 10;
int checkpoint_int = nsteps;
int start_time = 0;
// fluctuation_start in LBM_binary.H

// output parameters
int nvars_dump = 2;
int dump_start = 0;
int hydrovars_int = nsteps;
int SF_int = 0;
int distribution_int = 0;
int analysis_int = 0;
std::string analysis_filePath = "droplet_analysis.csv";
std::vector<std::string> col_headers;
int page_hold = 1000;

// default time stepping parameters

/* binary fluid parameter defaults in LBM_binary.H */

/* fluctuations parameter defaults in LBM_fluctuations.H */

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
  
  // Initial condition setup
  pp.query("init_cond", init_cond);
  pp.query("C1", C1); // Proportion of C1 in system. Only used for init_cond = 0 (mixed system)
  // Droplet properties when using init_cond = 2(droplet), 3(cylinder)
  pp.query("droplet_radius_prop", droplet_radius_prop);
  pp.query("phi_in", phi_in);
  pp.query("phi_out", phi_out);
  pp.query("rho_in", rho_in);
  pp.query("rho_out", rho_out);

  /* time stepping */
  pp.query("nsteps", nsteps);
  pp.query("n_checkpoint", checkpoint_int);
  pp.query("restore_string", start_time);
  pp.query("fluctuation_start", fluctuation_start);

  /* output parameters */
  pp.query("nvars_dump", nvars_dump);
  pp.query("dump_start", dump_start);
  pp.query("hydrovars_int", hydrovars_int);
  pp.query("SF_int", SF_int);
  pp.query("distribution_int", distribution_int);
  pp.query("analysis_int", analysis_int);
  pp.query("analysis_file", analysis_filePath);
  pp.query("page_hold", page_hold);

  /* binary fluid parameters */
  pp.query("chi", chi);
  pp.query("T", T);
  pp.query("kappa", kappa);
  pp.query("gamma", Gamma);

  /* fluctuations parameters*/
  pp.query("temperature", temperature);
  // pp.dumpTable()
  // pp.dumpTable(amrex::OutStream(), true);
}

inline void WriteHydrovars(int step, const Geometry& geom, const MultiFab& hydrovs){
  const int nvars = nvars_dump;
  const Vector<std::string> var_names = hydrovars_names(nvars);
  const std::string& pltfile = amrex::Concatenate("hydro_plt",step,9);
  WriteSingleLevelPlotfile(pltfile, hydrovs, var_names, geom, Real(step), step);
  // if (dump_hydro) {WriteSingleLevelPlotfile(pltfile, hydrovs, var_names, geom, Real(step), step);}
}

inline void WriteSF(int step, StructFact& structFact){
  // const int nvars = nvars_dump;
  // const Vector<std::string> var_names = hydrovars_names(nvars);
  // const std::string& pltfile = amrex::Concatenate("hydro_plt",step,9);
  // if (dump_hydro) {WriteSingleLevelPlotfile(pltfile, hydrovs, var_names, geom, Real(step), step);}
  const int zero_avg = 1;
  structFact.WritePlotFile(step, static_cast<Real>(step), "SF_plt", zero_avg);
  // if (dump_SF) {structFact.WritePlotFile(step, static_cast<Real>(step), "SF_plt", zero_avg);}
}

inline void WriteDists(int step,
                       const Geometry& geom,
                       const MultiFab& f,
                       const MultiFab& g) {
    std::string pltfile;
    Vector<std::string> var_names(nvel);

    // f distribution
    pltfile = amrex::Concatenate("dist_f",step,9);
    for (int i = 0; i < nvel; i++) {
      var_names[i] = "f" + std::to_string(i);
    }
    WriteSingleLevelPlotfile(pltfile, f, var_names, geom, Real(step), step);
    
    // g distribution
    pltfile = amrex::Concatenate("dist_g",step,9);
    for (int i = 0; i < nvel; i++) {
      var_names[i] = "g" + std::to_string(i);
    }
    WriteSingleLevelPlotfile(pltfile, g, var_names, geom, Real(step), step);
}

#ifndef AMREX_USE_CUDA
inline void write_csv(PrintToFile& AMReX_printObj, Array1D<Real, 0, 7> data_to_append){
  for (int i = 0; i < data_to_append.len(); i++){
    AMReX_printObj << data_to_append(i) << ",";
  }
  AMReX_printObj << std::endl;
}

inline void write_csv(PrintToFile& AMReX_printObj, std::vector<std::string> data_to_append){
  for (int i = 0; i < data_to_append.size(); i++){
    AMReX_printObj << data_to_append[i] << ",";
  }
  AMReX_printObj << std::endl;
}

inline void droplet_analysis(PrintToFile& AMReX_printObj, const int step, const MultiFab& hydrovs){
    BL_PROFILE_VAR("droplet_analysis()",droplet_analysis);
    Array1D<Real, 0, 7> droplet_data; //"Timestep, Radius, com_x, com_y, com_z, dx, dy, dz\n" 

    MultiFab droplet = binarize_droplet(calculate_C1(hydrovs), 0, 0.5);
    
    Real R = droplet_radius(droplet);
    // Real R = droplet_radius_profile_fit(droplet);
    // Print() << R << "\n";
    GpuArray<Real, 3> com = center_of_mass(droplet);
    GpuArray<Real, 3> dr = axial_radii(droplet);
    droplet_data(0) = step;
    droplet_data(1) = R;
    droplet_data(2) = com[0]; droplet_data(3) = com[1]; droplet_data(4) = com[2];
    droplet_data(5) = dr[0]; droplet_data(6) = dr[1]; droplet_data(7) = dr[2];
    write_csv(AMReX_printObj, droplet_data);
    // Print() << "step:" << step << ", R:" << R << ", com_x:" << com[0] << ", com_y:" << com[1] << ", com_z:" << com[2] << "\n";
}

// inline void add_to_string(std::string& output_str, std::vector<std::string> data_to_append){
//   for (int i = 0; i < data_to_append.size(); i++){
//     output_str += data_to_append[i] + ",";
//   }
//   output_str += "\n";
//   return output_str;
// }

// inline void add_to_string(std::string& output_str, Array1D<Real, 0, 7> data_to_append){
//   for (int i = 0; i < data_to_append.len(); i++){
//     output_str += static_cast<std::string>(data_to_append(i)) + ",";
//   }
//   output_str += "\n";
//   return output_str;
// }

// inline void droplet_analysis(std::string& output_str, const int step, const MultiFab& hydrovs){
//     BL_PROFILE_VAR("droplet_analysis()",droplet_analysis);
//     Array1D<Real, 0, 7> droplet_data; //"Timestep, Radius, com_x, com_y, com_z, dx, dy, dz\n" 

//     MultiFab droplet = binarize_droplet(calculate_C1(hydrovs), 0, 0.5);
    
//     Real R = droplet_radius(droplet);
//     // Real R = droplet_radius_profile_fit(droplet);
//     // Print() << R << "\n";
//     GpuArray<Real, 3> com = center_of_mass(droplet);
//     GpuArray<Real, 3> dr = axial_radii(droplet);
//     droplet_data(0) = step;
//     droplet_data(1) = R;
//     droplet_data(2) = com[0]; droplet_data(3) = com[1]; droplet_data(4) = com[2];
//     droplet_data(5) = dr[0]; droplet_data(6) = dr[1]; droplet_data(7) = dr[2];
//     write_csv(output_str, droplet_data);
//     // Print() << "step:" << step << ", R:" << R << ", com_x:" << com[0] << ", com_y:" << com[1] << ", com_z:" << com[2] << "\n";
// }
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
  #if FIXREF
  MultiFab reference(ba, dm, 2, nghost);
  reference.setVal(1.0, 0, 1, nghost); // set reference density to a constant value rho = 1.0
  reference.setVal(0., 1, 1, nghost); // set reference order parameter to a constant value phi = 0.0
  #endif

  // set up StructFact
  int nStructVars = 38;
  const Vector<std::string> var_names = hydrovars_names(nStructVars);
  Vector<int> pairA(nStructVars); std::iota(pairA.begin(), pairA.end(), 0); // idxs = [0, 1, ..., N-1]
  Vector<int> pairB(nStructVars); std::iota(pairB.begin(), pairB.end(), 0); // idxs = [0, 1, ..., N-1]
  const Vector<Real> var_scaling(pairA.size(), 1.0);
  StructFact structFact(ba, dm, var_names, var_scaling, pairA, pairB);

  static std::unique_ptr<PrintToFile> droplet_properties_output;
  if (!droplet_properties_output) {
    droplet_properties_output = std::make_unique<PrintToFile>(analysis_filePath, 0);
    droplet_properties_output->SetPrecision(14);
  }

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
      write_csv(*droplet_properties_output, col_headers);
      droplet_analysis(*droplet_properties_output, start_time, hydrovs);
      #endif
      break;
    case 3:
      LBM_init_cylinder(droplet_radius_prop, geom, fold, gold, hydrovs);
      start_time = 0;
      break;
    case 7:
      checkpointRestart(start_time, hydrovs, fold, gold, ba, dm); start_time--; //start_time is increased by 1 when checkpoint restart is done.
      #if FIXREF
      reference.Copy(reference, hydrovs, 0, 0, 2, nghost); // need 2 ghost cells
      #endif
      break;
    default:
      Print() << "Initial condition specified does not exist. Please enter a difference choice" << std::endl;
  }
  // if (hydrovars_int > 0) WriteOutput(start_time, geom, hydrovs, structFact);
  if (hydrovars_int > 0) WriteHydrovars(start_time, geom, hydrovs);
  if (checkpoint_int > 0) WriteCheckPoint(start_time, hydrovs); start_time++;
  Print() << "LB initialized lattice " << domain <<"\n" << ba << dm << std::endl;

  #ifndef AMREX_USE_CUDA
    #if AMREX_DEBUG
      unit_tests(geom, hydrovs);
    #endif
  #endif

  // TIMESTEP
  for (int step=start_time; step <= nsteps; ++step) {
    #if FIXREF
      LBM_timestep(geom, fold, gold, fnew, gnew, hydrovs, noise, reference, step);
    #else
      LBM_timestep(geom, fold, gold, fnew, gnew, hydrovs, noise, step);
    #endif
    structFact.FortStructure(hydrovs);

    if (hydrovars_int > 0 && step%hydrovars_int == 0 && step >= dump_start){
      WriteHydrovars(step, geom, hydrovs);
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
    #ifndef AMREX_USE_CUDA
    // if (analysis_int > 0 && step%analysis_int == 0){droplet_analysis(droplet_properties_output, step, hydrovs);}
    if (analysis_int > 0 && step%analysis_int == 0){droplet_analysis(*droplet_properties_output, step, hydrovs);}
    if (analysis_int > 0 && step%page_hold == 0) {
      droplet_properties_output.reset();
      droplet_properties_output = std::make_unique<PrintToFile>(analysis_filePath, 0);
      droplet_properties_output->SetPrecision(14);
      }
    #endif
  }

  Print() << "LB completed " << nsteps << " time steps" << std::endl;

  // Call the timer again and compute the maximum difference between the start time 
  // and stop time over all processors
  Real stop_time = ParallelDescriptor::second() - strt_time;
  ParallelDescriptor::ReduceRealMax(stop_time);
  amrex::Print() << "Run time = " << stop_time << " s (" << domain.numPts()*nsteps/stop_time << " LUP/s)" << std::endl;
}
