#!/bin/bash

#setup from the perspective of the folder it is run in
amrex_exec="../../../spatially_independent/main3d.gnu.MPI.ex"

# Defines the radius of the droplet as a proportion of system size. 0.5 is maximum
droplet_radii=("0.2" "0.25" "0.3" "0.35" "0.4")

# Define the base folder name
base_folder="R"

# Define property to be looped over and edited in input file
search_pattern1="droplet_radius_prop = 0.3"

# Loop through each radius in R
for droplet_radius in "${droplet_radii[@]}"; do
  folder="${base_folder}_${droplet_radius}"
  echo "Processing folder: $folder"
  
  # Check if the folder exists. Creates folder if it does not exist
  if [ ! -d "$folder" ]; then
    echo "$folder does not exist. Creating $folder"
    mkdir -p "$folder"
  fi

  cp inputs* $folder

  # Enters folder and executes commands before returning to parent directory
  cd "$folder"
  # Checks if run is complete and executes a run with appropriate modifications if it has not
  if [ -e "chk_hydro_0000100000" ]; then
    echo "Run complete"
    continue
  else
    echo "Executing commands in $folder"
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
    replacement_text="droplet_radius_prop = $droplet_radius"
    sed -i "s/$search_pattern1/$replacement_text/g" inputs_young-laplace # edits input file with appropriate setting

    mpirun -n 8 $amrex_exec inputs_young-laplace > young_laplace.txt
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
  fi
  cd ..

done
