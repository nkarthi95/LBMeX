#!/bin/bash

#setup from the perspective of the folder it is run in
amrex_exec="../../../../../spatially_independent/main3d.gnu.MPI.ex"
#amrex_exec="../../../../../main3d.gnu.MPI.ex"

# Defines the radius of the droplet as a proportion of system size. 0.5 is maximum
chi_s=("2.05" "2.1" "2.15" "2.2" "2.25" "2.3" "2.35" "2.45" "2.5") # "3.5" "3.75" "4" "4.25" "4.5" "4.75" "5")
T="0.3"
kappa="0.01"

# Define property to be looped over and edited in input file
#search_pattern1="droplet_radius_prop = 0.3"
search_pattern2="chi = 0.45"
search_pattern3="kappa = 0.01"
search_pattern4="T = 0.2"

for chi in "${chi_s[@]}"; do

folder="T_${T}/chi_${chi}/kappa_${kappa}"

if [ ! -d "$folder" ]; then
  echo "$folder does not exist. Creating $folder"
  mkdir -p "$folder"
fi

  cd $folder
  echo ${PWD}

  cp ../../../inputs_interface_width .

  # Enters folder and executes commands before returning to parent directory
  # Checks if run is complete and executes a run with appropriate modifications if it has not
  if [ -e "checkpoint_00030000" ]; then
    echo "Run complete"
    cd ../../../
    continue
  else
    echo "Executing commands in $folder"
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
    lambda=$(echo "$chi * $T" | bc -l)
    replacement_text="chi = ${lambda}"
    sed -i "s/$search_pattern2/$replacement_text/g" inputs_interface_width # edits input file with appropriate setting
    replacement_text="kappa = ${kappa}"
    sed -i "s/$search_pattern3/$replacement_text/g" inputs_interface_width # edits input file with appropriate setting
    replacement_text="T = ${T}"
    sed -i "s/$search_pattern4/$replacement_text/g" inputs_interface_width # edits input file with appropriate setting

    mpirun -n 1 $amrex_exec inputs_interface_width > outputs.txt
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
  fi
cd ../../../
done
