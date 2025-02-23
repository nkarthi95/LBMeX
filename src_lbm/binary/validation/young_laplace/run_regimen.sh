#!/bin/bash

#setup from the perspective of the folder it is run in
amrex_exec="../../../../../spatially_independent/main3d.gnu.MPI.ex"

# Defines the radius of the droplet as a proportion of system size. 0.5 is maximum
droplet_radii=("0.30" "0.35" "0.40")
chi_s=("2.15" "2.25" "2.5" "2.75" "3" "3.25")# "3.5" "3.75" "4" "4.25" "4.5" "4.75" "5")
T="0.2"
kappa="0.03"

# Define property to be looped over and edited in input file
search_pattern1="droplet_radius_prop = 0.3"
search_pattern2="chi = 0.45"
search_pattern3="kappa = 0.03"

for chi in "${chi_s[@]}"; do

base_folder="chi_${chi}/kappa_${kappa}"

if [ ! -d "$base_folder" ]; then
  echo "$base_folder does not exist. Creating $base_folder"
  mkdir -p "$base_folder"
fi

cd $base_folder
echo ${PWD}
# Loop through each radius in R
for droplet_radius in "${droplet_radii[@]}"; do
  folder="R_${droplet_radius}"
  echo "Processing folder: $folder"
  
  # Check if the folder exists. Creates folder if it does not exist
  if [ ! -d "$folder" ]; then
    echo "$folder does not exist. Creating $folder"
    mkdir -p "$folder"
  fi

  cp ../../inputs_young-laplace $folder

  # Enters folder and executes commands before returning to parent directory
  cd $folder
  # Checks if run is complete and executes a run with appropriate modifications if it has not
  if [ -e "checkpoint_000075000" ]; then
    echo "Run complete"
    cd ..
    continue
  else
    echo "Executing commands in $folder"
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
    replacement_text="droplet_radius_prop = $droplet_radius"
    sed -i "s/$search_pattern1/$replacement_text/g" inputs_young-laplace # edits input file with appropriate setting
    lambda=$(echo "$chi * $T" | bc -l)
    replacement_text="chi = ${lambda}"
    sed -i "s/$search_pattern2/$replacement_text/g" inputs_young-laplace # edits input file with appropriate setting
    replacement_text="kappa = ${kappa}"
    sed -i "s/$search_pattern2/$replacement_text/g" inputs_young-laplace # edits input file with appropriate setting

    mpirun -n 8 $amrex_exec inputs_young-laplace > young_laplace.txt
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
  fi
  cd ..
done
cd ../..
done
