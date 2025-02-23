#!/bin/bash

#setup from the perspective of the folder it is run in
amrex_exec="../../../../spatially_independent/main3d.gnu.MPI.ex"

# Defines the radius of the droplet as a proportion of system size. 0.5 is maximum
#("2.15" "2.25" "2.5" "2.75" "3" "3.25")
chi_s=("2.15" "2.25" "2.5" "2.75" "3" "3.25") #"3.5" "3.75" "4" "4.25" "4.5" "4.75" "5")
T="0.2"
kappa="0.03"

# Define property to be looped over and edited in input file
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

cp ../../inputs_interface_width .

if [ -e "checkpoint_000020000" ]; then
  echo "Run complete"
  cd ../../
  continue
else
  echo "Executing commands in $base_folder"
  # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
  lambda=$(echo "$chi * $T" | bc -l)
  replacement_text="chi = ${lambda}"
  sed -i "s/$search_pattern2/$replacement_text/g" inputs_interface_width # edits input file with appropriate setting
  replacement_text="kappa = ${kappa}"
  sed -i "s/$search_pattern2/$replacement_text/g" inputs_interface_width # edits input file with appropriate setting

  mpirun -n 4 $amrex_exec inputs_interface_width > output.txt
  # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
fi
cd ../..
done
