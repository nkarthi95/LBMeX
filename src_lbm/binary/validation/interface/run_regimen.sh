#!/bin/bash

# Defines whether to use the uncorrelated (0) or correlated noise (1)
noise_type=("0")
# Define the base folder name
base_folder="spatially_"

# Loop through each radius in R
for noise in "${noise_type[@]}"; do

  if [ "$noise" == "0" ]; then
    suffix="independent"
  else
    suffix="dependent"
  fi

  folder="${base_folder}${suffix}"
  amrex_exec="../../../${folder}/main3d.gnu.MPI.ex" #setup from the perspective of the folder it is run in

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
  if [ -e "chk_hydro_0000300000" ]; then
    echo "Run complete"
    continue
  else
    echo "Executing commands in $folder"
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
    mpirun -n 8 $amrex_exec inputs_production > run_output.txt
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
  fi
  cd ..

done
