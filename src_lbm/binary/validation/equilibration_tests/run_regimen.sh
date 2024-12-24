#!/bin/bash

amrex_exec="../../../main3d.gnu.TPROF.MPI.ex"

# Defines the radius of the droplet as a proportion of system size. 0.5 is maximum
noise_type=("0" "1")

# Define the base folder name
base_folder="spatially_"
search_pattern1="correlated_noise = 0"

# Loop through each radius in R
for noise in "${noise_type[@]}"; do

  if [ "$noise" == "0" ]; then
    suffix="independent"
  else
    suffix="dependent"
  fi

  folder="${base_folder}${suffix}"

  echo "Processing folder: $folder"
  
  # Check if the folder exists. Creates folder if it does not exist
  if [ ! -d "$folder" ]; then
    echo "$folder does not exist. Creating $folder"
    mkdir -p "$folder"
  fi

  cp inputs $folder
  cp inputs_restart $folder

  # Enters folder and executes commands before returning to parent directory
  cd "$folder"
  # Checks if run is complete and executes a run with appropriate modifications if it has not
  if [ -e "chk_hydro_0000050000" ]; then
    echo "Run complete"
    continue
  else
    echo "Executing commands in $folder"
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
    replacement_text="correlated_noise = $noise"
    sed -i "s/$search_pattern1/$replacement_text/g" inputs # edits input file with appropriate setting
    sed -i "s/$search_pattern1/$replacement_text/g" inputs_restart
    ./$amrex_exec inputs
    # EDIT COMMANDS HERE TO MAKE MODIFICATIONS TO RUNS #
  fi
  cd ..

done
