#!/bin/bash

# === CONFIGURATION ===
CODE_DIR="$(pwd)"                   # Current working directory
EXECUTABLE="a.out"    # Change to your actual executable name
INPUT_FILE="input.dat"              # If applicable, change or remove this
CODE_FILE="2d_mrlbm.f90"                # Your Fortran code file
CODE_FILE2="detect_vortex_cycle.f90"                # Your Fortran code file
BASE_DIR="$PWD/simulations"         # Base directory for all simulations
#BASE_DIR="$PWD/simulations"         # Base directory for all simulations

# === GET SIMULATION NAME ===
if [ -z "$1" ]; then
  echo "Error: Simulation directory name not provided."
  echo "Usage: ./run_case.sh <simulation_name>"
  exit 1
else
  SIM_NAME="$1"
fi

SIM_DIR="${BASE_DIR}/${SIM_NAME}"

# === CREATE SIMULATION DIRECTORY ===
mkdir -p "$SIM_DIR"
echo "Created simulation directory: $SIM_DIR"

# === COPY FILES ===
cp "$CODE_DIR/"*.[fF]* "$SIM_DIR/"       # Copy all Fortran source files
cp "$CODE_DIR/$CODE_FILE" "$SIM_DIR/"    # Copy the specific Fortran code file (code.f90)
cp "$CODE_DIR/$CODE_FILE2" "$SIM_DIR/"    # Copy the specific Fortran code file (code.f90)
cp "$CODE_DIR/$EXECUTABLE" "$SIM_DIR/"   # Copy the executable
[ -f "$CODE_DIR/$INPUT_FILE" ] && cp "$CODE_DIR/$INPUT_FILE" "$SIM_DIR/"  # Copy input file if it exists

# === OPEN NEW TERMINAL AND RUN SIMULATION ===
gnome-terminal -- bash -c "cd '$SIM_DIR' && gfortran -fopenmp -O3 2d_mrlbm.f90 && ./$EXECUTABLE; exec bash"

echo "Simulation is running in a new terminal. You can continue using this terminal."
