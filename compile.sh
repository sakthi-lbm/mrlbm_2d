#!/bin/bash

# === CONFIGURATION ===
CODE_DIR="$(pwd)"                   # Current working directory
INPUT_FILE="input.dat"
INPUT_FILE2="input_stat.dat"
CODE_FILE="2d_mrlbm.f90" 
CODE_FILE2="2d_statistics.f90"
MAKE_FILE="Makefile"
BASE_DIR="$PWD/simulations"

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
cp "$CODE_DIR/$CODE_FILE" "$SIM_DIR/"
cp "$CODE_DIR/$CODE_FILE2" "$SIM_DIR/" 
cp "$CODE_DIR/$INPUT_FILE" "$SIM_DIR/" 
cp "$CODE_DIR/$INPUT_FILE2" "$SIM_DIR/" 
cp "$CODE_DIR/$MAKE_FILE" "$SIM_DIR/" 

# === OPEN NEW TERMINAL AND RUN SIMULATION ===
gnome-terminal -- bash -c "cd '$SIM_DIR' && make clean && make; exec bash"

echo "Simulation is running in a new terminal. You can continue using this terminal."
