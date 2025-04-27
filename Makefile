# Compiler and flags
FC = gfortran
FFLAGS = -fopenmp -O3

# Files and output
SRC = 2d_mrlbm.f90
OUT = a.out

# Build and run
all:
	$(FC) $(FFLAGS) $(SRC) -o $(OUT)
	./$(OUT)

# Clean rule
.PHONY: clean
clean:
	rm -f $(OUT) *.o *.mod


















