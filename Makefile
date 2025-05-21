# Compiler and flags
FC = gfortran
FFLAGS = -fopenmp -O3

# Files and output
SRC = 2d_mrlbm.f90
OUT = a.out

DIRS = output restart data_probe data_geo data_mean

# Build and run
all: $(DIRS)
	$(FC) $(FFLAGS) $(SRC) -o $(OUT)
	./$(OUT)

# Rule to create directories if they don't exist
$(DIRS):
	mkdir -p $@

# Clean rule
.PHONY: clean
clean:
	rm -f $(OUT) *.o *.mod


















