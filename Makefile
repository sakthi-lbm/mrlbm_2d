# Compiler and flags
FC = gfortran
FFLAGS = -fopenmp -O3

# Files and output
SRC = 2d_mrlbm.f90
SRC2 = 2d_statistics.f90
OUT = a.out
OUT2 = stat.out
INPUT = input.dat
INPUT2 = input_stat.dat

run: $(OUT) $(OUT2)
	$(OUT): $(SRC)
	$(FC) $(FFLAGS) $(SRC) -o $(OUT)

	$(OUT2): $(SRC2)
	$(FC) $(FFLAGS) $(SRC2) -o $(OUT2)
	
	@read -p "Enter simulation directory name: " dir; \
	mkdir -p $$dir; \
	cp $(OUT) $(OUT2) $(INPUT) $(INPUT2) $$dir; \
	cd $$dir && ./$(OUT)

# Clean rule
.PHONY: clean
clean:
	rm -f $(OUT) $(OUT2) *.o *.mod
