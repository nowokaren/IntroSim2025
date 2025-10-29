# -------------------------------------------------
# Makefile para el proyecto Ising
# -------------------------------------------------
FC     = gfortran
FFLAGS = -O3
EXE    = ising
OBJS   = ziggurat.o ising.o

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXE) *.dat

ising: $(EXE)
