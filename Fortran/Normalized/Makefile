FC=mpifort
FFLAGS=
FFLAGSDEBUG=
FDEBUG=-g -Wall -Wextra -fcheck=all -fbacktrace
FOPTIM=-O3


APP=CEMFortranNormalizedUnits
SRC=$(APP).f90


FFLAGSDEBUG+=$(FDEBUG)
FFLAGS+=$(FOPTIM)

LIB=usr/lib/x86_64-linux-gnu/hdf5/serial
LIB2=usr/local/lib
LIB3=usr/lib/x86_64-linux-gnu/lapack
INC1=usr/include/hdf5/serial
INC2=usr/local/include
INC3=home/wilhelm/Desktop/lapack-3.11.0/SRC


%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -I/$(INC) -c $<

build: $(OBJ)
	$(FC) $(FFLAGS) -I/$(INC1) -I/$(INC2) -L/$(LIB) -L/$(LIB2) -L/$(LIB3) -o $(APP) $(SRC) -lhdf5 -lhdf5_fortran -lfftw3 -fopenmp -llapack -lblas

debug: $(OBJ)
	$(FC) $(FFLAGSDEBUG) -I/$(INC1) -I/$(INC2) -L/$(LIB) -L/$(LIB2) -o $(APP) $(SRC) -lhdf5 -lhdf5_fortran -lfftw3 -fopenmp -llapack -lblas

run:
	make build; ./$(APP)

debugrun:
	make debug; ./$(APP)

cleanrun:
	make cleanall; make run

clean:
	rm -rf $(OBJ)
	rm -rf *.mod
	rm -rf $(APP).x
	rm -rf ./$(APP)

cleanall:
	rm -rf $(OBJ)
	rm -rf *.mod
	rm -rf $(APP).x
	rm -rf ./$(APP)
	rm -rf *.h5