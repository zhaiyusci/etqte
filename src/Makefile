F90=gfortran
FFLAGS=-O3 -Wall
.PHONY: all clean

all: etqte.exe


etqte.exe: etqte.f90
	$(F90) $(FFLAGS) -llapack etqte.f90 -o etqte.exe

clean:
	rm -f *.exe
