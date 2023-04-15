F90=gfortran
FFLAGS=-O3 -Wall -fno-align-commons

CC=g++
CCFLAGS=--std=c++17 -Ofast -Wall -lfmt
.PHONY: all clean

all: etqte.cc.exe

etqte.cc.exe: etqte.cc
	$(CC) $(CCFLAGS) -llapack etqte.cc -o etqte.cc.exe

etqte.f90.exe: etqte.f90
	$(F90) $(FFLAGS) -llapack etqte.f90 -o etqte.f90.exe

clean:
	rm -f *.exe etqte.mod
