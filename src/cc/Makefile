F90=gfortran
FFLAGS=-O3 -Wall -fno-align-commons

CC=g++
# CCFLAGS=--std=c++17 -Ofast -W -Iinclude
CXXFLAGS=--std=c++17 -Ofast -W -isystem include
.PHONY: all clean

all: etqte.cc.exe

OBJECTS=etqte.o isosys.o ensemble.o main.o
etqte.cc.exe: $(OBJECTS)
	$(CC) $(OBJECTS) -o etqte.cc.exe


etqte.o: etqte.cc etqte.hh isosys.hh ensemble.hh
isosys.o: isosys.cc isosys.hh
ensemble.o: ensemble.cc ensemble.hh isosys.hh
main.o: main.cc etqte.hh isosys.hh ensemble.hh

clean:
	rm -f *.exe *.o
