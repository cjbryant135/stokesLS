#Makefile for the Chimera code
#################################################################
#Compiler options
#CXX = g++-8 #compiler for my mac
CXX = g++ #for the linux machines
OPT = -O3 #optimization options
DEBUG = -fbounds-check -Wall
CXXFLAGS = -std=c++11

#level set files
LD = -L./LevelSet/
ID = -I/usr/local/include -I./LevelSet/level/

LEVEL_SET = -llevel -llevelinit
LAPACK = -llapacke -llapack

#level set subroutines
levelSets:
	cd LevelSet; make lib; cd ..

#stokes solver
stokesGrid.o: stokesGrid.hpp stokesGrid.cpp
	$(CXX) -c $(ID) stokesGrid.cpp $(CXXFLAGS) $(OPT)

#fluid interface
fluidInterface.o: fluidInterface.hpp fluidInterface.cpp levelSets
	$(CXX) -c $(ID) -DUSE_INLINES fluidInterface.cpp $(CXXFLAGS) $(OPT)

#main file
shearDrop.o: shearDrop.cpp 
	$(CXX) -c $(ID) -DUSE_INLINES shearDrop.cpp $(CXXFLAGS) $(OPT)

shearDrop_exe: shearDrop.o fluidInterface.o stokesGrid.o levelSets
	$(CXX) -o shearDrop_exe shearDrop.o fluidInterface.o stokesGrid.o $(LAPACK) $(LD) $(LEVEL_SET) -lm $(CXXFLAGS) $(OPT)

clean:
	rm -r *_exe *.o *.out *.dSYM *.bin *.csv
