#!/bin/bash

# Steps come from:
# https://mprtmma.medium.com/c-shared-library-dynamic-linking-eps-1-bacf2c95d54f

echo "Running Solver"
# Cleanup
echo "--> Cleanup"
rm *.so
rm *.o

# Compile both files to .o files
echo "--> Build"
g++ -std=c++0x -c utils.cpp 
g++ -std=c++0x -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ -c grid.cpp 
g++ -std=c++0x -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ -fPIC -c numerical-methods.cpp -o numerical-methods.o
g++ -std=c++0x -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ -fPIC -c solver.cpp -o solver.o

# Link
echo "--> Link"
g++ -std=c++0x -shared -o libutils.so utils.o
g++ -std=c++0x -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ -shared -L. -lutils -o libgrid.so grid.o
g++ -std=c++0x -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ -shared -L. -lutils -lgrid -o libnumerical_methods.so numerical-methods.o
g++ -std=c++0x -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ solver.cpp -L. -lutils -lgrid -lnumerical_methods -o solver.o

# Run
echo "=== Running ==="
echo ""
LD_LIBRARY_PATH=. ./solver.o
