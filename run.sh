#!/bin/bash

# Steps come from:
# https://mprtmma.medium.com/c-shared-library-dynamic-linking-eps-1-bacf2c95d54f

echo "Run Laplace Operator 2"
# Cleanup
echo "--> Cleanup"
rm *.so
rm *.o

# Compile both files to .o files
echo "--> Build"
g++ -std=c++0x -c utils.cpp 
g++ -std=c++0x -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ -fPIC -c laplace-operator-2.cpp -o laplace-operator-2.o

# Link
echo "--> Link"
g++ -std=c++0x -shared -o libutils.so utils.o
g++ -std=c++0x -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ laplace-operator-2.cpp -L. -lutils -o laplace-operator-2.o

# Run
echo "=== Running ==="
echo ""
LD_LIBRARY_PATH=. ./laplace-operator-2.o
