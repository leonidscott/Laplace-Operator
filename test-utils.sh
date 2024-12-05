#!/bin/bash

# Build Utils + Test Files
echo " -> Building"
g++ -std=c++0x -c utils.cpp
g++ -std=c++0x -c utils-test.cpp

# Link
echo " -> Linking"
g++ -std=c++0x -shared -o libutils.so utils.o
g++ -std=c++0x utils-test.cpp -L. -lutils -o utils-test.o

# Run
echo "=== Running ==="
echo ""
LD_LIBRARY_PATH=. ./utils-test.o
