export CPATH=/opt/homebrew/include
export LIBRARY_PATH=/opt/homebrew/lib

#!/bin/bash
g++-12 main.cpp -o main -fopenmp -std=c++11 -O2 -larmadillo -o main.out
time ./main.out
