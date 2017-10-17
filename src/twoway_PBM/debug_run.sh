#!/bin/bash

make clean 
make -j 7
mpirun -np 1 xterm -e gdb ./model.out > output_debug.txt
