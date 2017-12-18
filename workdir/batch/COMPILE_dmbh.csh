#!/bin/bash

#initial conditions and power codes
cd ../utils/dist_init_dmbh
make clean
make
cd ../../batch
echo "compiled initial conditions"

#cubep3m integrator
cd ../source_threads
make clean
make -f Makefile
cd ../batch
echo "compiled cubep3m integrator"