#!/bin/bash 

#initial conditions
cd ../utils/dist_init_dmnu/
make clean
make
cd ../../batch

#cubep3m integrator
cd ../source_threads
make clean
make 
cd ../batch
echo "compiled cubep3m integrator"

#power spectra
cd ../utils/cic_power
make clean
make
cd ../../batch

#velocity spectra
cd ../utils/cic_velpower
make clean
make
cd ../../batch

echo ""
echo "finished compilations"
