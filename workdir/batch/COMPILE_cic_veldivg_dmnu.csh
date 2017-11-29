cd ../utils/cic_velpower

rm -f ngp_veldivg

mpif90 -shared-intel -fpp -g -O3 -qopenmp -mcmodel=medium -DNGP -DLOGBIN indexedsort.f90 cic_veldivg_dmnu.f90 -I$FFTW_INC -I$P3DFFT_INC -o ngp_veldivg -L$FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

cd ../../batch/

echo "Sourced COMPILE_cic_veldivg.csh" 

