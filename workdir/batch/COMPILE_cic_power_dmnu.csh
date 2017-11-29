cd ../utils/cic_power

rm -f ngp_power_dm
rm -f ngp_power_dmnu

mpif90 -shared-intel -qopenmp -mcmodel=medium -fpp -g -O3 -DNGP -DLOGBIN -DGROUPS cic_power_dmnu.f90 -I$FFTW_INC -I$P3DFFT_INC -o ngp_power_dm -L$FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

mpif90 -shared-intel -qopenmp -mcmodel=medium -fpp -g -O3 -DNGP -DLOGBIN -DGROUPS -DNEUTRINOS cic_power_dmnu.f90 -I$FFTW_INC -I$P3DFFT_INC -o ngp_power_dmnu -L$FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f


cd ../../batch/

echo "Sourced COMPILE_cic_power_dmnu.csh"

