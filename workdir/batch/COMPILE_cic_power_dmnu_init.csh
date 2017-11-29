cd ../utils/cic_power

rm -f ngp_power_dm_init
rm -f ngp_power_dmnu_init

mpif90 -shared-intel -qopenmp -mcmodel=medium -fpp -g -O3 -DNGP -DLOGBIN -DINITCONDITIONS cic_power_dmnu.f90 -I$FFTW_INC -I$P3DFFT_INC -o ngp_power_dm_init -L$FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

mpif90 -shared-intel -qopenmp -mcmodel=medium -fpp -g -O3 -DNGP -DLOGBIN -DINITCONDITIONS -DNEUTRINOS cic_power_dmnu.f90 -I$FFTW_INC -I$P3DFFT_INC -o ngp_power_dmnu_init -L$FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

cd ../../batch/

echo "Sourced COMPILE_cic_power_dmnu_init.csh"

