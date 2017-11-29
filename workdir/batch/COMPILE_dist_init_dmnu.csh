cd ../utils/dist_init

rm -f dist_init_dmnu_dm
rm -f dist_init_dmnu_nu

mpif90 -mcmodel=medium -shared-intel -fpp -g -O3 -xhost -qopenmp dist_init_dmnu.f90 -I$FFTW_INC -I$P3DFFT_INC -o dist_init_dmnu_dm -L$FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
mpif90 -mcmodel=medium -shared-intel -fpp -g -O3 -xhost -qopenmp -DNEUTRINOS -DVELTRANSFER dist_init_dmnu.f90 -I$FFTW_INC -I$P3DFFT_INC -o dist_init_dmnu_nu -L$FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f

cd ../../batch

echo "Sourced dist_init_p3dfft.csh"

