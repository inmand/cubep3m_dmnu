FC=mpif90

P3DFFT_LIB=${HOME}/lib/p3dfft_2.5.1/lib
P3DFFT_INC=${HOME}/lib/p3dfft_2.5.1/include

LDLIBS= -lp3dfft -lfftw3f 
FLAGS=-mcmodel=medium -shared-intel -fpp -g -O3 -xhost -qopenmp -DNGP -DLOGBIN

run: ngp_veldivg
	@echo "compiled cic veldivg"

ngp_veldivg: cic_veldivg_dmnu.f90
	$(FC) $(FLAGS) -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ indexedsort.f90 $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS) 

clean:
	rm -f ngp_veldivg *~
