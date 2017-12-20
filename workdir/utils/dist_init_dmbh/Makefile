FC=mpif90

P3DFFT_LIB=${HOME}/lib/p3dfft_2.5.1/lib
P3DFFT_INC=${HOME}/lib/p3dfft_2.5.1/include

LDLIBS= -lp3dfft -lfftw3f 
FLAGS=-mcmodel=medium -shared-intel -fpp -g -O3 -xhost -qopenmp

run: dist_init_dmbh ngp_power ngp_veldivg
	@echo "compiled dmbh"

dist_init_dmbh: dist_init_dmbh.f90
	$(FC) $(FLAGS) -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS)

ic_pbh: ic_pbh.f90
	$(FC) $(FLAGS) -o $@ $<

ngp_power: cic_power_dmnu.f90
	$(FC) $(FLAGS) -DNGP -DLOGBIN -DINITCONDITIONS -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS)

ngp_veldivg: cic_veldivg_dmnu.f90
	$(FC) $(FLAGS) -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ indexedsort.f90 $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS)

clean:
	rm -f dist_init_dmbh ic_pbh ngp_power ngp_veldivg *~
