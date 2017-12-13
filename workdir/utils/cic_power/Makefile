FC=mpif90

P3DFFT_LIB=${HOME}/lib/p3dfft_2.5.1/lib
P3DFFT_INC=${HOME}/lib/p3dfft_2.5.1/include

LDLIBS= -lp3dfft -lfftw3f 
FLAGS=-mcmodel=medium -shared-intel -fpp -g -O3 -xhost -qopenmp -DNGP -DLOGBIN -DGROUPS

run: ngp_power_dm ngp_power_dmnu ngp_power_dm_init ngp_power_dmnu_init
	@echo "compiled cic power"

ngp_power_dm: cic_power_dmnu.f90
	$(FC) $(FLAGS) -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS) 

ngp_power_dmnu: cic_power_dmnu.f90
	$(FC) $(FLAGS) -DNEUTRINOS -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS) 

ngp_power_dm_init: cic_power_dmnu.f90
	$(FC) $(FLAGS) -DINITCONDITIONS -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS) 

ngp_power_dmnu_init: cic_power_dmnu.f90
	$(FC) $(FLAGS) -DINITCONDITIONS -DNEUTRINOS -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS) 

clean:
	rm -f ngp_power_dm ngp_power_dmnu ngp_power_dm_init ngp_power_dmnu_init *~
