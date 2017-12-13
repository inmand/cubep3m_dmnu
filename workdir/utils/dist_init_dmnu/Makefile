FC=mpif90

P3DFFT_LIB=${HOME}/lib/p3dfft_2.5.1/lib
P3DFFT_INC=${HOME}/lib/p3dfft_2.5.1/include

LDLIBS= -lp3dfft -lfftw3f 
FLAGS=-mcmodel=medium -shared-intel -fpp -g -O3 -xhost -qopenmp

run: dist_init_dmnu_dm dist_init_dmnu_nu
	@echo "compiled dist init"

dist_init_dmnu_dm: dist_init_dmnu.f90
	$(FC) $(FLAGS) -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS) 

dist_init_dmnu_nu: dist_init_dmnu.f90
	$(FC) $(FLAGS) -DNEUTRINOS -DVELTRANSFER -I$(FFTW_INC) -I$(P3DFFT_INC) -o $@ $< -L$(FFTW_LIB) -L$(P3DFFT_LIB) $(LDLIBS) 

clean:
	rm -f dist_init_dmnu_dm dist_init_dmnu_nu *~
