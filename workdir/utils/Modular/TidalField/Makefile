fort=mpif90
xflags=-fpp -O3 -shared-intel -qopenmp -mcmodel=medium -traceback
oflags=-c $(xflags)
dflags= #-O0 -C -B -fp-stack-check -CB -fpe0
oflags+=$(dflags)
xflags+=$(dflags)

F90_FILES := $(wildcard *.f90)
OBJ_FILES := $(addprefix ,$(notdir $(F90_FILES:.f90=.o)))
EXC_FILES := $(addprefix ,$(notdir $(OBJ_FILES:.o=.x)))

MOD_FILES := $(wildcard ../Modules/*.f90)
LNK_FILES := $(wildcard ../Modules/*.o)

PENCIL_LIB=${HOME}/lib/p3dfft_2.5.1/lib
PENCIL_INC=${HOME}/lib/p3dfft_2.5.1/include

CUR_DIR=`basename "$(CURID)"`
MOD_DIR=../Modules

all: Modules $(EXC_FILES)
	@echo "Finished making modules: $(value MOD_FILES)"

Modules:
	cd $(MOD_DIR); make; cd $(CUR_DIR)

%.x: %.o 
	@echo "Files to be linked: $(value LNK_FILES)"
	$(fort) $(xflags)  -I$(FFTW_INC) -I$(PENCIL_INC) -o $@ $< $(value LNK_FILES) -L$(FFTW_LIB) -L$(PENCIL_LIB) -lp3dfft -lfftw3f

%.o: %.f90 
	$(fort) $(oflags) -I$(MOD_DIR) -I$(FFTW_INC) -I$(PENCIL_INC) -o $@ $< -L$(PENCIL_LIB) -L$(FFTW_LIB) -lp3dfft -lfftw3f

clean:
	rm *.mod *.o *.x *~ 
