EXEC = int2
#HOST = gcc
HOST = gcc-openmp
#HOST = intel
#HOST = intel-ompenmp

#
# For linux systems, it is assumed that the environment
# variable LD_LIBRARY_PATH contains the locations to libfmm3d.so
# and libsolvers3d.so, for Macosx, these .so files also need to be
# copied over /usr/local/lib


ifneq ($(OS),Windows_NT) 
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        LDF = /usr/local/lib
    endif
    ifeq ($(UNAME_S),Linux)
        LDF = ${HOME}/lib
    endif
endif


LIBS = -lfmm3d -lfmm3dbie 
ifeq ($(HOST),gcc)
    FC=gfortran -L${LDF} 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy 
endif

ifeq ($(HOST),intel)
    FC=ifort -L${LDF} 
    FFLAGS= -O3 -fPIC -march=native
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort -L${LDF} 
    FFLAGS= -O3 -fPIC -march=native -qopenmp
endif

SRC=../src

.PHONY: all clean 

OBJECTS =  comp_statj_lambdainf_paper.o \
    $(SRC)/staticj_gendeb_final_paper.o $(SRC)/surf_routs.o \
    $(SRC)/get_harm_vec_fields.o $(SRC)/lapbel_kernels.o \
    $(SRC)/lapbel2fast.o


#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f  
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.f90  
	$(FC) -c $(FFLAGS) $< -o $@

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) -L${LDF} $(LIBS)
	./$(EXEC)  

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



