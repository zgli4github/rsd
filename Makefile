F77 = mpif90
#F77 = ifort
Fcc = gcc
DUMMY_FOX= --enable-dummy
DEBUG = -g -O2 # -lmpich 

################################
DRIVER = radecz
################################

HEALPIX_DIR = /opt/Healpix_3.11
CFITS_DIR = /opt/cfitsio
LAPACK_DIR = /opt/intel/composer_xe_2013_sp1.0.080/mkl/lib/intel64
LAPACK_INCL = -I/opt/intel/composer_xe_2013_sp1.0.080/mkl/include/intel64/lp64
MPICH_DIR = /opt/mpich304

LAPACK_LIBS = -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread -lm
LAPACK_LIB = -L$(LAPACK_DIR) $(LAPACK_LIBS)

HEALPIX_INCL = -I$(HEALPIX_DIR)/include
HEALPIX_LIBS = -L$(HEALPIX_DIR)/lib
HEALPIX_LIB = $(HEALPIX_LIBS) -lhealpix

CFITS_INCL = -I$(CFITS_DIR)/include
CFITS_LIBS = -L$(CFITS_DIR)/lib
CFITS_LIB = $(CFITS_LIBS) -lcfitsio

MPICH_INCL = -I$(MPICH_DIR)/include
MPICH_LIB = -L$(MPICH_DIR)/lib -lmpichf90

INCLS = $(CFITS_INCL) $(HEALPIX_INCL) # $(LAPACK_INCL)
LIBS  = $(CFITS_LIB)  $(HEALPIX_LIB)  # $(LAPACK_LIB) 


run = $(DRIVER)
default: $(run)
all: $(run)
GOBJ =  integral.o utils_zg.o $(DRIVER).o healpix.o hr.o

$(DRIVER).o: utils_zg.o healpix.o hr.o
utils_zg.o: integral.o
hr.o: healpix.o

$(run):  $(GOBJ)
        $(F77) $(DEBUG) $(INCLS) $(GOBJ) $(LIBS) -o $@

%.o: %.f90
        $(F77) $(DEBUG)  $(INCLS) -c  $*.f90

clean:
        rm *.o *.mod $(run)
