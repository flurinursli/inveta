FC = mpif90

# flags GCC
#FFLAGS = -g -std=f2008 -fcheck=all -fbacktrace
FFLAGS = -O3 -std=f2008 -march=native

# flags INTEL
#FFLAGS = -O3 -xHost

# flags for lapack, fftw, etc.
#MKL       = -DMKL
#INCL_FFTW = /opt/intel/mkl/include/fftw
INCL_FFTW = $(FFTW_PATH)/include

#LINK_FLAGS = -lifcore -mkl
LINK_FLAGS = -llapack -lblas -lfftw3
