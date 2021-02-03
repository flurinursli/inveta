FC = mpif90

# flags GCC (add "-march=native" if compute node eqv. compile node)
#OPT = -g -std=f2008 -fcheck=all -fbacktrace
OPT = -O3 -std=f2008 -march=native

# flags INTEL (add "-xHost" if compute node eqv. compile node)
#OPT = -O3 -xHost

OMP = -fopenmp
#OMP = -qopenmp

# flags for lapack, fftw, etc.
#MKL = -DMKL
DEBUG = -DDEBUG
#INCL_FFTW = /opt/intel/mkl/include/fftw
INCL_FFTW = $(FFTW_PATH)/include

#LINK_FLAGS = -lifcore -mkl
LINK_FLAGS = -llapack -lblas -lfftw3
