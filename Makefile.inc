# mpi wrapper
FC = mpif90

# optimization flags for GCC (add "-march=native" if compute node eqv. compile node)
#OPT = -g -std=f2008 -fcheck=all -fbacktrace
OPT = -O3 -std=f2008 -march=native
OMP = -fopenmp

# optimization flags for INTEL (add "-xHost" if compute node eqv. compile node)
#OPT = -O3 -xHost
#OMP = -qopenmp

# flags for lapack, fftw, etc. libraries: these are system-dependent
# generic libraries (e.g. for GCC)
LINK_FLAGS = -llapack -lblas -lfftw3
INCL_FFTW = $(FFTW_PATH)/include

# comment line above and uncomment line below if you are using SED Bigstar machine
#INCL_FFTW = /usr/include

# make use of MKL libraries
#MKL = y
#INCL_FFTW = /opt/intel/mkl/include/fftw
#LINK_FLAGS = -lifcore -mkl

# Flag for error detection
ERROR_TRAP = y

# Flag for perfomance analysis
PERF = n
