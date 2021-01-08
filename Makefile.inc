FC = mpif90
#FFLAGS = -g -std=f2008 -fcheck=all -fbacktrace
FFLAGS = -O3 -std=f2008 -march=native
#FFLAGS = -O3 -xHost -ipo

#MKL       = -DMKL
#INCL_FFTW = /opt/intel/mkl/include/fftw
INCL_FFTW = $(FFTW_PATH)/include

#LINK_FLAGS = -lifcore -mkl
LINK_FLAGS = -llapack -lblas -lfftw3
