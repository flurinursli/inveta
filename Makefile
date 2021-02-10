include Makefile.inc

PP_FLAGS =

ifeq ($(ERROR_TRAP),y)
	PP_FLAGS += -DERROR_TRAP
endif

ifeq ($(PERF),y)
	PP_FLAGS += -DPERF
endif

ifeq ($(MKL),y)
	PP_FLAGS += -DMKL
endif

OBJ = main.o

DBG = debug.o

all : common llsq fft filter interpolation rtt na $(OBJ)
	$(FC) $(OPT) $(OMP) -o inveta.exe $(OBJ) -Lcommon/lib -Lllsq/lib -Lfft/lib -Linterpolation/lib -Lfilter/lib -Lrtt/lib -Lna/lib -L$(FFTW_PATH)/lib -L$(LAPACK_PATH) -L$(GSL_PATH) -lcommon -lllsqf -lfilterf -lfft -linterplf -lllsq -lrtt -lna $(LINK_FLAGS) -lgsl -lgslcblas -linterpl

debug : common llsq fft interpolation rtt $(DBG)
	$(FC) $(OPT) $(OMP) -o debug.exe $(DBG) -Lcommon/lib -Lllsq/lib -Lfft/lib -Linterpolation/lib -Lrtt/lib -L$(FFTW_PATH)/lib -L$(LAPACK_PATH) -L$(GSL_PATH) -lrtt -lfft -lllsq -lcommon $(LINK_FLAGS) -lgsl -lgslcblas -linterpl

# how to get main.o
$(OBJ) : %.o : %.f90
	$(FC) $(OPT) $(OMP) -cpp $(PP_FLAGS) -Icommon/include -Illsq/include32 -Iinterpolation/include32 -Ifft/include -Ifilter/include -Illsq/include32 -Irtt/include -Ina/include -c -o $@ $<

# how to get debug.o
$(DBG) : %.o : %.f90
	$(FC) $(OPT) $(OMP) -cpp -Icommon/include -Irtt/include -c -o $@ $<


common:
	cd $@; $(MAKE) $@ $<

filter:
	cd $@; $(MAKE) $@ $<

fft:
	cd $@; $(MAKE) $@ $<

llsq:
	cd $@; $(MAKE) $@ $<

interpolation:
	cd $@; $(MAKE) $@ $<

rtt:
	cd $@; $(MAKE) $@ $<

na:
	cd $@; $(MAKE) $@ $<

clean:
	cd common; $(MAKE) $@
	cd filter; $(MAKE) $@
	cd fft; $(MAKE) $@
	cd llsq; $(MAKE) $@
	cd interpolation; $(MAKE) $@
	cd rtt; $(MAKE) $@
	cd na; $(MAKE) $@
	rm -rf *.mod
	rm -rf *.o

.PHONY : common filter fft llsq interpolation rtt na clean debug
