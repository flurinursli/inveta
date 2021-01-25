include Makefile.inc

OBJ = main.o

all : common llsq fft filter interpolation rtt na $(OBJ)
	$(FC) $(FFLAGS) -fopenmp -o inveta.exe $(OBJ) -Lcommon/lib -Lllsq/lib -Lfft/lib -Linterpolation/lib -Lfilter/lib -Lrtt/lib -Lna/lib -L$(FFTW_PATH)/lib -L$(LAPACK_PATH) -L$(GSL_PATH) -lcommon -lllsqf -lfilterf -lfft -linterplf -lllsq -lrtt -lna $(LINK_FLAGS) -lgsl -lgslcblas -linterpl

# how to get main.o
$(OBJ) : %.o : %.f90
	$(FC) $(FFLAGS) -fopenmp -cpp $(DEBUG) -Icommon/include -Illsq/include32 -Iinterpolation/include32 -Ifft/include -Ifilter/include -Illsq/include32 -Irtt/include -Ina/include -c -o $@ $<

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

.PHONY : common filter fft llsq interpolation rtt na clean
