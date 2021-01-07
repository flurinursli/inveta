include ../Makefile.inc

F90_OBJ = m_na.o

LIB_NAME = libna.a

na: dirs $(LIB_NAME)

$(F90_OBJ) : %.o : %.f90 
	$(FC) $(FFLAGS) -cpp -I../common/include -c -o $@ $<

$(LIB_NAME) : $(F90_OBJ)
	ar qc $@ $(F90_OBJ)
	mv $(LIB_NAME) lib/
	cp *.mod include/

dirs:
	mkdir -p include lib

clean:
	rm -rf include
	rm -rf lib
	rm -rf *.o
	rm -rf *.mod

.PHONY: clean
