lib : common llsq fft filter interpolation rtt na 

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
