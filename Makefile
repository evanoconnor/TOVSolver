include make.inc

SOURCES=eosmodule.F90 readtable.F90 nuc_eos.F90 bisection.F90 findtemp.F90 \
	linterp_many.F90 constant_entropy.F90 constant_temperature.F90 \
	derivatives.F90 polytrope.F90 hybrid.F90 findrho.F90 beta_equil.F90
FSOURCES=linterp.f

CLEANSTUFF=rm -rf *.o *.mod *.a

OBJECTS=$(SOURCES:.F90=.o)
FOBJECTS=$(FSOURCES:.f=.o)

EXTRADEPS=

MODINC=$(HDF5INCS)

all: nuc_eos.a TOV

TOV: nuc_eos.a TOV.F90
	$(F90) $(F90FLAGS) -o TOV TOV.F90 nuc_eos.a $(HDF5LIBS)

nuc_eos.a: $(OBJECTS) $(FOBJECTS)
	ar r nuc_eos.a *.o

$(OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) -c $< -o $@

$(FOBJECTS): %.o: %.f $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) -c $< -o $@


clean: 
	$(CLEANSTUFF)
