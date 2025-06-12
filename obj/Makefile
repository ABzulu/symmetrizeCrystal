include make.inc

TARGET=symmetrizeCrystal

MAINSRC=$(TARGET).f90
MAINMOD=$(TARGET).o
MAINEXE=$(TARGET).x

SRCDIR=../src
OBJDIR=../obj
BINDIR=../bin

IFLAGS+= -I$(OBJDIR)
VPATH=${SRCDIR}:${SRCDIR}/fdf
LDFLAGS=-L$(OBJDIR) -lfdf

OBJECTS=inlineoptions_m.o input_m.o output_m.o m_getopts.o f2kcli.o \
	calculatereciprocallattice_m.o structureinit_m.o symmetry_m.o identifyCrystal_m.o \
	primitiveToConventional_m.o gcd_m.o spacegroup_db.o identifySpaceGroup_m.o
FDF=libfdf.a

all: $(MAINEXE)
$(MAINMOD): $(OBJECTS)
$(MAINEXE): $(OBJECTS) $(MAINMOD) 
	$(FC) -o $(BINDIR)/$(MAINEXE) \
		$(MAINMOD) $(IFLAGS) $(FCFLAGS) $(LDFLAGS) $(OBJECTS) $(LIBS) $(FDF)

inlineoptions_m.o: m_getopts.o
m_getopts: f2kcli.o
input_m.o: $(FDF) constants.o
$(FDF):
	(cd $(SRCDIR)/fdf ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)" module;cp fdf.mod $(FDF) ../$(OBJDIR))
constants.o: kind.o
kinds.o:
output_m.o: calculatereciprocallattice_m.o
calculatereciprocallattice_m.o:
structureinit_m.o: constants.o calculatereciprocallattice_m.o
symmetry_m.o:
identifyCrystal_m.o:
primitiveToConventional_m.o: constants.o gcd_m.o
spacegroup_db.o:
	$(FC) -O0 $(IFLAGS) $(FCFLAGS) $(MDFLAGS) -c $<
identifySpaceGroup_m.o:

.PHONY: clean

clean:
	cd $(BINDIR); rm -f *.x
	cd $(OBJDIR); rm -f *.o *.mod *.a
	cd $(SRCDIR)/fdf; make clean