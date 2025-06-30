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

OBJECTS=inlineoptions_m.o input_m.o f2kcli.o m_getopts.o  \
	calculatereciprocallattice_m.o cellReduction_m.o findRotationalSymmetry_m.o \
	identifyCrystal_m.o restricAtomicCoordinates_m.o \
	findSymmetryOperators_m.o spacegroup_db.o identifySpaceGroup_m.o \
	findDuplicateSymmetryOperators_m.o symmetrizeStructure_m.o writeOutput_m.o
FDF=libfdf.a

all: $(MAINEXE)
$(MAINMOD): $(OBJECTS)
$(MAINEXE): $(OBJECTS) $(MAINMOD) 
	$(FC) -o $(BINDIR)/$(MAINEXE) \
		$(MAINMOD) $(IFLAGS) $(FCFLAGS) $(LDFLAGS) $(OBJECTS) $(LIBS) $(FDF)

inlineoptions_m.o: m_getopts.o
m_getopts.o: f2kcli.o
f2kcli.o:
input_m.o: $(FDF) constants.o
$(FDF):
	(cd $(SRCDIR)/fdf ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)" module;cp fdf.mod $(FDF) ../$(OBJDIR))
constants.o: kind.o
kinds.o:
calculatereciprocallattice_m.o:
cellReduction_m.o:
findRotationalSymmetry_m.o:
identifyCrystal_m.o:
restricAtomicCoordinates_m.o:
findSymmetryOperators_m.o: findDuplicateSymmetryOperators_m.o constants.o
spacegroup_db.o:
identifySpaceGroup_m.o: constants.o
findDuplicateSymmetryOperators_m.o: constants.o
symmetrizeStructure_m.o: constants.o kinds.o
writeOutput_m.o: constants.o calculatereciprocallattice_m.o
spacegroup_db.o: spacegroup_db.f90
	$(FC) -O0 $(IFLAGS) $(MDFLAGS) -c $<

.PHONY: clean

clean:
	cd $(BINDIR); rm -f *.x
	cd $(OBJDIR); rm -f *.o *.mod *.a
	cd $(SRCDIR)/fdf; make clean