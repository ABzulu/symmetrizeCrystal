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

OBJECTS=input_m.o
FDF=libfdf.a

all: $(MAINEXE)
$(MAINMOD): $(OBJECTS)
$(MAINEXE): $(OBJECTS) $(MAINMOD) 
	$(FC) -o $(BINDIR)/$(MAINEXE) \
		$(MAINMOD) $(IFLAGS) $(FCFLAGS) $(LDFLAGS) $(OBJECTS) $(LIBS) $(FDF)

input_m.o: $(FDF) constants.o
$(FDF):
	(cd $(SRCDIR)/fdf ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)" module;cp fdf.mod $(FDF) ../$(OBJDIR))
constants.o: kind.o
kinds.o:

.PHONY: clean

clean:
	cd $(BINDIR); rm -f *.x
	cd $(OBJDIR); rm -f *.o *.mod *.a
	cd $(SRCDIR)/fdf; make clean