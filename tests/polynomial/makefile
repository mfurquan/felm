EXE = test_integration

SRCDIR = ../../src
MODDIR = ../../include
OBJDIR = ../../build

FC = gfortran
FFLAGS = -g -fbounds-check# -I$(MODDIR)

LD = $(FC)
LFLAGS = 

OBJ = global.o polynomial.o array_utils

$(EXE) : main.o $(OBJDIR)/global.o $(OBJDIR)/polynomial.o $(OBJDIR)/array_utils.o
	$(LD) -o $(EXE) main.o $(OBJDIR)/global.o $(OBJDIR)/polynomial.o \
		$(OBJDIR)/array_utils.o

main.o : main.f90 $(MODDIR)/global.mod $(MODDIR)/polynomial.mod
	$(FC) -c -o $@ $< $(FFLAGS)

$(OBJDIR)/global.o $(MODDIR)/global.mod : $(SRCDIR)/global.f90
	$(FC) -c -o $(OBJDIR)/global.o $< $(FFLAGS)

$(OBJDIR)/array_utils.o $(MODDIR)/array_utils.mod : \
	$(SRCDIR)/array_utils.f90 $(MODDIR)/global.mod
	$(FC) -c -o $(OBJDIR)/array_utils.o $< $(FFLAGS)

$(OBJDIR)/polynomial.o $(MODDIR)/polynomial.mod : \
	$(SRCDIR)/polynomial.f90 $(MODDIR)/global.mod $(MODDIR)/array_utils.mod
	$(FC) -c -o $(OBJDIR)/polynomial.o $< $(FFLAGS)

.PHONY:clean cleanall
	clean:
	rm *.o *.mod
cleanall:
	rm $(EXE) *.o *.mod
