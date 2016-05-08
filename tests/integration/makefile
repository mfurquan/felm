EXE = test_integration

SRCDIR = ../../src
MODDIR = ../../include
OBJDIR = ../../build

FC = gfortran
FFLAGS = -g -fbounds-check# -I$(MODDIR)

LD = $(FC)
LFLAGS = 

OBJ = global.o integration.o gauss_cubature.o

$(EXE) : main.o $(OBJDIR)/global.o $(OBJDIR)/gauss_cubature.o \
	$(OBJDIR)/integration.o
	$(LD) -o $(EXE) main.o $(OBJDIR)/global.o $(OBJDIR)/gauss_cubature.o \
	$(OBJDIR)/integration.o

main.o : main.f90 $(MODDIR)/global.mod $(MODDIR)/integration.mod
	$(FC) -c -o $@ $< $(FFLAGS)

$(OBJDIR)/global.o $(MODDIR)/global.mod : $(SRCDIR)/global.f90
	$(FC) -c -o $(OBJDIR)/global.o $< $(FFLAGS)

$(OBJDIR)/gauss_cubature.o $(MODDIR)/gauss_cubature.mod : \
	$(SRCDIR)/gauss_cubature.f90 $(MODDIR)/global.mod
	$(FC) -c -o $(OBJDIR)/gauss_cubature.o $< $(FFLAGS)

$(OBJDIR)/integration.o $(MODDIR)/integration.mod : \
	$(SRCDIR)/integration.f90 $(MODDIR)/global.mod \
	$(MODDIR)/gauss_cubature.mod
	$(FC) -c -o $(OBJDIR)/integration.o $< $(FFLAGS)

.PHONY:clean cleanall
	clean:
	rm *.o *.mod
cleanall:
	rm $(EXE) *.o *.mod
