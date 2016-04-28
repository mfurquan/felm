EXE = test_mesh

SRCDIR = ../../src
MODDIR = ../../include
OBJDIR = ../../build

FC = gfortran
FFLAGS = -g -fbounds-check -J$(MODDIR)

LD = $(FC)
LFLAGS = 

OBJ = main.o global.o mesh.o

$(EXE) : main.o $(OBJDIR)/global.o $(OBJDIR)/mesh.o
	$(LD) -o $(EXE) main.o $(OBJDIR)/global.o $(OBJDIR)/mesh.o

main.o : main.f90 $(MODDIR)/global.mod $(MODDIR)/mesh.mod
	$(FC) -c -o $@ $< $(FFLAGS)

$(OBJDIR)/global.o $(MODDIR)/global.mod : $(SRCDIR)/global.f90
	$(FC) -c -o $(OBJDIR)/global.o $< $(FFLAGS)

$(OBJDIR)/mesh.o $(MODDIR)/mesh.mod : \
	$(SRCDIR)/mesh.f90 $(MODDIR)/global.mod
	$(FC) -c -o $(OBJDIR)/mesh.o $< $(FFLAGS)

.PHONY:clean cleanall
clean:
	rm *.o *.mod
cleanall:
	rm $(EXE) *.o *.mod
