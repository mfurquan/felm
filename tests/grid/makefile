EXE = test_grid

SRCDIR = ../../src
MODDIR = ../../include
OBJDIR = ../../build

FC = gfortran
FFLAGS = -g -fbounds-check -J$(MODDIR)

LD = $(FC)
LFLAGS = 

OBJ = main.o global.o grid.o

$(EXE) : main.o $(OBJDIR)/global.o $(OBJDIR)/grid.o
	$(LD) -o $(EXE) main.o $(OBJDIR)/global.o $(OBJDIR)/grid.o

main.o : main.f90 $(MODDIR)/global.mod $(MODDIR)/grid.mod
	$(FC) -cpp -c -o $@ $< $(FFLAGS)

$(OBJDIR)/global.o $(MODDIR)/global.mod : $(SRCDIR)/global.f90
	$(FC) -cpp -c -o $(OBJDIR)/global.o $< $(FFLAGS)

$(OBJDIR)/grid.o $(MODDIR)/grid.mod : \
	$(SRCDIR)/grid.f90 $(MODDIR)/global.mod
	$(FC) -cpp -c -o $(OBJDIR)/grid.o $< $(FFLAGS)

.PHONY:clean cleanall
clean:
	rm *.o *.mod
cleanall:
	rm $(EXE) *.o *.mod
