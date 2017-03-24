#
# Which compiler to use
#
#COMPILER = INTEL
COMPILER = GNU
#
#
ifeq "$(COMPILER)" "INTEL"
   MPIFORTRAN = mpiifort
   FORTRAN = ifort
   FLAGS = -O3
endif
ifeq "$(COMPILER)" "GNU"
   MPIFORTRAN = mpif90
   FORTRAN = gfortran
#   FLAGS = -O3 -ffast-math
   FLAGS = -fbounds-check -Wunused
#   FLAGS = -g -pg 
endif
#
# Flag for Lapack:
LAPACKFLAGS = -llapack
#
# Build directories
#
BIN=../bin
 
###################################################################
#                                                                 #
# Generally no modifications are required after this.             #
#                                                                 #
###################################################################

all : $(BIN)/crfit $(BIN)/crfitone $(BIN)/createrandom

DIRECTORIES = bindir

MODULES = force_field.mod \
          constraints.mod \
          dihedrals.mod \
          linkedcells.mod \
          flashsortvars.mod 

COMOBJ = build_exclusions.o \
         compute_enbcell.o \
         compute_gnbcell.o \
         spg.o \
         cgnewton.o \
         optimize.o \
         common_modules.o \
         compute_eangle.o \
         compute_ebond.o \
         compute_edihed.o \
         compute_eimpr.o \
         compute_enbls.o \
         compute_enb.o \
         compute_fconst.o \
         computef.o \
         compute_gangle.o \
         compute_gbond.o \
         compute_gconst.o \
         compute_gdihed.o \
         computeg.o \
         compute_gimpr.o \
         compute_gnbls.o \
         compute_gnb.o \
         flashsort.o \
         ga_dihedrals.o \
         read_constraints.o \
         read_coor_pdb.o \
         readprm.o \
         readpsf.o \
         test_grad.o \
         addbest.o \
         random.o \
         newangles.o \
         align.o \
         solve_nb_cut.o \
         write_pdb.o 

MPIOBJ = stop_all_mpi.o 
SEROBJ = stop_all_serial.o

CRFIT_MAIN = crfit.o
CRFITONE_MAIN = crfitone.o
CREATERANDOM_MAIN = createrandom.o

$(BIN)/crfit: Makefile $(MODULES) $(COMOBJ) $(MPIOBJ) $(DIRECTORIES) $(CRFIT_MAIN)
	$(MPIFORTRAN) -o $(BIN)/crfit $(COMOBJ) $(MPIOBJ) $(CRFIT_MAIN) $(FLAGS) $(LAPACKFLAGS)

$(BIN)/crfitone: Makefile $(MODULES) $(COMOBJ) $(SEROBJ) $(DIRECTORIES) $(CRFITONE_MAIN)
	$(FORTRAN) -o $(BIN)/crfitone $(COMOBJ) $(SEROBJ) $(CRFITONE_MAIN) $(FLAGS) $(LAPACKFLAGS)

$(BIN)/createrandom: Makefile $(MODULES) $(COMOBJ) $(SEROBJ) $(DIRECTORIES) $(CREATERANDOM_MAIN)
	$(FORTRAN) -o $(BIN)/createrandom $(COMOBJ) $(SEROBJ) $(CREATERANDOM_MAIN) $(FLAGS) $(LAPACKFLAGS)

crfit.o : crfit.f90
	$(MPIFORTRAN) -c $(FLAGS) crfit.f90

crfitone.o : crfitone.f90
	$(FORTRAN) -c $(FLAGS) crfitone.f90

createrandom.o : createrandom.f90
	$(FORTRAN) -c $(FLAGS) createrandom.f90

bindir : 
	mkdir -p ../bin

common_modules.o force_field.mod constraints.mod dihedrals.mod : common_modules.f90
	@rm -f force_field.mod constraints.mod dihedrals.mod
	$(FORTRAN) -c $(FLAGS) common_modules.f90

write_pdb.o : write_pdb.f90
	$(FORTRAN) -c $(FLAGS) write_pdb.f90

flashsort.o flashsortvars.mod : flashsort.f90
	@rm -f flashsortvars.mod
	$(FORTRAN) -c $(FLAGS) flashsort.f90

compute_enbcell.o : compute_enbcell.f90
	$(FORTRAN) -c $(FLAGS) compute_enbcell.f90

compute_gnbcell.o : compute_gnbcell.f90
	$(FORTRAN) -c $(FLAGS) compute_gnbcell.f90

build_exclusions.o : build_exclusions.f90
	$(FORTRAN) -c $(FLAGS) build_exclusions.f90

solve_nb_cut.o : solve_nb_cut.f90
	$(FORTRAN) -c $(FLAGS) solve_nb_cut.f90

spg.o : spg.f90
	$(FORTRAN) -c $(FLAGS) spg.f90

cgnewton.o : cgnewton.f90
	$(FORTRAN) -c $(FLAGS) cgnewton.f90

optimize.o : optimize.f90
	$(FORTRAN) -c $(FLAGS) optimize.f90

compute_eangle.o : compute_eangle.f90
	$(FORTRAN) -c $(FLAGS) compute_eangle.f90

compute_ebond.o : compute_ebond.f90
	$(FORTRAN) -c $(FLAGS) compute_ebond.f90

compute_edihed.o : compute_edihed.f90
	$(FORTRAN) -c $(FLAGS) compute_edihed.f90

compute_eimpr.o : compute_eimpr.f90
	$(FORTRAN) -c $(FLAGS) compute_eimpr.f90

compute_enbls.o : compute_enbls.f90
	$(FORTRAN) -c $(FLAGS) compute_enbls.f90

compute_enb.o : compute_enb.f90
	$(FORTRAN) -c $(FLAGS) compute_enb.f90

compute_fconst.o : compute_fconst.f90
	$(FORTRAN) -c $(FLAGS) compute_fconst.f90

computef.o : computef.f90
	$(FORTRAN) -c $(FLAGS) computef.f90

compute_gangle.o : compute_gangle.f90
	$(FORTRAN) -c $(FLAGS) compute_gangle.f90

compute_gbond.o : compute_gbond.f90
	$(FORTRAN) -c $(FLAGS) compute_gbond.f90

compute_gconst.o : compute_gconst.f90
	$(FORTRAN) -c $(FLAGS) compute_gconst.f90

compute_gdihed.o : compute_gdihed.f90
	$(FORTRAN) -c $(FLAGS) compute_gdihed.f90

computeg.o : computeg.f90
	$(FORTRAN) -c $(FLAGS) computeg.f90

compute_gimpr.o : compute_gimpr.f90
	$(FORTRAN) -c $(FLAGS) compute_gimpr.f90

compute_gnbls.o : compute_gnbls.f90
	$(FORTRAN) -c $(FLAGS) compute_gnbls.f90

compute_gnb.o : compute_gnb.f90
	$(FORTRAN) -c $(FLAGS) compute_gnb.f90

ga_dihedrals.o : ga_dihedrals.f90
	$(FORTRAN) -c $(FLAGS) ga_dihedrals.f90

read_constraints.o : read_constraints.f90
	$(FORTRAN) -c $(FLAGS) read_constraints.f90

read_coor_pdb.o : read_coor_pdb.f90
	$(FORTRAN) -c $(FLAGS) read_coor_pdb.f90

readprm.o : readprm.f90
	$(FORTRAN) -c $(FLAGS) readprm.f90

readpsf.o : readpsf.f90
	$(FORTRAN) -c $(FLAGS) readpsf.f90

stop_all_serial.o : stop_all_serial.f90
	$(FORTRAN) -c $(FLAGS) stop_all_serial.f90

stop_all_mpi.o : stop_all_mpi.f90
	$(MPIFORTRAN) -c $(FLAGS) stop_all_mpi.f90

addbest.o : addbest.f90
	$(FORTRAN) -c $(FLAGS) addbest.f90

newangles.o : newangles.f90
	$(FORTRAN) -c $(FLAGS) newangles.f90

align.o : align.f90
	$(FORTRAN) -c $(FLAGS) align.f90 $(LAPACKFLAGS)

random.o : random.f90
	$(FORTRAN) -c $(FLAGS) random.f90

test_grad.o : test_grad.f90
	$(FORTRAN) -c $(FLAGS) test_grad.f90

clean : 
	rm -f ./*.o ./*.mod $(BIN)/*

