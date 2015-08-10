# Linux GNU compiler
FC = gfortran -fconvert=big-endian
#FFLAGS = -g -fbounds-check -ffpe-trap=invalid,zero,overflow
FFLAGS = -O3 -funroll-loops

# Linux Intel compiler
#FC = ifort -convert big_endian
#FFLAGS = -g -check bounds -fpe0 -ftrapuv -debug semantic_stepping -debug variable_locations -fpp
#FFLAGS = -O3 -ipo -openmp -no-prec-div

# Linux PGI compiler
#F90 = pgf90 -Mfree -byteswapio 
#FFLAGS = -g -C -Ktrap=fp
#FFLAGS = -fastsse -O4 -Mipa=fast,safe

# IBM compiler and its options
#F90 = xlf90 -qsuffix=f=90
#FFLAGS = -g -C -qflttrap=ov:zero:inv:en -q64 -pg
#FFLAGS = -q64 -qipa -O5 -qhot -qtune=pwr5 -qarch=pwr5 -qfree=f90 -qstrict

# Names and option to fortran compiler
# Dec compiler and its options
#FC = f90 -convert big_endian
#compiler options for debugging
#FFLAGS = -ladebug -g -free -check bounds -check overflow -fpe
# compiler options for fast code
#FFLAGS = -fast -O5

# SGI compiler and its options
#FC = f90
#FFLAGS = -g -C -freeform -mips4 -DEBUG:trap_uninitialized=ON:verbose_ru
#FFLAGS = -Ofast -freeform -mips4 -64 -API

# Sun compiler and its options
#FC = f90
#FFLAGS = -g -free  -dalign -ftrap=%all,no%inexact,no%underflow
#FFLAGS = -free -dalign -native -fast -O5

#names of some MPI directories

# precedence rules
%.o:%.f90
	$(FC) -c $(FFLAGS) $<

SOURCE = cgridops.f90 grid.f90 shallow.f90 shparams.f90 shtstep.f90 \
	trossby2.f90 kinds.f90 shinit.f90 shrhs.f90 shinit.f90


shallowOBJS = kinds.o cgridops.o grid.o shallow.o shtstep.o shinit.o \
	shrhs.o shparams.o
shallow: $(shallowOBJS)
	$(FC) $(FFLAGS) $(shallowOBJS) -o shallow

clean:
	rm *.o *.mod

depend:
	sfmakedepend --depend=obj $(SOURCE)

# DO NOT DELETE THIS LINE - used by make depend
grid.o: kinds.o
shallow.o: grid.o kinds.o shparams.o shrhs.o shtstep.o shinit.o
shinit.o: kinds.o
shparams.o: kinds.o
shrhs.o: kinds.o cgridops.o shparams.o
shtstep.o: kinds.o shparams.o shrhs.o
tlinwave.o: kinds.o shparams.o
trossby2.o: kinds.o
shinit.mod: tlinwave.o
