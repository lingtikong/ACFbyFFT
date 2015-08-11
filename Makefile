.SUFFIXES : .o .c .f .f90
# compiler and flags
FC     = gfortran
FFLAGS = -O3 -Wall
OFLAGS = -O3
FREE   = -ffree-form
#
# fftw 3 library
INC    = -I/opt/libs/fftw/3.2.2/include
LIB    = -L/opt/libs/fftw/3.2.2/lib -lfftw3
#
#====================================================================
# executable name
BASE   = correlation
EXE    = ${BASE}
#====================================================================
# source and rules
SOURCE = vardef.o real.o real_cross.o complex.o complex_cross.o main.o

#====================================================================
all:  ${EXE}

${EXE}:  $(SOURCE)
	$(FC) $(SOURCE) $(OFLAGS) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod ${EXE}

tar:
	rm -f ${BASE}.tar; tar -czvf ${BASE}.tar.gz *.f90 Makefile README


.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
