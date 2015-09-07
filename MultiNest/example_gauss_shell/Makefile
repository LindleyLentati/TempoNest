#FC = mpif90
#CC = mpicc
#CXX = mpicc


NESTLIBDIR = ../
LIBS = -L$(NESTLIBDIR) -lnest3 $(LIBS_LAPACK)
 
OBJFILES = params.o like.o nestwrap.o main.o

all: gauss_shell 

%.o: %.f90
	$(FC) $(FFLAGS) -fPIC -I$(NESTLIBDIR) -c $*.f90
 
gauss_shell: $(OBJFILES)
	$(FC) -o ../gauss_shell $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod ../gauss_shell
