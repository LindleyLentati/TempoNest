CC = mpicc -cc=gcc
FC = mpif90 -fc=gfortran

CHORDLIBDIR = ../src/
LIBS =  -L$(CHORDLIBDIR) -lchord -lstdc++
 
OBJFILES = eggbox.o

all: eggbox 

%.o: %.cpp
	$(CC) $(CFLAGS) -c $*.cpp
 
eggbox: $(OBJFILES)
	$(FC) $(FFLAGS) -o eggboxC $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod eggboxC
