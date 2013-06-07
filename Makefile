FC = gfortran
CC = gcc
CXX = g++
FFLAGS += -O3 
CFLAGS += -O3 -DHAVE_CONFIG_H  

LAPACKLIB =  -llapack -lblas

NESTLIBDIR = /data/ltl21/Code/Tempo2Fit/MultiNest
TEMPO2LIB = /data/ltl21/Code/Tempo2Fit/Tempo2CVS/lib
TEMPO2INC = /data/ltl21/Code/Tempo2Fit/Tempo2CVS/include
TEMPO2SRC = /data/ltl21/Code/Tempo2Fit/Tempo2CVS/Tempo2Src

LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)  -L$(TEMPO2LIB) -lsofa -ltempo2 -lgsl -lgslcblas

OBJFILES = dpotrs.o dgesvd.o dgemm.o dgemv.o dpotri.o dpotrf.o MultiNestParams.o TempoNestParams.o TempoNestTextOutput.o TempoNestUpdateLinear.o TempoNestFindMax.o TempoNestUtilities.o TempoNestSim.o TempoNestLikeFuncs.o TempoNest.o

all: TempoNest

%.o: %.c
	$(CXX) $(CLAGS) -I$(NESTLIBDIR) -I$(TEMPO2INC) -I$(TEMPO2SRC) -c $*.c

 
TempoNest : $(OBJFILES)
	$(FC) -o TempoNest  $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o  TempoNest
