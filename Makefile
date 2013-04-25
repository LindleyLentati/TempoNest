FC = gfortran
CC = gcc
CXX = g++
FFLAGS += -O3 
CFLAGS += -O3 -DHAVE_CONFIG_H  

LAPACKLIB =  -llapack -lblas

NESTLIBDIR = /data/ltl21/Code/Tempo2Fit/TempoNestRelease/MultiNest
TEMPO2LIB = /data/ltl21/Code/Tempo2Fit/Tempo2-2013.3.1Install/lib
TEMPO2INC = /data/ltl21/Code/Tempo2Fit/Tempo2-2013.3.1Install/include
TEMPO2SRC = /data/ltl21/Code/Tempo2Fit/Tempo2-3.1/tempo2-2013.3.1

LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)  -L$(TEMPO2LIB) -lsofa -ltempo2 -lgsl -lgslcblas

OBJFILES = dpotrs.o dgesvd.o dgemm.o dgemv.o dpotri.o dpotrf.o MultiNestParams.o TempoNestParams.o TempoNestTextOutput.o TempoNestUpdateLinear.o TempoNestFindMax.o TempoNestUtilities.o TempoNestSim.o TempoNestLinearLikeFuncs.o TempoNestLikeFuncs.o TempoNestNoGPU.o

all: TempoNest 

%.o: %.c
	$(CXX) $(CLAGS) -I$(NESTLIBDIR) -I$(TEMPO2INC) -I$(TEMPO2SRC) -c $*.c

 
TempoNest : $(OBJFILES)
	$(FC) -o TempoNest  $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o  TempoNest
