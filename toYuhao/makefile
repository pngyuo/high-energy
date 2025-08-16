OBJS = main.o 
EXE = main.exe

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)

PYTHIAINC=${PYTHIA8}/include/
PYTHIALIBS=${PYTHIA8}/lib/


INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++  -g
GFC = gfortran -g
LD = g++ -g

all: $(EXE)

$(EXE): $(OBJS)
	$(LD)  -I./  -o $(EXE) $(OBJS) -L$(PYTHIALIBS)/ -lpythia8 $(ROOTLIBS) 

main.o: main.C 
	$(CXX)  -I./ -I$(PYTHIAINC)  $(ROOTCFLAGS) -c main.C 


.PHONY : clean
clean :
	-rm -f $(OBJS) $(EXE)
