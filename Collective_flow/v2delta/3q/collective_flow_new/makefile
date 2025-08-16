#OBJS = analyse.o AMPT.o
#EXE = analyse
OBJS = analyse_flow_longRange.o AMPT.o
EXE = analyse_flow_longRange


ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)

FASTJETINC=${HOME}/install/fastjet-3.3.2/include/
FASTJETLIBS=${HOME}/install/fastjet-3.3.2/lib/


INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++
FLAGS = -Wall -g $(ROOTCFLAGS) -I$(FASTJET)/include -I$(PYTHIA8)/include/

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -I$(FASTJETINC) $(ROOTCFLAGS) -o $(EXE) $(OBJS) $(ROOTLIBS) -L$(FASTJET)/lib/ -lfastjet -L$(PYTHIA8)/lib -lpythia8

%.o: %.cxx
	$(COMPILE)  $<

.PHONY : clean
clean :
	-rm -f $(OBJS) $(EXE)
