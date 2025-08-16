OBJS = analyse.o AMPT.o
EXE = analyse
OBJS_FLOW = analyse_decayer.o AMPT.o
EXE_FLOW = analyse_decayer



ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)

FASTJETINC=${HOME}/install/fastjet-3.3.2/include/
FASTJETLIBS=${HOME}/install/fastjet-3.3.2/lib/


INCFLAGS = -I$(ROOTSYS)/include 
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++
FLAGS = -Wall -g $(ROOTCFLAGS) -I$(FASTJETINC) -I$(PYTHIA8)/include/

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE_FLOW)

$(EXE): $(OBJS)
	$(CXX) -I$(FASTJETINC) $(ROOTCFLAGS) -o $(EXE) $(OBJS) $(ROOTLIBS) $(FASTJETLIBS)/libfastjet.a -L$(PYTHIA8)/lib -lpythia8
#	$(CXX) $(ROOTCFLAGS) -I$(HOME)/Macro/CODE/ -o $(EXE) $(OBJS) $(ROOTLIBS) 

$(EXE_FLOW): $(OBJS_FLOW)
	$(CXX) -I$(FASTJETINC) $(ROOTCFLAGS) -o $(EXE_FLOW) $(OBJS_FLOW) $(ROOTLIBS) $(FASTJETLIBS)/libfastjet.a -L$(PYTHIA8)/lib -lpythia8
#	$(CXX) $(ROOTCFLAGS) -I$(HOME)/Macro/CODE/ -o $(EXE_FLOW) $(OBJS_FLOW) $(ROOTLIBS) 


%.o: %.cxx
	$(COMPILE) -I$(HOME)/Macro/CODE/ $<

.PHONY : clean
clean :
	-rm -f $(OBJS) $(EXE) $(EXE_FLOW) 
