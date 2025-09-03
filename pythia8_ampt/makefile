#*************************************************************************
#   > File Name: makefile
#   > Author: Liang Zheng
#   > Mail: liangzhphy@gmail.com 
#   > Created Time: Wed Dec 23 13:19:30 2020
#************************************************************************

OBJS = main.o run_ampt.o linana.o amptsub.o zpc.o art1f.o hipyset1.35.o hijing1.383_ampt.o
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
	$(LD)  -I./  -o $(EXE) $(OBJS) -L$(PYTHIALIBS)/ -lpythia8 $(ROOTLIBS) -lgfortran  -fPIC;
	rsync -a makefile main.C *.f *.h ./bak/

main.o: main.C amptWrapper.h
	$(CXX)  -I./ -I$(PYTHIAINC)  $(ROOTCFLAGS) -c main.C 

linana.o: linana.f
	$(GFC)  -c linana.f

hijing1.383_ampt.o: hijing1.383_ampt.f
	$(GFC)  -c hijing1.383_ampt.f

amptsub.o: amptsub.f
	$(GFC)  -c amptsub.f

zpc.o: zpc.f
	$(GFC)  -c zpc.f

art1f.o: art1f.f
	$(GFC)  -c art1f.f

hipyset1.35.o: hipyset1.35.f
	$(GFC)  -c hipyset1.35.f

run_ampt.o: run_ampt.f
	$(GFC)  -c run_ampt.f



.PHONY : clean
clean :
	-rm -f $(OBJS) $(EXE)
