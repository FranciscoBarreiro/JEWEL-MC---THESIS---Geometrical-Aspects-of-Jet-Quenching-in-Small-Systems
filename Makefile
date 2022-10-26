all: jewel-gajos_do_glauber-vac jewel-gajos_do_glauber-simple
#all: jewel-gajos_do_glauber-simple

CC := g++
CCFLAGS := -std=c++14 -g
#INC := -I src

FC := gfortran
FFLAGS := -g -static

# path to LHAPDF library
LHAPDF_PATH := /home/barreiro/local/lib

ROOTINC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --glibs)
LDLIB = -lMathMore

jewel-gajos_do_glauber-vac: jewel-2.3.0.o glauber-medium-vac.o pythia6425mod-lhapdf6.o meix.o
	@$(FC) -o $@ -L$(LHAPDF_PATH) $^ -lLHAPDF

jewel-gajos_do_glauber-simple: jewel-2.3.0.o temp_glauber_medium_simple.o glauber-medium-simple.o pythia6425mod-lhapdf6.o meix.o
	@$(FC) -o $@ $^ $(ROOTINC) $(ROOTLIB) $(LDLIB) -lstdc++ -L$(LHAPDF_PATH) -lLHAPDF

%.o : %.C
	@echo "compiling $<... [$@]"
	@$(CC) $(CCFLAGS) -c -o $@ $< $(ROOTINC) 

%.o : %.f
	@echo "compiling $<... [$@]"
	@$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -fv *.o *.exe *.out

test:
	@echo "Test 1"
	@echo $(ROOTINC)
	@echo $(ROOTLIB)
	@echo $(LDLIB)

.PHONY: all clean test
