CXX=mpic++
CCFLAGS=-Wall -O3
LDFLAGS= -lpthread -lm 

GTEST_DIR = ./gtest
CCFLAGS += -I.

test: CARMA_test
	@echo "usage: mpirun -np <#processors> ./CARMA_test <m> <k> <n>"
	@echo "### TESTING WITH 4 PROCESSES, [m,k,n] = [2000,2000,2000]  ###" 
	mpirun -np 4 ./CARMA_test 2000 2000 2000

CARMA_test: CARMA_test.o CARMA.o 
	$(CXX) -g -o $@ $^ $(LDFLAGS) -mkl

#gtest-all.o : $(GTEST_DIR)/gtest-all.cc $(GTEST_DIR)/gtest.h
#	$(CXX) $(CCFLAGS) -c $(GTEST_DIR)/gtest-all.cc

%.o: %.cpp %.h
	$(CXX) -g $(CCFLAGS) -c $<

%.o: %.cpp
	$(CXX) -g $(CCFLAGS) -c $<
