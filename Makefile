CXX=mpic++
CCFLAGS=-Wall -O3
LDFLAGS= -lpthread -lm

GTEST_DIR = ./gtest
CCFLAGS += -I.

test: CARMA_test
	echo "### TESTING WITH 4 PROCESSES ###"; mpirun -np 4 ./CARMA_test

CARMA_test: CARMA_test.o CARMA.o gtest-all.o
	$(CXX) -o $@ $^ $(LDFLAGS)

gtest-all.o : $(GTEST_DIR)/gtest-all.cc $(GTEST_DIR)/gtest.h
	$(CXX) $(CCFLAGS) -c $(GTEST_DIR)/gtest-all.cc

%.o: %.cpp %.h
	$(CXX) $(CCFLAGS) -c $<

%.o: %.cpp
	$(CXX) $(CCFLAGS) -c $<