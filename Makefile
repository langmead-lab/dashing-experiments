CXX?=g++

all: bin/pairwise_benchmark bin/testhash bin/dsexp
BONSAI_DIR=bonsai
FLAGS=-O3 -fopenmp -march=native -std=c++14 -Wno-ignored-attributes
INC= -I $(BONSAI_DIR) -I $(BONSAI_DIR)/hll -I $(BONSAI_DIR)/hll/vec -I $(BONSAI_DIR)/hll/libpopcnt -I$(BONSAI_DIR)/circularqueue -I $(BONSAI_DIR)/clhash/include

bin/dsexp: dsexp/dsexp.cpp
	$(CXX) $(FLAGS) -o bin/dsexp dsexp/dsexp.cpp $(INC)
bin/pairwise_benchmark: accuracy/pairwise_benchmark.cpp
	$(CXX) $(FLAGS) -o bin/pairwise_benchmark accuracy/pairwise_benchmark.cpp  bonsai/clhash/src/clhash.c $(INC) -lz #I $(BONSAI_DIR) -I $(BONSAI_DIR)/hll -I $(BONSAI_DIR)/hll/vec
bin/testhash: hash/testhash.cpp
	$(CXX) $(FLAGS) -o bin/testhash hash/testhash.cpp bonsai/clhash/src/clhash.c  $(INC)
