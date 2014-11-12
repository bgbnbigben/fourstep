ifeq ($(MAKECMDGOALS),mpi)
	CXX=mpic++
endif

all: varPoint.o varRegion.o particle.o point.o mesh_search.o main.o
	$(CXX) $^ -Lclp/build/lib/ -lCoinUtils -llapack -lbz2 -lz -lm libspatialindex/build/lib/libspatialindex.a -o main -ggdb -std=c++1y -Wall -Wextra $(CPPFLAGS) -fopenmp

mpi: all

main.o: main.cpp
	$(CXX) $< -c -ggdb -std=c++1y -Wall -Wextra $(CPPFLAGS) -fopenmp

%.o: %.cpp %.h
	$(CXX) $< -c -ggdb -std=c++1y -Wall -Wextra $(CPPFLAGS) -fopenmp
