ifeq ($(MAKECMDGOALS),mpi)
	CXX=mpic++
endif

all: varPoint.o varRegion.o particle.o point.o mesh_search.o main.o
	$(CXX) $^ -Lclp/build/lib/ -lCoinUtils -llapack -lbz2 -lz -lm libspatialindex/build/lib/libspatialindex.a -o main -g -std=c++1y -Wall -Wextra $(CPPFLAGS)

mpi: all

main.o: main.cpp
	$(CXX) $< -c -g -std=c++1y -Wall -Wextra $(CPPFLAGS)

%.o: %.cpp %.h
	$(CXX) $< -c -g -std=c++1y -Wall -Wextra $(CPPFLAGS)
