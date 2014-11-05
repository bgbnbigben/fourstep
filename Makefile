all: varPoint.o varRegion.o particle.o point.o mesh_search.o
	g++ main.cpp $^ libspatialindex/libspatialindex.a -o main -g -std=c++1y -Wall -Wextra

%.o: %.cpp %.h
	g++ $< -c -g -std=c++1y -Wall -Wextra
