all: varPoint.o varRegion.o
	g++ main.cpp varRegion.o varPoint.o libspatialindex/libspatialindex.a -o main -g -std=c++11 -Wall -Wextra

varPoint.o: varPoint.cpp varPoint.h
	g++ varPoint.cpp -c -g -std=c++11 -Wall -Wextra

varRegion.o: varRegion.cpp varRegion.h
	g++ varRegion.cpp -c -g -std=c++11 -Wall -Wextra
