all: varPoint.o varRegion.o particle.o point.o
	g++ main.cpp varRegion.o varPoint.o particle.o point.o libspatialindex/libspatialindex.a -o main -g -std=c++11 -Wall -Wextra

varPoint.o: varPoint.cpp varPoint.h
	g++ varPoint.cpp -c -g -std=c++11 -Wall -Wextra

varRegion.o: varRegion.cpp varRegion.h
	g++ varRegion.cpp -c -g -std=c++11 -Wall -Wextra

particle.o: particle.cpp particle.h
	g++ particle.cpp -c -g -std=c++11 -Wall -Wextra

point.o: point.cpp point.h
	g++ point.cpp -c -g -std=c++11 -Wall -Wextra
