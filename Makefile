all:
	g++ varPoint.cpp -c -g -std=c++11 -Wall -Wextra
	g++ varRegion.cpp -c -g -std=c++11 -Wall -Wextra
	g++ main.cpp -o main -g -std=c++11 -Wall -Wextra
