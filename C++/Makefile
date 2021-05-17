CXX = g++
SOP = $(WINDOWS)
WINDOWS = -lpsapi
LINUX = 
INCLUDE_EIGEN_W = -I C:/Users/andrea1/Desktop/eigen-3.3.9
INCLUDE_EIGEN_L = -I ~/Scrivania/eigen-3.3.9
main.exe: main.o
	$(CXX) main.o $(SOP) -fopenmp -o main.exe

main.o: main.cpp
	$(CXX) -c main.cpp $(INCLUDE_EIGEN_W)  -o main.o

.PHONY: clean
clean:
	rm -r *.o *.exe
