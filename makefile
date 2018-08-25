# ML Lab 1 makefile
# Xhanti Singatha
CC=g++


pca: PCAAlgorithm.o 
	$(CC)  PCAAlgorithm.o -std=c++11 -o pca

PCAAlgorithm.o: PCAAlgorithm.cpp
	$(CC) PCAAlgorithm.cpp -std=c++11 -c

clean:
	@rm -f *.o
