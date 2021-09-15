CC=g++
CXX=g++
CXXFLAGS=-I. -fopenmp -g -static-libstdc++ -O3 -fno-exceptions -fno-rtti -march=native -funroll-loops -pedantic -Wall
DEPS = sha512.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: kseqalign.o
	$(CC) -o kseqalign kseqalign.o $(CXXFLAGS)

.PHONY: clean

clean:
	rm *.o kseqalign
