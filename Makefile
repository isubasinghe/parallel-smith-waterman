CC=g++
CXXFLAGS=-I. -O3 -fopenmp -D_GLIBCXX_PARALLEL -march=native -flto -pedantic -Wall
DEPS = sha512.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

kseqalign: kseqalign.o
	$(CC) -o kseqalign kseqalign.o $(CXXFLAGS)

.PHONY: clean

clean:
	rm *.o kseqalign