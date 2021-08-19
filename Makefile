CC=g++
CXXFLAGS=-I. -g -march=native -flto -pedantic -Wall
DEPS = sha512.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

kseqalign: kseqalign.o
	$(CC) -o kseqalign kseqalign.o $(CXXFLAGS)

.PHONY: clean

clean:
	rm *.o kseqalign