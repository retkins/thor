
CC = gcc
CFLAGS = -O3 -march=native 
CFILES = src/*.c

all:
	$(CC) $(CFLAGS) $(CFILES) -lm -o main.o
	./main.o

macos:
	gcc-15 $(CFLAGS) $(CFILES) -g -fopenmp -fopt-info-vec-all=vec.log -lm -o main.o
	./main.o