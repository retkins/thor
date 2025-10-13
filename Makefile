
CC = gcc
CFLAGS = -O3 -march=native 
CFILES = src/*.c

all:
	$(CC) $(CFLAGS) $(CFILES) -lm -o main.o
	./main.o