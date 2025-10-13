
CC = GCC 
CFLAGS = -O3 -march=native 
CFILES = src/main.c 

all:
	$(CC) $(CFLAGS) $(CFILES) -lm -o main.o
	./main.o