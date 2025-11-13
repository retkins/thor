
CFLAGS = -O3 -march=native 
CFILES = ../src/biotsavart.c ../src/octree.c ../src/utils.c

wire: 
	gcc $(CFLAGS) test_wire.c ../src/biotsavart.c ../src/octree.c ../src/utils.c -o test_wire.o -lm
	./test_wire.o

loop:
	gcc $(CFLAGS) test_loop.c $(CFILES) -o test_loop.o -lm
	./test_loop.o 10000 10000 1
