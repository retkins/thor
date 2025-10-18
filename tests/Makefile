
CFLAGS = -O3 -march=native 

wire: 
	gcc $(CFLAGS) test_wire.c ../src/biotsavart.c ../src/octree.c ../src/utils.c -o test_wire.o -lm
	./test_wire.o