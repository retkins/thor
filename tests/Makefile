
CFLAGS = -O3 -march=native 
CFILES = ../src/biotsavart.c ../src/octree.c ../src/utils.c

wire: 
	gcc $(CFLAGS) test_wire.c ../src/biotsavart.c ../src/octree.c ../src/utils.c -o test_wire.o -lm
	./test_wire.o

loop:
	gcc $(CFLAGS) test_loop.c $(CFILES) -o test_loop.o -lm
	# Crossover point at n=1000
	./test_loop.o 100000 100000 1 0.1

perftest:
	gcc -O3 -g -march=native -mno-avx512f \
	../src/biotsavart.c ../src/octree.c ../src/utils.c \
	test_loop.c -o test_loop_perf.o -lm 
	valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=cg.out ./test_loop_perf.o 1000 1000 1 0.1
	valgrind --leak-check=yes ./test_loop_perf.o 1000 1000 1 0.1