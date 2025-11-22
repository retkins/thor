# Make file for thor tests
# (c) 2025 ryan@freestatelabs
#

# --- 
# Compiler options
# ---

CFLAGS = -O3 -march=native 
CFILES = ../src/*.c

# ---
# Make targets 
# ---

# Test current-carrying loop and/or solenoid
loop:
	gcc $(CFLAGS) test_loop.c $(CFILES) -o test_loop.o -lm
	./test_loop.o 10000 10000 1 0.005

# Test simple solenoid
solenoid: 
	gcc $(CFLAGS) test_solenoid.c $(CFILES) -o test_solenoid.o -lm
	./test_solenoid.o 10000 10000 1 0.1


# Use valgrind to check for memory leaks and cache misses
perftest:
	gcc -O3 -g -march=native -mno-avx512f \
	../src/biotsavart.c ../src/octree.c ../src/utils.c \
	test_loop.c -o test_loop_perf.o -lm 
	valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=cg.out ./test_loop_perf.o 1000 1000 1 0.1
	valgrind --leak-check=yes ./test_loop_perf.o 1000 1000 1 0.1

# (not used)
wire: 
	gcc $(CFLAGS) test_wire.c ../src/biotsavart.c ../src/octree.c ../src/utils.c -o test_wire.o -lm
	./test_wire.o