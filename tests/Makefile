# Make file for thor tests
# (c) 2025 ryan@freestatelabs
#

# --- 
# Compiler options
# ---

CFLAGS = -O3 -march=native -Wall
CFILES = ../src/*.c

# ---
# Make targets 
# ---

# Test current-carrying loop and/or solenoid
loop:
	gcc $(CFLAGS) test_loop.c $(CFILES) -o build/test_loop -lm
	./build/test_loop 10000 10000 1 0.05

# Test simple solenoid
solenoid: 
	gcc $(CFLAGS) test_solenoid.c $(CFILES) -o build/test_solenoid -lm
	./build/test_solenoid 100000 100000 1 0.1

# Test dense solenoid 
dense_solenoid: 
	gcc $(CFLAGS) test_dense_solenoid.c $(CFILES) -o build/test_dense_solenoid -lm 
	./build/test_dense_solenoid

# Test self fieldsdense solenoid 
dense_self: 
	gcc $(CFLAGS) test_self_dense_solenoid.c $(CFILES) -o build/test_self_dense_solenoid -lm 
	./build/test_self_dense_solenoid

# Test direct simd 
simd: 
	gcc $(CFLAGS) test_simd.c $(CFILES) -o build/test_simd -lm 
	./build/test_simd


# Use valgrind to check for memory leaks and cache misses
perftest:
	gcc -O3 -g -march=native -mno-avx512f \
	../src/biotsavart.c ../src/octree.c ../src/utils.c \
	test_loop.c -o build/test_loop_perf -lm 
	valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file=cgut ./build/test_loop_perf 1000 1000 1 0.1
	valgrind --leak-check=yes ./build/test_loop_perf 1000 1000 1 0.1

# (not used)
wire: 
	gcc $(CFLAGS) test_wire.c ../src/biotsavart.c ../src/octree.c ../src/utils.c -o build/test_wire -lm
	./build/test_wire

mesh_solenoid: $(CFILES) test_mesh_solenoid.c 
	mkdir -p build
	$(CC) $(CFLAGS) $(CFILES) test_mesh_solenoid.c -o build/mesh_solenoid -lm
	./build/mesh_solenoid
