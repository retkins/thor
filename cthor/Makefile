# Makefile for thor
# (c) 2025 ryan@freestatelabs
# 
# Tested on MacOS (Apple silicon, M1 Pro) and Ubuntu (Ryzen 9 9950X)
#

# ---
# MODIFY THESE AS NEEDED
# ---

CC := gcc
CFLAGS := -O3 -march=native 
CFILES := src/*.c
BUILD_DIR := build

# --- 
# Make targets
# ---

# Build a shared library to call from Julia 
# build with `make lib`
# Creates a .so file: thor/build/thorlib.so
BUILD_DIR := build
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/thorlib.so: $(CFILES) | $(BUILD_DIR)
	$(CC) -shared $(CFLAGS) -fPIC -o $@ $(CFILES)

lib: $(BUILD_DIR)/thorlib.so


# For running the main program (currently not functional)
all:
	$(CC) $(CFLAGS) $(CFILES) -lm -o main
	./main
	
# Experimenting with some vectorization reports on GCC via homebrew, macos
macos:
	gcc-15 $(CFLAGS) $(CFILES) -g -fopenmp -fopt-info-vec-all=vec.log -lm -o main
	./main