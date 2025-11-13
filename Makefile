
CC := gcc-15
CFLAGS := -O3 -march=native 
CFILES := src/*.c
BUILD_DIR := build

all:
	$(CC) $(CFLAGS) $(CFILES) -lm -o main.o
	./main.o

macos:
	gcc-15 $(CFLAGS) $(CFILES) -g -fopenmp -fopt-info-vec-all=vec.log -lm -o main.o
	./main.o

# Build a shared library to call from Julia 
# build with `make lib`
# Creates a .so file: thor/build/thorlib.so
BUILD_DIR := build
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/thorlib.so: $(CFILES) | $(BUILD_DIR)
	$(CC) -shared $(CFLAGS) -fPIC -o $@ $(CFILES)

lib: $(BUILD_DIR)/thorlib.so