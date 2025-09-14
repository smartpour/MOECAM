CXX = g++
CC = gcc
CXXFLAGS = -std=c++17 -O3 -Wall -fPIC
CFLAGS = -std=c99 -O3 -Wall -fPIC

# Paths
WFG_DIR = ../../../csources/wfg/WFG_1.15
WFG_SOURCES = $(WFG_DIR)/wfg.c $(WFG_DIR)/read.c
WFG_HEADERS = $(WFG_DIR)/wfg.h

# Targets
all: test_wfg_simple libwfg_simple.so

# Compile WFG C sources (core without main function)
wfg_core.o: wfg_core.c $(WFG_HEADERS)
	$(CC) $(CFLAGS) -I$(WFG_DIR) -c wfg_core.c -o wfg_core.o

read.o: $(WFG_DIR)/read.c $(WFG_HEADERS)
	$(CC) $(CFLAGS) -c $(WFG_DIR)/read.c -o read.o

# Compile simple C interface
wfg_simple_interface.o: wfg_simple_interface.c wfg_simple_interface.h $(WFG_HEADERS)
	$(CC) $(CFLAGS) -I$(WFG_DIR) -c wfg_simple_interface.c -o wfg_simple_interface.o

# Create shared library for Python
libwfg_simple.so: wfg_simple_interface.o wfg_core.o read.o
	$(CXX) -shared -fPIC $^ -o $@

# Test executable
test_wfg_simple: test_wfg_simple.cpp wfg_simple_interface.o wfg_core.o read.o
	$(CXX) $(CXXFLAGS) -I$(WFG_DIR) $^ -o $@

# Clean
clean:
	rm -f *.o test_wfg_cpp libwfg_interface.so

# Test run
test: test_wfg_cpp
	./test_wfg_cpp

.PHONY: all clean test
