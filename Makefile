# Useful "automatic variables" for makefiles.
# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.

# Compiler and flags needed for calling it.
CC = clang++
CFLAGS = -c -std=c++11 -stdlib=libc++

# Main executables.
TARGETS = gintsieve gintsieve_test

# Variables with some relevant files.
CORE = src/BaseSieve.cpp src/OctantSieve.cpp include/BaseSieve.hpp include/OctantSieve.hpp

EXTENDED = src/BaseSieve.cpp src/OctantSieve.cpp src/OctantDonutSieve.cpp \
		   include/BaseSieve.hpp include/OctantSieve.hpp include/OctantDonutSieve.hpp

EVERYTHING = src/BaseSieve.cpp src/OctantSieve.cpp src/OctantDonutSieve.cpp \
	   	     src/BlockSieve.cpp src/BlockDonutSieve.cpp src/SectorSieve.cpp \
		     include/BaseSieve.hpp include/OctantSieve.hpp include/OctantDonutSieve.hpp \
		     include/BlockSieve.hpp include/BlockDonutSieve.hpp include/SectorSieve.hpp

# All object files but main.o and tests.o
OBJECTS = opt/BaseSieve.o opt/OctantSieve.o opt/OctantDonutSieve.o \
          opt/BlockSieve.o opt/BlockDonutSieve.o opt/SectorSieve.o


# Telling make to compile every object.
.PHONY: all
all: $(OBJECTS) opt $(TARGETS)

# Make opt/ directory if it doesn't exist.
$(shell mkdir -p opt/)

# Doing all the compiling
opt/BaseSieve.o: src/BaseSieve.cpp include/BaseSieve.hpp
	$(CC) $(CFLAGS) src/BaseSieve.cpp -o $@

opt/OctantSieve.o: $(CORE)
	$(CC) $(CFLAGS) src/OctantSieve.cpp -o $@

opt/OctantDonutSieve.o: $(EXTENDED)
	$(CC) $(CFLAGS) src/OctantDonutSieve.cpp -o $@

opt/BlockSieve.o: $(CORE) src/BlockSieve.cpp include/BlockSieve.hpp
	$(CC) $(CFLAGS) src/BlockSieve.cpp -o $@

opt/BlockDonutSieve.o: $(EXTENDED) src/BlockDonutSieve.cpp include/BlockDonutSieve.hpp
	$(CC) $(CFLAGS) src/BlockDonutSieve.cpp -o $@

opt/SectorSieve.o: $(EXTENDED) src/SectorSieve.cpp include/SectorSieve.hpp
	$(CC) $(CFLAGS) src/SectorSieve.cpp -o $@

opt/main.o: $(EVERYTHING) src/main.cpp
	$(CC) $(CFLAGS) src/main.cpp -o $@

opt/test.o: $(EVERYTHING) src/test.cpp
	$(CC) $(CFLAGS) src/test.cpp -o $@

gintsieve: $(OBJECTS) opt/main.o
	$(CC) -o $@ $^

gintsieve_test: $(OBJECTS) opt/test.o
	$(CC) -o $@ $^


.PHONY: clean
clean:
	rm -r opt/ gintsieve gintsieve_test

.PHONY: test
test:
	./gintsieve_test