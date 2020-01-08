# Useful "automatic variables" for makefiles.
# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.

CC = g++
CFLAGS = -c -std=c++11 -stdlib=libc++

# Main executable.
target = gintsieve
# Not including main.cpp in these three variables
sources = src/BaseSieve.cpp src/QuadrantSieve.cpp src/OctantSieve.cpp src/DonutSieve.cpp src/SegmentedSieve.cpp
headers = include/BaseSieve.hpp include/QuadrantSieve.hpp include/OctantSieve.hpp include/DonutSieve.hpp include/SegmentedSieve.hpp
objects = opt/BaseSieve.o opt/QuadrantSieve.o opt/OctantSieve.o opt/DonutSieve.o opt/SegmentedSieve.o

.PHONY: all
all: $(objects) opt $(target)

opt/BaseSieve.o: src/BaseSieve.cpp include/BaseSieve.hpp
	$(CC) $(CFLAGS) src/BaseSieve.cpp -o $@

opt/QuadrantSieve.o: src/BaseSieve.cpp include/BaseSieve.hpp src/QuadrantSieve.cpp include/QuadrantSieve.hpp
	$(CC) $(CFLAGS) src/QuadrantSieve.cpp -o $@

opt/OctantSieve.o: src/BaseSieve.cpp include/BaseSieve.hpp src/QuadrantSieve.cpp include/QuadrantSieve.hpp src/OctantSieve.cpp include/OctantSieve.hpp
	$(CC) $(CFLAGS) src/OctantSieve.cpp -o $@

opt/DonutSieve.o: $(sources) $(headers)
	$(CC) $(CFLAGS) src/DonutSieve.cpp -o $@

opt/SegmentedSieve.o: $(sources) $(headers)
	$(CC) $(CFLAGS) src/SegmentedSieve.cpp -o $@

opt/main.o: $(sources) $(headers) src/main.cpp
	$(CC) $(CFLAGS) src/main.cpp -o $@

$(target): $(objects) opt/main.o
	$(CC) -o $@ $^

.PHONY: clean
clean:
	rm opt/*.o gintsieve