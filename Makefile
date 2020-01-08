# Useful "automatic variables":
# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.

CC = g++
CFLAGS = -c -std=c++11 -stdlib=libc++

TARGET = gintsieve
SOURCES  := src/*.cpp
INCLUDES := include/*.hpp
OBJECTS  := bin/BaseSieve.o bin/QuadrantSieve.o bin/OctantSieve.o bin/DonutSieve.o bin/SegmentedSieve.o bin/main.o

.PHONY: all
all: $(OBJECTS) $(TARGET)

bin/BaseSieve.o: src/BaseSieve.cpp include/BaseSieve.hpp
	$(CC) $(CFLAGS) src/BaseSieve.cpp -o $@

bin/QuadrantSieve.o: src/BaseSieve.cpp include/BaseSieve.hpp src/QuadrantSieve.cpp include/QuadrantSieve.hpp
	$(CC) $(CFLAGS) src/QuadrantSieve.cpp -o $@

bin/OctantSieve.o: src/BaseSieve.cpp include/BaseSieve.hpp src/QuadrantSieve.cpp include/QuadrantSieve.hpp src/OctantSieve.cpp include/OctantSieve.hpp
	$(CC) $(CFLAGS) src/OctantSieve.cpp -o $@

bin/DonutSieve.o: $(INCLUDES) $(SOURCES)
	$(CC) $(CFLAGS) src/DonutSieve.cpp -o $@

bin/SegmentedSieve.o: $(INCLUDES) $(SOURCES)
	$(CC) $(CFLAGS) src/SegmentedSieve.cpp -o $@

bin/main.o: $(INCLUDES) $(SOURCES)
	$(CC) $(CFLAGS) src/main.cpp -o $@

$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGET) $(OBJECTS)

.PHONY: clean
clean:
	rm bin/*.o gintsieve