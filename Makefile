# Useful "automatic variables" for makefiles.
# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.

# Compiler and flags needed for calling it.
CC = clang++
CFLAGS = -std=c++11 -stdlib=libc++ -I include/

# Main executables.
TARGETS = gintsieve ginttest gintmoat

# Variables with some relevant files.
CORE = src/BaseSieve.cpp src/OctantSieve.cpp include/BaseSieve.hpp include/OctantSieve.hpp

EXTENDED = src/BaseSieve.cpp src/OctantSieve.cpp src/OctantDonutSieve.cpp \
		   include/BaseSieve.hpp include/OctantSieve.hpp include/OctantDonutSieve.hpp

EVERYTHING = src/BaseSieve.cpp src/OctantSieve.cpp src/OctantDonutSieve.cpp \
	   	     src/BlockSieve.cpp src/BlockDonutSieve.cpp src/SectorSieve.cpp \
		     include/BaseSieve.hpp include/OctantSieve.hpp include/OctantDonutSieve.hpp \
		     include/BlockSieve.hpp include/BlockDonutSieve.hpp include/SectorSieve.hpp

MOAT = src/OctantMoat.cpp src/SegmentedMoat.cpp src/VerticalMoat.cpp include/Moat.hpp

# All object files from sources in EVERYTHING
OBJECTS = obj/BaseSieve.o obj/OctantSieve.o obj/OctantDonutSieve.o \
          obj/BlockSieve.o obj/BlockDonutSieve.o obj/SectorSieve.o \
          obj/OctantMoat.o obj/SegmentedMoat.o obj/VerticalMoat.o


# Telling make to compile every object.
.PHONY: all
all: $(OBJECTS) obj $(TARGETS)

# Make obj/ directory if it doesn't exist.
$(shell mkdir -p obj/)

# Doing all the compiling
obj/BaseSieve.o: src/BaseSieve.cpp include/BaseSieve.hpp
	$(CC) $(CFLAGS) -c src/BaseSieve.cpp -o $@

obj/OctantSieve.o: $(CORE)
	$(CC) $(CFLAGS) -c src/OctantSieve.cpp -o $@

obj/OctantDonutSieve.o: $(EXTENDED)
	$(CC) $(CFLAGS) -c src/OctantDonutSieve.cpp -o $@

obj/BlockSieve.o: $(CORE) src/BlockSieve.cpp include/BlockSieve.hpp
	$(CC) $(CFLAGS) -c src/BlockSieve.cpp -o $@

obj/BlockDonutSieve.o: $(EXTENDED) src/BlockDonutSieve.cpp include/BlockDonutSieve.hpp
	$(CC) $(CFLAGS) -c src/BlockDonutSieve.cpp -o $@

obj/SectorSieve.o: $(EXTENDED) src/SectorSieve.cpp include/SectorSieve.hpp
	$(CC) $(CFLAGS) -c src/SectorSieve.cpp -o $@

obj/OctantMoat.o: $(CORE) src/OctantMoat.cpp include/Moat.hpp
	$(CC) $(CFLAGS) -c src/OctantMoat.cpp -o $@

obj/VerticalMoat.o: $(CORE) src/VerticalMoat.cpp include/Moat.hpp
	$(CC) $(CFLAGS) -c src/VerticalMoat.cpp -o $@

obj/SegmentedMoat.o: $(CORE) src/SegmentedMoat.cpp include/Moat.hpp
	$(CC) $(CFLAGS) -c src/SegmentedMoat.cpp -o $@

# Building to-be executables into objects
obj/gintsieve.o: $(EVERYTHING) src/gintsieve.cpp
	$(CC) $(CFLAGS) -c src/gintsieve.cpp -o $@

obj/gintmoat.o: $(CORE) $(MOAT) src/gintmoat.cpp
	$(CC) $(CFLAGS) -c src/gintmoat.cpp -o $@

obj/ginttest.o: $(EVERYTHING) $(MOAT) src/ginttest.cpp
	$(CC) $(CFLAGS) -c src/ginttest.cpp -o $@

# Linking objects to executables
gintsieve: $(OBJECTS) obj/gintsieve.o
	$(CC) $(CFLAGS) -o $@ $^

ginttest: $(OBJECTS) obj/ginttest.o
	$(CC) $(CFLAGS) -o $@ $^

gintmoat: $(OBJECTS) obj/gintmoat.o
	$(CC) $(CFLAGS) -o $@ $^


clean:
	rm -r obj/ gintsieve ginttest gintmoat

test:
	./ginttest


pybuild:
	python setup.py sdist bdist_wheel clean --all
	# showing wheel contents
	for f in dist/*.whl; do unzip -l $$f; done
	# showing tar contents
	tar -tvf dist/*.tar.gz
	twine check dist/*

pyinstall:
	python -m pip install .

pytest:
	python -m pytest
