# ----- Make Macros -----

CXX     = CC
CXXFLAGS = -std=c++11 -march=native -O3 -fopenmp -Wall -Wextra -Wstrict-aliasing -Wshadow -Wpedantic -Wno-unused-result -ffast-math

LD_FLAGS = -lgomp 

TARGETS = mom
OBJECTS = mom.o solve.o integrate.o setup_mlfma.o bicgs.o farfield.o
INCLUDES = -Iinclude

# ----- Make Rules -----
all:    $(TARGETS)

%.o : %.cpp
	${CXX} ${CXXFLAGS} ${INCLUDES} $^ -c -o $@

mom: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LD_FLAGS)

clean:
	rm -f $(TARGETS) *.o *.txt *.bin core
