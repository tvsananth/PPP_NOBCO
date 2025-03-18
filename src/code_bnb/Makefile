SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

# ---------------------------------------------------------------------
# Compiler selection
# ---------------------------------------------------------------------

CCC = g++ -O3

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
CCFLAGS = $(CCOPT)

all: leblanc_solver

leblanc_solver: solve_leblanc.o selfish_routing.o fr_solver.o
	$(CCC) $(CCFLAGS) -o leblanc_solver solve_leblanc.o selfish_routing.o fr_solver.o

solve_leblanc.o: ./solve_leblanc.cpp
	$(CCC) -c $(CCFLAGS) ./solve_leblanc.cpp

selfish_routing.o: ./selfish_routing.cpp ./selfish_routing.h
	$(CCC) -c $(CCFLAGS) ./selfish_routing.cpp

fr_solver.o: ./fr_solver.cpp ./fr_solver.h
	$(CCC) -c $(CCFLAGS) ./fr_solver.cpp

clean:
	rm -f *.o leblanc_solver
