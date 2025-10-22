# ==========================================================
# Makefile for MPI 1D Euler Solver with Gatherv
# ==========================================================

CXX      = mpicxx
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra
TARGET   = euler_mpi
SRC      = euler_mpi.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

run:
	mpirun -np 4 ./$(TARGET)

clean:
	rm -rf $(TARGET) *.o *.dat *.png logs*

.PHONY: all clean run
