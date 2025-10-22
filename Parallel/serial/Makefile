# ==========================================================
# Makefile for 1D Euler Serial Solver
# ==========================================================

# Compiler and flags
CXX       = g++
CXXFLAGS  = -O2 -std=c++17 -Wall -Wextra

# Executable name
TARGET    = euler_serial

# Source files
SRC       = euler_serial.cpp

# ==========================================================
# Build rules
# ==========================================================

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o solution_*.dat *.png

# Phony targets (not real files)
.PHONY: all clean run
