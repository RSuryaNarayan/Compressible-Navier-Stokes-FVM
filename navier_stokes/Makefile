# Makefile for gridGen

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -O3

# macOS-specific settings
BREW_PREFIX := $(shell brew --prefix)
OPENMP_FLAGS = -Xpreprocessor -fopenmp
OPENMP_INCLUDES = -I$(BREW_PREFIX)/opt/libomp/include
OPENMP_LIBS = -L$(BREW_PREFIX)/opt/libomp/lib -lomp

# Source files
SRC_RAMP = gridGenRamp.cpp
SRC_CYLINDER = gridGenCylinder.cpp
SRC_READER = readGrid.cpp
SRC_TEST = testSparseMatrixInv.cpp

# Executable names
TARGET_RAMP = gridGenRamp
TARGET_CYLINDER = gridGenCylinder
TARGET_READER = readGrid
# TARGET_TEST = testSparseMatrixInv

# Default target
all: $(TARGET_RAMP) $(TARGET_CYLINDER) $(TARGET_READER) $(TARGET_TEST)

# Build the executables
$(TARGET_RAMP): $(SRC_RAMP)
	$(CXX) $(CXXFLAGS) -I/opt/homebrew/opt/eigen/include/eigen3 -o $(TARGET_RAMP) $(SRC_RAMP)

$(TARGET_CYLINDER): $(SRC_CYLINDER)
	$(CXX) $(CXXFLAGS) -I/opt/homebrew/opt/eigen/include/eigen3 -o $(TARGET_CYLINDER) $(SRC_CYLINDER)

$(TARGET_READER): $(SRC_READER)
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(OPENMP_INCLUDES) $(OPENMP_LIBS) -I/opt/homebrew/opt/eigen/include/eigen3 -o $(TARGET_READER) $(SRC_READER)

# $(TARGET_TEST): $(SRC_TEST)
# 	$(CXX) $(CXXFLAGS) -I/opt/homebrew/opt/eigen/include/eigen3 -o $(TARGET_TEST) $(SRC_TEST)

# Clean up build files
clean:
	rm -f $(TARGET_RAMP) $(TARGET_CYLINDER) $(TARGET_READER) $(TARGET_TEST) *.csv