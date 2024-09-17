# Compiler
CXX = g++
CXXFLAGS = -Wall -O2 -std=c++11

# Paths
INCLUDE_DIR = include
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Libraries (Add PETSc or other libraries here)
# Use environment variables set by the module
PETSC_DIR ?= $(shell echo $$PETSC_DIR)
PETSC_ARCH ?= $(shell echo $$PETSC_ARCH)
LIBS = -L$(PETSC_DIR)/lib -lpetsc

# Include directories
CXXFLAGS += -I$(INCLUDE_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include/$(PETSC_ARCH)

# Files to compile (only include the necessary files)
SRCS = $(SRC_DIR)/main.cpp  \
	   $(SRC_DIR)/Grid.cpp

# Corresponding object files
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Executable name
EXEC = $(BIN_DIR)/main

# Target: Build executable
$(EXEC): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

# Rule to create the object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Directories for the binaries and objects
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean rule to remove compiled files
clean:
	rm -rf $(OBJ_DIR)/*.o $(EXEC)

.PHONY: clean

