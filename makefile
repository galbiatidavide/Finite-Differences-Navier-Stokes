# Compiler
#CXX = g++
CXX = mpicxx
CXXFLAGS = -Wall -O2 -std=c++14

# Paths
INCLUDE_DIR = include
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin


# Libraries
PETSC_DIR ?= $(shell echo $$PETSC_DIR)
PETSC_ARCH ?= $(shell echo $$PETSC_ARCH)
VTK_DIR ?= $(shell echo $$mkVtkPrefix)
LIBS = -L$(PETSC_DIR)/lib -lpetsc -L$(VTK_DIR)/lib 
LIBS = -L/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/petsc/3.15.1/lib -lpetsc \
       -L/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/vtk/9.0.3/lib \
       -lvtkCommonCore-9.0 -lvtkIOCore-9.0 -lvtkFiltersCore-9.0 \
       -lvtkCommonDataModel-9.0 -lvtkIOXML-9.0 -lvtkRenderingCore-9.0 \
       -lvtkCommonExecutionModel-9.0 -lvtkIOGeometry-9.0


# Include directories
CXXFLAGS += -I$(INCLUDE_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include/$(PETSC_ARCH) -I$(VTK_DIR)/include/vtk-9.0

# Files to compile
SRCS = $(SRC_DIR)/main.cpp $(SRC_DIR)/Parabolic.cpp 

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
