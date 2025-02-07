MPICXX ?= $(shell which mpicxx)
CXX = $(MPICXX)
CXXFLAGS = -Wall -O3 -std=c++20

# Paths
INCLUDE_DIR = include
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
RESULTS_DIR = results

# Libraries
PETSC_DIR ?= $(shell echo $$PETSC_DIR)
PETSC_ARCH ?= $(shell echo $$PETSC_ARCH)
VTK_DIR ?= $(shell echo $$mkVtkPrefix)
LIBS = -L/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/petsc/3.15.1/lib -lpetsc \
       -L/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/vtk/9.0.3/lib \
       -lvtkCommonCore-9.0 -lvtkIOCore-9.0 -lvtkFiltersCore-9.0 \
       -lvtkCommonDataModel-9.0 -lvtkIOXML-9.0 -lvtkRenderingCore-9.0 \
       -lvtkCommonExecutionModel-9.0 -lvtkIOGeometry-9.0

# Include directories
CXXFLAGS += -I$(INCLUDE_DIR) -I$(PETSC_DIR)/include -I$(PETSC_DIR)/include/$(PETSC_ARCH) -I$(VTK_DIR)/include/vtk-9.0 -I.

# Source and object files
SRC_ALL = $(wildcard $(SRC_DIR)/*.cpp)
SRC_MANDATORY = $(SRC_DIR)/main.cpp $(SRC_DIR)/utils.cpp
SRC_SELECTED = $(if $(FILES),$(addprefix $(SRC_DIR)/, $(FILES)),$(SRC_ALL))

ifneq ($(filter $(SRC_DIR)/navier_stokes.cpp, $(SRC_SELECTED)),)
    SRC_SELECTED += $(SRC_DIR)/parabolic.cpp $(SRC_DIR)/transport.cpp $(SRC_DIR)/poisson.cpp
    CXXFLAGS += -DCOMPILE_NAVIER_STOKES
endif

ifneq ($(filter $(SRC_DIR)/inviscid_euler.cpp, $(SRC_SELECTED)),)
    SRC_SELECTED += $(SRC_DIR)/poisson.cpp $(SRC_DIR)/transport.cpp
    CXXFLAGS += -DCOMPILE_EULER
endif

ifneq ($(filter $(SRC_DIR)/stokes.cpp, $(SRC_SELECTED)),)
    SRC_SELECTED += $(SRC_DIR)/poisson.cpp $(SRC_DIR)/parabolic.cpp
    CXXFLAGS += -DCOMPILE_STOKES
endif

ifneq ($(filter $(SRC_DIR)/advection_diffusion.cpp, $(SRC_SELECTED)),)
    SRC_SELECTED += $(SRC_DIR)/transport.cpp $(SRC_DIR)/parabolic.cpp
    CXXFLAGS += -DCOMPILE_ADVECTION_DIFFUSION
endif

ifneq ($(filter $(SRC_DIR)/parabolic.cpp, $(SRC_SELECTED)),)
    SRC_SELECTED += $(SRC_DIR)/parabolic.cpp
    CXXFLAGS += -DCOMPILE_PARABOLIC
endif


SRC_FINAL = $(sort $(SRC_MANDATORY) $(SRC_SELECTED))
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FINAL))

# Executable name
EXEC = $(BIN_DIR)/main

# Build executable
$(EXEC): $(OBJS) | $(BIN_DIR) $(RESULTS_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

# Build object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create directories if not exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(RESULTS_DIR):
	mkdir -p $(RESULTS_DIR)

# Clean rule
clean:
	rm -rf $(OBJ_DIR)/*.o $(EXEC)

.PHONY: clean
