# Compiler and flags
CXX = mpicxx
CXXFLAGS = -O3 -Wall -std=c++17

# PETSc configuration (using pkg-config)
PETSC_CFLAGS = $(shell pkg-config --cflags PETSc)
PETSC_LIBS = $(shell pkg-config --libs PETSc)

# VTK configuration
VTK_INC = $(mkVtkInc)
VTK_LIB = $(mkVtkLib)
VTK_CFLAGS = -I$(VTK_INC)
VTK_LIBS = /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/vtk/9.0.3/lib/libvtkCommonCore-9.0.so \
           /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/vtk/9.0.3/lib/libvtkCommonDataModel-9.0.so \
           /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/vtk/9.0.3/lib/libvtkIOCore-9.0.so \
           /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/vtk/9.0.3/lib/libvtkIOLegacy-9.0.so \
           /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/vtk/9.0.3/lib/libvtkCommonExecutionModel-9.0.so \
    	   /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/vtk/9.0.3/lib/libvtkIOGeometry-9.0.so


# Source files and target
MAIN = main_problem.cpp
SRCS = $(filter-out $(MAIN), $(wildcard *.cpp))  # All other .cpp files
OBJS = $(SRCS:.cpp=.o)
TARGET = my_program

# Rules
all: $(TARGET)

$(TARGET): $(OBJS) $(MAIN:.cpp=.o)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(PETSC_LIBS) $(VTK_LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(PETSC_CFLAGS) $(VTK_CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(MAIN:.cpp=.o) $(TARGET)
