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

# Directories
SRCDIR = src
INCDIR = include
BINDIR = bin
OBJDIR = obj

# Source files and target
MAIN = $(SRCDIR)/main_problem.cpp
SRCS = $(filter-out $(MAIN), $(wildcard $(SRCDIR)/*.cpp))  # All other .cpp files
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS)) $(OBJDIR)/main_problem.o
TARGET = $(BINDIR)/my_program

# Rules
all: $(TARGET)

$(TARGET): $(OBJS)
	mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(PETSC_LIBS) $(VTK_LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(PETSC_CFLAGS) $(VTK_CFLAGS) -I$(INCDIR) -c $< -o $@

clean:
	rm -rf $(OBJDIR) $(BINDIR)

.PHONY: all clean
