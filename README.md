# Finite Differences Navier-Stokes Solver

This project implements a solver based on finite differences for the Navier-Stokes equations.
The solver can be used on a default cubic geometry as well as on an inputed .stl orginal mesh for the brinkmann solver.

## **Requirements**

To compile and run the project, make sure you have the following tools and libraries installed and configured:

- **Compiler**: GCC with C++20 support
- **Libraries**:
  - PETSc
  - VTK
  - Eigen

---

## **Installation**

1. **Load the required modules**:  
   Before compiling or running the program, load the necessary modules. Ensure that your environment is correctly set up by running:
   ```bash
   source /u/sw/etc/bash.bashrc
   module load gcc-glibc
   module load Eigen
   module load vtk
   module load petsc

2. **Compile the project**:
    After loading the modules, compile the project by running:
    ```bash
    make 

3. **Usage**
    To run the solver, execute the following command:
    ```bash
    ./bin/main 
    ```
    The program will compute the solution and store the results in the results directory

