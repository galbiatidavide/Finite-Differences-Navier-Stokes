# Finite Differences Navier-Stokes Solver

This project implements a solver based on finite differences for the Navier-Stokes equations. The solver is designed to handle multiple physics problems, with the final goal of solving Brinkman flow for a generic geometry. It allows computations on a default cubic domain as well as on an input mesh using the Brinkman penalization method.

Our multiphysics solver is capable of solving a variety of problems, including Navier-Stokes and parabolic problems. This flexibility enables a broad range of simulations in fluid dynamics.

## **Requirements**

To compile and run the project, make sure you have the following tools and libraries installed and configured:

- **Compiler**: GCC with C++20 support
- **Libraries**:
  - PETSc
  - VTK
  - Eigen
- **Environment modules:** https://github.com/pcafrica/mk

---

## **Installation**

1. **Load the required modules**:  
   Before compiling or running the program, load the necessary modules. Ensure that your environment is correctly set up by running:
   ```bash
   source /u/sw/etc/bash.bashrc
   module load gcc-glibc
   module load eigen
   module load vtk
   module load petsc
   ```
   
2. **Set config_problem.hpp**: This file contains the configuration for the problem to solve. User can set their parameters and set more flow variables if needed. This file is well documented and should receive careful consideration by the user. Current configuration is set for an analytical solution of of 3D Navier-Stokes Equations with no Brinkman Penalization (cfr. Doxygen)
    After loading the modules, compile the project according to the chosen by running:


3. **Compile the project**:
    After loading the modules, compile the project according to the chosen by running:
    ```bash
    make clean
    make FILES="navier_stokes.cpp" #for Navier-Stokes and Brinkman problem
    make FILES="stokes.cpp" #for Navier-Stokes and Brinkman problem
    make FILES="inviscid_euler.cpp" #for Navier-Stokes and Brinkman problem
    make FILES="advection_diffusion.cpp" #for Navier-Stokes and Brinkman problem
    make FILES="parabolic.cpp" #for parabolic problems

    ```

4. **Usage**
    To run the solver, execute the following command:
    ```bash
    ./bin/main
    ```
    For the Brinkman solver, the name of the `.stl` file must be specified. An example mesh is provided. To run the solver on the example mesh, execute the following command:
    ```bash
    ./bin/main caroitd.stl
    ```
    The program will compute the solution and store the results in the `results` directory.

## Simulation Results

<table>
  <tr>
    <td align="center">
      <strong>Velocity Magnitude for Navier-Stokes incompressible flow at Re=1</strong><br>
      <img src="graphic_examples/magnitude_Re_1.gif" width="300">
    </td>
    <td align="center">
      <strong>Velocity Magnitude for Navier-Stokes incompressible flow at Re=2000</strong><br>
      <img src="graphic_examples/magnitude_Re_2000.gif" width="300">
    </td>
  </tr>
  <tr>
    <td align="center">
      <strong>Parabolic Flow with μ=10</strong><br>
      <img src="graphic_examples/parabolic_mu_10.gif" width="300">
    </td>
    <td align="center">
      <strong>Brinkman Flow Simulation Re=200</strong><br>
      <img src="graphic_examples/brinkman_Re_200_dt_1e-3.gif" width="300">
    </td>
  </tr>
</table>

These examples showcase the solver's ability to handle various flow regimes and conditions.

---

Developed by **Davide Galbiati** and **Alessandra Gotti**
