/******************************************************************************
 *                                                                            *
 *  Project:    Brinkman - Navier-Stokes multiple solver                      *                  
 *  Author:     Dave & Ale                                                    *
 *  Created:    26 October 2023                                               *
 *                                                                            *
 *  Copyright © 2023 Dave & Ale. All rights reserved.                         *
 *                                                                            *
 *  Redistribution and use in source and binary forms, with or without        *
 *  modification, are permitted provided that the following conditions are    *
 *  met:                                                                      *
 *                                                                            *
 *  1. Redistributions of source code must retain the above copyright notice, *
 *     this list of conditions, and the following disclaimer.                 *
 *  2. Redistributions in binary form must reproduce the above copyright      *
 *     notice, this list of conditions, and the following disclaimer in the   *
 *     documentation and/or other materials provided with the distribution.   *
 *                                                                            *
 *  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY EXPRESS OR IMPLIED          *
 *  WARRANTIES, INCLUDING BUT NOT LIMITED TO MERCHANTABILITY OR FITNESS FOR   *
 *  A PARTICULAR PURPOSE. THE AUTHOR SHALL NOT BE HELD LIABLE FOR ANY CLAIMS  *
 *  OR DAMAGES ARISING FROM THE USE OF THIS SOFTWARE.                         *
 *                                                                            *
 ******************************************************************************/
#include "macros.hpp"
#include "utils.hpp"
#include "config_problem.hpp"

using namespace problem_setting;

#ifndef POISSON_PROBLEM_HPP
#define POISSON_PROBLEM_HPP
/**
 * @class poisson_problem
 * @brief Represents a Poisson equation solver for pressure correction in fluid simulations.
 *
 * This class solves the Poisson equation for a variabble like pressure, with homogeneous Neumann boundary conditions.
 * It manages the assembly of matrices, calculation of divergence. In our framework, it is used to solve the pressure correction in the Navier-Stokes equations.
 * For our purposes, enforcing compatibility condition was not required. Beware that for a stand-alone problem, compatible-to-null-bc's source must be provided.
 */
class poisson_problem
{
private:

DM dmGrid_staggered_x; ///< Discretized grid for staggered x-direction.
DM dmGrid_staggered_y; ///< Discretized grid for staggered y-direction.
DM dmGrid_staggered_z; ///< Discretized grid for staggered z-direction.
DM dmGrid_centered;    ///< Discretized grid for centered formulation.
DM dmGrid_cent_rich;   ///< Refined centered grid for pressure computation.

Vec P;   ///< Pressure field.
Vec P_x; ///< Pressure derivative in x-direction.
Vec P_y; ///< Pressure derivative in y-direction.
Vec P_z; ///< Pressure derivative in z-direction.

Mat A;   ///< Matrix for the discretized Poisson equation.

/**
 * @brief Assembles the left-hand side (LHS) matrix for the Poisson problem.
 */
PetscErrorCode const assemble_lhs();
/**
 * @brief Computes the divergence term for velocity fields and outputs on dmGrid_centered_rich grid.
 * @param div Output divergence vector.
 * @param U Velocity component in the x-direction.
 * @param V Velocity component in the y-direction.
 * @param W Velocity component in the z-direction.
 */
PetscErrorCode const assemble_divergence(Vec & div, Vec const & U, Vec const &  V, Vec const & W);
/**
 * @brief Assemble routines for a previously found divergence of the given velocity field and migrates it on the correct dmGrid_centered grid.
 * @param div Output divergence vector.
 * @param U Velocity component in the x-direction.
 * @param V Velocity component in the y-direction.
 * @param W Velocity component in the z-direction.
 */
PetscErrorCode const compute_divergence(Vec & div, Vec const & U_n, Vec const & V_n, Vec const & W_n);
/**
 * @brief Computes the derivative of the pressure field in the x-direction and outputes on dmGrid_centered_rich grid.
 * @param P_x_shifted Output shifted pressure derivative.
 * @param vec Input pressure field.
 */
PetscErrorCode const derive_x_P(Vec & P_x_shifted, Vec const & vec);
/**
 * @brief Computes the derivative of the pressure field in the y-direction and outputes on dmGrid_centered_rich grid.
 * @param P_y_shifted Output shifted pressure derivative.
 * @param vec Input pressure field.
 */
PetscErrorCode const derive_y_P(Vec & P_y_shifted, Vec const & vec);
/**
 * @brief Computes the derivative of the pressure field in the z-direction and outputes on dmGrid_centered_rich grid.
 * @param P_z_shifted Output shifted pressure derivative.
 * @param vec Input vector.
 */
PetscErrorCode const derive_z_P(Vec & P_z_shifted, Vec const & vec);

public:

/**
 * @brief Constructor to initialize the Poisson problem with pre-allocated grids.
 * @param dmGrid_staggered_x Grid for staggered x-direction.
 * @param dmGrid_staggered_y Grid for staggered y-direction.
 * @param dmGrid_staggered_z Grid for staggered z-direction.
 * @param dmGrid_centered Grid for centered formulation.
 * @param dmGrid_cent_rich Grid for refined pressure computation.
 */
poisson_problem(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, DM const & dmGrid_centered, DM const & dmGrid_cent_rich)
    : dmGrid_staggered_x(dmGrid_staggered_x), dmGrid_staggered_y(dmGrid_staggered_y), dmGrid_staggered_z(dmGrid_staggered_z), dmGrid_centered(dmGrid_centered), dmGrid_cent_rich(dmGrid_cent_rich)

{
    DMCreateGlobalVector(dmGrid_centered, &P);
    DMCreateGlobalVector(dmGrid_staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_staggered_z, &P_z);
    DMCreateMatrix(dmGrid_centered, &A);
    assemble_lhs();
}

//sistemo questo costruttore
/*
navier_stokes_problem()
{
    //Allocate the grids
    CreateGrid(&dmGrid_staggered_x, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_staggered_y, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_staggered_z, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_cent_rich, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    //Create parallel vectors
    DMCreateGlobalVector(dmGrid_staggered_x, &U_up);
    DMCreateGlobalVector(dmGrid_staggered_y, &V_up);
    DMCreateGlobalVector(dmGrid_staggered_z, &W_up);
    CreateAnalyticalU(dmGrid_staggered_x, U_up, 0);
    CreateAnalyticalV(dmGrid_staggered_y, V_up, 0);
    CreateAnalyticalW(dmGrid_staggered_z, W_up, 0);
    PetscReal norm_U;
    VecNorm(U_up, NORM_INFINITY, &norm_U);
    std::cout << "Norm U: " << norm_U << std::endl;

    DMCreateGlobalVector(dmGrid_centered, &P);
    DMCreateGlobalVector(dmGrid_staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_staggered_z, &P_z);
    DMCreateMatrix(dmGrid_centered, &A);

    DMCreateGlobalVector(dmGrid_staggered_x, &mask_U);
    DMCreateGlobalVector(dmGrid_staggered_y, &mask_V);
    DMCreateGlobalVector(dmGrid_staggered_z, &mask_W);

    if(brinkmann)
    {
        createMaskU(dmGrid_staggered_x, mask_U, vertices, faces);
        createMaskV(dmGrid_staggered_y, mask_V, vertices, faces);
        createMaskW(dmGrid_staggered_z, mask_W, vertices, faces);
    }
    else {
        VecSet(mask_U, 0.0);
        VecSet(mask_V, 0.0);
        VecSet(mask_W, 0.0);
    }

};
*/

/**
 * @brief Exports simulation results in .vtk format for visualization.
 */
PetscErrorCode const exodus(size_t i);
/**
 * @brief Solves the elliptic Poisson equation with GMRES and multigrid preconditioner.
 * @param U_up Velocity field in the x-direction.
 * @param V_up Velocity field in the y-direction.
 * @param W_up Velocity field in the z-direction.
 * @param P Pressure field.
 */
PetscErrorCode const manage_pressure(Vec const & U_up, Vec const & V_up, Vec const & W_up, Vec const & P);
/**
 * @brief Assembles the pressure derivative in the x-direction on the staggered grid.
 * @param P Pressure field.
 * @param P_x Updated pressure derivative in the x-direction.
 */
PetscErrorCode const manage_pressure_x(Vec const & P, Vec const & P_x);
/**
 * @brief Assembles the pressure derivative in the y-direction on the staggered grid.
 * @param P Pressure field.
 * @param P_y Updated pressure derivative in the y-direction.
 */
PetscErrorCode const manage_pressure_y(Vec const & P, Vec const & P_y);
/**
 * @brief Assembles the pressure derivative in the z-direction on the staggered grid.
 * @param P Pressure field.
 * @param P_z Updated pressure derivative in the z-direction.
 */
PetscErrorCode const manage_pressure_z(Vec const & P, Vec const & P_z);
/**
 * @brief Destructor to clean up allocated resources. required by PETSc management of native objects.
 */
~poisson_problem()
{
    VecDestroy(&P);
    VecDestroy(&P_x);
    VecDestroy(&P_y);
    VecDestroy(&P_z);
    MatDestroy(&A);
    DMDestroy(&dmGrid_staggered_x);
    DMDestroy(&dmGrid_staggered_y);
    DMDestroy(&dmGrid_staggered_z);
    DMDestroy(&dmGrid_centered);
    DMDestroy(&dmGrid_cent_rich);
    std::cout << "Poisson Destructor Called" << std::endl;
}

};

#endif // POISSON_PROBLEM_HPP


