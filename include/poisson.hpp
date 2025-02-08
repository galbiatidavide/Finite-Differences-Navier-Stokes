/// \file
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
#include "utils.hpp"

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
Vec U_up; ///< Velocity component in the x-direction.
Vec V_up; ///< Velocity component in the y-direction.
Vec W_up; ///< Velocity component in the z-direction.

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
poisson_problem(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, DM const & dmGrid_centered, DM const & dmGrid_cent_rich);

/**
 * @brief Constructor to solve stand-alone problem
 */
poisson_problem();

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
PetscErrorCode const manage_pressure(std::optional<std::reference_wrapper<Vec>> U_opt = std::nullopt,
std::optional<std::reference_wrapper<Vec>> V_up_opt = std::nullopt,
std::optional<std::reference_wrapper<Vec>> W_up_opt = std::nullopt, 
std::optional<std::reference_wrapper<Vec>> P_opt = std::nullopt);
/**
 * @brief Assembles the pressure derivative in the x-direction on the staggered grid.
 * @param P Pressure field.
 * @param P_x Updated pressure derivative in the x-direction.
 */
PetscErrorCode const manage_pressure_x(std::optional<std::reference_wrapper<Vec>> P_opt = std::nullopt, std::optional<std::reference_wrapper<Vec>> P_x_opt = std::nullopt);
/**
 * @brief Assembles the pressure derivative in the y-direction on the staggered grid.
 * @param P Pressure field.
 * @param P_y Updated pressure derivative in the y-direction.
 */
PetscErrorCode const manage_pressure_y(std::optional<std::reference_wrapper<Vec>> P_opt = std::nullopt, std::optional<std::reference_wrapper<Vec>> P_y_opt = std::nullopt);
/**
 * @brief Assembles the pressure derivative in the z-direction on the staggered grid.
 * @param P Pressure field.
 * @param P_z Updated pressure derivative in the z-direction.
 */
PetscErrorCode const manage_pressure_z(std::optional<std::reference_wrapper<Vec>> P_opt = std::nullopt, std::optional<std::reference_wrapper<Vec>> P_z_opt = std::nullopt);
/**
 * @brief Destructor to clean up allocated resources. required by PETSc management of native objects.
 */
~poisson_problem();

};

#endif // POISSON_PROBLEM_HPP


