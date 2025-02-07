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

#include "parabolic.hpp"
#include "transport.hpp"

using namespace problem_setting;

/**
 * @class advection diffusion_problem
 * @brief Solves the Advection Diffusion evolutionary incompressible equations.
 * It consists of a transport problem for each velocity component and a parabolic problem for each velocity component.
 * Current implementation allows Dirichlet bc's only. Boundary conditions are accessed from reference solution, which is a function of space and time.
 * Current implementation allows for Brinkamn penalty method.
 */

#ifndef ADVECTION_DIFFUSION_PROBLEM_HPP
#define ADVECTION_DIFFUSION_PROBLEM_HPP

class advection_diffusion_problem
{
private:

DM dmGrid_staggered_x; ///< Discretized grid for the staggered x-direction.
DM dmGrid_staggered_y; ///< Discretized grid for the staggered y-direction.
DM dmGrid_staggered_z; ///< Discretized grid for the staggered z-direction.
DM dmGrid_centered;       ///< Discretized grid for magnitude
DM dmGrid_cent_rich;      ///< Discretized grid for shifted quantities.
DM dmGrid_stag_transp; ///< Staggered grid for transport equations.
DM dmGrid_shift_transp; ///< Shifted grid for transport computations.

Vec Magnitude; ///< Magnitude of velocity vectors.

Vec U_up, V_up, W_up; ///< Velocity fields in the x, y, and z directions.

Vec mask_U, mask_V, mask_W; ///< Mask vectors for boundary conditions.

//PetscErrorCode const update_bc_U(PetscReal const & theta);

//PetscErrorCode const update_bc_V(PetscReal const & theta);

//PetscErrorCode const update_bc_W(PetscReal const & theta);

/**
 * @brief Updates the velocity field, opearating final pressure correction for CT scheme.
 */
PetscErrorCode const update_velocity(PetscReal const & theta);

/**
 * @brief Post-processing function: assembles the magnitude of velocity vectors, locating a new vector in the cell-centers.
 * @param Magnitude_Shifted Output shifted magnitude vector on dmGrid_cent_rich (dofs on faces and cell=centers.)
 * @param U Velocity component in the x-direction.
 * @param V Velocity component in the y-direction.
 * @param W Velocity component in the z-direction.
 */
PetscErrorCode const assemble_magnitude(Vec & Magnitude_Shifted, Vec const & U, Vec const & V, Vec const & W);
/**
 * @brief Final assembly of magnitude vector on cell-centers.
 */
PetscErrorCode const compute_magnitude();

public:
/**
 * @brief Constructor that initializes the Advection-Diffusion problem with given grids and velocity fields.
 * This costructor has been defined for consistency and if future implementations will require solving Advection-Diffusion in a broader context.
 * For our applications only a stand-alone constructor is used.
 */
advection_diffusion_problem(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, DM const & dmGrid_centered, DM const & dmGrid_cent_rich, DM const & dmGrid_stag_transp, DM const dmGrid_shift_transp, Vec const & U_up, Vec const & V_up, Vec const W_up);

/**
 * @brief Default constructor that initializes stand-alone Advection-Diffusion problem with automatically created grids.
 */
advection_diffusion_problem();

/**
 * @brief Exports simulation in .vtk format format for visualization of x,y,z-componets, pressure and magnitude.
 */
PetscErrorCode exodus(size_t i);
/**
 * @brief Solves the Advection-Diffusion equations leveraging parabolic_problem_x, y, z, transport_problem_x, y, z, and poisson_problem classes.
 */
PetscErrorCode const solve();
/**
 * @brief Destructor to clean up allocated resources. Automatically calls sub-problems destructors. After destruction, a message is printed.
 */
~advection_diffusion_problem();

};

#endif // ADVECTION_DIFFUSION_PROBLEM_HPP


