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
#include "poisson.hpp"

using namespace problem_setting;

/**
 * @class stokes_problem
 * @brief Solves the Stokes evolutionary incompressible Stokes equations, applying a first order Chorin-Temam algorithm.
 * It consists of a parabolic problem for each velocity component, and a Poisson problem for the pressure.
 * Current implementation allows Dirichlet bc's only. Boundary conditions are accessed from reference solution, which is a function of space and time.
 * Current implementation allows for Brinkamn penalty method.
 */

#ifndef STOKES_PROBLEM_HPP
#define STOKES_PROBLEM_HPP

class stokes_problem
{
private:

DM dmGrid_staggered_x; ///< Discretized grid for the staggered x-direction.
DM dmGrid_staggered_y; ///< Discretized grid for the staggered y-direction.
DM dmGrid_staggered_z; ///< Discretized grid for the staggered z-direction.
DM dmGrid_centered;    ///< Discretized grid for the centered formulation.
DM dmGrid_cent_rich;   ///< Refined centered grid for pressure computation.

Vec P;       ///< Pressure field.
Vec P_x;     ///< Pressure derivative in the x-direction.
Vec P_y;     ///< Pressure derivative in the y-direction.
Vec P_z;     ///< Pressure derivative in the z-direction.
Vec Magnitude; ///< Magnitude of velocity vectors.

Vec U_up, V_up, W_up; ///< Velocity fields in the x, y, and z directions.

Vec mask_U, mask_V, mask_W; ///< Mask vectors for boundary conditions.


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
 * @brief Constructor that initializes the Stokes problem with given grids and velocity fields.
 * This costructor has been defined for consistency and if future implementations will require solving Stokes in a broader context.
 * For our applications only a stand-alone constructor is used.
 */
stokes_problem(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, DM const & dmGrid_centered, DM const & dmGrid_cent_rich, Vec const & U_up, Vec const & V_up, Vec const W_up)
    : dmGrid_staggered_x(dmGrid_staggered_x), dmGrid_staggered_y(dmGrid_staggered_y), dmGrid_staggered_z(dmGrid_staggered_z), dmGrid_centered(dmGrid_centered), dmGrid_cent_rich(dmGrid_cent_rich), U_up(U_up), V_up(V_up), W_up(W_up)

{
    DMCreateGlobalVector(dmGrid_centered, &P);
    DMCreateGlobalVector(dmGrid_staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_staggered_z, &P_z);
    DMCreateGlobalVector(dmGrid_centered, &Magnitude);

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
}
/**
 * @brief Default constructor that initializes stand-alone Stokes problem with automatically created grids.
 */
stokes_problem()
{
    //Allocate the grids
    CreateGrid(&dmGrid_staggered_x, 0, 1, 0);
    CreateGrid(&dmGrid_staggered_y, 0, 1, 0);
    CreateGrid(&dmGrid_staggered_z, 0, 1, 0);
    CreateGrid(&dmGrid_centered, 0, 0, 1);
    CreateGrid(&dmGrid_cent_rich, 0, 1, 1);

    //Create parallel vectors
    DMCreateGlobalVector(dmGrid_staggered_x, &U_up);
    DMCreateGlobalVector(dmGrid_staggered_y, &V_up);
    DMCreateGlobalVector(dmGrid_staggered_z, &W_up);
    CreateAnalyticalU(dmGrid_staggered_x, U_up, 0);
    CreateAnalyticalV(dmGrid_staggered_y, V_up, 0);
    CreateAnalyticalW(dmGrid_staggered_z, W_up, 0);

    DMCreateGlobalVector(dmGrid_centered, &P);
    DMCreateGlobalVector(dmGrid_staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_staggered_z, &P_z);
    DMCreateGlobalVector(dmGrid_centered, &Magnitude);

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
/**
 * @brief Exports simulation in .vtk format format for visualization of x,y,z-componets, pressure and magnitude.
 */
PetscErrorCode exodus(size_t i);
/**
 * @brief Solves the Stokes equations leveraging parabolic_problem_x, y, z, transport_problem_x, y, z, and poisson_problem classes.
 */
PetscErrorCode const solve();
/**
 * @brief Destructor to clean up allocated resources. Automatically calls sub-problems destructors. After destruction, a message is printed.
 */
~stokes_problem()
{
    VecDestroy(&P);
    VecDestroy(&P_x);
    VecDestroy(&P_y);
    VecDestroy(&P_z);
    VecDestroy(&Magnitude);
    VecDestroy(&U_up);
    VecDestroy(&V_up);
    VecDestroy(&W_up);
    DMDestroy(&dmGrid_staggered_x);
    DMDestroy(&dmGrid_staggered_y);
    DMDestroy(&dmGrid_staggered_z);
    DMDestroy(&dmGrid_centered);
    DMDestroy(&dmGrid_cent_rich);
    VecDestroy(&mask_U);
    VecDestroy(&mask_V);
    VecDestroy(&mask_W);
    std::cout << "Stokes Destructor Called" << std::endl;

}

};

#endif // STOKES_PROBLEM_HPP


