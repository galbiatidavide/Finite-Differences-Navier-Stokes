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
advection_diffusion_problem(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, DM const & dmGrid_centered, DM const & dmGrid_cent_rich, DM const & dmGrid_stag_transp, DM const dmGrid_shift_transp, Vec const & U_up, Vec const & V_up, Vec const W_up)
    : dmGrid_staggered_x(dmGrid_staggered_x), dmGrid_staggered_y(dmGrid_staggered_y), dmGrid_staggered_z(dmGrid_staggered_z), dmGrid_centered(dmGrid_centered), dmGrid_cent_rich(dmGrid_cent_rich), dmGrid_stag_transp(dmGrid_stag_transp), dmGrid_shift_transp(dmGrid_shift_transp), U_up(U_up), V_up(V_up), W_up(W_up)

{
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
 * @brief Default constructor that initializes stand-alone Advection-Diffusion problem with automatically created grids.
 */
advection_diffusion_problem()
{
    //Allocate the grids
    CreateGrid(&dmGrid_staggered_x, 0, 1, 0);
    CreateGrid(&dmGrid_staggered_y, 0, 1, 0);
    CreateGrid(&dmGrid_staggered_z, 0, 1, 0);
    CreateGrid(&dmGrid_centered, 0, 0, 1);
    CreateGrid(&dmGrid_cent_rich, 0, 1, 1);
    CreateGrid(&dmGrid_shift_transp, 1, 1, 0);
    CreateGrid(&dmGrid_stag_transp, 0, 1, 0);

    //Create parallel vectors
    DMCreateGlobalVector(dmGrid_staggered_x, &U_up);
    DMCreateGlobalVector(dmGrid_staggered_y, &V_up);
    DMCreateGlobalVector(dmGrid_staggered_z, &W_up);
    CreateAnalyticalU(dmGrid_staggered_x, U_up, 0);
    CreateAnalyticalV(dmGrid_staggered_y, V_up, 0);
    CreateAnalyticalW(dmGrid_staggered_z, W_up, 0);

    DMCreateGlobalVector(dmGrid_centered, &Magnitude);

    /*DMCreateGlobalVector(dmGrid_staggered_x, &U_prova);
    DMCreateGlobalVector(dmGrid_staggered_y, &V_prova);
    DMCreateGlobalVector(dmGrid_staggered_z, &W_prova);*/

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
 * @brief Solves the Advection-Diffusion equations leveraging parabolic_problem_x, y, z, transport_problem_x, y, z, and poisson_problem classes.
 */
PetscErrorCode const solve();
/**
 * @brief Destructor to clean up allocated resources. Automatically calls sub-problems destructors. After destruction, a message is printed.
 */
~advection_diffusion_problem()
{
    VecDestroy(&Magnitude);
    VecDestroy(&U_up);
    VecDestroy(&V_up);
    VecDestroy(&W_up);
    DMDestroy(&dmGrid_staggered_x);
    DMDestroy(&dmGrid_staggered_y);
    DMDestroy(&dmGrid_staggered_z);
    DMDestroy(&dmGrid_centered);
    DMDestroy(&dmGrid_cent_rich);
    DMDestroy(&dmGrid_shift_transp);
    DMDestroy(&dmGrid_stag_transp);
    VecDestroy(&mask_U);
    VecDestroy(&mask_V);
    VecDestroy(&mask_W);


    /*VecDestroy(&U_prova);
    VecDestroy(&V_prova);
    VecDestroy(&W_prova);*/
    // NO NEED TO CALL: HO STAMPATO QUANDO ENTRA NEL DISTRUTTORE DI PARABOLIC_PROBLEM E LO FA DOPO cout << "Destructor called" << std::endl;
    /*pb_x.~parabolic_problem_x();
    pb_y.~parabolic_problem_y();
    pb_z.~parabolic_problem_z();*/
    std::cout << "Advection-Diffusion Destructor Called" << std::endl;

}

};

#endif // ADVECTION_DIFFUSION_PROBLEM_HPP


