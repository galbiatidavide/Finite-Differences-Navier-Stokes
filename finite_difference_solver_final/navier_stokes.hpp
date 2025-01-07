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

//penalizzazione
//rivedo slide formaggia
//templare args solve_step

using namespace problem_setting;

#ifndef NAVIER_STOKES_PROBLEM_HPP
#define NAVIER_STOKES_PROBLEM_HPP

class navier_stokes_problem
{
private:

//Cannot be declared as constant
DM dmGrid_staggered_x;
DM dmGrid_staggered_y;
DM dmGrid_staggered_z;
DM dmGrid_centered;
DM dmGrid_cent_rich;
DM dmGrid_shift_transp;
DM dmGrid_stag_transp;

Vec P;
Vec P_x;
Vec P_y;
Vec P_z;
Vec Magnitude;

Mat A;

Vec U_up, V_up, W_up;

Vec mask_U, mask_V, mask_W;

PetscErrorCode const assemble_lhs();

PetscErrorCode const assemble_divergence(Vec & div, Vec const & U, Vec const &  V, Vec const & W);

PetscErrorCode const compute_divergence(Vec & div, Vec const & U_n, Vec const & V_n, Vec const & W_n);

PetscErrorCode const derive_x_P(Vec & P_x_shifted, Vec const & vec);

PetscErrorCode const derive_y_P(Vec & P_y_shifted, Vec const & vec);

PetscErrorCode const derive_z_P(Vec & P_z_shifted, Vec const & vec);

PetscErrorCode const manage_pressure_x();

PetscErrorCode const manage_pressure_y();

PetscErrorCode const manage_pressure_z();

PetscErrorCode const manage_pressure();

PetscErrorCode const update_bc_U(PetscReal const & theta);

PetscErrorCode const update_bc_V(PetscReal const & theta);

PetscErrorCode const update_bc_W(PetscReal const & theta);

PetscErrorCode const update_velocity(PetscReal const & theta);

PetscErrorCode const assemble_magnitude(Vec & Magnitude_Shifted, Vec const & U, Vec const & V, Vec const & W);

PetscErrorCode const compute_magnitude();

public:

navier_stokes_problem(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, DM const & dmGrid_centered, DM const & dmGrid_cent_rich, DM const & dmGrid_stag_transp, DM const dmGrid_shift_transp, Vec const & U_up, Vec const & V_up, Vec const W_up)
    : dmGrid_staggered_x(dmGrid_staggered_x), dmGrid_staggered_y(dmGrid_staggered_y), dmGrid_staggered_z(dmGrid_staggered_z), dmGrid_centered(dmGrid_centered), dmGrid_cent_rich(dmGrid_cent_rich), dmGrid_stag_transp(dmGrid_stag_transp), dmGrid_shift_transp(dmGrid_shift_transp), U_up(U_up), V_up(V_up), W_up(W_up)

{
    DMCreateGlobalVector(dmGrid_centered, &P);
    DMCreateGlobalVector(dmGrid_staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_staggered_z, &P_z);
    DMCreateGlobalVector(dmGrid_centered, &Magnitude);
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
}

navier_stokes_problem()
{
    //Allocate the grids
    CreateGrid(&dmGrid_staggered_x, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_staggered_y, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_staggered_z, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_cent_rich, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_shift_transp, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid_stag_transp, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

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
    DMCreateMatrix(dmGrid_centered, &A);

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

PetscErrorCode exodus(size_t i);

PetscErrorCode const solve();

~navier_stokes_problem()
{
    VecDestroy(&P);
    VecDestroy(&P_x);
    VecDestroy(&P_y);
    VecDestroy(&P_z);
    VecDestroy(&Magnitude);
    MatDestroy(&A);
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
    std::cout << "Destructor called" << std::endl;

}

};

#endif // NAVIER_STOKES_PROBLEM_HPP


