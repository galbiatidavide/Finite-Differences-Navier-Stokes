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

#include "problem_setting.hpp"
using namespace problem_setting;

#ifndef PARABOLIC_PROBLEM_X_HPP
#define PARABOLIC_PROBLEM_X_HPP

class parabolic_problem_x {
private:

    DM dmGrid;
    Mat A;
    Vec rhs;
    Vec U_up;

    Vec mask_U;
    
    const PetscErrorCode assemble_rhs(PetscReal const & theta);
   
    const PetscErrorCode assemble_rhs(PetscReal const & theta, Vec const & U_up);

    PetscErrorCode exodus(size_t const & i);

public:

    parabolic_problem_x(DM const & dmGrid) :
    dmGrid(dmGrid)
    {   
        DMCreateGlobalVector(dmGrid, &U_up);
        DMCreateMatrix(dmGrid, &A);
        DMCreateGlobalVector(dmGrid, &rhs);

        DMCreateGlobalVector(dmGrid, &mask_U);
        if(brinkmann){
            createMaskU(dmGrid, mask_U, vertices, faces);
        }
        else {
            VecSet(mask_U, 0.0);
        }
    }

    parabolic_problem_x()
    {   
        CreateGrid(&dmGrid, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMCreateMatrix(dmGrid, &A);
        DMCreateGlobalVector(dmGrid, &rhs);
        DMCreateGlobalVector(dmGrid, &U_up);
        CreateAnalyticalU(dmGrid, U_up, 0);

        DMCreateGlobalVector(dmGrid, &mask_U);
        if(brinkmann){
            createMaskU(dmGrid, mask_U, vertices, faces);
        }
        else {
            VecSet(mask_U, 0.0);
        }
        std::cout << "Creating parabolic_problem_x" << std::endl;
    }

    Vec get_U();

    void set_U(Vec const & U); 

    PetscErrorCode assemble_lhs();

    PetscErrorCode solve_step(PetscReal const & theta);

    PetscErrorCode solve_step(PetscReal const & theta, Vec const & U_up);
   
    PetscErrorCode solve();

    ~parabolic_problem_x()
    {
        MatDestroy(&A);
        VecDestroy(&rhs);
        VecDestroy(&U_up);
        DMDestroy(&dmGrid);
        VecDestroy(&mask_U);
    }

};

#endif // PARABOLIC_PROBLEM_X_HPP

#ifndef PARABOLIC_PROBLEM_Y_HPP
#define PARABOLIC_PROBLEM_Y_HPP

class parabolic_problem_y {
private:
    DM  dmGrid;
    Mat A;
    Vec rhs;
    Vec V_up;

    Vec mask_V;

    const PetscErrorCode assemble_rhs(PetscReal const & theta);

    const PetscErrorCode assemble_rhs(PetscReal const & theta, Vec const & V_up);

    PetscErrorCode exodus(size_t const & i);


public:

    parabolic_problem_y(DM const & dmGrid) :
    dmGrid(dmGrid)
    {
        DMCreateGlobalVector(dmGrid, &V_up);
        DMCreateMatrix(dmGrid, &A);
        DMCreateGlobalVector(dmGrid, &rhs);
        DMCreateGlobalVector(dmGrid, &mask_V);
        if(brinkmann){
            createMaskV(dmGrid, mask_V, vertices, faces);
        }
        else {
            VecSet(mask_V, 0.0);
        }
    }

    parabolic_problem_y()
    {   CreateGrid(&dmGrid, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMCreateMatrix(dmGrid, &A);
        DMCreateGlobalVector(dmGrid, &rhs);
        DMCreateGlobalVector(dmGrid, &V_up);
        CreateAnalyticalV(dmGrid, V_up, 0);
        DMCreateGlobalVector(dmGrid, &mask_V);
        if(brinkmann){
            createMaskV(dmGrid, mask_V, vertices, faces);
        }
        else {
            VecSet(mask_V, 0.0);
        }
    }

    Vec get_V(); 

    void set_V(Vec const & V);

    PetscErrorCode assemble_lhs();
    
    PetscErrorCode solve_step(PetscReal const & theta);
  
    PetscErrorCode solve_step(PetscReal const & theta, Vec const & V_up);

    PetscErrorCode solve();

    ~parabolic_problem_y()
    {
        MatDestroy(&A);
        VecDestroy(&rhs);
        VecDestroy(&V_up);
        DMDestroy(&dmGrid);
        VecDestroy(&mask_V);
    }

};

#endif // PARABOLIC_PROBLEM_Y_HPP

#ifndef PARABOLIC_PROBLEM_Z_HPP
#define PARABOLIC_PROBLEM_Z_HPP

class parabolic_problem_z {
private:

    DM dmGrid;
    Mat A;
    Vec rhs;
    Vec W_up;

    Vec mask_W;

    PetscErrorCode const assemble_rhs(PetscReal const & theta);

    PetscErrorCode const assemble_rhs(PetscReal const & theta, Vec const & W_up);

    PetscErrorCode exodus(size_t const & i);


public:

    parabolic_problem_z(DM const & dmGrid) :
    dmGrid(dmGrid)
    {
        DMCreateGlobalVector(dmGrid, &W_up);
        DMCreateMatrix(dmGrid, &A);
        DMCreateGlobalVector(dmGrid, &rhs);
        DMCreateGlobalVector(dmGrid, &mask_W);
        if(brinkmann){
            createMaskW(dmGrid, mask_W, vertices, faces);
        }
        else {
            VecSet(mask_W, 0.0);
        }
    }

    parabolic_problem_z()
    {   CreateGrid(&dmGrid, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMCreateMatrix(dmGrid, &A);
        DMCreateGlobalVector(dmGrid, &rhs);
        DMCreateGlobalVector(dmGrid, &W_up);
        CreateAnalyticalW(dmGrid, W_up, 0);
        DMCreateGlobalVector(dmGrid, &mask_W);
        if(brinkmann){
            createMaskW(dmGrid, mask_W, vertices, faces);
        }
        else {
            VecSet(mask_W, 0.0);
        }
    }

    Vec get_W();

    void set_W(Vec const & W);

    PetscErrorCode assemble_lhs();

    PetscErrorCode solve_step(PetscReal const & theta);

    PetscErrorCode solve_step(PetscReal const & theta, Vec const & W_up);

    PetscErrorCode solve();
    
    ~parabolic_problem_z()
    {
        MatDestroy(&A);
        VecDestroy(&rhs);
        VecDestroy(&W_up);
        DMDestroy(&dmGrid);
        VecDestroy(&mask_W);
    }

};

#endif // PARABOLIC_PROBLEM_Z_HPP
