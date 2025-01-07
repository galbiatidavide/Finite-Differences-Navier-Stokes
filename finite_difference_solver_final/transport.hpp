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

#ifndef TRANSPORT_PROBLEM_X_HPP
#define TRANSPORT_PROBLEM_X_HPP

class transport_problem_x {

protected:

    DM dmGrid_Shifted;
    DM dmGrid_Centered;
    DM dmGrid_Staggered;


    Vec U_n, V_n, W_n;

    Vec mask_U, mask_V, mask_W;
    //Vec U_int;

    PetscErrorCode FirstShiftU_y(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok

    PetscErrorCode FirstShiftU_z(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok

    PetscErrorCode FirstShiftV_y(Vec & VShifted, Vec const & vec, PetscScalar const & theta); //ok

    PetscErrorCode FirstShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta); //ok

    PetscErrorCode FirstDerive_y(Vec & AB_y, Vec const & AB); //ok

    PetscErrorCode FirstDerive_z(Vec & AB_z, Vec const & AB); //ok

    PetscErrorCode CenterU(Vec & UCenter, Vec const & vec, PetscReal const & theta);

    PetscErrorCode Derive_x(Vec & U2_x, Vec const & vec, PetscReal const & theta);

public:

    transport_problem_x(DM const & dmGrid_Shifted, DM const & dmGrid_Staggered, DM const & dmGrid_Centered) :
    dmGrid_Shifted(dmGrid_Shifted), dmGrid_Staggered(dmGrid_Staggered), dmGrid_Centered(dmGrid_Centered)
    {
        //DMCreateGlobalVector(dmGrid_Staggered, &U_int);
        DMCreateGlobalVector(dmGrid_Staggered, &U_n);
        DMCreateGlobalVector(dmGrid_Staggered, &V_n);
        DMCreateGlobalVector(dmGrid_Staggered, &W_n);
        /*VecCopy(U_0_, U_0);
        VecCopy(V_0_, V_0);
        VecCopy(W_0_, W_0);*/
        DMCreateGlobalVector(dmGrid_Staggered, &mask_U);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_V);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_W);

        if(brinkmann)
        {
            createMaskU(dmGrid_Staggered, mask_U, vertices, faces);
            createMaskV(dmGrid_Staggered, mask_V, vertices, faces);
            createMaskW(dmGrid_Staggered, mask_W, vertices, faces);
        }
        else {
            VecSet(mask_U, 0.0);
            VecSet(mask_V, 0.0);
            VecSet(mask_W, 0.0);
        }
    }

    transport_problem_x()
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMCreateGlobalVector(dmGrid_Staggered, &U_n);
        DMCreateGlobalVector(dmGrid_Staggered, &V_n);
        DMCreateGlobalVector(dmGrid_Staggered, &W_n);
        CreateAnalyticalU(dmGrid_Staggered, U_n, 0);
        CreateAnalyticalV(dmGrid_Staggered, V_n, 0);
        CreateAnalyticalW(dmGrid_Staggered, W_n, 0);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_U);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_V);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_W);

        if(brinkmann)
        {
            createMaskU(dmGrid_Staggered, mask_U, vertices, faces);
            createMaskV(dmGrid_Staggered, mask_V, vertices, faces);
            createMaskW(dmGrid_Staggered, mask_W, vertices, faces);
        }
        else {
            VecSet(mask_U, 0.0);
            VecSet(mask_V, 0.0);
            VecSet(mask_W, 0.0);
        }
    }

    PetscErrorCode const solve_step_x(PetscScalar const & theta);

    PetscErrorCode const solve_step_x(PetscScalar const & theta, Vec const & U_n, Vec const & V_n, Vec const & W_n);

    /*PetscErrorCode const solve_x()
    {
        PetscInt i;
        for (i = 0; i < 16; ++i) {
            solve_step_x(0.0);
        }
        PetscFunctionReturn(0);
    }*/

    //Vec get_U() const { return U_int; }

    /*void set_U(Vec const & U) { VecCopy(U, U_0); }
    void set_V(Vec const & V) { VecCopy(V, V_0); }
    void set_W(Vec const & W) { VecCopy(W, W_0); }*/


    ~transport_problem_x()
    {
        /*VecDestroy(&U_0);
        VecDestroy(&V_0);
        VecDestroy(&W_0);*/
        //VecDestroy(&U_int);
        VecDestroy(&U_n);
        VecDestroy(&V_n);
        VecDestroy(&W_n);
        DMDestroy(&dmGrid_Shifted);
        DMDestroy(&dmGrid_Staggered);
        DMDestroy(&dmGrid_Centered);
        VecDestroy(&mask_U);
        VecDestroy(&mask_V);
        VecDestroy(&mask_W);
   }

};

#endif // TRANSPORT_PROBLEM_X_HPP


#ifndef TRANSPORT_PROBLEM_Y_HPP
#define TRANSPORT_PROBLEM_Y_HPP

class transport_problem_y {
protected:

    DM dmGrid_Shifted;
    DM dmGrid_Centered;
    DM dmGrid_Staggered;

    Vec U_n, V_n, W_n;
    Vec mask_U, mask_V, mask_W;
    //Vec V_int;

    PetscErrorCode FirstShiftU_y(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok

    PetscErrorCode FirstShiftU_z(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok

    PetscErrorCode FirstShiftV_y(Vec & VShifted, Vec const & vec, PetscScalar const & theta); //ok

    PetscErrorCode FirstShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta); //ok

    PetscErrorCode SecondShiftV_z(Vec & VShifted, Vec const & vec, PetscScalar const & theta); //ok

    PetscErrorCode SecondShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta); //ok

    PetscErrorCode SecondDerive_x(Vec & AB_x, Vec const & AB); //ok

    PetscErrorCode SecondDerive_z(Vec & AB_z, Vec const & AB); //ok

    PetscErrorCode CenterV(Vec & VCenter, Vec const & vec, PetscReal const & theta); 

    PetscErrorCode Derive_y(Vec & V2_y, Vec const & vec, PetscReal const & theta);

public:

    transport_problem_y(DM const & dmGrid_Shifted, DM const & dmGrid_Staggered, DM const & dmGrid_Centered) :
    dmGrid_Shifted(dmGrid_Shifted), dmGrid_Staggered(dmGrid_Staggered), dmGrid_Centered(dmGrid_Centered)
    {
        DMCreateGlobalVector(dmGrid_Staggered, &U_n);
        DMCreateGlobalVector(dmGrid_Staggered, &V_n);
        DMCreateGlobalVector(dmGrid_Staggered, &W_n);
        /*DMCreateGlobalVector(dmGrid_Staggered, &V_0);
        DMCreateGlobalVector(dmGrid_Staggered, &W_0);
        DMCreateGlobalVector(dmGrid_Staggered, &V_int);
        VecCopy(U_0_, U_0);
        VecCopy(V_0_, V_0);
        VecCopy(W_0_, W_0);*/
        DMCreateGlobalVector(dmGrid_Staggered, &mask_U);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_V);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_W);

        if(brinkmann)
        {
            createMaskU(dmGrid_Staggered, mask_U, vertices, faces);
            createMaskV(dmGrid_Staggered, mask_V, vertices, faces);
            createMaskW(dmGrid_Staggered, mask_W, vertices, faces);
        }
        else {
            VecSet(mask_U, 0.0);
            VecSet(mask_V, 0.0);
            VecSet(mask_W, 0.0);
        }
    }

    transport_problem_y()
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMCreateGlobalVector(dmGrid_Staggered, &U_n);
        DMCreateGlobalVector(dmGrid_Staggered, &V_n);
        DMCreateGlobalVector(dmGrid_Staggered, &W_n);
        CreateAnalyticalU(dmGrid_Staggered, U_n, 0);
        CreateAnalyticalV(dmGrid_Staggered, V_n, 0);
        CreateAnalyticalW(dmGrid_Staggered, W_n, 0);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_U);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_V);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_W);

        if(brinkmann)
        {
            createMaskU(dmGrid_Staggered, mask_U, vertices, faces);
            createMaskV(dmGrid_Staggered, mask_V, vertices, faces);
            createMaskW(dmGrid_Staggered, mask_W, vertices, faces);
        }
        else {
            VecSet(mask_U, 0.0);
            VecSet(mask_V, 0.0);
            VecSet(mask_W, 0.0);
        }
    }

    PetscErrorCode const solve_step_y(PetscScalar const & theta);

    PetscErrorCode const solve_step_y(PetscScalar const & theta, Vec const & U_n, Vec const & V_n, Vec const & W_n);

    /*Vec get_V() const { return V_int; }

    void set_U(Vec const & U) { VecCopy(U, U_0); }
    void set_V(Vec const & V) { VecCopy(V, V_0); }
    void set_W(Vec const & W) { VecCopy(W, W_0); }*/

    ~transport_problem_y()
    {
        VecDestroy(&U_n);
        VecDestroy(&V_n);
        VecDestroy(&W_n);
        DMDestroy(&dmGrid_Shifted);
        DMDestroy(&dmGrid_Staggered);
        DMDestroy(&dmGrid_Centered);
        //VecDestroy(&V_int);
        VecDestroy(&mask_U);
        VecDestroy(&mask_V);
        VecDestroy(&mask_W);
    }

};

#endif // TRANSPORT_PROBLEM__Y_HPP

#ifndef TRANSPORT_PROBLEM_Z_HPP
#define TRANSPORT_PROBLEM_Z_HPP

class transport_problem_z {

protected:

    DM dmGrid_Shifted;
    DM dmGrid_Centered;
    DM dmGrid_Staggered;

    Vec U_n, V_n, W_n;
    Vec mask_U, mask_V, mask_W;
    //Vec W_int;

    PetscErrorCode FirstShiftU_y(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok

    PetscErrorCode FirstShiftU_z(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok

    PetscErrorCode FirstShiftV_y(Vec & VShifted, Vec const & vec, PetscScalar const & theta); //ok

    PetscErrorCode FirstShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta); //ok

    PetscErrorCode SecondShiftV_z(Vec & VShifted, Vec const & vec, PetscScalar const & theta); //ok
   
    PetscErrorCode SecondShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta); //ok
 
    PetscErrorCode ThirdDerive_x(Vec & AB_x, Vec const & AB); //ok

    PetscErrorCode ThirdDerive_y(Vec & AB_y, Vec const & AB); //ok

    PetscErrorCode CenterW(Vec & WCenter, Vec const & vec, PetscReal const & theta); 

    PetscErrorCode Derive_z(Vec & W2_z, Vec const & vec, PetscReal const & theta);

public:

    transport_problem_z(DM const & dmGrid_Shifted, DM const & dmGrid_Staggered, DM const & dmGrid_Centered) :
    dmGrid_Shifted(dmGrid_Shifted), dmGrid_Staggered(dmGrid_Staggered), dmGrid_Centered(dmGrid_Centered)
    {
        DMCreateGlobalVector(dmGrid_Staggered, &U_n);
        DMCreateGlobalVector(dmGrid_Staggered, &V_n);
        DMCreateGlobalVector(dmGrid_Staggered, &W_n);
        /*DMCreateGlobalVector(dmGrid_Staggered, &V_0);
        DMCreateGlobalVector(dmGrid_Staggered, &W_0);
        DMCreateGlobalVector(dmGrid_Staggered, &W_int);
        VecCopy(U_0_, U_0);
        VecCopy(V_0_, V_0);
        VecCopy(W_0_, W_0);*/
        DMCreateGlobalVector(dmGrid_Staggered, &mask_U);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_V);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_W);

        if(brinkmann)
        {
            createMaskU(dmGrid_Staggered, mask_U, vertices, faces);
            createMaskV(dmGrid_Staggered, mask_V, vertices, faces);
            createMaskW(dmGrid_Staggered, mask_W, vertices, faces);
        }
        else {
            VecSet(mask_U, 0.0);
            VecSet(mask_V, 0.0);
            VecSet(mask_W, 0.0);
        }
    }

    transport_problem_z()
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMCreateGlobalVector(dmGrid_Staggered, &U_n);
        DMCreateGlobalVector(dmGrid_Staggered, &V_n);
        DMCreateGlobalVector(dmGrid_Staggered, &W_n);
        CreateAnalyticalU(dmGrid_Staggered, U_n, 0);
        CreateAnalyticalV(dmGrid_Staggered, V_n, 0);
        CreateAnalyticalW(dmGrid_Staggered, W_n, 0);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_U);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_V);
        DMCreateGlobalVector(dmGrid_Staggered, &mask_W);

        if(brinkmann)
        {
            createMaskU(dmGrid_Staggered, mask_U, vertices, faces);
            createMaskV(dmGrid_Staggered, mask_V, vertices, faces);
            createMaskW(dmGrid_Staggered, mask_W, vertices, faces);
        }
        else {
            VecSet(mask_U, 0.0);
            VecSet(mask_V, 0.0);
            VecSet(mask_W, 0.0);
        }
    }

    PetscErrorCode const solve_step_z(PetscScalar const & theta);

    PetscErrorCode const solve_step_z(PetscScalar const & theta, Vec const & U_n, Vec const & V_n, Vec const & W_n);

    /*Vec get_W() const { return W_int; }

    void set_U(Vec const & U) { VecCopy(U, U_0); }
    void set_V(Vec const & V) { VecCopy(V, V_0); }
    void set_W(Vec const & W) { VecCopy(W, W_0); }*/

    ~transport_problem_z()
    {
        VecDestroy(&U_n);
        VecDestroy(&V_n);
        VecDestroy(&W_n);
        DMDestroy(&dmGrid_Shifted);
        DMDestroy(&dmGrid_Staggered);
        DMDestroy(&dmGrid_Centered);
        VecDestroy(&mask_U);
        VecDestroy(&mask_V);
        VecDestroy(&mask_W);
    }

};

#endif // TRANSPORT_PROBLEM_Z_HPP
