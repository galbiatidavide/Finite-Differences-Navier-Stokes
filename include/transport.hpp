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

/**
 * @class transport_problem_x
 * @brief Represents a transport problem in the x-direction.
 *
 * This class solves an advection-transport problem in the x-direction using 
 * 2nd order centered differences and a fully explicit approach.
 */

#ifndef TRANSPORT_PROBLEM_X_HPP
#define TRANSPORT_PROBLEM_X_HPP

class transport_problem_x {

private:

    DM dmGrid_Shifted;   ///< DMGrid with dofs on faces and edges.
    DM dmGrid_Staggered; ///< DMGrid with dofs of faces.
    DM dmGrid_Centered;  ///< DMGrid with dofs of faces and centers.

    Vec U_n, V_n, W_n;   ///< Velocity fields at the current time step.

    Vec mask_U, mask_V, mask_W; ///< Mask vectors for boundary conditions.

    /**
     * @brief Performs the shift for x-component in a 3D transport non-linear problem in the y-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode FirstShiftU_y(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok
    /**
     * @brief Performs the shift for x-component in a 3D transport non-linear problem in the z-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode FirstShiftU_z(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok
    /**
     * @brief Performs the shift for y-component in a 3D transport non-linear problem in the y-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode FirstShiftV_y(Vec & VShifted, Vec const & vec, PetscScalar const & theta); //ok
    /**
     * @brief Performs the shift for z-component in a 3D transport non-linear problem in the z-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode FirstShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta); //ok
    /**
     * @brief Computes the first derivative in the y-direction, placing the result on faces.
     */
    PetscErrorCode FirstDerive_y(Vec & AB_y, Vec const & AB); //ok
    /**
     * @brief Computes the first derivative in the z-direction, placing the result on faces.
     */
    PetscErrorCode FirstDerive_z(Vec & AB_z, Vec const & AB); //ok
    /**
     * @brief Centers the x-velocity component in cell-centers.
     */
    PetscErrorCode CenterU(Vec & UCenter, Vec const & vec, PetscReal const & theta);
    /**
     * @brief Computes the first derivative in the x-direction, placing the result on faces.
     */
    PetscErrorCode Derive_x(Vec & U2_x, Vec const & vec, PetscReal const & theta);

public:
    /**
     * @brief Constructor initializing the transport problem with given grids.
     */
    transport_problem_x(DM const & dmGrid_Shifted, DM const & dmGrid_Staggered, DM const & dmGrid_Centered); 

    /**
     * @brief Default constructor for stand-alone transport problem.
     */
    transport_problem_x();
  

    /**
     * @brief Performs a single time step for solving the transport equation.
     */

    PetscErrorCode const solve_step_x(PetscScalar const & theta, 
    std::optional<std::reference_wrapper<Vec>> U_n_opt = std::nullopt, 
    std::optional<std::reference_wrapper<Vec>> V_n_opt = std::nullopt, 
    std::optional<std::reference_wrapper<Vec>> W_n_opt = std::nullopt);

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


    /**
     * @brief Destructor to clean up allocated resources.
     */
    ~transport_problem_x();


};

#endif // TRANSPORT_PROBLEM_X_HPP

/**
 * @class transport_problem_y
 * @brief Represents a transport problem in the y-direction.
 *
 * This class solves an advection-transport problem in the y-direction using 
 * 2nd order centered differences and a fully explicit approach.
 */

#ifndef TRANSPORT_PROBLEM_Y_HPP
#define TRANSPORT_PROBLEM_Y_HPP

class transport_problem_y {
private:

    DM dmGrid_Shifted;   ///< DMGrid with dofs on faces and edges.
    DM dmGrid_Staggered; ///< DMGrid with dofs of faces.
    DM dmGrid_Centered;  ///< DMGrid with dofs of faces and centers.

    Vec U_n, V_n, W_n;   ///< Velocity fields at the current time step.

    Vec mask_U, mask_V, mask_W; ///< Mask vectors for boundary conditions.

    /**
     * @brief Performs the shift for x-component in a 3D transport non-linear problem in the y-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode FirstShiftU_y(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok
    /**
     * @brief Performs the shift for y-component in a 3D transport non-linear problem in the y-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode FirstShiftV_y(Vec & VShifted, Vec const & vec, PetscScalar const & theta); //ok
    /**
     * @brief Performs the shift for y-component in a 3D transport non-linear problem in the z-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode SecondShiftV_z(Vec & VShifted, Vec const & vec, PetscScalar const & theta); //ok
    /**
     * @brief Performs the shift for z-component in a 3D transport non-linear problem in the z-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode SecondShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta); //ok
    /**
     * @brief Computes the first derivative in the x-direction, placing the result on faces.
     */
    PetscErrorCode SecondDerive_x(Vec & AB_x, Vec const & AB); //ok
    /**
     * @brief Computes the first derivative in the z-direction, placing the result on faces.
     */
    PetscErrorCode SecondDerive_z(Vec & AB_z, Vec const & AB); //ok
    /**
     * @brief Centers the y-velocity component in cell-centers.
     */
    PetscErrorCode CenterV(Vec & VCenter, Vec const & vec, PetscReal const & theta); 
    /**
     * @brief Computes the first derivative in the y-direction, placing the result on faces.
     */
    PetscErrorCode Derive_y(Vec & V2_y, Vec const & vec, PetscReal const & theta);

public:
    /**
     * @brief Constructor initializing the transport problem with given grids.
     */
    transport_problem_y(DM const & dmGrid_Shifted, DM const & dmGrid_Staggered, DM const & dmGrid_Centered); 
   
    /**
     * @brief Default constructor for stand-alone transport problem.
     */
    transport_problem_y();

    /**
     * @brief Performs a single time step for solving the transport equation.
     */

    PetscErrorCode const solve_step_y(PetscScalar const & theta,     
    std::optional<std::reference_wrapper<Vec>> U_n_opt = std::nullopt, 
    std::optional<std::reference_wrapper<Vec>> V_n_opt = std::nullopt, 
    std::optional<std::reference_wrapper<Vec>> W_n_opt = std::nullopt);

    /*Vec get_V() const { return V_int; }

    void set_U(Vec const & U) { VecCopy(U, U_0); }
    void set_V(Vec const & V) { VecCopy(V, V_0); }
    void set_W(Vec const & W) { VecCopy(W, W_0); }*/
    /**
     * @brief Destructor to clean up allocated resources.
     */
    ~transport_problem_y();

};

#endif // TRANSPORT_PROBLEM__Y_HPP

/**
 * @class transport_problem_z
 * @brief Represents a transport problem in the z-direction.
 *
 * This class solves an advection-transport problem in the z-direction using 
 * 2nd order centered differences and a fully explicit approach.
 */

#ifndef TRANSPORT_PROBLEM_Z_HPP
#define TRANSPORT_PROBLEM_Z_HPP

class transport_problem_z {

private:

    DM dmGrid_Shifted;   ///< DMGrid with dofs on faces and edges.
    DM dmGrid_Staggered; ///< DMGrid with dofs of faces.
    DM dmGrid_Centered;  ///< DMGrid with dofs of faces and centers.

    Vec U_n, V_n, W_n;   ///< Velocity fields at the current time step.

    Vec mask_U, mask_V, mask_W; ///< Mask vectors for boundary conditions.

    /**
     * @brief Performs the shift for x-component in a 3D transport non-linear problem in the z-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode FirstShiftU_z(Vec & UShifted, Vec const & vec, PetscScalar const & theta); //okok
    /**
     * @brief Performs the shift for z-component in a 3D transport non-linear problem in the z-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode FirstShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta); //ok
    /**
     * @brief Performs the shift for y-component in a 3D transport non-linear problem in the z-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode SecondShiftV_z(Vec & VShifted, Vec const & vec, PetscScalar const & theta); //ok
    /**
     * @brief Performs the shift for z-component in a 3D transport non-linear problem in the z-direction. Interpolates velocity component on edges.
     */
    PetscErrorCode SecondShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta); //ok
    /**
     * @brief Computes the first derivative in the x-direction, placing the result on faces.
    */
    PetscErrorCode ThirdDerive_x(Vec & AB_x, Vec const & AB); //ok
    /**
     * @brief Computes the first derivative in the y-direction, placing the result on faces.
    */
    PetscErrorCode ThirdDerive_y(Vec & AB_y, Vec const & AB); //ok
    /**
     * @brief Centers the z-velocity component in cell-centers.
     */
    PetscErrorCode CenterW(Vec & WCenter, Vec const & vec, PetscReal const & theta); 
    /**
     * @brief Computes the first derivative in the z-direction, placing the result on faces.
     */
    PetscErrorCode Derive_z(Vec & W2_z, Vec const & vec, PetscReal const & theta);

public:
    /**
     * @brief Constructor initializing the transport problem with given grids.
     */
    transport_problem_z(DM const & dmGrid_Shifted, DM const & dmGrid_Staggered, DM const & dmGrid_Centered); 
    /**
     * @brief Default constructor for stand-alone transport problem.
     */
    transport_problem_z();

    /**
     * @brief Performs a single time step for solving the transport equation.
     */

    PetscErrorCode const solve_step_z(PetscScalar const & theta,     
    std::optional<std::reference_wrapper<Vec>> U_n_opt = std::nullopt, 
    std::optional<std::reference_wrapper<Vec>> V_n_opt = std::nullopt, 
    std::optional<std::reference_wrapper<Vec>> W_n_opt = std::nullopt);

    /*Vec get_W() const { return W_int; }

    void set_U(Vec const & U) { VecCopy(U, U_0); }
    void set_V(Vec const & V) { VecCopy(V, V_0); }
    void set_W(Vec const & W) { VecCopy(W, W_0); }*/
    /**
     * @brief Destructor to clean up allocated resources.
     */
    ~transport_problem_z();

};

#endif // TRANSPORT_PROBLEM_Z_HPP
