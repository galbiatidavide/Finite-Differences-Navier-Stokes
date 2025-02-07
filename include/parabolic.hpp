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
 * @class parabolic_problem_x
 * @brief Represents a parabolic problem in the x-direction.
 *
 * This class solves a generic evolutionary diffusive problem in the x-direction with fully Dirichlet bc's. 
 * Boundary conditions are imposed by means of a reference solution.
 * Discretization is performed on a staggered grid and in x-direction only, meaning that variables are located in position LEFT and RIGHT.
 * Linear problem is solved by means of a PETSc KSP solver, with GMRES and Jacobi preconditioner.
 */

#ifndef PARABOLIC_PROBLEM_X_HPP
#define PARABOLIC_PROBLEM_X_HPP

class parabolic_problem_x {
private:

    DM dmGrid;    ///< Discretized grid for the problem.
    Mat A;        ///< Matrix representing the left-hand side of the system.
    Vec rhs;      ///< Right-hand side vector.
    Vec U_up;     ///< Solution vector for the x-direction. If not provided, it must be passed as reference in solve_step(...)

    Vec mask_U;   ///< Mask vector used for Brinkman flow.

    /**
     * @brief Assembles the right-hand side (RHS) vector.
     */    
    PetscErrorCode assemble_rhs(PetscReal const & theta, Vec const & U_up);

    /**
     * @brief Assembles the left-hand side (LHS) matrix.
     */
    PetscErrorCode assemble_lhs();

    /**
     * @brief Exports results in .vtk format for post-processing.
     */
    PetscErrorCode exodus(size_t const & i);

public:
    /**
     * @brief Constructor that initializes the problem with a given grid. Use when parabolic problem is just a step of a bigger problem, like in Chorin-Temam method.
     * @param dmGrid staggered grid petsc-object must be already be defined.
     */
    parabolic_problem_x(DM const & dmGrid);

    /**
     * @brief Default constructor for stand-alone problem.
    */
    parabolic_problem_x();

    /**
     * @brief Performs a single time step of the numerical solution.
     */
    PetscErrorCode solve_step(PetscReal const & theta, std::optional<std::reference_wrapper<Vec>> U_up_opt = std::nullopt);   
    /**
     * @brief Solves the entire parabolic problem in the x-direction.
     */
    PetscErrorCode solve();

    /**
     * @brief Destructor to clean up allocated resources. Fundamental in PETSc implementation lo avoid leaks and unexpected RAM overhead.
     */
    ~parabolic_problem_x();

};

#endif // PARABOLIC_PROBLEM_X_HPP

/**
 * @class parabolic_problem_y
 * @brief Represents a parabolic problem in the y-direction.
 *
 * This class solves a generic evolutionary diffusive problem in the y-direction with fully Dirichlet bc's. 
 * Boundary conditions are imposed by means of a reference solution.
 * Discretization is performed on a staggered grid and in x-direction only, meaning that variables are located in position DOWN and UP.
 * Linear problem is solved by means of a PETSc KSP solver, with GMRES and Jacobi preconditioner.
 */

#ifndef PARABOLIC_PROBLEM_Y_HPP
#define PARABOLIC_PROBLEM_Y_HPP

class parabolic_problem_y {
private:
    DM dmGrid;    ///< Discretized grid for the problem.
    Mat A;        ///< Matrix representing the left-hand side of the system.
    Vec rhs;      ///< Right-hand side vector.
    Vec V_up;     ///< Solution vector for the y-direction. If not provided, it must be passed as reference in solve_step(...)

    Vec mask_V;   ///< Mask vector used for Brinkman flow.

    /**
     * @brief Assembles the right-hand side (RHS) vector.
     */ 
    PetscErrorCode assemble_rhs(PetscReal const & theta, Vec const & V_up);

    /**
     * @brief Assembles the left-hand side (LHS) matrix.
     */
    PetscErrorCode assemble_lhs();

    /**
     * @brief Exports results in .vtk format for post-processing.
     */
    PetscErrorCode exodus(size_t const & i);


public:

    /**
     * @brief Constructor that initializes the problem with a given grid. Use when parabolic problem is just a step of a bigger problem, like in Chorin-Temam method.
     * @param dmGrid staggered grid petsc-object must be already be defined.
     */

    parabolic_problem_y(DM const & dmGrid);

    /**
     * @brief Default constructor for stand-alone problem.
    */

    parabolic_problem_y();

    /**
     * @brief Performs a single time step of the numerical solution.
     */  
    PetscErrorCode solve_step(PetscReal const & theta, std::optional<std::reference_wrapper<Vec>> V_up_opt = std::nullopt);    
    /**
     * @brief Solves the entire parabolic problem in the Y-direction.
     */
    PetscErrorCode solve();
    /**
     * @brief Destructor to clean up allocated resources. Fundamental in PETSc implementation lo avoid leaks and unexpected RAM overhead.
     */
    ~parabolic_problem_y();

};

#endif // PARABOLIC_PROBLEM_Y_HPP

/**
 * @class parabolic_problem_z
 * @brief Represents a parabolic problem in the z-direction.
 *
 * This class solves a generic evolutionary diffusive problem in the z-direction with fully Dirichlet bc's. 
 * Boundary conditions are imposed by means of a reference solution.
 * Discretization is performed on a staggered grid and in x-direction only, meaning that variables are located in position BACK and FRONT.
 * Linear problem is solved by means of a PETSc KSP solver, with GMRES and Jacobi preconditioner.
 */

#ifndef PARABOLIC_PROBLEM_Z_HPP
#define PARABOLIC_PROBLEM_Z_HPP

class parabolic_problem_z {
private:

    DM dmGrid;    ///< Discretized grid for the problem.
    Mat A;        ///< Matrix representing the left-hand side of the system.
    Vec rhs;      ///< Right-hand side vector.
    Vec W_up;     ///< Solution vector for the x-direction. If not provided, it must be passed as reference in solve_step(...)

    Vec mask_W;   ///< Mask vector used for Brinkman flow.
    /**
     * @brief Assembles the right-hand side (RHS) vector.
     */ 
    PetscErrorCode assemble_rhs(PetscReal const & theta, Vec const & W_up);
    /**
     * @brief Assembles the left-hand side (LHS) matrix.
     */
    PetscErrorCode assemble_lhs();
    /**
     * @brief Exports results in .vtk format for post-processing.
     */
    PetscErrorCode exodus(size_t const & i);


public:
    /**
     * @brief Constructor that initializes the problem with a given grid. Use when parabolic problem is just a step of a bigger problem, like in Chorin-Temam method.
     * @param dmGrid staggered grid petsc-object must be already be defined.
     */
    parabolic_problem_z(DM const & dmGrid);
    /**
     * @brief Default constructor for stand-alone problem.
     */
    parabolic_problem_z();

    /**
     * @brief Performs a single time step of the numerical solution.
     */
    PetscErrorCode solve_step(PetscReal const & theta, std::optional<std::reference_wrapper<Vec>> W_up_opt = std::nullopt);
    /**
     * @brief Solves the entire parabolic problem in the Z-direction.
     */
    PetscErrorCode solve();
    /**
     * @brief Destructor to clean up allocated resources. Fundamental in PETSc implementation lo avoid leaks and unexpected RAM overhead.
     */
    ~parabolic_problem_z();

};

#endif // PARABOLIC_PROBLEM_Z_HPP
