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

#include "config_problem.hpp"

using namespace problem_setting;

/**
 * @file utils.hpp
 * @brief Utility functions for grid creation, analytical solutions, and geometric operations.
 *
 * This file contains various utility functions used in numerical simulations, including:
 * - Functions for generating analytical velocity solutions.
 * - Grid creation and setup.
 * - Mesh operations such as point-in-mesh checks and ray-triangle intersections.
 * - Functions for creating Brinkman masks.
 */

#ifndef UTILS_HPP
#define UTILS_HPP

/**
 * @brief Domain length in the x-direction.
 */
PetscReal constexpr D_x{Lx - Lx_0};
/**
 * @brief Domain length in the y-direction.
 */
PetscReal constexpr D_y{Ly - Ly_0};
/**
 * @brief Domain length in the z-direction.
 */
PetscReal constexpr D_z{Lz - Lz_0};

/**
 * @brief List of mesh vertices.
 */
inline std::vector<std::array<double, 3>> vertices;

/**
 * @brief List of mesh faces.
 */
inline std::vector<std::array<int, 3>> faces;
/**
 * @brief Filename of the geometry file.
 */
inline std::string filename;

/**
 * @brief Checks the difference between a computed solution and a reference solution and evaulates L2-norm.
 * @param sol Computed solution vector.
 * @param solRef Reference solution vector.
 */
PetscErrorCode CheckSolution(Vec const & sol, Vec const & solRef, std::string const & comp);
/**
 * @brief Creates an analytical solution on a staggered grid in the x-direction (dofs in position LEFT and RIGHT)
 * @param dmGrid Discretized grid.
 * @param vec Output velocity field.
 */
PetscErrorCode CreateAnalyticalU(DM const & dmGrid, Vec & vec, PetscReal const & theta);
/**
 * @brief Creates an analytical solution on a staggered grid in the y-direction (dofs in position DOWN and UP)
 * @param dmGrid Discretized grid.
 * @param vec Output velocity field.
 */
PetscErrorCode CreateAnalyticalV(DM const & dmGrid, Vec & vec, PetscReal const & theta);
/**
 * @brief Creates an analytical solution on a staggered grid in the z-direction (dofs in position BACK and FRONT)
 * @param dmGrid Discretized grid.
 * @param vec Output velocity field.
 */
PetscErrorCode CreateAnalyticalW(DM const & dmGrid, Vec & vec, PetscReal const & theta);

//PetscErrorCode CreateAnalyticalP(DM const & dmGrid, Vec & vec, PetscReal const & theta);
/**
 * @brief Creates a computational grid for simulation. Allows to specify dofs ONLY on edges, faces and cell-cenetrs. Pass 0 to not allow dofs in some position, 1 allocate dofs.
 * @param dmGrid Pointer to the grid.
 * @param dof1 Degrees of freedom on edges.
 * @param dof2 Degrees of freedom of faces.
 * @param dof3 Degrees of freedom at cell-ceneters.
 */
PetscErrorCode CreateGrid(DM * const dmGrid, PetscInt const & dof1, PetscInt const & dof2, PetscInt const & dof3);

PetscErrorCode CreateAnalyticalP(DM const & dmGrid, Vec & vec, PetscReal const & theta);


/**
 * @brief Prints the simulation parameters to the console.
 */
void PrintSimulationParameters();

/**
 * @brief Checks if a ray intersects with a triangle in 3D space (part of ray-casting algorithm to create Brinkman masks)
 * @param rayOrigin Starting point of the ray.
 * @param rayVector Direction of the ray.
 * @param v0 First vertex of the triangle.
 * @param v1 Second vertex of the triangle.
 * @param v2 Third vertex of the triangle.
 * @return True if the ray intersects the triangle, false otherwise.
 */
bool rayIntersectsTriangle(const std::array<double, 3>& rayOrigin, const std::array<double, 3>& rayVector, const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2);
/**
 * @brief Determines if a point is inside a 3D mesh.
 * @param point The point to be checked.
 * @param vertices List of vertices of the mesh.
 * @param faces List of faces of the mesh.
 * @return True if the point is inside the mesh, false otherwise.
 */
bool isPointInsideMesh(const std::array<double, 3>& point, const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<int, 3>>& faces);
/**
 * @brief Reads a mesh file and extracts vertices and faces.
 * @param filename The name of the mesh file.
 * @param vertices Output list of vertices.
 * @param faces Output list of faces.
 */
void reader(const std::string& filename, std::vector<std::array<double, 3>>& vertices, std::vector<std::array<int, 3>>& faces);
/**
 * @brief Creates a mask for the x-component component.
 * @param dmGrid Discretized grid.
 * @param vec_stag Output mask vector.
 * @param vertices List of mesh vertices.
 * @param faces List of mesh faces.
 */
PetscErrorCode createMaskU(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces);
/**
 * @brief Creates a mask for the y-component component.
 * @param dmGrid Discretized grid.
 * @param vec_stag Output mask vector.
 * @param vertices List of mesh vertices.
 * @param faces List of mesh faces.
 */
PetscErrorCode createMaskV(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces);
/**
 * @brief Creates a mask for the z-component component.
 * @param dmGrid Discretized grid.
 * @param vec_stag Output mask vector.
 * @param vertices List of mesh vertices.
 * @param faces List of mesh faces.
 */
PetscErrorCode createMaskW(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces);


#endif // UTILS_HPP