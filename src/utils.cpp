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

PetscErrorCode CheckSolution(Vec const & sol, Vec const & solRef, std::string const & comp)
{
    Vec       diff;
    PetscReal normsolRef, errAbs, errRel;
    PetscFunctionBegin;

    PetscFunctionBegin;
    VecDuplicate(sol, &diff);
    VecCopy(sol, diff);
    VecAXPY(diff, -1.0, solRef);
    VecNorm(diff, NORM_2, &errAbs);
    VecNorm(solRef, NORM_2, &normsolRef);
    errRel = errAbs / normsolRef;
    PetscPrintf(PETSC_COMM_WORLD, "Error %s: Absolute = %.6g, Relative = %.6g\n", comp.c_str(), (double)errAbs, (double)errRel);
    VecDestroy(&diff);
    PetscFunctionReturn(0);


}

PetscErrorCode CreateAnalyticalU(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iux, icux[3];
    DM              dmCoord;
    Vec             vecLocal, coord, coordLocal;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, LEFT, 0, &icux[0]);
    DMStagGetLocationSlot(dmCoord, LEFT, 1, &icux[1]);
    DMStagGetLocationSlot(dmCoord, LEFT, 2, &icux[2]);  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iux] = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode CreateAnalyticalV(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iuy, icuy[3];
    Vec             vecLocal, coord, coordLocal;
    DM              dmCoord;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, DOWN, 0, &icuy[0]);
    DMStagGetLocationSlot(dmCoord, DOWN, 1, &icuy[1]);
    DMStagGetLocationSlot(dmCoord, DOWN, 2, &icuy[2]);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iuy] = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode CreateAnalyticalW(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iuz, icuz[3];
    Vec             vecLocal, coord, coordLocal;
    DM              dmCoord;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, BACK, 0, &icuz[0]);
    DMStagGetLocationSlot(dmCoord, BACK, 1, &icuz[1]);
    DMStagGetLocationSlot(dmCoord, BACK, 2, &icuz[2]);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iuz] = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}
/*
PetscErrorCode CreateAnalyticalP(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iux, icux[3], iuy, icuy[3], iuz, icuz[3], iue, icue[3];
    DM              dmCoord;
    Vec             vecLocal, coord, coordLocal;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, LEFT, 0, &icux[0]);
    DMStagGetLocationSlot(dmCoord, LEFT, 1, &icux[1]);
    DMStagGetLocationSlot(dmCoord, LEFT, 2, &icux[2]); 
    DMStagGetLocationSlot(dmCoord, DOWN, 0, &icuy[0]);
    DMStagGetLocationSlot(dmCoord, DOWN, 1, &icuy[1]);
    DMStagGetLocationSlot(dmCoord, DOWN, 2, &icuy[2]);
    DMStagGetLocationSlot(dmCoord, BACK, 0, &icuz[0]);
    DMStagGetLocationSlot(dmCoord, BACK, 1, &icuz[1]);
    DMStagGetLocationSlot(dmCoord, BACK, 2, &icuz[2]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 0, &icue[0]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 1, &icue[1]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 2, &icue[2]);     
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iue);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iue] = pRef(arrCoord[ez][ey][ex][icue[0]], arrCoord[ez][ey][ex][icue[1]], arrCoord[ez][ey][ex][icue[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}
*/

PetscErrorCode CreateGrid(DM * const dmGrid, PetscInt const & dof1, PetscInt const & dof2, PetscInt const & dof3, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz)
{
    const PetscInt dof0 = 0;
    const PetscInt stencilWidth = 1;

    PetscFunctionBegin;

    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, dmGrid);
    DMSetFromOptions(*dmGrid);
    DMSetUp(*dmGrid);
    DMStagSetUniformCoordinatesExplicit(*dmGrid, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    PetscFunctionReturn(0);
}

void PrintSimulationParameters()
{
    PetscPrintf(PETSC_COMM_WORLD, "Simulation Parameters:\n");
    PetscPrintf(PETSC_COMM_WORLD, "------------------------\n");

    // Grid Parameters
    PetscPrintf(PETSC_COMM_WORLD, "Grid Dimensions:\n");
    PetscPrintf(PETSC_COMM_WORLD, "  nx = %d (Number of grid points in x-direction)\n", nx);
    PetscPrintf(PETSC_COMM_WORLD, "  ny = %d (Number of grid points in y-direction)\n", ny);
    PetscPrintf(PETSC_COMM_WORLD, "  nz = %d (Number of grid points in z-direction)\n", nz);

    // Domain Boundaries
    PetscPrintf(PETSC_COMM_WORLD, "\nDomain Extents:\n");
    PetscPrintf(PETSC_COMM_WORLD, "  x-direction: [%.2f, %.2f]\n", Lx_0, Lx);
    PetscPrintf(PETSC_COMM_WORLD, "  y-direction: [%.2f, %.2f]\n", Ly_0, Ly);
    PetscPrintf(PETSC_COMM_WORLD, "  z-direction: [%.2f, %.2f]\n", Lz_0, Lz);

    // Time Parameters
    PetscPrintf(PETSC_COMM_WORLD, "\nTime Parameters:\n");
    PetscPrintf(PETSC_COMM_WORLD, "  dt = %.5f (Time step size)\n", dt);
    PetscPrintf(PETSC_COMM_WORLD, "  iter = %.0f (Number of iterations)\n", iter);
    PetscPrintf(PETSC_COMM_WORLD, "  theta = %.2f (Starting time)\n", theta);


    // Physical Parameters
    PetscPrintf(PETSC_COMM_WORLD, "\nPhysical Parameters:\n");
    PetscPrintf(PETSC_COMM_WORLD, "  Re = %.2f (Reynolds number)\n", Re);

    PetscPrintf(PETSC_COMM_WORLD, "------------------------\n");
    PetscPrintf(PETSC_COMM_WORLD, "End of Parameter List.\n");
}


bool rayIntersectsTriangle(const std::array<double, 3>& rayOrigin, const std::array<double, 3>& rayVector, const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2) {
    const double EPSILON = 1e-8;
    std::array<double, 3> edge1, edge2, h, s, q;
    double a, f, u, v;

    // Calculate edges
    for (int i = 0; i < 3; ++i) {
        edge1[i] = v1[i] - v0[i];
        edge2[i] = v2[i] - v0[i];
    }

    // Calculate determinant
    h[0] = rayVector[1] * edge2[2] - rayVector[2] * edge2[1];
    h[1] = rayVector[2] * edge2[0] - rayVector[0] * edge2[2];
    h[2] = rayVector[0] * edge2[1] - rayVector[1] * edge2[0];

    a = edge1[0] * h[0] + edge1[1] * h[1] + edge1[2] * h[2];

    if (a > -EPSILON && a < EPSILON)
        return false; // This means the ray is parallel to the triangle.

    f = 1.0 / a;

    // Calculate u parameter and test bound
    for (int i = 0; i < 3; ++i) {
        s[i] = rayOrigin[i] - v0[i];
    }

    u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
    if (u < 0.0 || u > 1.0)
        return false;

    // Calculate v parameter and test bound
    q[0] = s[1] * edge1[2] - s[2] * edge1[1];
    q[1] = s[2] * edge1[0] - s[0] * edge1[2];
    q[2] = s[0] * edge1[1] - s[1] * edge1[0];

    v = f * (rayVector[0] * q[0] + rayVector[1] * q[1] + rayVector[2] * q[2]);
    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage, we can compute t to find out where the intersection point is on the line.
    double t = f * (edge2[0] * q[0] + edge2[1] * q[1] + edge2[2] * q[2]);
    if (t > EPSILON) // ray intersection
        return true;

    return false; // No intersection
}

bool isPointInsideMesh(const std::array<double, 3>& point, const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<int, 3>>& faces) {
    // Define a ray direction. We'll use the positive x-direction.
    std::array<double, 3> rayDirection = {1.0, 0.0, 0.0};

    int intersectionCount = 0;

    // Iterate over all faces (triangles)
    for (const auto& face : faces) {
        const std::array<double, 3>& v0 = vertices[face[0]];
        const std::array<double, 3>& v1 = vertices[face[1]];
        const std::array<double, 3>& v2 = vertices[face[2]];

        if (rayIntersectsTriangle(point, rayDirection, v0, v1, v2)) {
            intersectionCount++;
        }
    }

    // If the intersection count is odd, the point is inside; otherwise, it's outside.
    return (intersectionCount % 2 == 1);
}

void reader(const std::string& filename, std::vector<std::array<double, 3>>& vertices, std::vector<std::array<int, 3>>& faces) {
    // Create a reader for the STL file
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // Get the polydata object
    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

    // Extract vertices
    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
    vertices.resize(points->GetNumberOfPoints());

    double x_min = std::numeric_limits<double>::max();
    double y_min = std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();

    double x_max = std::numeric_limits<double>::lowest();
    double y_max = std::numeric_limits<double>::lowest();
    double z_max = std::numeric_limits<double>::lowest();

    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        double p[3];
        points->GetPoint(i, p);
        vertices[i] = {{ p[0], p[1], p[2] }};

        x_min = std::min(x_min, p[0]);
        y_min = std::min(y_min, p[1]);
        z_min = std::min(z_min, p[2]);

        x_max = std::max(x_max, p[0]);
        y_max = std::max(y_max, p[1]);
        z_max = std::max(z_max, p[2]);
    }

    std::cout << "Min X: " << x_min << "\n";
    std::cout << "Max X: " << x_max << "\n";
    std::cout << "Min Y: " << y_min << "\n";
    std::cout << "Max Y: " << y_max << "\n";
    std::cout << "Min Z: " << z_min << "\n";
    std::cout << "Max Z: " << z_max << "\n";

    // Calculate the maximum dimension of the bounding box
    /*double max_dimension = std::max({x_max - x_min, y_max - y_min, z_max - z_min});

    // Normalize vertices using the maximum dimension
    for (auto& v : vertices) {
        v[0] = (v[0] - x_min) / max_dimension;
        v[1] = (v[1] - y_min) / max_dimension;
        v[2] = (v[2] - z_min) / max_dimension;
    }

    // Initialize normalized min and max values
    double norm_min_x = std::numeric_limits<double>::max();
    double norm_max_x = std::numeric_limits<double>::lowest();
    double norm_min_y = std::numeric_limits<double>::max();
    double norm_max_y = std::numeric_limits<double>::lowest();
    double norm_min_z = std::numeric_limits<double>::max();
    double norm_max_z = std::numeric_limits<double>::lowest();

    // Find the normalized min and max values
    for (const auto& v : vertices) {
        norm_min_x = std::min(norm_min_x, v[0]);
        norm_max_x = std::max(norm_max_x, v[0]);

        norm_min_y = std::min(norm_min_y, v[1]);
        norm_max_y = std::max(norm_max_y, v[1]);

        norm_min_z = std::min(norm_min_z, v[2]);
        norm_max_z = std::max(norm_max_z, v[2]);
    }

    // Output the results (replace with desired output method)
    std::cout << "Normalized Min X: " << norm_min_x << "\n";
    std::cout << "Normalized Max X: " << norm_max_x << "\n";
    std::cout << "Normalized Min Y: " << norm_min_y << "\n";
    std::cout << "Normalized Max Y: " << norm_max_y << "\n";
    std::cout << "Normalized Min Z: " << norm_min_z << "\n";
    std::cout << "Normalized Max Z: " << norm_max_z << "\n";*/

    vtkSmartPointer<vtkCellArray> triangles = polyData->GetPolys();

    triangles->InitTraversal();
    vtkIdType npts;
    const vtkIdType *pts;
    while (triangles->GetNextCell(npts, pts)) {
        if (npts == 3) {
            faces.push_back({{ static_cast<int>(pts[0]), static_cast<int>(pts[1]), static_cast<int>(pts[2]) }});
        }
    }
}


PetscErrorCode createMaskU(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces) {

    PetscInt icux_right[3], icux_left[3], iux_right, iux_left;
    PetscInt startx, starty, startz, N[3], ex, ey, ez, nx, ny, nz, d;
    DM dmCoord;
    Vec coord, coordLocal, vec_stag_local;
    PetscReal ****arrCoord, ****arrVecStag;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
    }  

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);

    DMCreateLocalVector(dmGrid, &vec_stag_local);
    DMGlobalToLocalBegin(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMGlobalToLocalEnd(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMStagVecGetArray(dmGrid, vec_stag_local, &arrVecStag);

    std::vector<std::vector<PetscInt>> icux;
    icux.push_back({icux_right[0], icux_right[1], icux_right[2]});
    icux.push_back({icux_left[0], icux_left[1], icux_left[2]});

    std::vector<PetscInt> iux = {iux_right, iux_left};

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                for(auto i : icux) {
                    for (auto j : iux){
                    std::array<double, 3> point = {arrCoord[ez][ey][ex][i[0]], arrCoord[ez][ey][ex][i[1]], arrCoord[ez][ey][ex][i[2]]};
                    if (isPointInsideMesh(point, vertices, faces)) {
                        arrVecStag[ez][ey][ex][j] = 0.0;
                    }
                    else {
                        arrVecStag[ez][ey][ex][j] = eps;
                        }
                    }
                }

            }
        }
    }
    
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArray(dmGrid, vec_stag_local, &arrVecStag);
    DMLocalToGlobal(dmGrid, vec_stag_local, INSERT_VALUES, vec_stag);
    DMRestoreLocalVector(dmGrid, &vec_stag_local);

    PetscFunctionReturn(0);
}

PetscErrorCode createMaskV(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces) {

    PetscInt icuy_up[3], icuy_down[3], iuy_up, iuy_down;
    PetscInt startx, starty, startz, N[3], ex, ey, ez, nx, ny, nz, d;
    DM dmCoord;
    Vec coord, coordLocal, vec_stag_local;
    PetscReal ****arrCoord, ****arrVecStag;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy_down[d]);
    }  

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);

    DMCreateLocalVector(dmGrid, &vec_stag_local);
    DMGlobalToLocalBegin(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMGlobalToLocalEnd(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMStagVecGetArray(dmGrid, vec_stag_local, &arrVecStag);

    std::vector<std::vector<PetscInt>> icuy;
    icuy.push_back({icuy_up[0], icuy_up[1], icuy_up[2]});
    icuy.push_back({icuy_down[0], icuy_down[1], icuy_down[2]});

    std::vector<PetscInt> iuy = {iuy_up, iuy_down};

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                for(auto i : icuy) {
                    for (auto j : iuy){
                    std::array<double, 3> point = {arrCoord[ez][ey][ex][i[0]], arrCoord[ez][ey][ex][i[1]], arrCoord[ez][ey][ex][i[2]]};
                    if (isPointInsideMesh(point, vertices, faces)) {
                        arrVecStag[ez][ey][ex][j] = 0.0;
                    }
                    else {
                        arrVecStag[ez][ey][ex][j] = eps;
                        }
                    }
                }

            }
        }
    }
    
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArray(dmGrid, vec_stag_local, &arrVecStag);
    DMLocalToGlobal(dmGrid, vec_stag_local, INSERT_VALUES, vec_stag);
    DMRestoreLocalVector(dmGrid, &vec_stag_local);

    PetscFunctionReturn(0);
}

PetscErrorCode createMaskW(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces) {

    PetscInt icuz_front[3], icuz_back[3], iuz_front, iuz_back;
    PetscInt startx, starty, startz, N[3], ex, ey, ez, nx, ny, nz, d;
    DM dmCoord;
    Vec coord, coordLocal, vec_stag_local;
    PetscReal ****arrCoord, ****arrVecStag;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz_back[d]);
    }  

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);

    DMCreateLocalVector(dmGrid, &vec_stag_local);
    DMGlobalToLocalBegin(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMGlobalToLocalEnd(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMStagVecGetArray(dmGrid, vec_stag_local, &arrVecStag);

    std::vector<std::vector<PetscInt>> icuz;
    icuz.push_back({icuz_front[0], icuz_front[1], icuz_front[2]});
    icuz.push_back({icuz_back[0], icuz_back[1], icuz_back[2]});

    std::vector<PetscInt> iuz = {iuz_front, iuz_back};

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                for(auto i : icuz) {
                    for (auto j : iuz){
                    std::array<double, 3> point = {arrCoord[ez][ey][ex][i[0]], arrCoord[ez][ey][ex][i[1]], arrCoord[ez][ey][ex][i[2]]};
                    if (isPointInsideMesh(point, vertices, faces)) {
                        arrVecStag[ez][ey][ex][j] = 0.0;
                    }
                    else {
                        arrVecStag[ez][ey][ex][j] = eps;
                        }
                    }
                }

            }
        }
    }
    
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArray(dmGrid, vec_stag_local, &arrVecStag);
    DMLocalToGlobal(dmGrid, vec_stag_local, INSERT_VALUES, vec_stag);
    DMRestoreLocalVector(dmGrid, &vec_stag_local);

    PetscFunctionReturn(0);
}
