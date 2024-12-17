//
// Created by dave on 26/10/23.
//

#include "parameters.hpp"
#include "parabolic.hpp"
#include "petscvec.h" 
#include "navier_stokes.hpp"
#include "Reader.cpp"

PetscErrorCode CreateAnalyticalU(DM const & dmGrid, Vec & vec, PetscReal const & theta, Vec const & penal)
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

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid, penal, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid, penal, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid, penalizationLocal, &arrPenalization);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                if(arrPenalization[ez][ey][ex][iux] == 0.0) {
                    arrVec[ez][ey][ex][iux] = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
                } else {
                    arrVec[ez][ey][ex][iux] = 0.0;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMStagVecRestoreArrayRead(dmGrid, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode CreateAnalyticalV(DM const & dmGrid, Vec & vec, PetscReal const & theta, Vec const & penal)
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

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid, penal, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid, penal, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid, penalizationLocal, &arrPenalization);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                if(arrPenalization[ez][ey][ex][iuy] == 0.0) {
                    arrVec[ez][ey][ex][iuy] = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                } else {
                    arrVec[ez][ey][ex][iuy] = 0.0;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMStagVecRestoreArrayRead(dmGrid, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode CreateAnalyticalW(DM const & dmGrid, Vec & vec, PetscReal const & theta, Vec const & penal)
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

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid, penal, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid, penal, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid, penalizationLocal, &arrPenalization);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                if(arrPenalization[ez][ey][ex][iuz] == 0.0) {
                    arrVec[ez][ey][ex][iuz] = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                } else {
                    arrVec[ez][ey][ex][iuz] = 0.0;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMStagVecRestoreArrayRead(dmGrid, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalizationLocal);

    PetscFunctionReturn(0);
}
// Create Domain routines
PetscErrorCode CreateGrid(DM * dmGrid, PetscInt const & dof1, PetscInt const & dof2, PetscInt const & dof3, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz)
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



int main(int argc, char **argv)
{   

    auto start = std::chrono::high_resolution_clock::now();

    PetscInitialize(&argc, &argv, NULL, NULL);

        if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <surface_cut.stl>" << std::endl;
        return EXIT_FAILURE;
    }

        // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z;//need to declare to due to solving linear system laplacian messing things up
    {
        CreateGrid(&dmGrid_Shifted, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered_x, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMClone(dmGrid_Staggered_x, &dmGrid_Staggered_y);
        DMClone(dmGrid_Staggered_x, &dmGrid_Staggered_z);
    }



    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> faces;
    std::string filename = argv[1];
    reader(filename, vertices, faces);

    /*DM dmGrid, dm_cent;
    CreateGrid(&dm_cent, 0, 0, 1, nx, nx, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);*/

    /*Vec vec_stag, vec_cent;
    DMCreateGlobalVector(dm_cent, &vec_cent);
    DMCreateGlobalVector(dmGrid, &vec_stag);*/

    Vec penalization_U, penalization_V, penalization_W;
    DMCreateGlobalVector(dmGrid_Staggered_x, &penalization_U);
    DMCreateGlobalVector(dmGrid_Staggered_y, &penalization_V);
    DMCreateGlobalVector(dmGrid_Staggered_z, &penalization_W);
    createMaskU(dmGrid_Staggered_x, penalization_U, vertices, faces);
    std::cout << "Mask U created" << std::endl;
    createMaskV(dmGrid_Staggered_y, penalization_V, vertices, faces);
    std::cout << "Mask V created" << std::endl;
    createMaskW(dmGrid_Staggered_z, penalization_W, vertices, faces);
    std::cout << "Mask W created" << std::endl;


   

    Vec U_0, V_0, W_0;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_0);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_0);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_0);
    CreateAnalyticalU(dmGrid_Staggered_x, U_0, 0, penalization_U);
    CreateAnalyticalV(dmGrid_Staggered_y, V_0, 0, penalization_V);
    CreateAnalyticalW(dmGrid_Staggered_z, W_0, 0, penalization_W);



    //staggered


    {  

    //ProblemSetting<Stokes> setting(nx, ny, nz, Lx_0, Ly_0, Lz_0, Lx, Ly, Lz, dt, iter, 1.0/Re);    

    /*stokes_problem stokes(U_0, V_0, W_0, dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, dt, 16, Re);
    stokes.solve();*/
    //setting.~ProblemSetting();
    //stokes.~stokes_problem();


    }

    {
        ProblemSetting<Transport> setting(nx, ny, nz, Lx_0, Ly_0, Lz_0, Lx,  Ly,  Lz, dt, 1);

        navier_stokes_problem navier_stokes(setting, U_0, V_0, W_0, penalization_U, penalization_V, penalization_W, dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, dt, 16, Re);
        navier_stokes.solve();
    }



    VecDestroy(&U_0);
    VecDestroy(&V_0);
    VecDestroy(&W_0);
    
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_x);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_y);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_z); 
    PetscObjectDestroy((PetscObject*)&dmGrid_Shifted);
    PetscObjectDestroy((PetscObject*)&dmGrid_Centered);



    PetscFinalize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Trinity test successfully completed. Ad maiora!" << std::endl;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    PetscFunctionReturn(0); 
}



