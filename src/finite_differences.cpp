//
// Created by dave on 26/10/23.
//

#include "parameters.hpp"
#include "parabolic.hpp"
#include "stokes.hpp"
#include "euler.hpp"

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

// First non-linear members: mixed

PetscErrorCode FirstShiftU_y(DM const & dmGrid, Vec & UShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_up_left[3], icux_up_right[3], icux_down_left[3], icux_down_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pUShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icux_up_left[d]);
        DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icux_up_right[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icux_down_left[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icux_down_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid,&l);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,l);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,l);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey != N[1] - 1 and ex != N[0] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = RIGHT;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey + 1;
                    row_next.k = ez;
                    row_next.loc = RIGHT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_RIGHT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);


                }

                if (ey == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN_LEFT;
                    row.c = 0;

                    inter = uxRef(arrCoord[ez][ey][ex][icux_down_left[0]], arrCoord[ez][ey][ex][icux_down_left[1]],
                                  arrCoord[ez][ey][ex][icux_down_left[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_LEFT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = LEFT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey + 1;
                    
                    inter = uxRef(arrCoord[ez][ey][ex][icux_up_left[0]], arrCoord[ez][ey][ex][icux_up_left[1]],
                                  arrCoord[ez][ey][ex][icux_up_left[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == N[1] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_RIGHT;
                    row.c = 0;
                    inter = uxRef(arrCoord[ez][ey][ex][icux_up_right[0]], arrCoord[ez][ey][ex][icux_up_right[1]],
                                  arrCoord[ez][ey][ex][icux_up_right[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN_RIGHT;
                    row.c = 0;
                    inter = uxRef(arrCoord[ez][ey][ex][icux_down_right[0]], arrCoord[ez][ey][ex][icux_down_right[1]],
                                  arrCoord[ez][ey][ex][icux_down_right[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);


                }
            }
        }
    }



    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    VecAssemblyBegin(UShifted);
    VecAssemblyEnd(UShifted);


    PetscObjectDestroy((PetscObject*)&l);

    return 0;


}

PetscErrorCode FirstShiftU_z(DM const & dmGrid, Vec & UShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pUShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
        DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid,&l);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,l);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez != N[2] - 1 and ex != N[0] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = RIGHT;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey;
                    row_next.k = ez + 1;
                    row_next.loc = RIGHT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_RIGHT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_LEFT;
                    row.c = 0;
                    inter = uxRef(arrCoord[ez][ey][ex][icux_back_left[0]], arrCoord[ez][ey][ex][icux_back_left[1]],
                                  arrCoord[ez][ey][ex][icux_back_left[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_LEFT;
                    row.c = 0;
                    inter = uxRef(arrCoord[ez][ey][ex][icux_front_left[0]], arrCoord[ez][ey][ex][icux_front_left[1]],
                                  arrCoord[ez][ey][ex][icux_front_left[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == N[2] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_RIGHT;
                    row.c = 0;
                    inter = uxRef(arrCoord[ez][ey][ex][icux_front_right[0]], arrCoord[ez][ey][ex][icux_front_right[1]],
                                  arrCoord[ez][ey][ex][icux_front_right[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_RIGHT;
                    row.c = 0;
                    inter = uxRef(arrCoord[ez][ey][ex][icux_back_right[0]], arrCoord[ez][ey][ex][icux_back_right[1]],
                                  arrCoord[ez][ey][ex][icux_back_right[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(UShifted);
    VecAssemblyEnd(UShifted);

    PetscObjectDestroy((PetscObject*)&l);

    return 0;
}

PetscErrorCode FirstShiftV_y(DM const & dmGrid, Vec & VShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up_left[3], icuy_up_right[3], icuy_down_left[3], icuy_down_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pVShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icuy_up_left[d]);
        DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icuy_up_right[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icuy_down_left[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icuy_down_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid, &l);
    DMGlobalToLocalBegin(dmGrid, solRef, INSERT_VALUES, l);
    DMGlobalToLocalEnd(dmGrid, solRef, INSERT_VALUES, l);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ex != N[0] - 1 and ey != N[1] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = UP;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex + 1;
                    row_next.j = ey;
                    row_next.k = ez;
                    row_next.loc = UP;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_RIGHT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN_LEFT;
                    row.c = 0;
                    inter = uyRef(arrCoord[ez][ey][ex][icuy_down_left[0]], arrCoord[ez][ey][ex][icuy_down_left[1]],
                                  arrCoord[ez][ey][ex][icuy_down_left[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_LEFT;
                    row.c = 0;
                    inter = uyRef(arrCoord[ez][ey][ex][icuy_up_left[0]], arrCoord[ez][ey][ex][icuy_up_left[1]],
                                  arrCoord[ez][ey][ex][icuy_up_left[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == N[1] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_RIGHT;
                    row.c = 0;
                    inter = uyRef(arrCoord[ez][ey][ex][icuy_up_right[0]], arrCoord[ez][ey][ex][icuy_up_right[1]],
                                  arrCoord[ez][ey][ex][icuy_up_right[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN_RIGHT;
                    row.c = 0;
                    inter = uyRef(arrCoord[ez][ey][ex][icuy_down_right[0]], arrCoord[ez][ey][ex][icuy_down_right[1]],
                                  arrCoord[ez][ey][ex][icuy_down_right[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VShifted);
    VecAssemblyEnd(VShifted);

    PetscObjectDestroy((PetscObject*)&l);



    return 0;
}

PetscErrorCode FirstShiftW_z(DM const & dmGrid, Vec & WShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pWShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
        DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid,&l);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,l);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez != N[2] - 1 and ex != N[0] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = FRONT;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex + 1;
                    row_next.j = ey;
                    row_next.k = ez;
                    row_next.loc = FRONT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_RIGHT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_LEFT;
                    row.c = 0;
                    inter = uzRef(arrCoord[ez][ey][ex][icux_back_left[0]], arrCoord[ez][ey][ex][icux_back_left[1]],
                                  arrCoord[ez][ey][ex][icux_back_left[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_LEFT;
                    row.c = 0;
                    inter = uzRef(arrCoord[ez][ey][ex][icux_front_left[0]], arrCoord[ez][ey][ex][icux_front_left[1]],
                                  arrCoord[ez][ey][ex][icux_front_left[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == N[2] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_RIGHT;
                    row.c = 0;
                    inter = uzRef(arrCoord[ez][ey][ex][icux_front_right[0]], arrCoord[ez][ey][ex][icux_front_right[1]],
                                  arrCoord[ez][ey][ex][icux_front_right[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_RIGHT;
                    row.c = 0;
                    inter = uzRef(arrCoord[ez][ey][ex][icux_back_right[0]], arrCoord[ez][ey][ex][icux_back_right[1]],
                                  arrCoord[ez][ey][ex][icux_back_right[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WShifted);
    VecAssemblyEnd(WShifted);

    PetscObjectDestroy((PetscObject*)&l);


    return 0;
}


/*
PetscErrorCode FirstShiftU_y(DM const & dmGrid, Vec & UShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_up_left[3], icux_up_right[3], icux_down_left[3], icux_down_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pUShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icux_up_left[d]);
        DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icux_up_right[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icux_down_left[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icux_down_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid,&l);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,l);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,l);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey != N[1] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = LEFT;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey + 1;
                    row_next.k = ez;
                    row_next.loc = LEFT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_LEFT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN_LEFT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = LEFT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey + 1;
                    row_second.k = ez;
                    row_second.loc = LEFT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey + 2;
                    row_third.k = ez;
                    row_third.loc = LEFT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }



                if (ey == N[1] - 1) {

                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_LEFT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = LEFT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey - 1;
                    row_second.k = ez;
                    row_second.loc = LEFT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey - 2;
                    row_third.k = ez;
                    row_third.loc = LEFT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }


                if (ex == N[0] - 1) {

                    if(ey != N[1] - 1)
                    {
                        PetscScalar inter, next, current;
                        DMStagStencil row_current;
                        row_current.i = ex;
                        row_current.j = ey;
                        row_current.k = ez;
                        row_current.loc = RIGHT;
                        row_current.c = 0;
                        DMStagStencil row_next;
                        row_next.i = ex;
                        row_next.j = ey + 1;
                        row_next.k = ez;
                        row_next.loc = RIGHT;
                        row_next.c = 0;
                        DMStagStencil row;
                        row.i = ex;
                        row.j = ey;
                        row.k = ez;
                        row.loc = UP_RIGHT;
                        row.c = 0;

                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                        inter = (next + current) / 2.0;
                        DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                    
                    }

                    if(ey == 0){
                        DMStagStencil row;
                        PetscScalar inter;
                        row.i = ex;
                        row.j = ey;
                        row.k = ez;
                        row.loc = DOWN_RIGHT;
                        row.c = 0;

                        DMStagStencil row_first;
                        PetscScalar first;
                        row_first.i = ex;
                        row_first.j = ey;
                        row_first.k = ez;
                        row_first.loc = RIGHT;
                        row_first.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                        DMStagStencil row_second;
                        PetscScalar second;
                        row_second.i = ex;
                        row_second.j = ey + 1;
                        row_second.k = ez;
                        row_second.loc = RIGHT;
                        row_second.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                        DMStagStencil row_third;
                        PetscScalar third;
                        row_third.i = ex;
                        row_third.j = ey + 2;
                        row_third.k = ez;
                        row_third.loc = RIGHT;
                        row_third.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                        inter = (first + 0.5*second - 0.5*third);
                        DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                    }

                    if(ey == N[1] - 1){
                        DMStagStencil row;
                        PetscScalar inter;
                        row.i = ex;
                        row.j = ey;
                        row.k = ez;
                        row.loc = UP_RIGHT;
                        row.c = 0;

                        DMStagStencil row_first;
                        PetscScalar first;
                        row_first.i = ex;
                        row_first.j = ey;
                        row_first.k = ez;
                        row_first.loc = RIGHT;
                        row_first.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                        DMStagStencil row_second;
                        PetscScalar second;
                        row_second.i = ex;
                        row_second.j = ey - 1;
                        row_second.k = ez;
                        row_second.loc = RIGHT;
                        row_second.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                        DMStagStencil row_third;
                        PetscScalar third;
                        row_third.i = ex;
                        row_third.j = ey - 2;
                        row_third.k = ez;
                        row_third.loc = RIGHT;
                        row_third.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                        inter = (first + 0.5*second - 0.5*third);
                        DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                    }
                }
            }
        }
    }



    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    VecAssemblyBegin(UShifted);
    VecAssemblyEnd(UShifted);


    PetscObjectDestroy((PetscObject*)&l);

    return 0;


}

PetscErrorCode FirstShiftU_z(DM const & dmGrid, Vec & UShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pUShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
        DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid,&l);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,l);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez != N[2] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = LEFT;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey;
                    row_next.k = ez + 1;
                    row_next.loc = LEFT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_LEFT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_LEFT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = LEFT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey;
                    row_second.k = ez + 1;
                    row_second.loc = LEFT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey;
                    row_third.k = ez + 2;
                    row_third.loc = LEFT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == N[2] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_LEFT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = LEFT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey;
                    row_second.k = ez - 1;
                    row_second.loc = LEFT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey;
                    row_third.k = ez - 2;
                    row_third.loc = LEFT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == N[0] - 1)
                {
                    if(ez != N[2] - 1)
                    {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = RIGHT;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey;
                    row_next.k = ez + 1;
                    row_next.loc = RIGHT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_RIGHT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                    }
                    if(ez == 0)
                    {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_RIGHT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = RIGHT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey;
                    row_second.k = ez + 1;
                    row_second.loc = RIGHT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey;
                    row_third.k = ez + 2;
                    row_third.loc = RIGHT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                    }
                    if(ez == N[2] - 1)
                    {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_RIGHT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = RIGHT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey;
                    row_second.k = ez - 1;
                    row_second.loc = RIGHT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey;
                    row_third.k = ez - 2;
                    row_third.loc = RIGHT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, UShifted, 1, &row, &inter, INSERT_VALUES);
                    }
                }
            }
        }
    }           

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(UShifted);
    VecAssemblyEnd(UShifted);

    PetscObjectDestroy((PetscObject*)&l);

    return 0;
}

PetscErrorCode FirstShiftV_y(DM const & dmGrid, Vec & VShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up_left[3], icuy_up_right[3], icuy_down_left[3], icuy_down_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pVShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icuy_up_left[d]);
        DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icuy_up_right[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icuy_down_left[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icuy_down_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid, &l);
    DMGlobalToLocalBegin(dmGrid, solRef, INSERT_VALUES, l);
    DMGlobalToLocalEnd(dmGrid, solRef, INSERT_VALUES, l);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ex != N[0] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = DOWN;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex + 1;
                    row_next.j = ey;
                    row_next.k = ez;
                    row_next.loc = DOWN;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN_RIGHT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }



                if (ex == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN_LEFT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = DOWN;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex + 1;
                    row_second.j = ey;
                    row_second.k = ez;
                    row_second.loc = DOWN;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex + 2;
                    row_third.j = ey;
                    row_third.k = ez;
                    row_third.loc = DOWN;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }



                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN_RIGHT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = DOWN;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex - 1;
                    row_second.j = ey;
                    row_second.k = ez;
                    row_second.loc = DOWN;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex - 2;
                    row_third.j = ey;
                    row_third.k = ez;
                    row_third.loc = DOWN;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);                 
                }

                if (ey == N[1] - 1) {
                if (ex != N[0] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = UP;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex + 1;
                    row_next.j = ey;
                    row_next.k = ez;
                    row_next.loc = UP;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_RIGHT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_LEFT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = UP;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex + 1;
                    row_second.j = ey;
                    row_second.k = ez;
                    row_second.loc = UP;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex + 2;
                    row_third.j = ey;
                    row_third.k = ez;
                    row_third.loc = UP;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP_RIGHT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = UP;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex - 1;
                    row_second.j = ey;
                    row_second.k = ez;
                    row_second.loc = UP;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex - 2;
                    row_third.j = ey;
                    row_third.k = ez;
                    row_third.loc = UP;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);                 
                }

                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VShifted);
    VecAssemblyEnd(VShifted);

    PetscObjectDestroy((PetscObject*)&l);



    return 0;
}

PetscErrorCode FirstShiftW_z(DM const & dmGrid, Vec & WShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pWShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
        DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid,&l);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,l);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ex != N[0] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = FRONT;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex + 1;
                    row_next.j = ey;
                    row_next.k = ez;
                    row_next.loc = FRONT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_RIGHT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_LEFT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = FRONT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex + 1;
                    row_second.j = ey;
                    row_second.k = ez;
                    row_second.loc = FRONT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex + 2;
                    row_third.j = ey;
                    row_third.k = ez;
                    row_third.loc = FRONT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_RIGHT;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = FRONT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex - 1;
                    row_second.j = ey;
                    row_second.k = ez;
                    row_second.loc = FRONT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex - 2;
                    row_third.j = ey;
                    row_third.k = ez;
                    row_third.loc = FRONT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == 0) {
                    if (ex != N[0] - 1) {
                        PetscScalar inter, next, current;
                        DMStagStencil row_current;
                        row_current.i = ex;
                        row_current.j = ey;
                        row_current.k = ez;
                        row_current.loc = BACK;
                        row_current.c = 0;
                        DMStagStencil row_next;
                        row_next.i = ex + 1;
                        row_next.j = ey;
                        row_next.k = ez;
                        row_next.loc = BACK;
                        row_next.c = 0;
                        DMStagStencil row;
                        row.i = ex;
                        row.j = ey;
                        row.k = ez;
                        row.loc = BACK_RIGHT;
                        row.c = 0;

                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                        inter = (next + current) / 2.0;
                        DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                    }

                    if (ex == 0) {
                        DMStagStencil row;
                        PetscScalar inter;
                        row.i = ex;
                        row.j = ey;
                        row.k = ez;
                        row.loc = BACK_LEFT;
                        row.c = 0;

                        DMStagStencil row_first;
                        PetscScalar first;
                        row_first.i = ex;
                        row_first.j = ey;
                        row_first.k = ez;
                        row_first.loc = BACK;
                        row_first.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                        DMStagStencil row_second;
                        PetscScalar second;
                        row_second.i = ex + 1;
                        row_second.j = ey;
                        row_second.k = ez;
                        row_second.loc = BACK;
                        row_second.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                        DMStagStencil row_third;
                        PetscScalar third;
                        row_third.i = ex + 2;
                        row_third.j = ey;
                        row_third.k = ez;
                        row_third.loc = BACK;
                        row_third.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                        inter = (first + 0.5*second - 0.5*third);
                        DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                    }

                    if (ex == N[0] - 1) {
                        DMStagStencil row;
                        PetscScalar inter;
                        row.i = ex;
                        row.j = ey;
                        row.k = ez;
                        row.loc = BACK_RIGHT;
                        row.c = 0;

                        DMStagStencil row_first;
                        PetscScalar first;
                        row_first.i = ex;
                        row_first.j = ey;
                        row_first.k = ez;
                        row_first.loc = BACK;
                        row_first.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                        DMStagStencil row_second;
                        PetscScalar second;
                        row_second.i = ex - 1;
                        row_second.j = ey;
                        row_second.k = ez;
                        row_second.loc = BACK;
                        row_second.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                        DMStagStencil row_third;
                        PetscScalar third;
                        row_third.i = ex - 2;
                        row_third.j = ey;
                        row_third.k = ez;
                        row_third.loc = BACK;
                        row_third.c = 0;
                        DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                        inter = (first + 0.5*second - 0.5*third);
                        DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                    }
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WShifted);
    VecAssemblyEnd(WShifted);

    PetscObjectDestroy((PetscObject*)&l);


    return 0;
}
*/

static PetscErrorCode FirstDerive_y(DM const & dmSol, Vec & AB_y, Vec const & AB) //ok
{
    Vec             AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hy;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pAB_y);

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hy = 1.0 / N[1];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                DMStagStencil row_down;
                PetscScalar val_down;
                row_down.i = ex;
                row_down.j = ey;
                row_down.k = ez;
                row_down.loc = DOWN_LEFT;
                row_down.c = 0;

                DMStagStencil row_up;
                PetscScalar val_up;
                row_up.i = ex;
                row_up.j = ey;
                row_up.k = ez;
                row_up.loc = UP_LEFT;
                row_up.c = 0;

                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_down, &val_down);
                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_up, &val_up);

                DMStagStencil row;
                PetscScalar der;
                row.i = ex;
                row.j = ey;
                row.k = ez;
                row.loc = LEFT;
                row.c = 0;
                der = (val_up - val_down) / hy;

                DMStagVecSetValuesStencil(dmSol, AB_y, 1, &row, &der, INSERT_VALUES);

                if (ex == N[0] - 1) {
                    DMStagStencil row_down;
                    PetscScalar val_down;
                    row_down.i = ex;
                    row_down.j = ey;
                    row_down.k = ez;
                    row_down.loc = DOWN_RIGHT;
                    row_down.c = 0;

                    DMStagStencil row_up;
                    PetscScalar val_up;
                    row_up.i = ex;
                    row_up.j = ey;
                    row_up.k = ez;
                    row_up.loc = UP_RIGHT;
                    row_up.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_down, &val_down);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_up, &val_up);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    der = (val_up - val_down) / hy;
                    DMStagVecSetValuesStencil(dmSol, AB_y, 1, &row, &der, INSERT_VALUES);
                }


                
                
                /*if (ex == N[0] - 1) {
                    DMStagStencil row_down;
                    PetscScalar val_down;
                    row_down.i = ex;
                    row_down.j = ey;
                    row_down.k = ez;
                    row_down.loc = DOWN_RIGHT;
                    row_down.c = 0;

                    DMStagStencil row_up;
                    PetscScalar val_up;
                    row_up.i = ex;
                    row_up.j = ey;
                    row_up.k = ez;
                    row_up.loc = UP_RIGHT;
                    row_up.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_down, &val_down);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_up, &val_up);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    der = (val_up - val_down) / hy;
                    DMStagVecSetValuesStencil(dmSol, AB_y, 1, &row, &der, INSERT_VALUES);

                } else {
                    DMStagStencil row_down;
                    PetscScalar val_down;
                    row_down.i = ex;
                    row_down.j = ey;
                    row_down.k = ez;
                    row_down.loc = DOWN_LEFT;
                    row_down.c = 0;

                    DMStagStencil row_up;
                    PetscScalar val_up;
                    row_up.i = ex;
                    row_up.j = ey;
                    row_up.k = ez;
                    row_up.loc = UP_LEFT;
                    row_up.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_down, &val_down);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_up, &val_up);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = LEFT;
                    row.c = 0;
                    der = (val_up - val_down) / hy;

                    DMStagVecSetValuesStencil(dmSol, AB_y, 1, &row, &der, INSERT_VALUES);
                }*/
            }
        }
    }

    VecAssemblyBegin(AB_y);
    VecAssemblyEnd(AB_y);


    PetscObjectDestroy((PetscObject*)&AB_local);


    return 0;
}

static PetscErrorCode FirstDerive_z(DM const & dmSol, Vec & AB_z, Vec const & AB) //ok
{
    Vec             AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hz;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pAB_z);

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hz = 1.0 / N[2];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                DMStagStencil row_back;
                PetscScalar val_back;
                row_back.i = ex;
                row_back.j = ey;
                row_back.k = ez;
                row_back.loc = BACK_LEFT;
                row_back.c = 0;

                DMStagStencil row_front;
                PetscScalar val_front;
                row_front.i = ex;
                row_front.j = ey;
                row_front.k = ez;
                row_front.loc = FRONT_LEFT;
                row_front.c = 0;

                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_back, &val_back);
                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_front, &val_front);

                DMStagStencil row;
                PetscScalar der;
                row.i = ex;
                row.j = ey;
                row.k = ez;
                row.loc = LEFT;
                row.c = 0;
                der = (val_front - val_back) / hz;
                DMStagVecSetValuesStencil(dmSol, AB_z, 1, &row, &der, INSERT_VALUES);

                if (ex == N[0] - 1) {
                    DMStagStencil row_back;
                    PetscScalar val_back;
                    row_back.i = ex;
                    row_back.j = ey;
                    row_back.k = ez;
                    row_back.loc = BACK_RIGHT;
                    row_back.c = 0;

                    DMStagStencil row_front;
                    PetscScalar val_front;
                    row_front.i = ex;
                    row_front.j = ey;
                    row_front.k = ez;
                    row_front.loc = FRONT_RIGHT;
                    row_front.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_back, &val_back);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_front, &val_front);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    der = (val_front - val_back) / hz;
                    DMStagVecSetValuesStencil(dmSol, AB_z, 1, &row, &der, INSERT_VALUES);
                }
                
            }
        }
    }

    VecAssemblyBegin(AB_z);
    VecAssemblyEnd(AB_z);


    PetscObjectDestroy((PetscObject*)&AB_local);


    return 0;
}


// Second non-linear members: mixed
PetscErrorCode SecondShiftV_z(DM const & dmGrid, Vec & VShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_front_up[3], icuy_front_down[3], icuy_back_up[3], icuy_back_down[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pVShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icuy_front_up[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icuy_front_down[d]);
        DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icuy_back_up[d]);
        DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icuy_back_down[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid, &l);
    DMGlobalToLocalBegin(dmGrid, solRef, INSERT_VALUES, l);
    DMGlobalToLocalEnd(dmGrid, solRef, INSERT_VALUES, l);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey != N[1] - 1 and ez != N[2] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = UP;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey;
                    row_next.k = ez + 1;
                    row_next.loc = UP;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_UP;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_DOWN;
                    row.c = 0;
                    inter = uyRef(arrCoord[ez][ey][ex][icuy_front_down[0]], arrCoord[ez][ey][ex][icuy_front_down[1]],
                                  arrCoord[ez][ey][ex][icuy_front_down[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == N[2] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_UP;
                    row.c = 0;
                    inter = uyRef(arrCoord[ez][ey][ex][icuy_front_up[0]], arrCoord[ez][ey][ex][icuy_front_up[1]],
                                  arrCoord[ez][ey][ex][icuy_front_up[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == N[1] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_UP;
                    row.c = 0;
                    inter = uyRef(arrCoord[ez][ey][ex][icuy_back_up[0]], arrCoord[ez][ey][ex][icuy_back_up[1]],
                                  arrCoord[ez][ey][ex][icuy_back_up[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_DOWN;
                    row.c = 0;
                    inter = uyRef(arrCoord[ez][ey][ex][icuy_back_down[0]], arrCoord[ez][ey][ex][icuy_back_down[1]],
                                  arrCoord[ez][ey][ex][icuy_back_down[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VShifted);
    VecAssemblyEnd(VShifted);
    PetscObjectDestroy((PetscObject*)&l);
    


    return 0;
}

PetscErrorCode SecondShiftW_z(DM const & dmGrid, Vec & WShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_up[3], icux_back_down[3], icux_front_up[3], icux_front_down[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pWShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icux_front_up[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icux_front_down[d]);
        DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icux_back_up[d]);
        DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icux_back_down[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid,&l);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,l);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez != N[2] - 1 and ey != N[1] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = FRONT;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey + 1;
                    row_next.k = ez;
                    row_next.loc = FRONT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_UP;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_DOWN;
                    row.c = 0;
                    inter = uzRef(arrCoord[ez][ey][ex][icux_back_down[0]], arrCoord[ez][ey][ex][icux_back_down[1]],
                                  arrCoord[ez][ey][ex][icux_back_down[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_DOWN;
                    row.c = 0;
                    inter = uzRef(arrCoord[ez][ey][ex][icux_front_down[0]], arrCoord[ez][ey][ex][icux_front_down[1]],
                                  arrCoord[ez][ey][ex][icux_front_down[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == N[2] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_UP;
                    row.c = 0;
                    inter = uzRef(arrCoord[ez][ey][ex][icux_front_up[0]], arrCoord[ez][ey][ex][icux_front_up[1]],
                                  arrCoord[ez][ey][ex][icux_front_up[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == N[1] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_UP;
                    row.c = 0;
                    inter = uzRef(arrCoord[ez][ey][ex][icux_back_up[0]], arrCoord[ez][ey][ex][icux_back_up[1]],
                                  arrCoord[ez][ey][ex][icux_back_up[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WShifted);
    VecAssemblyEnd(WShifted);

    PetscObjectDestroy((PetscObject*)&l);


    return 0;
}

/*
// Second non-linear members: mixed
PetscErrorCode SecondShiftV_z(DM const & dmGrid, Vec & VShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_front_up[3], icuy_front_down[3], icuy_back_up[3], icuy_back_down[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pVShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icuy_front_up[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icuy_front_down[d]);
        DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icuy_back_up[d]);
        DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icuy_back_down[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid, &l);
    DMGlobalToLocalBegin(dmGrid, solRef, INSERT_VALUES, l);
    DMGlobalToLocalEnd(dmGrid, solRef, INSERT_VALUES, l);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez != N[2] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = UP;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey;
                    row_next.k = ez + 1;
                    row_next.loc = UP;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_UP;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == N[2] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_UP;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = UP;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey;
                    row_second.k = ez - 1;
                    row_second.loc = UP;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey;
                    row_third.k = ez - 2;
                    row_third.loc = UP;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);

                }

                if (ez == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_UP;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = UP;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey;
                    row_second.k = ez + 1;
                    row_second.loc = UP;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey;
                    row_third.k = ez + 2;
                    row_third.loc = UP;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);

                }

                if(ey == 0)
                {
                if (ez != N[2] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = DOWN;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey;
                    row_next.k = ez + 1;
                    row_next.loc = DOWN;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_DOWN;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ez == N[2] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_DOWN;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = DOWN;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey;
                    row_second.k = ez - 1;
                    row_second.loc = DOWN;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey;
                    row_third.k = ez - 2;
                    row_third.loc = DOWN;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);

                }

                if (ez == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_DOWN;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = DOWN;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey;
                    row_second.k = ez + 1;
                    row_second.loc = DOWN;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey;
                    row_third.k = ez + 2;
                    row_third.loc = DOWN;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, VShifted, 1, &row, &inter, INSERT_VALUES);

                }
                }            
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VShifted);
    VecAssemblyEnd(VShifted);
    PetscObjectDestroy((PetscObject*)&l);
    


    return 0;
}

PetscErrorCode SecondShiftW_z(DM const & dmGrid, Vec & WShifted, Vec const & solRef, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_up[3], icux_back_down[3], icux_front_up[3], icux_front_down[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pWShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icux_front_up[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icux_front_down[d]);
        DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icux_back_up[d]);
        DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icux_back_down[d]);
    }

    Vec l;
    DMCreateLocalVector(dmGrid,&l);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,l);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey != N[1] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = FRONT;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey + 1;
                    row_next.k = ez;
                    row_next.loc = FRONT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_UP;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_DOWN;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = FRONT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey + 1;
                    row_second.k = ez;
                    row_second.loc = FRONT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey + 2;
                    row_third.k = ez;
                    row_third.loc = FRONT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == N[1] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT_UP;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = FRONT;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey - 1;
                    row_second.k = ez;
                    row_second.loc = FRONT;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey - 2;
                    row_third.k = ez;
                    row_third.loc = FRONT;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if(ez == 0)
                {
                if (ey != N[1] - 1) {
                    PetscScalar inter, next, current;
                    DMStagStencil row_current;
                    row_current.i = ex;
                    row_current.j = ey;
                    row_current.k = ez;
                    row_current.loc = BACK;
                    row_current.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey + 1;
                    row_next.k = ez;
                    row_next.loc = BACK;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_UP;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == 0) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_DOWN;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = BACK;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey + 1;
                    row_second.k = ez;
                    row_second.loc = BACK;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey + 2;
                    row_third.k = ez;
                    row_third.loc = BACK;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }

                if (ey == N[1] - 1) {
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK_UP;
                    row.c = 0;

                    DMStagStencil row_first;
                    PetscScalar first;
                    row_first.i = ex;
                    row_first.j = ey;
                    row_first.k = ez;
                    row_first.loc = BACK;
                    row_first.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_first, &first);

                    DMStagStencil row_second;
                    PetscScalar second;
                    row_second.i = ex;
                    row_second.j = ey - 1;
                    row_second.k = ez;
                    row_second.loc = BACK;
                    row_second.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_second, &second);

                    DMStagStencil row_third;
                    PetscScalar third;
                    row_third.i = ex;
                    row_third.j = ey - 2;
                    row_third.k = ez;
                    row_third.loc = BACK;
                    row_third.c = 0;
                    DMStagVecGetValuesStencil(dmGrid, l, 1, &row_third, &third);

                    inter = (first + 0.5*second - 0.5*third);
                    DMStagVecSetValuesStencil(dmGrid, WShifted, 1, &row, &inter, INSERT_VALUES);
                }
                
                }









            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WShifted);
    VecAssemblyEnd(WShifted);

    PetscObjectDestroy((PetscObject*)&l);


    return 0;
}
*/
/*
PetscErrorCode SecondDerive_x(DM const & dmGrid, Vec & AB_x, Vec const & vec)
{
    PetscInt iux_up, iux_down, iux_down_left, iux_down_right, iux_up_left, iux_up_right;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    Vec vecLocal, vecOutLocal;
    PetscReal ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0]; 

    DMStagGetLocationSlot(dmGrid, UP, 0, &iux_up);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iux_down);
    DMStagGetLocationSlot(dmGrid, DOWN_LEFT, 0, &iux_down_left);
    DMStagGetLocationSlot(dmGrid, DOWN_RIGHT, 0, &iux_down_right);
    DMStagGetLocationSlot(dmGrid, UP_LEFT, 0, &iux_up_left);
    DMStagGetLocationSlot(dmGrid, UP_RIGHT, 0, &iux_up_right);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                PetscReal inter, next, prev;
                next = arrVec[ez][ey][ex][iux_down_right];
                prev = arrVec[ez][ey][ex][iux_down_left];
                inter = (next - prev) / hx;
                arrOut[ez][ey][ex][iux_down] = inter;

                if (ey == N[1] - 1) {
                    PetscReal inter, next, prev;
                    next = arrVec[ez][ey][ex][iux_up_right];
                    prev = arrVec[ez][ey][ex][iux_up_left];
                    inter = (next - prev) / hx;
                    arrOut[ez][ey][ex][iux_up] = inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, AB_x);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    PetscFunctionReturn(0);
}

PetscErrorCode SecondDerive_z(DM const & dmGrid, Vec & AB_z, Vec const & vec)
{
    PetscInt iuz_up, iuz_down, iuz_front_up, iuz_front_down, iuz_back_up, iuz_back_down;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    Vec vecLocal, vecOutLocal;
    PetscReal ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hz = 1.0 / N[2]; 

    DMStagGetLocationSlot(dmGrid, UP, 0, &iuz_up);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuz_down);
    DMStagGetLocationSlot(dmGrid, FRONT_UP, 0, &iuz_front_up);
    DMStagGetLocationSlot(dmGrid, FRONT_DOWN, 0, &iuz_front_down);
    DMStagGetLocationSlot(dmGrid, BACK_UP, 0, &iuz_back_up);
    DMStagGetLocationSlot(dmGrid, BACK_DOWN, 0, &iuz_back_down);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal inter, next, prev;
                next = arrVec[ez][ey][ex][iuz_front_down];
                prev = arrVec[ez][ey][ex][iuz_back_down];
                inter = (next - prev) / hz;
                arrOut[ez][ey][ex][iuz_down] = inter;

                if (ey == N[1] - 1) {
                    PetscReal inter, next, prev;
                    next = arrVec[ez][ey][ex][iuz_front_up];
                    prev = arrVec[ez][ey][ex][iuz_back_up];
                    inter = (next - prev) / hz;
                    arrOut[ez][ey][ex][iuz_up] = inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, AB_z);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    PetscFunctionReturn(0);
}
*/

static PetscErrorCode SecondDerive_x(DM const & dmSol, Vec & AB_x, Vec const & AB) //ok
{
    Vec            AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hx;
    PetscFunctionBeginUser;
   // DMCreateGlobalVector(dmSol, pAB_x);

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hx = 1.0 / N[0];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                DMStagStencil row_left;
                PetscScalar val_left;
                row_left.i = ex;
                row_left.j = ey;
                row_left.k = ez;
                row_left.loc = DOWN_LEFT;
                row_left.c = 0;

                DMStagStencil row_right;
                PetscScalar val_right;
                row_right.i = ex;
                row_right.j = ey;
                row_right.k = ez;
                row_right.loc = DOWN_RIGHT;
                row_right.c = 0;

                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_left, &val_left);
                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_right, &val_right);

                DMStagStencil row;
                PetscScalar der;
                row.i = ex;
                row.j = ey;
                row.k = ez;
                row.loc = DOWN;
                row.c = 0;
                der = (val_right - val_left) / hx;

                DMStagVecSetValuesStencil(dmSol, AB_x, 1, &row, &der, INSERT_VALUES);

                if (ey == N[1] - 1) {
                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = UP_LEFT;
                    row_left.c = 0;

                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = UP_RIGHT;
                    row_right.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_left, &val_left);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_right, &val_right);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP;
                    row.c = 0;
                    der = (val_right - val_left) / hx;

                    DMStagVecSetValuesStencil(dmSol, AB_x, 1, &row, &der, INSERT_VALUES);
                }

                /*
                if (ey == N[1] - 1) {
                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = UP_LEFT;
                    row_left.c = 0;

                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = UP_RIGHT;
                    row_right.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_left, &val_left);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_right, &val_right);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP;
                    row.c = 0;
                    der = (val_right - val_left) / hx;

                    DMStagVecSetValuesStencil(dmSol, AB_x, 1, &row, &der, INSERT_VALUES);

                } else {
                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = DOWN_LEFT;
                    row_left.c = 0;

                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = DOWN_RIGHT;
                    row_right.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_left, &val_left);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_right, &val_right);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN;
                    row.c = 0;
                    der = (val_right - val_left) / hx;

                    DMStagVecSetValuesStencil(dmSol, AB_x, 1, &row, &der, INSERT_VALUES);
                }*/
            }
        }
    }

    VecAssemblyBegin(AB_x);
    VecAssemblyEnd(AB_x);
    PetscObjectDestroy((PetscObject*)&AB_local);

    return 0;
}

static PetscErrorCode SecondDerive_z(DM const & dmSol, Vec & AB_z, Vec const & AB) //ok
{
    Vec             AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hz;
    PetscFunctionBeginUser;
   // DMCreateGlobalVector(dmSol, pAB_z);

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hz = 1.0 / N[2];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                DMStagStencil row_back;
                PetscScalar val_back;
                row_back.i = ex;
                row_back.j = ey;
                row_back.k = ez;
                row_back.loc = BACK_DOWN;
                row_back.c = 0;

                DMStagStencil row_front;
                PetscScalar val_front;
                row_front.i = ex;
                row_front.j = ey;
                row_front.k = ez;
                row_front.loc = FRONT_DOWN;
                row_front.c = 0;

                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_back, &val_back);
                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_front, &val_front);

                DMStagStencil row;
                PetscScalar der;
                row.i = ex;
                row.j = ey;
                row.k = ez;
                row.loc = DOWN;
                row.c = 0;
                der = (val_front - val_back) / hz;
                DMStagVecSetValuesStencil(dmSol, AB_z, 1, &row, &der, INSERT_VALUES);

                if (ey == N[1] - 1) {
                    DMStagStencil row_back;
                    PetscScalar val_back;
                    row_back.i = ex;
                    row_back.j = ey;
                    row_back.k = ez;
                    row_back.loc = BACK_UP;
                    row_back.c = 0;

                    DMStagStencil row_front;
                    PetscScalar val_front;
                    row_front.i = ex;
                    row_front.j = ey;
                    row_front.k = ez;
                    row_front.loc = FRONT_UP;
                    row_front.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_back, &val_back);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_front, &val_front);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP;
                    row.c = 0;
                    der = (val_front - val_back) / hz;
                    DMStagVecSetValuesStencil(dmSol, AB_z, 1, &row, &der, INSERT_VALUES);

                }
            }
        }
    }

    VecAssemblyBegin(AB_z);
    VecAssemblyEnd(AB_z);
    PetscObjectDestroy((PetscObject*)&AB_local);


    return 0;
}

// Third non-linear members: mixed
/*
PetscErrorCode ThirdDerive_x(DM const & dmGrid, Vec & AB_x, Vec const & vec)
{
    PetscInt iux_front, iux_back, iux_front_left, iux_front_right, iux_back_left, iux_back_right;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    Vec vecLocal, vecOutLocal;
    PetscReal ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0]; 

    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iux_front);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iux_back);
    DMStagGetLocationSlot(dmGrid, FRONT_LEFT, 0, &iux_front_left);
    DMStagGetLocationSlot(dmGrid, FRONT_RIGHT, 0, &iux_front_right);
    DMStagGetLocationSlot(dmGrid, BACK_LEFT, 0, &iux_back_left);
    DMStagGetLocationSlot(dmGrid, BACK_RIGHT, 0, &iux_back_right);  

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal inter, next, prev;
                next = arrVec[ez][ey][ex][iux_back_right];
                prev = arrVec[ez][ey][ex][iux_back_left];
                inter = (next - prev) / hx;
                arrOut[ez][ey][ex][iux_back] = inter;

                if (ez == N[2] - 1) {
                    PetscReal inter, next, prev;
                    next = arrVec[ez][ey][ex][iux_front_right];
                    prev = arrVec[ez][ey][ex][iux_front_left];
                    inter = (next - prev) / hx;
                    arrOut[ez][ey][ex][iux_front] = inter;
                } 
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, AB_x);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    PetscFunctionReturn(0);
}

PetscErrorCode ThirdDerive_y(DM const & dmGrid, Vec & AB_y, Vec const & vec)
{
    PetscInt iuy_front, iuy_back, iuy_front_up, iuy_front_down, iuy_back_up, iuy_back_down;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    Vec vecLocal, vecOutLocal;
    PetscReal ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hy = 1.0 / N[1]; 

    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuy_front);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuy_back);
    DMStagGetLocationSlot(dmGrid, FRONT_UP, 0, &iuy_front_up);
    DMStagGetLocationSlot(dmGrid, FRONT_DOWN, 0, &iuy_front_down);
    DMStagGetLocationSlot(dmGrid, BACK_UP, 0, &iuy_back_up);
    DMStagGetLocationSlot(dmGrid, BACK_DOWN, 0, &iuy_back_down);    

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal inter, next, prev;
                next = arrVec[ez][ey][ex][iuy_back_up];
                prev = arrVec[ez][ey][ex][iuy_back_down];
                inter = (next - prev) / hy;
                arrOut[ez][ey][ex][iuy_back] = inter;

                if (ez == N[2] - 1) {
                    PetscReal inter, next, prev;
                    next = arrVec[ez][ey][ex][iuy_front_up];
                    prev = arrVec[ez][ey][ex][iuy_front_down];
                    inter = (next - prev) / hy;
                    arrOut[ez][ey][ex][iuy_front] = inter;
                } 
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, AB_y);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    PetscFunctionReturn(0);
}
*/

static PetscErrorCode ThirdDerive_x(DM const & dmSol, Vec & AB_x, Vec const & AB) //ok
{
    Vec             AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hx;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pAB_x);

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hx = 1.0 / N[0];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                DMStagStencil row_left;
                PetscScalar val_left;
                row_left.i = ex;
                row_left.j = ey;
                row_left.k = ez;
                row_left.loc = BACK_LEFT;
                row_left.c = 0;

                DMStagStencil row_right;
                PetscScalar val_right;
                row_right.i = ex;
                row_right.j = ey;
                row_right.k = ez;
                row_right.loc = BACK_RIGHT;
                row_right.c = 0;

                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_left, &val_left);
                DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_right, &val_right);

                DMStagStencil row;
                PetscScalar der;
                row.i = ex;
                row.j = ey;
                row.k = ez;
                row.loc = BACK;
                row.c = 0;
                der = (val_right - val_left) / hx;

                DMStagVecSetValuesStencil(dmSol, AB_x, 1, &row, &der, INSERT_VALUES);


                if (ez == N[2] - 1) {
                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = FRONT_LEFT;
                    row_left.c = 0;

                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = FRONT_RIGHT;
                    row_right.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_left, &val_left);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_right, &val_right);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT;
                    row.c = 0;
                    der = (val_right - val_left) / hx;

                    DMStagVecSetValuesStencil(dmSol, AB_x, 1, &row, &der, INSERT_VALUES);

                } 
            }
        }
    }

    VecAssemblyBegin(AB_x);
    VecAssemblyEnd(AB_x);
    PetscObjectDestroy((PetscObject*)&AB_local);

    return 0;
}

static PetscErrorCode ThirdDerive_y(DM const & dmSol, Vec & AB_y, Vec const & AB) //ok
{
    Vec             AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hy;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pAB_y);

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hy = 1.0 / N[1];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                    DMStagStencil row_down;
                    PetscScalar val_down;
                    row_down.i = ex;
                    row_down.j = ey;
                    row_down.k = ez;
                    row_down.loc = BACK_DOWN;
                    row_down.c = 0;

                    DMStagStencil row_up;
                    PetscScalar val_up;
                    row_up.i = ex;
                    row_up.j = ey;
                    row_up.k = ez;
                    row_up.loc = BACK_UP;
                    row_up.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_down, &val_down);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_up, &val_up);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK;
                    row.c = 0;
                    der = (val_up - val_down) / hy;

                    DMStagVecSetValuesStencil(dmSol, AB_y, 1, &row, &der, INSERT_VALUES);
                if (ez == N[2] - 1) {
                    
                    DMStagStencil row_down;
                    PetscScalar val_down;
                    row_down.i = ex;
                    row_down.j = ey;
                    row_down.k = ez;
                    row_down.loc = FRONT_DOWN;
                    row_down.c = 0;

                    DMStagStencil row_up;
                    PetscScalar val_up;
                    row_up.i = ex;
                    row_up.j = ey;
                    row_up.k = ez;
                    row_up.loc = FRONT_UP;
                    row_up.c = 0;

                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_down, &val_down);
                    DMStagVecGetValuesStencil(dmSol, AB_local, 1, &row_up, &val_up);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT;
                    row.c = 0;
                    der = (val_up - val_down) / hy;
                    DMStagVecSetValuesStencil(dmSol, AB_y, 1, &row, &der, INSERT_VALUES);

                } 
            }
        }
    }

    VecAssemblyBegin(AB_y);
    VecAssemblyEnd(AB_y);
    PetscObjectDestroy((PetscObject*)&AB_local);

    return 0;
}

PetscErrorCode CenterU(DM const & dmGrid, Vec & UCenter, Vec const & vec, PetscReal const & theta)
{
    PetscInt icux_left[3], icux_right[3], iux_left, iux_right, iux_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, vecULocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrU;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iux_element);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecULocal);
    DMStagVecGetArray(dmGrid, vecULocal, &arrU);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ex == N[0] - 1) {
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex][iux_left];
                    //next = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]], theta);
                    next = arrVec[ez][ey][ex][iux_right];
                    inter = ((prev + next) / 2.0) * ((prev + next) / 2.0);
                    arrU[ez][ey][ex][iux_element] = inter;
                } else if(ex == 0) {
                    PetscReal inter, prev, next;
                    //prev = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]], arrCoord[ez][ey][ex][icux_left[2]], theta);
                    prev = arrVec[ez][ey][ex][iux_left];
                    next = arrVec[ez][ey][ex][iux_right];
                    inter = ((prev + next) / 2.0) * ((prev + next) / 2.0);
                    arrU[ez][ey][ex][iux_element] = inter;
                } else {
                    PetscReal inter, next, prev;
                    next = arrVec[ez][ey][ex][iux_right];
                    prev = arrVec[ez][ey][ex][iux_left];
                    inter = ((next + prev) / 2.0) * ((next + prev) / 2.0);
                    arrU[ez][ey][ex][iux_element] = inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecULocal, &arrU);
    DMLocalToGlobal(dmGrid, vecULocal, INSERT_VALUES, UCenter);
    DMRestoreLocalVector(dmGrid, &vecULocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode CenterV(DM const & dmGrid, Vec & VCenter, Vec const & vec, PetscReal const & theta) 
{
    PetscInt icuy_down[3], icuy_up[3], iuy_up, iuy_down, iuy_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, vecVLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrV;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy_down[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
    } 
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iuy_element);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecVLocal);
    DMStagVecGetArray(dmGrid, vecVLocal, &arrV); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey == N[1] - 1) {
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex][iuy_down];
                    //next = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]], theta);
                    next = arrVec[ez][ey][ex][iuy_up];
                    inter = ((prev + next) / 2.0) * ((prev + next) / 2.0);
                    arrV[ez][ey][ex][iuy_element] = inter;
                 } else if(ey == 0) {
                    PetscReal inter, prev, next;
                    //prev = uyRef(arrCoord[ez][ey][ex][icuy_down[0]], arrCoord[ez][ey][ex][icuy_down[1]], arrCoord[ez][ey][ex][icuy_down[2]], theta);
                    prev = arrVec[ez][ey][ex][iuy_down];
                    next = arrVec[ez][ey][ex][iuy_up];
                    inter = ((prev + next) / 2.0) * ((prev + next) / 2.0);
                    arrV[ez][ey][ex][iuy_element] = inter;
                } else {
                    PetscReal inter, next, prev;
                    next = arrVec[ez][ey][ex][iuy_up];
                    prev = arrVec[ez][ey][ex][iuy_down];
                    inter = ((next + prev) / 2.0) * ((next + prev) / 2.0);
                    arrV[ez][ey][ex][iuy_element] = inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecVLocal, &arrV);
    DMLocalToGlobal(dmGrid, vecVLocal, INSERT_VALUES, VCenter);
    DMRestoreLocalVector(dmGrid, &vecVLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode CenterW(DM const & dmGrid, Vec & WCenter, Vec const & vec, PetscReal const & theta) 
{
    PetscInt icuz_back[3], icuz_front[3], iuz_back, iuz_front, iuz_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, vecWLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrW;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz_back[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front); 
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iuz_element);
    
    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecWLocal);
    DMStagVecGetArray(dmGrid, vecWLocal, &arrW); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez == N[2] - 1) {
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex][iuz_back];
                    //next = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]], theta);
                    next = arrVec[ez][ey][ex][iuz_front];
                    inter = ((prev + next) / 2.0) * ((prev + next) / 2.0);
                    arrW[ez][ey][ex][iuz_element] = inter;
                } else if(ez == 0) {
                    PetscReal inter, prev, next;
                    //prev = uzRef(arrCoord[ez][ey][ex][icuz_back[0]], arrCoord[ez][ey][ex][icuz_back[1]], arrCoord[ez][ey][ex][icuz_back[2]], theta);
                    prev = arrVec[ez][ey][ex][iuz_back];
                    next = arrVec[ez][ey][ex][iuz_front];
                    inter = ((prev + next) / 2.0) * ((prev + next) / 2.0);
                    arrW[ez][ey][ex][iuz_element] = inter;
                } else {
                    PetscReal inter, next, prev;
                    next = arrVec[ez][ey][ex][iuz_front];
                    prev = arrVec[ez][ey][ex][iuz_back];
                    inter = ((next + prev) / 2.0) * ((next + prev) / 2.0);
                    arrW[ez][ey][ex][iuz_element] = inter;   
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecWLocal, &arrW);
    DMLocalToGlobal(dmGrid, vecWLocal, INSERT_VALUES, WCenter);
    DMRestoreLocalVector(dmGrid, &vecWLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode Derive_x(DM const & dmGrid, Vec & U2_x, Vec const & vec, PetscReal const & theta)
{
    PetscInt icux_element[3], iux_left, iux_right, iux_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, vecULocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrU;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icux_element[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iux_element);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecULocal);
    DMStagVecGetArray(dmGrid, vecULocal, &arrU);
  
    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
               
                if (ex != 0) {
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex - 1][iux_element];
                    next = arrVec[ez][ey][ex][iux_element];
                    inter = (next - prev) / hx;
                    arrU[ez][ey][ex][iux_left] = inter;                  
                }
                if(ex == 0) {
                    /*PetscReal inter, prev, next;
                    next = arrVec[ez][ey][ex][iux_element];
                    prev = uxRef(arrCoord[ez][ey][ex][icux_element[0]] - hx, arrCoord[ez][ey][ex][icux_element[1]], arrCoord[ez][ey][ex][icux_element[2]], theta);
                    inter = (next - prev*prev) / hx;
                    arrU[ez][ey][ex][iux_left] = inter;*/
                    PetscReal first, second, third, der;
                    first = arrVec[ez][ey][ex][iux_element];
                    second = arrVec[ez][ey][ex + 1][iux_element];
                    third = arrVec[ez][ey][ex + 2][iux_element];
                    der = (-2*first + 3*second - third)/hx;
                    arrU[ez][ey][ex][iux_left] = der;
                }        
                if(ex == N[0] - 1){
                    /*PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex][iux_element];
                    next = uxRef(arrCoord[ez][ey][ex][icux_element[0]] + hx, arrCoord[ez][ey][ex][icux_element[1]], arrCoord[ez][ey][ex][icux_element[2]], theta);
                    inter = (next*next - prev) / hx;
                    arrU[ez][ey][ex][iux_right] = inter;*/
                    PetscReal first, second, third, der;
                    first = arrVec[ez][ey][ex][iux_element];
                    second = arrVec[ez][ey][ex - 1][iux_element];
                    third = arrVec[ez][ey][ex - 2][iux_element];
                    der = (-2*first + 3*second - third)/hx;
                    arrU[ez][ey][ex][iux_right] = der;

                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecULocal, &arrU);
    DMLocalToGlobal(dmGrid, vecULocal, INSERT_VALUES, U2_x);
    DMRestoreLocalVector(dmGrid, &vecULocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode Derive_y(DM const & dmGrid, Vec & V2_y, Vec const & vec, PetscReal const & theta)
{
    PetscInt icuy_element[3], iuy_up, iuy_down, iuy_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, vecVLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrV;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hy = 1.0 / N[1];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icuy_element[d]);
    } 
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iuy_element);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecVLocal);
    DMStagVecGetArray(dmGrid, vecVLocal, &arrV); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if(ey != 0){
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey - 1][ex][iuy_element];
                    next = arrVec[ez][ey][ex][iuy_element];
                    inter = (next - prev) / hy;
                    arrV[ez][ey][ex][iuy_down] = inter;
                }

                if(ey == 0) {
                    /*PetscReal inter, prev, next;
                    next = arrVec[ez][ey][ex][iuy_element];
                    prev = uyRef(arrCoord[ez][ey][ex][icuy_element[0]], arrCoord[ez][ey][ex][icuy_element[1]] - hy, arrCoord[ez][ey][ex][icuy_element[2]], theta);
                    inter = (next - prev*prev) / hy;
                    arrV[ez][ey][ex][iuy_down] = inter;*/
                    PetscReal first, second, third, der;
                    first = arrVec[ez][ey][ex][iuy_element];
                    second = arrVec[ez][ey + 1][ex][iuy_element];
                    third = arrVec[ez][ey + 2][ex][iuy_element];
                    der = (-2*first + 3*second - third)/hy;
                    arrV[ez][ey][ex][iuy_down] = der;
                }

                if(ey == N[1] - 1){
                    /*PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex][iuy_element];
                    next = uyRef(arrCoord[ez][ey][ex][icuy_element[0]], arrCoord[ez][ey][ex][icuy_element[1]] + hy, arrCoord[ez][ey][ex][icuy_element[2]], theta);
                    inter = (next*next - prev) / hy;
                    arrV[ez][ey][ex][iuy_up] = inter;*/
                    PetscReal first, second, third, der;
                    first = arrVec[ez][ey][ex][iuy_element];
                    second = arrVec[ez][ey - 1][ex][iuy_element];
                    third = arrVec[ez][ey - 2][ex][iuy_element];
                    der = (-2*first + 3*second - third)/hy;
                    arrV[ez][ey][ex][iuy_up] = der;

                }


            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecVLocal, &arrV);
    DMLocalToGlobal(dmGrid, vecVLocal, INSERT_VALUES, V2_y);
    DMRestoreLocalVector(dmGrid, &vecVLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode Derive_z(DM const & dmGrid, Vec & W2_z, Vec const & vec, PetscReal const & theta)
{
    PetscInt icuz_element[3], iuz_back, iuz_front, iuz_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, vecWLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrW;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icuz_element[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front); 
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iuz_element);
    
    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecWLocal);
    DMStagVecGetArray(dmGrid, vecWLocal, &arrW); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ez != 0) {
                    PetscReal inter, prev, next;
                    prev = arrVec[ez - 1][ey][ex][iuz_element];
                    next = arrVec[ez][ey][ex][iuz_element];
                    inter = (next - prev) / hz;
                    arrW[ez][ey][ex][iuz_back] = inter;
                }
                if(ez == 0) {
                    /*PetscReal inter, prev, next;
                    next = arrVec[ez][ey][ex][iuz_element];
                    prev = uzRef(arrCoord[ez][ey][ex][icuz_element[0]], arrCoord[ez][ey][ex][icuz_element[1]], arrCoord[ez][ey][ex][icuz_element[2]] - hz, theta);
                    inter = (next - prev*prev) / hz;
                    arrW[ez][ey][ex][iuz_back] = inter;*/
                    PetscReal first, second, third, der;
                    first = arrVec[ez][ey][ex][iuz_element];
                    second = arrVec[ez + 1][ey][ex][iuz_element];
                    third = arrVec[ez + 2][ey][ex][iuz_element];
                    der = (-2*first + 3*second - third)/hz;
                    arrW[ez][ey][ex][iuz_back] = der;
                }
                if(ez == N[2] - 1){
                    /*PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex][iuz_element];
                    next = uzRef(arrCoord[ez][ey][ex][icuz_element[0]], arrCoord[ez][ey][ex][icuz_element[1]], arrCoord[ez][ey][ex][icuz_element[2]] + hz, theta);
                    inter = (next*next - prev) / hz;
                    arrW[ez][ey][ex][iuz_front] = inter;*/
                    PetscReal first, second, third, der;
                    first = arrVec[ez][ey][ex][iuz_element];
                    second = arrVec[ez - 1][ey][ex][iuz_element];
                    third = arrVec[ez - 2][ey][ex][iuz_element];
                    der = (-2*first + 3*second - third)/hz;
                    arrW[ez][ey][ex][iuz_front] = der;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecWLocal, &arrW);
    DMLocalToGlobal(dmGrid, vecWLocal, INSERT_VALUES, W2_z);
    DMRestoreLocalVector(dmGrid, &vecWLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

// Assembling advection term
PetscErrorCode ManageAdvection_x(PetscReal const & dt, Vec & U_int, Vec const & U_n, Vec const & V_n, Vec const & W_n, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz, PetscReal const & theta)
{
    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered;
    PetscFunctionBegin;
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    }

    Vec U_shift;
    DMCreateGlobalVector(dmGrid_Shifted, &U_shift);
    DMStagMigrateVec(dmGrid_Staggered, U_n, dmGrid_Shifted, U_shift);
    Vec V_shift;
    DMCreateGlobalVector(dmGrid_Shifted, &V_shift);
    DMStagMigrateVec(dmGrid_Staggered, V_n, dmGrid_Shifted, V_shift);
    Vec W_shift;
    DMCreateGlobalVector(dmGrid_Shifted, &W_shift);
    DMStagMigrateVec(dmGrid_Staggered, W_n, dmGrid_Shifted, W_shift);

    // Managing first mixed non-linear term, please DO NOT remove any comment
    Vec UV_y, UW_z;
    DMCreateGlobalVector(dmGrid_Shifted, &UV_y);
    DMCreateGlobalVector(dmGrid_Shifted, &UW_z);

    {
        Vec U_y, U_z, V_y, W_z, UV, UW;
        DMCreateGlobalVector(dmGrid_Shifted, &U_y);
        DMCreateGlobalVector(dmGrid_Shifted, &U_z);
        DMCreateGlobalVector(dmGrid_Shifted, &V_y);
        DMCreateGlobalVector(dmGrid_Shifted, &W_z);
        DMCreateGlobalVector(dmGrid_Shifted, &UV);
        DMCreateGlobalVector(dmGrid_Shifted, &UW);
        FirstShiftU_y(dmGrid_Shifted, U_y, U_shift, theta);
        FirstShiftU_z(dmGrid_Shifted, U_z, U_shift, theta);
        FirstShiftV_y(dmGrid_Shifted, V_y, V_shift, theta);          
        FirstShiftW_z(dmGrid_Shifted, W_z, W_shift, theta); 
        VecPointwiseMult(UV, U_y, V_y);
        VecPointwiseMult(UW, U_z, W_z);
        FirstDerive_y(dmGrid_Shifted, UV_y, UV);
        FirstDerive_z(dmGrid_Shifted, UW_z, UW);
        VecAXPY(UV_y, 1.0, UW_z);
        VecDestroy(&UV);
        VecDestroy(&UW);
        VecDestroy(&U_y);
        VecDestroy(&U_z);
        VecDestroy(&V_y);
        VecDestroy(&W_z);
        VecDestroy(&U_shift);
        VecDestroy(&V_shift);
        VecDestroy(&W_shift);
    }

    Vec mixedFirst;
    DMCreateGlobalVector(dmGrid_Staggered, &mixedFirst);
    DMStagMigrateVec(dmGrid_Shifted, UV_y, dmGrid_Staggered, mixedFirst); 



    Vec U_center;
    DMCreateGlobalVector(dmGrid_Centered, &U_center);
    DMStagMigrateVec(dmGrid_Staggered, U_n, dmGrid_Centered, U_center);
    Vec U_c;
    DMCreateGlobalVector(dmGrid_Centered, &U_c);
    CenterU(dmGrid_Centered, U_c, U_center, theta);

    /*Vec benchmark;
    DMCreateGlobalVector(dmGrid_prova, &benchmark);
    CreateReferenceSolutionTry(dmGrid_prova, benchmark, 0);
    Vec center;
    DMCreateGlobalVector(dmGrid_prova, &center);
    DMStagMigrateVec(dmGrid_Centered, U_c, dmGrid_prova, center);
    CheckSolution(center, benchmark);*/

                /*PetscViewer viewer_u;
            DM DM_u;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmGrid_Staggered, homoFirst, LEFT, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "ciao");
            char filename_u[50]; 
            sprintf(filename_u, "results/ciao.vtr");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u); */

    Vec U2_x;
    DMCreateGlobalVector(dmGrid_Centered, &U2_x);
    Derive_x(dmGrid_Centered, U2_x, U_c, theta);
    Vec homoFirst;
    DMCreateGlobalVector(dmGrid_Staggered, &homoFirst);
    DMStagMigrateVec(dmGrid_Centered, U2_x, dmGrid_Staggered, homoFirst);

    VecAXPBYPCZ(U_n, -dt, -dt, 1.0, homoFirst, mixedFirst);
    VecCopy(U_n, U_int);

    VecDestroy(&UV_y);
    VecDestroy(&UW_z);
    VecDestroy(&mixedFirst);
    VecDestroy(&U_c);
    VecDestroy(&U2_x);
    VecDestroy(&U_center);
    VecDestroy(&homoFirst);
    DMDestroy(&dmGrid_Shifted);
    DMDestroy(&dmGrid_Centered);
    DMDestroy(&dmGrid_Staggered);

    PetscFunctionReturn(0);
}

PetscErrorCode ManageAdvection_y(PetscReal const & dt, Vec & V_int, Vec const & U_n, Vec const & V_n, Vec const & W_n, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz, PetscReal const & theta)
{
    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered;
    PetscFunctionBegin;
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    }

    Vec U_shift;
    DMCreateGlobalVector(dmGrid_Shifted, &U_shift);
    DMStagMigrateVec(dmGrid_Staggered, U_n, dmGrid_Shifted, U_shift);
    Vec V_shift;
    DMCreateGlobalVector(dmGrid_Shifted, &V_shift);
    DMStagMigrateVec(dmGrid_Staggered, V_n, dmGrid_Shifted, V_shift);
    Vec W_shift;
    DMCreateGlobalVector(dmGrid_Shifted, &W_shift);
    DMStagMigrateVec(dmGrid_Staggered, W_n, dmGrid_Shifted, W_shift);
    
    // Managing second mixed non-linear term, please DO NOT remove any comment
    Vec VU_x, VW_z;
    DMCreateGlobalVector(dmGrid_Shifted, &VU_x);
    DMCreateGlobalVector(dmGrid_Shifted, &VW_z);

    {
        Vec V_x, V_z, U_x, W_z, VU, VW;
        DMCreateGlobalVector(dmGrid_Shifted, &V_x);
        DMCreateGlobalVector(dmGrid_Shifted, &V_z);
        DMCreateGlobalVector(dmGrid_Shifted, &U_x);
        DMCreateGlobalVector(dmGrid_Shifted, &W_z);
        DMCreateGlobalVector(dmGrid_Shifted, &VU);
        DMCreateGlobalVector(dmGrid_Shifted, &VW);
        FirstShiftV_y(dmGrid_Shifted, V_x, V_shift, theta);
        SecondShiftV_z(dmGrid_Shifted, V_z, V_shift, theta);
        FirstShiftU_y(dmGrid_Shifted, U_x, U_shift, theta);            
        SecondShiftW_z(dmGrid_Shifted, W_z, W_shift, theta);

            /*PetscViewer viewer_u;
            DM DM_u;
            DMStagCreateCompatibleDMStag(dmGrid_Shifted, 0, 1, 0, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmGrid_Shifted, W_z, FRONT_UP, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "ciao");
            char filename_u[50]; 
            sprintf(filename_u, "results/ciao.vtr");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Shifted), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u);*/ 
            
        VecPointwiseMult(VU, V_x, U_x);
        VecPointwiseMult(VW, V_z, W_z);
        SecondDerive_x(dmGrid_Shifted, VU_x, VU);
        SecondDerive_z(dmGrid_Shifted, VW_z, VW);
        VecAXPY(VU_x, 1.0, VW_z);
        PetscObjectDestroy((PetscObject*)&VU);
        PetscObjectDestroy((PetscObject*)&VW);
        PetscObjectDestroy((PetscObject*)&V_x);
        PetscObjectDestroy((PetscObject*)&V_z);
        PetscObjectDestroy((PetscObject*)&U_x);
        PetscObjectDestroy((PetscObject*)&W_z);
        PetscObjectDestroy((PetscObject*)&U_shift);
        PetscObjectDestroy((PetscObject*)&V_shift);
        PetscObjectDestroy((PetscObject*)&W_shift);
    }

    Vec mixedSecond;
    DMCreateGlobalVector(dmGrid_Staggered, &mixedSecond);
    DMStagMigrateVec(dmGrid_Shifted, VU_x, dmGrid_Staggered, mixedSecond);

    Vec V_center;
    DMCreateGlobalVector(dmGrid_Centered, &V_center);
    DMStagMigrateVec(dmGrid_Staggered, V_n, dmGrid_Centered, V_center);
    Vec V_c;
    DMCreateGlobalVector(dmGrid_Centered, &V_c);
    CenterV(dmGrid_Centered, V_c, V_center, theta);
    Vec V2_y;
    DMCreateGlobalVector(dmGrid_Centered, &V2_y);
    Derive_y(dmGrid_Centered, V2_y, V_c, theta);
    Vec homoSecond;
    DMCreateGlobalVector(dmGrid_Staggered, &homoSecond);
    DMStagMigrateVec(dmGrid_Centered, V2_y, dmGrid_Staggered, homoSecond);



    VecAXPBYPCZ(V_n, -dt, -dt, 1.0, homoSecond, mixedSecond);
    VecCopy(V_n, V_int);

    VecDestroy(&VU_x);
    VecDestroy(&VW_z);    
    VecDestroy(&mixedSecond);
    VecDestroy(&V_c);
    VecDestroy(&V2_y);
    VecDestroy(&V_center);
    VecDestroy(&homoSecond);
    DMDestroy(&dmGrid_Shifted);
    DMDestroy(&dmGrid_Centered);
    DMDestroy(&dmGrid_Staggered);
    
    PetscFunctionReturn(0);
}

PetscErrorCode ManageAdvection_z(PetscReal const & dt, Vec & W_int, Vec const & U_n, Vec const & V_n, Vec const & W_n, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz, PetscReal const & theta)
{
    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered;
    PetscFunctionBegin;
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    }

    Vec U_shift;
    DMCreateGlobalVector(dmGrid_Shifted, &U_shift);
    DMStagMigrateVec(dmGrid_Staggered, U_n, dmGrid_Shifted, U_shift);
    Vec V_shift;
    DMCreateGlobalVector(dmGrid_Shifted, &V_shift);
    DMStagMigrateVec(dmGrid_Staggered, V_n, dmGrid_Shifted, V_shift);
    Vec W_shift;
    DMCreateGlobalVector(dmGrid_Shifted, &W_shift);
    DMStagMigrateVec(dmGrid_Staggered, W_n, dmGrid_Shifted, W_shift);
    
    // Managing third mixed non-linear term, please DO NOT remove any comment, ESPECCIALLY the "==[] ones"
    Vec WU_x, WV_y;
    DMCreateGlobalVector(dmGrid_Shifted, &WU_x);
    DMCreateGlobalVector(dmGrid_Shifted, &WV_y);

    {
        Vec W_x, W_y, U_x, V_y, WU, WV;
        DMCreateGlobalVector(dmGrid_Shifted, &W_x);
        DMCreateGlobalVector(dmGrid_Shifted, &W_y);
        DMCreateGlobalVector(dmGrid_Shifted, &U_x);
        DMCreateGlobalVector(dmGrid_Shifted, &V_y);
        DMCreateGlobalVector(dmGrid_Shifted, &WU);
        DMCreateGlobalVector(dmGrid_Shifted, &WV);
        FirstShiftW_z(dmGrid_Shifted, W_x, W_shift, theta);// ==FirsShiftW_z
        SecondShiftW_z(dmGrid_Shifted, W_y, W_shift, theta);// ==SecondShiftW_z
        FirstShiftU_z(dmGrid_Shifted, U_x, U_shift, theta);// ==FirsShiftU_z
        SecondShiftV_z(dmGrid_Shifted, V_y, V_shift, theta);// ==SecondShiftV_
        VecPointwiseMult(WU, W_x, U_x);
        VecPointwiseMult(WV, W_y, V_y);
        ThirdDerive_x(dmGrid_Shifted, WU_x, WU);
        ThirdDerive_y(dmGrid_Shifted, WV_y, WV);
        VecAXPY(WU_x, 1.0, WV_y);
        VecDestroy(&W_x);
        VecDestroy(&W_y);
        VecDestroy(&U_x);
        VecDestroy(&V_y);
        VecDestroy(&WU);
        VecDestroy(&WV);
        VecDestroy(&U_shift);
        VecDestroy(&V_shift);
        VecDestroy(&W_shift);
    }

    Vec mixedThird;
    DMCreateGlobalVector(dmGrid_Staggered, &mixedThird);
    DMStagMigrateVec(dmGrid_Shifted, WU_x, dmGrid_Staggered, mixedThird);

    Vec W_center;
    DMCreateGlobalVector(dmGrid_Centered, &W_center);
    DMStagMigrateVec(dmGrid_Staggered, W_n, dmGrid_Centered, W_center);
    Vec W_c;
    DMCreateGlobalVector(dmGrid_Centered, &W_c);
    CenterW(dmGrid_Centered, W_c, W_center, theta);
    Vec W2_z;
    DMCreateGlobalVector(dmGrid_Centered, &W2_z);
    Derive_z(dmGrid_Centered, W2_z, W_c, theta);
    Vec homoThird;
    DMCreateGlobalVector(dmGrid_Staggered, &homoThird);
    DMStagMigrateVec(dmGrid_Centered, W2_z, dmGrid_Staggered, homoThird);

    /*Vec benchmark;
    DMCreateGlobalVector(dmGrid_Staggered, &benchmark);
    CreateReferenceSolutionTry(dmGrid_Staggered, benchmark, 0);

    CheckSolution(homoThird, benchmark);*/

    VecAXPBYPCZ(W_n, -dt, -dt, 1.0, homoThird, mixedThird);
    VecCopy(W_n, W_int);

    Vec benchmark;
    DMCreateGlobalVector(dmGrid_Staggered, &benchmark);
    CreateReferenceSolutionTry(dmGrid_Staggered, benchmark, 0);
    CheckSolution(W_n, benchmark);

    VecDestroy(&WU_x);
    VecDestroy(&WV_y);
    VecDestroy(&mixedThird);
    VecDestroy(&W_c);
    VecDestroy(&W2_z);
    VecDestroy(&W_center);
    VecDestroy(&homoThird);
    DMDestroy(&dmGrid_Shifted);
    DMDestroy(&dmGrid_Centered);
    DMDestroy(&dmGrid_Staggered);

    PetscFunctionReturn(0);
}

// Second step: viscosity members to solve implicit laplacian
PetscErrorCode Assemble_x(DM const & dmGrid, Mat & A, Vec & rhs, Vec const & vec, PetscReal const & dt, PetscReal const & Re, PetscReal const & theta) {

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux[3], icux_right[3];
    PetscReal const Ret = Re/dt;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }

    Vec local;
    DMCreateLocalVector(dmGrid, &local);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, local);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                /* Right Boundary velocity Dirichlet */
                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    PetscReal valRhs;
                    const PetscReal valA = 1.0;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }      
                /* Equation on left face of this element */
                if (ex == 0) {
                    /* Left velocity Dirichlet */
                    DMStagStencil row;
                    PetscReal valRhs;
                    const PetscReal valA = 1.0;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = LEFT;
                    row.c = 0;

                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                } else {
                    /* X-momentum interior equation : (u_xx + u_yy + u_zz) - p_x = f^x */
                    DMStagStencil row, col[7];
                    PetscReal valA[7], valRhs;
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = LEFT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]-hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz, theta);
                            valRhs = -Ret*valRhs - bc_1/(hy*hy) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                              
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] - hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]+hz, theta);
                            valRhs = -Ret*valRhs - bc_1/(hy*hy) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = LEFT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_2;
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] - hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy);                           
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]+hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz, theta);
                            valRhs = -Ret*valRhs -bc_1/(hy*hy) - bc_2/(hz*hz);                         
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                          
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] + hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]] + hz, theta);                            
                            valRhs = -Ret*valRhs -bc_1/(hy*hy) - bc_2/(hz*hz);       
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = LEFT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_2;
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]+hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        }
                    } else if (ez == 0) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = LEFT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = LEFT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = LEFT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = LEFT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = LEFT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        PetscReal bc_1;
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz, theta);
                        valRhs = -Ret*valRhs - bc_1/(hz*hz);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    } else if (ez == N[2] - 1) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = LEFT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = LEFT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = LEFT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = LEFT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = LEFT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        PetscReal bc_1;
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]+hz, theta);
                        valRhs = -Ret*valRhs - bc_1/(hz*hz);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    } else {
                        nEntries = 7;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = LEFT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = LEFT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = LEFT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = LEFT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = LEFT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        col[6].i = ex;
                        col[6].j = ey;
                        col[6].k = ez + 1;
                        col[6].loc = LEFT;
                        col[6].c = 0;
                        valA[6] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        valRhs = -Ret*valRhs;
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    }
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }
                
            }
        }
    }
    VecDestroy(&local);

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);

    PetscFunctionReturn(0); 
}

PetscErrorCode Assemble_y(DM const & dmGrid, Mat & A, Vec & rhs, Vec const & solRef, PetscReal const & dt, PetscReal const & Re, PetscReal const & theta) {

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up[3];
    PetscReal const Ret = Re/dt;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);        
    }

    Vec local;
    DMCreateLocalVector(dmGrid,&local);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,local);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {  
                if (ey == N[1] - 1) {
                    /* Top boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscReal       valRhs;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = UP;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }                
                /* Equation on bottom face of this element */
                if (ey == 0) {
                    /* Bottom boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscReal       valRhs;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = DOWN;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                } else {
                    /* Y-momentum equation, (v_xx + v_yy + v_zz) - p_y = f^y */
                    DMStagStencil row, col[7];
                    PetscReal   valA[7], valRhs;
                    PetscInt      nEntries;
                    row.i   = ex;
                    row.j   = ey;
                    row.k   = ez;
                    row.loc = DOWN;
                    row.c   = 0;
                    if (ex == 0) {
                        if (ez == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            PetscReal bc_1, bc_2;
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz, theta);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else if (ez == N[2] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            PetscReal bc_1, bc_2;
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz, theta);
                            valRhs = -Ret*valRhs -bc_1/(hx*hx) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = DOWN;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz);
                            PetscReal bc_1;
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        }
                    } else if (ex == N[0] - 1) {
                        if (ez == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex - 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            PetscReal bc_1, bc_2;
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz, theta);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                         
                        } else if (ez == N[2] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex - 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            PetscReal bc_1, bc_2;
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz, theta);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex - 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = DOWN;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz);
                            PetscReal bc_1;
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        }
                    } else if (ez == 0) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = DOWN;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = DOWN;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = DOWN;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = DOWN;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = DOWN;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = DOWN;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        PetscReal bc_2;
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                        
                        bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz, theta);
                        valRhs = -Ret*valRhs - bc_2/(hz*hz);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);

                    } else if (ez == N[2] - 1) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = DOWN;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = DOWN;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = DOWN;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = DOWN;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = DOWN;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez - 1;
                        col[5].loc = DOWN;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        PetscReal bc_1;
                        bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz, theta);
                        valRhs = -Ret*valRhs - bc_1/(hz*hz);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                    
                    } else {
                        nEntries   = 7;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = DOWN;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = DOWN;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = DOWN;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = DOWN;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = DOWN;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez - 1;
                        col[5].loc = DOWN;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        col[6].i   = ex;
                        col[6].j   = ey;
                        col[6].k   = ez + 1;
                        col[6].loc = DOWN;
                        col[6].c   = 0;
                        valA[6]    = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        valRhs = -Ret*valRhs;
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    }
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }
                
            }
        }
    }

    VecDestroy(&local);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);

    PetscFunctionReturn(0);  

}

PetscErrorCode Assemble_z(DM const & dmGrid, Mat & A, Vec & rhs, Vec const & solRef, PetscReal const & dt, PetscReal const & Re, PetscReal const & theta) {
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuz[3], icuz_front[3];
    PetscReal const Ret = Re/dt;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;

    PetscFunctionBegin;


    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal hx = 1.0 / N[0];
    PetscReal hy = 1.0 / N[1];
    PetscReal hz = 1.0 / N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }

    Vec local;
    DMCreateLocalVector(dmGrid,&local);
    DMGlobalToLocalBegin(dmGrid,solRef,INSERT_VALUES,local);
    DMGlobalToLocalEnd(dmGrid,solRef,INSERT_VALUES,local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {              
                if (ez == N[2] - 1) {
                    /* Front boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscReal       valRhs;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = FRONT;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }                
                /* Equation on back face of this element */
                if (ez == 0) {
                    /* Back boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscReal       valRhs;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = BACK;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                    DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                } else {
                    /* Z-momentum equation, (w_xx + w_yy + w_zz) - p_z = f^z */
                    DMStagStencil row, col[7];
                    PetscReal   valA[7], valRhs;
                    PetscInt      nEntries;
                    row.i   = ex;
                    row.j   = ey;
                    row.k   = ez;
                    row.loc = BACK;
                    row.c   = 0;
                    if (ex == 0) {
                        if (ey == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) - 2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);                            
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES); 
                        } else if (ey == N[1] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] + hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                             
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = BACK;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES); 
                        }
                    } else if (ex == N[0] - 1) {
                        if (ey == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] + hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);       
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES); 
                        } else if (ey == N[1] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] + hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] + hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES); 
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex - 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = BACK;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            PetscReal bc_1;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] + hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES); 
                        }
                    } else if (ey == 0) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey + 1;
                        col[1].k   = ez;
                        col[1].loc = BACK;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex - 1;
                        col[2].j   = ey;
                        col[2].k   = ez;
                        col[2].loc = BACK;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hx * hx);
                        col[3].i   = ex + 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = BACK;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez - 1;
                        col[4].loc = BACK;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = BACK;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        PetscReal bc_2;
                        bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                        valRhs = -Ret*valRhs - bc_2/(hy*hy);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                         
                    } else if (ey == N[1] - 1) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) - 2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = BACK;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex - 1;
                        col[2].j   = ey;
                        col[2].k   = ez;
                        col[2].loc = BACK;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hx * hx);
                        col[3].i   = ex + 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = BACK;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez - 1;
                        col[4].loc = BACK;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = BACK;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        PetscReal bc_2;
                        bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] + hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                        valRhs = -Ret*valRhs - bc_2/(hy*hy);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES); 
                    } else {
                        nEntries   = 7;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = BACK;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = BACK;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = BACK;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = BACK;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez - 1;
                        col[5].loc = BACK;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        col[6].i   = ex;
                        col[6].j   = ey;
                        col[6].k   = ez + 1;
                        col[6].loc = BACK;
                        col[6].c   = 0;
                        valA[6]    = 1.0 / (hz * hz);

                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        valRhs = -Ret*valRhs;
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                         
                    }
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecDestroy(&local);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);

    PetscFunctionReturn(0); 
}

// Assembling Viscosity term
PetscErrorCode ManageViscosity(DM const & dmGrid_x, DM const & dmGrid_y, DM const & dmGrid_z, PetscReal const & dt, PetscReal const & Re, Vec & U_pre, Vec & V_pre, Vec & W_pre, Vec const & U_int, Vec const & V_int, Vec const & W_int, PetscReal const & theta){

    Mat A_x, A_y, A_z;
    Vec rhs_x, rhs_y, rhs_z;
    KSP       ksp_x, ksp_y, ksp_z;
    PC        pc_x, pc_y, pc_z;

    PetscFunctionBegin;

    
    DMCreateGlobalVector(dmGrid_x, &rhs_x);
    DMCreateGlobalVector(dmGrid_y, &rhs_y);
    DMCreateGlobalVector(dmGrid_z, &rhs_z);
    DMCreateMatrix(dmGrid_x, &A_x);
    DMCreateMatrix(dmGrid_y, &A_y);
    DMCreateMatrix(dmGrid_z, &A_z);

    {
    Assemble_x(dmGrid_x, A_x, rhs_x, U_int, dt, Re, theta);

    KSPCreate(PETSC_COMM_WORLD, &ksp_x);
    KSPSetType(ksp_x, KSPCG);
    KSPSetOperators(ksp_x, A_x, A_x);
    KSPGetPC(ksp_x, &pc_x);
    PCSetType(pc_x, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc_x, PETSC_TRUE);
    //PCSetType(pc_x, PCJACOBI);
    KSPSetFromOptions(ksp_x);
    KSPSolve(ksp_x, rhs_x, U_pre);

    KSPDestroy(&ksp_x);
    MatDestroy(&A_x);
    VecDestroy(&rhs_x);
    }


    {
    Assemble_y(dmGrid_y, A_y, rhs_y, V_int, dt, Re, theta);

    KSPCreate(PETSC_COMM_WORLD, &ksp_y);
    KSPSetType(ksp_y, KSPCG);
    KSPSetOperators(ksp_y, A_y, A_y);
    KSPGetPC(ksp_y, &pc_y);



    PCSetType(pc_y, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc_y, PETSC_TRUE);
    //PCSetType(pc_x, PCJACOBI);

    KSPSetFromOptions(ksp_y);
    KSPSolve(ksp_y, rhs_y, V_pre);

    KSPDestroy(&ksp_y);
    MatDestroy(&A_y);
    VecDestroy(&rhs_y);
    }


    {
    Assemble_z(dmGrid_z, A_z, rhs_z, W_int, dt, Re, theta);

    KSPCreate(PETSC_COMM_WORLD, &ksp_z);
    KSPSetType(ksp_z, KSPCG);
    KSPSetOperators(ksp_z, A_z, A_z);
    KSPGetPC(ksp_z, &pc_z);
    PCSetType(pc_z, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc_z, PETSC_TRUE);
    //PCSetType(pc_x, PCJACOBI);
    KSPSetFromOptions(ksp_z);
    KSPSolve(ksp_z, rhs_z, W_pre);

    KSPDestroy(&ksp_z);
    MatDestroy(&A_z);
    VecDestroy(&rhs_z);
    }

    PetscFunctionReturn(0);
}


PetscErrorCode AttachNullspace(DM const & dmGrid, Mat & A)
{
  Vec          basis;
  MatNullSpace matNullSpace;
  PetscFunctionBegin;
  DMCreateGlobalVector(dmGrid, &basis);
  VecSet(basis, 1.0);
  MatNullSpaceCreate(PetscObjectComm((PetscObject)dmGrid), PETSC_FALSE, 1, &basis, &matNullSpace);
  VecDestroy(&basis);
  MatSetNullSpace(A, matNullSpace);
  MatNullSpaceDestroy(&matNullSpace);
  PetscFunctionReturn(0);
}

PetscErrorCode Assemble_P(DM const & dmGrid, Mat & A, Vec & rhs, Vec const & source) 
{
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];

    Vec local;
    DMCreateLocalVector(dmGrid,&local);
    DMGlobalToLocalBegin(dmGrid,source,INSERT_VALUES,local);
    DMGlobalToLocalEnd(dmGrid,source,INSERT_VALUES,local);

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {

                    DMStagStencil row, col[6];
                    PetscReal valA[6], valRhs;
                    PetscInt nEntries;

                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {                  
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx); 
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                           
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);                            
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                          
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                      
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        }
                    } else if (ez == 0) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    } else if (ez == N[2] - 1) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    } else {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = LEFT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    }
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                } 
                else if (ex == 0) {
                    DMStagStencil row, col[6];
                    PetscReal valA[6], valRhs;
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                               
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                          
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                      
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        }
                    } else if (ez == 0) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    } else if (ez == N[2] - 1) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    } else {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = LEFT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    }
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                } else {
                    /* X-momentum interior equation : (u_xx + u_yy + u_zz) - p_x = f^x */
                    DMStagStencil row, col[7];
                    PetscReal valA[7], valRhs;
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                               
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                          
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                      
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            /* Missing up term */
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        }
                    } else if (ez == 0) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    } else if (ez == N[2] - 1) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    } else {
                        nEntries = 7;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = LEFT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        col[6].i = ex;
                        col[6].j = ey;
                        col[6].k = ez + 1;
                        col[6].loc = ELEMENT;
                        col[6].c = 0;
                        valA[6] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    }
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }                
            }
        }
    }
    VecDestroy(&local);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);

    PetscFunctionReturn(0); 
}

PetscErrorCode AssembleDivergence(DM const & dmGrid, Vec & div, Vec const & U, Vec const &  V, Vec const & W, PetscReal const & dt) 
{
    PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iu_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iu_right);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iu_up);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iu_down);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iu_front);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iu_back);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iu_element);

    DMCreateLocalVector(dmGrid, &vecULocal);
    DMGlobalToLocalBegin(dmGrid, U, INSERT_VALUES, vecULocal);
    DMGlobalToLocalEnd(dmGrid, U, INSERT_VALUES, vecULocal);
    DMStagVecGetArrayRead(dmGrid, vecULocal, &arrU);

    DMCreateLocalVector(dmGrid, &vecVLocal);
    DMGlobalToLocalBegin(dmGrid, V, INSERT_VALUES, vecVLocal);
    DMGlobalToLocalEnd(dmGrid, V, INSERT_VALUES, vecVLocal);
    DMStagVecGetArrayRead(dmGrid, vecVLocal, &arrV);

    DMCreateLocalVector(dmGrid, &vecWLocal);
    DMGlobalToLocalBegin(dmGrid, W, INSERT_VALUES, vecWLocal);
    DMGlobalToLocalEnd(dmGrid, W, INSERT_VALUES, vecWLocal);
    DMStagVecGetArrayRead(dmGrid, vecWLocal, &arrW);    

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal inter, left, right, up, down, front, back;
                left = arrU[ez][ey][ex][iu_left];
                right = arrU[ez][ey][ex][iu_right];
                up = arrV[ez][ey][ex][iu_up];
                down = arrV[ez][ey][ex][iu_down];
                front = arrW[ez][ey][ex][iu_front];
                back = arrW[ez][ey][ex][iu_back];
                inter = ((up - down) / hy + (right - left) / hx + (front - back) / hz)/dt;
                arrOut[ez][ey][ex][iu_element] = inter;
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecULocal, &arrU);
    DMStagVecRestoreArrayRead(dmGrid, vecVLocal, &arrV);
    DMStagVecRestoreArrayRead(dmGrid, vecWLocal, &arrW);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, div);    
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecULocal);
    DMRestoreLocalVector(dmGrid, &vecVLocal);
    DMRestoreLocalVector(dmGrid, &vecWLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode ComputeDivergence(DM const & dmGrid_centered, DM const & dmGrid_shifted, DM const & dmGrid_staggered, Vec & div, Vec const & U_n, Vec const & V_n, Vec const & W_n, PetscReal const & dt) 
{
        PetscFunctionBegin;

        Vec U_shifted, V_shifted, W_shifted;
        DMCreateGlobalVector(dmGrid_shifted, &U_shifted);
        DMCreateGlobalVector(dmGrid_shifted, &V_shifted);
        DMCreateGlobalVector(dmGrid_shifted, &W_shifted);
        DMStagMigrateVec(dmGrid_staggered, U_n, dmGrid_shifted, U_shifted);
        DMStagMigrateVec(dmGrid_staggered, V_n, dmGrid_shifted, V_shifted);
        DMStagMigrateVec(dmGrid_staggered, W_n, dmGrid_shifted, W_shifted);

        Vec div_shifted;
        DMCreateGlobalVector(dmGrid_shifted, &div_shifted);
        AssembleDivergence(dmGrid_shifted, div_shifted, U_shifted, V_shifted, W_shifted, dt);
        DMStagMigrateVec(dmGrid_shifted, div_shifted, dmGrid_centered, div);

        VecDestroy(&U_shifted);
        VecDestroy(&V_shifted);
        VecDestroy(&W_shifted);
        VecDestroy(&div_shifted);

        PetscFunctionReturn(0); 
}

PetscErrorCode Derive_x_P(DM const & dmGrid, Vec & P_x, Vec const & vec)
{
    PetscInt iux_left, iux_right, iux_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iux_element);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);
  
    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
               
                if (ex != 0) {
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex - 1][iux_element];
                    next = arrVec[ez][ey][ex][iux_element];
                    inter = (next - prev) / hx;
                    arrOut[ez][ey][ex][iux_left] = inter;                  
                }
                if(ex == 0) {
                    arrOut[ez][ey][ex][iux_left] = 0.0;
                }        
                if(ex == N[0] - 1){
                    arrOut[ez][ey][ex][iux_right] = 0.0;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, P_x);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode Derive_y_P(DM const & dmGrid, Vec & P_y, Vec const & vec)
{
    PetscInt iuy_up, iuy_down, iuy_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hy = 1.0 / N[1];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iuy_element);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if(ey != 0){
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey - 1][ex][iuy_element];
                    next = arrVec[ez][ey][ex][iuy_element];
                    inter = (next - prev) / hy;
                    arrOut[ez][ey][ex][iuy_down] = inter;
                }
                if(ey == 0) {
                    arrOut[ez][ey][ex][iuy_down] = 0.0;
                }
                if(ey == N[1] - 1){
                    arrOut[ez][ey][ex][iuy_up] = 0.0;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, P_y);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode Derive_z_P(DM const & dmGrid, Vec & P_z, Vec const & vec)
{
    PetscInt iuz_back, iuz_front, iuz_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front); 
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iuz_element);
    
    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ez != 0) {
                    PetscReal inter, prev, next;
                    prev = arrVec[ez - 1][ey][ex][iuz_element];
                    next = arrVec[ez][ey][ex][iuz_element];
                    inter = (next - prev) / hz;
                    arrOut[ez][ey][ex][iuz_back] = inter;
                }
                if(ez == 0) {
                    arrOut[ez][ey][ex][iuz_back] = 0.0;
                }
                if(ez == N[2] - 1){
                    arrOut[ez][ey][ex][iuz_front] = 0.0;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, P_z);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0); 
}

//capire se preconditioner va bene ref articolo islamico
PetscErrorCode ManagePressure(DM const & dmGrid_centered, DM const & dmGrid_shifted, DM const & dmGrid_staggered, PetscReal const & dt, Vec & P, Vec const & U_n, Vec const & V_n, Vec const & W_n)
{
    Mat A;
    Vec rhs;
    KSP ksp;
    PC  pc;

    PetscFunctionBegin;

    DMCreateGlobalVector(dmGrid_centered, &rhs);
    DMCreateMatrix(dmGrid_centered, &A);

    Vec div;
    DMCreateGlobalVector(dmGrid_centered, &div);
    ComputeDivergence(dmGrid_centered, dmGrid_shifted, dmGrid_staggered, div, U_n, V_n, W_n, dt); 

    Assemble_P(dmGrid_centered, A, rhs, div);
    VecDestroy(&div);

    AttachNullspace(dmGrid_centered, A);

    /*//questo sistema non e' precondizionato: il preconditioner e' nullo. Bisogna farlo ed e' importante: punto piu' lento del codice
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    //KSPGetPC(ksp, &pc);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, P);*/

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCHYPRE);
    PCHYPRESetType(pc, "pilut");
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, P);

    MatDestroy(&A);
    VecDestroy(&rhs);
    KSPDestroy(&ksp);
       
    PetscFunctionReturn(0); 
}

//sistemare il create compatible dmstag (evito di passare tutti quei parametri) ma non so perche non riesco a farlo
PetscErrorCode ManagePressure_x(DM const & dmGrid_staggered, DM const & dmGrid_centered, DM const & dmGrid_shifted, Vec & P_x, Vec const & P)
{
    PetscFunctionBegin;

    Vec P_shifted;
    DMCreateGlobalVector(dmGrid_shifted, &P_shifted);
    DMStagMigrateVec(dmGrid_centered, P, dmGrid_shifted, P_shifted);
    
    Vec P_x_shifted;
    DMCreateGlobalVector(dmGrid_shifted, &P_x_shifted);
    Derive_x_P(dmGrid_shifted, P_x_shifted, P_shifted);
    DMStagMigrateVec(dmGrid_shifted, P_x_shifted, dmGrid_staggered, P_x);

    VecDestroy(&P_x_shifted);
    VecDestroy(&P_shifted);
    /*DMDestroy(&dmGrid_centered);
    DMDestroy(&dmGrid_shifted);*/

    PetscFunctionReturn(0); 
}

PetscErrorCode ManagePressure_y(DM const & dmGrid_staggered, DM const & dmGrid_centered, DM const & dmGrid_shifted, Vec & P_y, Vec const & P)
{
    PetscFunctionBegin

    Vec P_shifted;
    DMCreateGlobalVector(dmGrid_shifted, &P_shifted);
    DMStagMigrateVec(dmGrid_centered, P, dmGrid_shifted, P_shifted);

    Vec P_y_shifted;
    DMCreateGlobalVector(dmGrid_shifted, &P_y_shifted);
    Derive_y_P(dmGrid_shifted, P_y_shifted, P_shifted);
    DMStagMigrateVec(dmGrid_shifted, P_y_shifted, dmGrid_staggered, P_y);

    VecDestroy(&P_y_shifted);
    VecDestroy(&P_shifted);

    PetscFunctionReturn(0); 
}

PetscErrorCode ManagePressure_z(DM const & dmGrid_staggered, DM const & dmGrid_centered, DM const & dmGrid_shifted, Vec & P_z, Vec const & P)
{
    PetscFunctionBegin

    Vec P_shifted;
    DMCreateGlobalVector(dmGrid_shifted, &P_shifted);
    DMStagMigrateVec(dmGrid_centered, P, dmGrid_shifted, P_shifted);

    Vec P_z_shifted;
    DMCreateGlobalVector(dmGrid_shifted, &P_z_shifted);
    Derive_z_P(dmGrid_shifted, P_z_shifted, P_shifted);
    DMStagMigrateVec(dmGrid_shifted, P_z_shifted, dmGrid_staggered, P_z);

    VecDestroy(&P_z_shifted);
    VecDestroy(&P_shifted);

    PetscFunctionReturn(0); 
}


PetscErrorCode UpdatebcU(DM const & dmGrid, Vec & U_up, PetscReal const & theta)
{
    PetscInt icux_left[3], icux_right[3], iux_left, iux_right;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec; 

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, U_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, U_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ex == N[0] - 1) {
                    PetscReal val;
                    val = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]], theta);
                    arrVec[ez][ey][ex][iux_right] = val;
                    
                } else if(ex == 0) {
                    PetscReal val;
                    val = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]], arrCoord[ez][ey][ex][icux_left[2]], theta);
                    arrVec[ez][ey][ex][iux_left] = val;
                }

            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, U_up);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode UpdatebcV(DM const & dmGrid, Vec & V_up, PetscReal const & theta) 
{

    PetscInt icuy_down[3], icuy_up[3], iuy_down, iuy_up;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy_down[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
    } 
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, V_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, V_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey == N[1] - 1) {
                    PetscReal val;
                    val = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]], theta);
                    arrVec[ez][ey][ex][iuy_up] = val;
                } else if(ey == 0) {
                    PetscReal val;
                    val = uyRef(arrCoord[ez][ey][ex][icuy_down[0]], arrCoord[ez][ey][ex][icuy_down[1]], arrCoord[ez][ey][ex][icuy_down[2]], theta);
                    arrVec[ez][ey][ex][iuy_down] = val;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, V_up);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode UpdatebcW(DM const & dmGrid, Vec & W_up, PetscReal const & theta) 
{

    PetscInt icuz_back[3], icuz_front[3], iuz_back, iuz_front;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz_back[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front); 
    
    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, W_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, W_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez == N[2] - 1) {
                    PetscReal val;
                    val = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]], theta);
                    arrVec[ez][ey][ex][iuz_front] = val;
                } else if(ez == 0) {
                    PetscReal val;
                    val = uzRef(arrCoord[ez][ey][ex][icuz_back[0]], arrCoord[ez][ey][ex][icuz_back[1]], arrCoord[ez][ey][ex][icuz_back[2]], theta);
                    arrVec[ez][ey][ex][iuz_back] = val;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, W_up);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}


PetscErrorCode UpdateVelocity(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, PetscReal const & dt, Vec & U_up, Vec & V_up, Vec & W_up, Vec const & P_x, Vec const & P_y, Vec const & P_z, Vec const & U_pre, Vec const & V_pre, Vec const & W_pre, PetscReal const & theta)
{
    PetscFunctionBegin;

    VecAXPY(U_pre, -dt, P_x);
    VecAXPY(V_pre, -dt, P_y);
    VecAXPY(W_pre, -dt, P_z);

    VecCopy(U_pre, U_up);
    VecCopy(V_pre, V_up);
    VecCopy(W_pre, W_up);

    UpdatebcU(dmGrid_staggered_x, U_up, theta);
    UpdatebcV(dmGrid_staggered_y, V_up, theta);
    UpdatebcW(dmGrid_staggered_z, W_up, theta);


    PetscFunctionReturn(0);
}



PetscErrorCode AssembleMagnitude(DM const & dmGrid, Vec & magnitude, Vec const & U, Vec const & V, Vec const & W) 
{
    PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
    PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iu_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iu_right);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iu_up);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iu_down);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iu_front);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iu_back);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iu_element);

    DMCreateLocalVector(dmGrid, &vecULocal);
    DMGlobalToLocalBegin(dmGrid, U, INSERT_VALUES, vecULocal);
    DMGlobalToLocalEnd(dmGrid, U, INSERT_VALUES, vecULocal);
    DMStagVecGetArrayRead(dmGrid, vecULocal, &arrU);

    DMCreateLocalVector(dmGrid, &vecVLocal);
    DMGlobalToLocalBegin(dmGrid, V, INSERT_VALUES, vecVLocal);
    DMGlobalToLocalEnd(dmGrid, V, INSERT_VALUES, vecVLocal);
    DMStagVecGetArrayRead(dmGrid, vecVLocal, &arrV);

    DMCreateLocalVector(dmGrid, &vecWLocal);
    DMGlobalToLocalBegin(dmGrid, W, INSERT_VALUES, vecWLocal);
    DMGlobalToLocalEnd(dmGrid, W, INSERT_VALUES, vecWLocal);
    DMStagVecGetArrayRead(dmGrid, vecWLocal, &arrW);    

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal inter, left, right, up, down, front, back;
                left = arrU[ez][ey][ex][iu_left];
                right = arrU[ez][ey][ex][iu_right];
                up = arrV[ez][ey][ex][iu_up];
                down = arrV[ez][ey][ex][iu_down];
                front = arrW[ez][ey][ex][iu_front];
                back = arrW[ez][ey][ex][iu_back];
                inter = sqrt(((right + left)*(right + left)) / 4 + ((up + down)*(up + down)) / 4 + ((front + back)*(front + back)) / 4);
                arrOut[ez][ey][ex][iu_element] = inter;
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecULocal, &arrU);
    DMStagVecRestoreArrayRead(dmGrid, vecVLocal, &arrV);
    DMStagVecRestoreArrayRead(dmGrid, vecWLocal, &arrW);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, magnitude);    
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecULocal);
    DMRestoreLocalVector(dmGrid, &vecVLocal);
    DMRestoreLocalVector(dmGrid, &vecWLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);    
}

PetscErrorCode ComputeMagnitude(DM const & dmGrid_Staggered_x, DM const & dmGrid_Staggered_y, DM const & dmGrid_Staggered_z, DM const & dmGrid_Centered, DM const & dmGrid_Shifted, Vec & magnitude, Vec const & U, Vec const & V, Vec const & W)
{
    
    PetscFunctionBegin;

    Vec U_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &U_shifted);
    DMStagMigrateVec(dmGrid_Staggered_x, U, dmGrid_Shifted, U_shifted);
    Vec V_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &V_shifted);
    DMStagMigrateVec(dmGrid_Staggered_y, V, dmGrid_Shifted, V_shifted);
    Vec W_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &W_shifted);
    DMStagMigrateVec(dmGrid_Staggered_z, W, dmGrid_Shifted, W_shifted);


    Vec magnitude_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &magnitude_shifted);
    AssembleMagnitude(dmGrid_Shifted, magnitude_shifted, U_shifted, V_shifted, W_shifted);
    DMStagMigrateVec(dmGrid_Shifted, magnitude_shifted, dmGrid_Centered, magnitude);

    VecDestroy(&magnitude_shifted);
    VecDestroy(&U_shifted);  
    VecDestroy(&V_shifted);
    VecDestroy(&W_shifted);
            
    PetscFunctionReturn(0); 

}

/*
annotazione importante: la condizione di compatibilita' per il toerema di Neumann in teoria non e' rispettata in Chorin Teman
*/

PetscErrorCode laplacian_x_explicit(DM const & dmGrid, Vec & laplacian, Vec const & vec) 
{
    PetscInt iux_left, iux_right;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {

                    if (ey == 0) {
                        if (ez == 0) {   
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_right];
                            second_x = arrVec[ez][ey][ex-1][iux_right];
                            third_x = arrVec[ez][ey][ex-2][iux_right];
                            fourth_x = arrVec[ez][ey][ex-3][iux_right];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_right];
                            second_y = arrVec[ez][ey+1][ex][iux_right];
                            third_y = arrVec[ez][ey+2][ex][iux_right];
                            fourth_y = arrVec[ez][ey+3][ex][iux_right];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_right];    
                            second_z = arrVec[ez+1][ey][ex][iux_right];
                            third_z = arrVec[ez+2][ey][ex][iux_right];
                            fourth_z = arrVec[ez+3][ey][ex][iux_right];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_right] = laplacian;

                            PetscReal prev_x_left, next_x_left, inter_x_left;
                            prev_x_left = arrVec[ez][ey][ex-1][iux_left];
                            next_x_left = arrVec[ez][ey][ex+1][iux_left];
                            inter_x_left = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y_left, second_y_left, third_y_left, fourth_y_left;
                            first_y_left = arrVec[ez][ey][ex][iux_left];
                            second_y_left = arrVec[ez][ey+1][ex][iux_left];
                            third_y_left = arrVec[ez][ey+2][ex][iux_left];
                            fourth_y_left = arrVec[ez][ey+3][ex][iux_left];

                            PetscReal first_z_left, second_z_left, third_z_left, fourth_z_left;
                            first_z_left = arrVec[ez][ey][ex][iux_left];
                            second_z_left = arrVec[ez+1][ey][ex][iux_left];
                            third_z_left = arrVec[ez+2][ey][ex][iux_left];
                            fourth_z_left = arrVec[ez+3][ey][ex][iux_left];

                            PetscReal laplacian_left = (prev_x_left - 2*inter_x_left + next_x_left) / (hx*hx) + (2*first_y_left - 5*second_y_left + 4*third_y_left - fourth_y_left) / (hy*hy) + (2*first_z_left - 5*second_z_left + 4*third_z_left - fourth_z_left) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian_left;                   
                                 
                        } else if (ez == N[2] - 1) {
                            
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_right];
                            second_x = arrVec[ez][ey][ex-1][iux_right];
                            third_x = arrVec[ez][ey][ex-2][iux_right];
                            fourth_x = arrVec[ez][ey][ex-3][iux_right];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_right];
                            second_y = arrVec[ez][ey+1][ex][iux_right];
                            third_y = arrVec[ez][ey+2][ex][iux_right];
                            fourth_y = arrVec[ez][ey+3][ex][iux_right];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_right];    
                            second_z = arrVec[ez-1][ey][ex][iux_right];
                            third_z = arrVec[ez-2][ey][ex][iux_right];
                            fourth_z = arrVec[ez-3][ey][ex][iux_right];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_right] = laplacian; 

                            PetscReal prev_x_left, next_x_left, inter_x_left;
                            prev_x_left = arrVec[ez][ey][ex-1][iux_left];
                            next_x_left = arrVec[ez][ey][ex+1][iux_left];
                            inter_x_left = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y_left, second_y_left, third_y_left, fourth_y_left;
                            first_y_left = arrVec[ez][ey][ex][iux_left];
                            second_y_left = arrVec[ez][ey+1][ex][iux_left];
                            third_y_left = arrVec[ez][ey+2][ex][iux_left];
                            fourth_y_left = arrVec[ez][ey+3][ex][iux_left];

                            PetscReal first_z_left, second_z_left, third_z_left, fourth_z_left;
                            first_z_left = arrVec[ez][ey][ex][iux_left];
                            second_z_left = arrVec[ez-1][ey][ex][iux_left];
                            third_z_left = arrVec[ez-2][ey][ex][iux_left];
                            fourth_z_left = arrVec[ez-3][ey][ex][iux_left];

                            PetscReal laplacian_left = (prev_x_left - 2*inter_x_left + next_x_left) / (hx*hx) + (2*first_y_left - 5*second_y_left + 4*third_y_left - fourth_y_left) / (hy*hy) + (2*first_z_left - 5*second_z_left + 4*third_z_left - fourth_z_left) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian_left; 
                        
                        } else {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_right];
                            second_x = arrVec[ez][ey][ex-1][iux_right];
                            third_x = arrVec[ez][ey][ex-2][iux_right];
                            fourth_x = arrVec[ez][ey][ex-3][iux_right];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_right];
                            second_y = arrVec[ez][ey+1][ex][iux_right];
                            third_y = arrVec[ez][ey+2][ex][iux_right];
                            fourth_y = arrVec[ez][ey+3][ex][iux_right];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iux_right];
                            prev_z = arrVec[ez-1][ey][ex][iux_right];
                            inter_z = arrVec[ez][ey][ex][iux_right];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_right] = laplacian;

                            PetscReal prev_x_left, next_x_left, inter_x_left;
                            prev_x_left = arrVec[ez][ey][ex-1][iux_left];
                            next_x_left = arrVec[ez][ey][ex+1][iux_left];
                            inter_x_left = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y_left, second_y_left, third_y_left, fourth_y_left;
                            first_y_left = arrVec[ez][ey][ex][iux_left];
                            second_y_left = arrVec[ez][ey+1][ex][iux_left];
                            third_y_left = arrVec[ez][ey+2][ex][iux_left];
                            fourth_y_left = arrVec[ez][ey+3][ex][iux_left];

                            PetscReal prev_z_left, next_z_left, inter_z_left;
                            prev_z_left = arrVec[ez-1][ey][ex][iux_left];
                            next_z_left = arrVec[ez+1][ey][ex][iux_left];
                            inter_z_left = arrVec[ez][ey][ex][iux_left];

                            PetscReal laplacian_left = (prev_x_left - 2*inter_x_left + next_x_left) / (hx*hx) + (2*first_y_left - 5*second_y_left + 4*third_y_left - fourth_y_left) / (hy*hy) + (prev_z_left - 2*inter_z_left + next_z_left) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian_left;
                                                 
                        }
                        





                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                                
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_right];
                            second_x = arrVec[ez][ey][ex-1][iux_right];
                            third_x = arrVec[ez][ey][ex-2][iux_right];
                            fourth_x = arrVec[ez][ey][ex-3][iux_right];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_right];
                            second_y = arrVec[ez][ey-1][ex][iux_right];
                            third_y = arrVec[ez][ey-2][ex][iux_right];
                            fourth_y = arrVec[ez][ey-3][ex][iux_right];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_right];    
                            second_z = arrVec[ez+1][ey][ex][iux_right];
                            third_z = arrVec[ez+2][ey][ex][iux_right];
                            fourth_z = arrVec[ez+3][ey][ex][iux_right];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_right] = laplacian; 

                            PetscReal prev_x_left, next_x_left, inter_x_left;
                            prev_x_left = arrVec[ez][ey][ex-1][iux_left];
                            next_x_left = arrVec[ez][ey][ex+1][iux_left];
                            inter_x_left = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y_left, second_y_left, third_y_left, fourth_y_left;
                            first_y_left = arrVec[ez][ey][ex][iux_left];
                            second_y_left = arrVec[ez][ey-1][ex][iux_left];
                            third_y_left = arrVec[ez][ey-2][ex][iux_left];
                            fourth_y_left = arrVec[ez][ey-3][ex][iux_left];

                            PetscReal first_z_left, second_z_left, third_z_left, fourth_z_left;
                            first_z_left = arrVec[ez][ey][ex][iux_left];
                            second_z_left = arrVec[ez+1][ey][ex][iux_left];
                            third_z_left = arrVec[ez+2][ey][ex][iux_left];
                            fourth_z_left = arrVec[ez+3][ey][ex][iux_left];

                            PetscReal laplacian_left = (prev_x_left - 2*inter_x_left + next_x_left) / (hx*hx) + (2*first_y_left - 5*second_y_left + 4*third_y_left - fourth_y_left) / (hy*hy) + (2*first_z_left - 5*second_z_left + 4*third_z_left - fourth_z_left) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian_left; 
                                
                                                
                        } else if (ez == N[2] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_right];
                            second_x = arrVec[ez][ey][ex-1][iux_right];
                            third_x = arrVec[ez][ey][ex-2][iux_right];
                            fourth_x = arrVec[ez][ey][ex-3][iux_right];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_right];
                            second_y = arrVec[ez][ey-1][ex][iux_right];
                            third_y = arrVec[ez][ey-2][ex][iux_right];
                            fourth_y = arrVec[ez][ey-3][ex][iux_right];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_right];    
                            second_z = arrVec[ez-1][ey][ex][iux_right];
                            third_z = arrVec[ez-2][ey][ex][iux_right];
                            fourth_z = arrVec[ez-3][ey][ex][iux_right];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_right] = laplacian; 

                            PetscReal prev_x_left, next_x_left, inter_x_left;
                            prev_x_left = arrVec[ez][ey][ex-1][iux_left];
                            next_x_left = arrVec[ez][ey][ex+1][iux_left];
                            inter_x_left = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y_left, second_y_left, third_y_left, fourth_y_left;
                            first_y_left = arrVec[ez][ey][ex][iux_left];
                            second_y_left = arrVec[ez][ey-1][ex][iux_left];
                            third_y_left = arrVec[ez][ey-2][ex][iux_left];
                            fourth_y_left = arrVec[ez][ey-3][ex][iux_left];

                            PetscReal first_z_left, second_z_left, third_z_left, fourth_z_left;
                            first_z_left = arrVec[ez][ey][ex][iux_left];
                            second_z_left = arrVec[ez-1][ey][ex][iux_left];
                            third_z_left = arrVec[ez-2][ey][ex][iux_left];
                            fourth_z_left = arrVec[ez-3][ey][ex][iux_left];

                            PetscReal laplacian_left = (prev_x_left - 2*inter_x_left + next_x_left) / (hx*hx) + (2*first_y_left - 5*second_y_left + 4*third_y_left - fourth_y_left) / (hy*hy) + (2*first_z_left - 5*second_z_left + 4*third_z_left - fourth_z_left) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian_left; 
                                
                        
                        } else {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_right];
                            second_x = arrVec[ez][ey][ex-1][iux_right];
                            third_x = arrVec[ez][ey][ex-2][iux_right];
                            fourth_x = arrVec[ez][ey][ex-3][iux_right];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_right];
                            second_y = arrVec[ez][ey-1][ex][iux_right];
                            third_y = arrVec[ez][ey-2][ex][iux_right];
                            fourth_y = arrVec[ez][ey-3][ex][iux_right];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iux_right];
                            prev_z = arrVec[ez-1][ey][ex][iux_right];
                            inter_z = arrVec[ez][ey][ex][iux_right];


                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);

                            arrOut[ez][ey][ex][iux_right] = laplacian;

                            PetscReal prev_x_left, next_x_left, inter_x_left;
                            prev_x_left = arrVec[ez][ey][ex-1][iux_left];
                            next_x_left = arrVec[ez][ey][ex+1][iux_left];
                            inter_x_left = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y_left, second_y_left, third_y_left, fourth_y_left;
                            first_y_left = arrVec[ez][ey][ex][iux_left];
                            second_y_left = arrVec[ez][ey-1][ex][iux_left];
                            third_y_left = arrVec[ez][ey-2][ex][iux_left];
                            fourth_y_left = arrVec[ez][ey-3][ex][iux_left];

                            PetscReal prev_z_left, next_z_left, inter_z_left;
                            prev_z_left = arrVec[ez-1][ey][ex][iux_left];
                            next_z_left = arrVec[ez+1][ey][ex][iux_left];
                            inter_z_left = arrVec[ez][ey][ex][iux_left];

                            PetscReal laplacian_left = (prev_x_left - 2*inter_x_left + next_x_left) / (hx*hx) + (2*first_y_left - 5*second_y_left + 4*third_y_left - fourth_y_left) / (hy*hy) + (prev_z_left - 2*inter_z_left + next_z_left) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian_left;                            
                        }
                    } else if (ez == 0) {

                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iux_right];
                        second_x = arrVec[ez][ey][ex-1][iux_right];
                        third_x = arrVec[ez][ey][ex-2][iux_right];
                        fourth_x = arrVec[ez][ey][ex-3][iux_right];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iux_right];
                        next_y = arrVec[ez][ey+1][ex][iux_right];
                        inter_y = arrVec[ez][ey][ex][iux_right];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iux_right];
                        second_z = arrVec[ez+1][ey][ex][iux_right];
                        third_z = arrVec[ez+2][ey][ex][iux_right];
                        fourth_z = arrVec[ez+3][ey][ex][iux_right];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (next_y - 2*inter_y + prev_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);

                        arrOut[ez][ey][ex][iux_right] = laplacian;

                        PetscReal prev_x_left, next_x_left, inter_x_left;
                        prev_x_left = arrVec[ez][ey][ex-1][iux_left];
                        next_x_left = arrVec[ez][ey][ex+1][iux_left];
                        inter_x_left = arrVec[ez][ey][ex][iux_left];

                        PetscReal prev_y_left, next_y_left, inter_y_left;
                        prev_y_left = arrVec[ez][ey-1][ex][iux_left];
                        next_y_left = arrVec[ez][ey+1][ex][iux_left];
                        inter_y_left = arrVec[ez][ey][ex][iux_left];

                        PetscReal first_z_left, second_z_left, third_z_left, fourth_z_left;
                        first_z_left = arrVec[ez][ey][ex][iux_left];
                        second_z_left = arrVec[ez+1][ey][ex][iux_left];
                        third_z_left = arrVec[ez+2][ey][ex][iux_left];
                        fourth_z_left = arrVec[ez+3][ey][ex][iux_left];

                        PetscReal laplacian_left = (prev_x_left - 2*inter_x_left + next_x_left) / (hx*hx) + (next_y_left - 2*inter_y_left + prev_y_left) / (hy*hy) + (2*first_z_left - 5*second_z_left + 4*third_z_left - fourth_z_left) / (hz*hz);
                        arrOut[ez][ey][ex][iux_left] = laplacian_left;  




                                             
                    } else if (ez == N[2] - 1) {

                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iux_right];
                        second_x = arrVec[ez][ey][ex-1][iux_right];
                        third_x = arrVec[ez][ey][ex-2][iux_right];
                        fourth_x = arrVec[ez][ey][ex-3][iux_right];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iux_right];
                        next_y = arrVec[ez][ey+1][ex][iux_right];
                        inter_y = arrVec[ez][ey][ex][iux_right];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iux_right];
                        second_z = arrVec[ez-1][ey][ex][iux_right];
                        third_z = arrVec[ez-2][ey][ex][iux_right];
                        fourth_z = arrVec[ez-3][ey][ex][iux_right];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (next_y - 2*inter_y + prev_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);

                        arrOut[ez][ey][ex][iux_right] = laplacian;


                        PetscReal prev_x_left, next_x_left, inter_x_left;
                        prev_x_left = arrVec[ez][ey][ex-1][iux_left];
                        next_x_left = arrVec[ez][ey][ex+1][iux_left];
                        inter_x_left = arrVec[ez][ey][ex][iux_left];

                        PetscReal prev_y_left, next_y_left, inter_y_left;
                        prev_y_left = arrVec[ez][ey-1][ex][iux_left];
                        next_y_left = arrVec[ez][ey+1][ex][iux_left];
                        inter_y_left = arrVec[ez][ey][ex][iux_left];

                        PetscReal first_z_left, second_z_left, third_z_left, fourth_z_left;
                        first_z_left = arrVec[ez][ey][ex][iux_left];
                        second_z_left = arrVec[ez-1][ey][ex][iux_left];
                        third_z_left = arrVec[ez-2][ey][ex][iux_left];
                        fourth_z_left = arrVec[ez-3][ey][ex][iux_left];

                        PetscReal laplacian_left = (prev_x_left - 2*inter_x_left + next_x_left) / (hx*hx) + (next_y_left - 2*inter_y_left + prev_y_left) / (hy*hy) + (2*first_z_left - 5*second_z_left + 4*third_z_left - fourth_z_left) / (hz*hz);
                        arrOut[ez][ey][ex][iux_left] = laplacian_left;  
                                              
                    } else {

                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iux_right];
                        second_x = arrVec[ez][ey][ex-1][iux_right];
                        third_x = arrVec[ez][ey][ex-2][iux_right];
                        fourth_x = arrVec[ez][ey][ex-3][iux_right];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iux_right];
                        next_y = arrVec[ez][ey+1][ex][iux_right];
                        inter_y = arrVec[ez][ey][ex][iux_right];

                        PetscReal next_z, prev_z, inter_z;
                        next_z = arrVec[ez+1][ey][ex][iux_right];
                        prev_z = arrVec[ez-1][ey][ex][iux_right];
                        inter_z = arrVec[ez][ey][ex][iux_right];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (next_y - 2*inter_y + prev_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);

                        arrOut[ez][ey][ex][iux_right] = laplacian;

                        PetscReal prev_x_left, next_x_left, inter_x_left;
                        prev_x_left = arrVec[ez][ey][ex-1][iux_left];
                        next_x_left = arrVec[ez][ey][ex+1][iux_left];
                        inter_x_left = arrVec[ez][ey][ex][iux_left];

                        PetscReal prev_y_left, next_y_left, inter_y_left;
                        prev_y_left = arrVec[ez][ey-1][ex][iux_left];
                        next_y_left = arrVec[ez][ey+1][ex][iux_left];
                        inter_y_left = arrVec[ez][ey][ex][iux_left];

                        PetscReal prev_z_left, next_z_left, inter_z_left;
                        prev_z_left = arrVec[ez-1][ey][ex][iux_left];
                        next_z_left = arrVec[ez+1][ey][ex][iux_left];
                        inter_z_left = arrVec[ez][ey][ex][iux_left];

                        PetscReal laplacian_left = (next_x_left - 2*inter_x_left + prev_x_left) / (hx*hx) + (next_y_left - 2*inter_y_left + prev_y_left) / (hy*hy) + (next_z_left - 2*inter_z_left + prev_z_left) / (hz*hz);


                        arrOut[ez][ey][ex][iux_left] = laplacian_left;

                    }




                } 
                else if (ex == 0) {

                    if (ey == 0) {
                        if (ez == 0) {
                            
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_left];
                            second_x = arrVec[ez][ey][ex+1][iux_left];
                            third_x = arrVec[ez][ey][ex+2][iux_left];
                            fourth_x = arrVec[ez][ey][ex+3][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey+1][ex][iux_left];
                            third_y = arrVec[ez][ey+2][ex][iux_left];
                            fourth_y = arrVec[ez][ey+3][ex][iux_left];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_left];    
                            second_z = arrVec[ez+1][ey][ex][iux_left];
                            third_z = arrVec[ez+2][ey][ex][iux_left];
                            fourth_z = arrVec[ez+3][ey][ex][iux_left];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 
                            
                                                 
                        } else if (ez == N[2] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_left];
                            second_x = arrVec[ez][ey][ex+1][iux_left];
                            third_x = arrVec[ez][ey][ex+2][iux_left];
                            fourth_x = arrVec[ez][ey][ex+3][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey+1][ex][iux_left];
                            third_y = arrVec[ez][ey+2][ex][iux_left];
                            fourth_y = arrVec[ez][ey+3][ex][iux_left];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_left];    
                            second_z = arrVec[ez-1][ey][ex][iux_left];
                            third_z = arrVec[ez-2][ey][ex][iux_left];
                            fourth_z = arrVec[ez-3][ey][ex][iux_left];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian;                           


                        } else {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_left];
                            second_x = arrVec[ez][ey][ex+1][iux_left];
                            third_x = arrVec[ez][ey][ex+2][iux_left];
                            fourth_x = arrVec[ez][ey][ex+3][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey+1][ex][iux_left];
                            third_y = arrVec[ez][ey+2][ex][iux_left];
                            fourth_y = arrVec[ez][ey+3][ex][iux_left];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iux_left];
                            prev_z = arrVec[ez-1][ey][ex][iux_left];
                            inter_z = arrVec[ez][ey][ex][iux_left];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 
                           
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                                
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_left];
                            second_x = arrVec[ez][ey][ex+1][iux_left];
                            third_x = arrVec[ez][ey][ex+2][iux_left];
                            fourth_x = arrVec[ez][ey][ex+3][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey-1][ex][iux_left];
                            third_y = arrVec[ez][ey-2][ex][iux_left];
                            fourth_y = arrVec[ez][ey-3][ex][iux_left];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_left];    
                            second_z = arrVec[ez+1][ey][ex][iux_left];
                            third_z = arrVec[ez+2][ey][ex][iux_left];
                            fourth_z = arrVec[ez+3][ey][ex][iux_left];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 
                                

                          
                        } else if (ez == N[2] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_left];
                            second_x = arrVec[ez][ey][ex+1][iux_left];
                            third_x = arrVec[ez][ey][ex+2][iux_left];
                            fourth_x = arrVec[ez][ey][ex+3][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey-1][ex][iux_left];
                            third_y = arrVec[ez][ey-2][ex][iux_left];
                            fourth_y = arrVec[ez][ey-3][ex][iux_left];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_left];    
                            second_z = arrVec[ez-1][ey][ex][iux_left];
                            third_z = arrVec[ez-2][ey][ex][iux_left];
                            fourth_z = arrVec[ez-3][ey][ex][iux_left];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 
                                                    

                        } else {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iux_left];
                            second_x = arrVec[ez][ey][ex+1][iux_left];
                            third_x = arrVec[ez][ey][ex+2][iux_left];
                            fourth_x = arrVec[ez][ey][ex+3][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey-1][ex][iux_left];
                            third_y = arrVec[ez][ey-2][ex][iux_left];
                            fourth_y = arrVec[ez][ey-3][ex][iux_left];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iux_left];
                            prev_z = arrVec[ez-1][ey][ex][iux_left];
                            inter_z = arrVec[ez][ey][ex][iux_left];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 
                           
                        }
                    } else if (ez == 0) {

                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iux_left];
                        second_x = arrVec[ez][ey][ex+1][iux_left];
                        third_x = arrVec[ez][ey][ex+2][iux_left];
                        fourth_x = arrVec[ez][ey][ex+3][iux_left];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iux_left];
                        next_y = arrVec[ez][ey+1][ex][iux_left];
                        inter_y = arrVec[ez][ey][ex][iux_left];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iux_left];
                        second_z = arrVec[ez+1][ey][ex][iux_left];
                        third_z = arrVec[ez+2][ey][ex][iux_left];
                        fourth_z = arrVec[ez+3][ey][ex][iux_left];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (next_y - 2*inter_y + prev_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);

                        arrOut[ez][ey][ex][iux_left] = laplacian;
                        
                      
                    } else if (ez == N[2] - 1) {

                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iux_left];
                        second_x = arrVec[ez][ey][ex+1][iux_left];
                        third_x = arrVec[ez][ey][ex+2][iux_left];
                        fourth_x = arrVec[ez][ey][ex+3][iux_left];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iux_left];
                        next_y = arrVec[ez][ey+1][ex][iux_left];
                        inter_y = arrVec[ez][ey][ex][iux_left];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iux_left];
                        second_z = arrVec[ez-1][ey][ex][iux_left];
                        third_z = arrVec[ez-2][ey][ex][iux_left];
                        fourth_z = arrVec[ez-3][ey][ex][iux_left];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (next_y - 2*inter_y + prev_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);

                        arrOut[ez][ey][ex][iux_left] = laplacian;
                     
                    } else {

                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iux_left];
                        second_x = arrVec[ez][ey][ex+1][iux_left];
                        third_x = arrVec[ez][ey][ex+2][iux_left];
                        fourth_x = arrVec[ez][ey][ex+3][iux_left];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iux_left];
                        next_y = arrVec[ez][ey+1][ex][iux_left];
                        inter_y = arrVec[ez][ey][ex][iux_left];

                        PetscReal next_z, prev_z, inter_z;
                        next_z = arrVec[ez+1][ey][ex][iux_left];
                        prev_z = arrVec[ez-1][ey][ex][iux_left];
                        inter_z = arrVec[ez][ey][ex][iux_left];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (next_y - 2*inter_y + prev_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);

                        arrOut[ez][ey][ex][iux_left] = laplacian;

                    }
                } else {
                   
                    if (ey == 0) {
                        if (ez == 0) {

                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iux_left];
                            next_x = arrVec[ez][ey][ex+1][iux_left];
                            inter_x = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey+1][ex][iux_left];
                            third_y = arrVec[ez][ey+2][ex][iux_left];
                            fourth_y = arrVec[ez][ey+3][ex][iux_left];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_left];    
                            second_z = arrVec[ez+1][ey][ex][iux_left];
                            third_z = arrVec[ez+2][ey][ex][iux_left];
                            fourth_z = arrVec[ez+3][ey][ex][iux_left];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 
                            
                                                      
                        } else if (ez == N[2] - 1) {

                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iux_left];
                            next_x = arrVec[ez][ey][ex+1][iux_left];
                            inter_x = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey+1][ex][iux_left];
                            third_y = arrVec[ez][ey+2][ex][iux_left];
                            fourth_y = arrVec[ez][ey+3][ex][iux_left];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_left];    
                            second_z = arrVec[ez-1][ey][ex][iux_left];
                            third_z = arrVec[ez-2][ey][ex][iux_left];
                            fourth_z = arrVec[ez-3][ey][ex][iux_left];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 
                            
                        } else {

                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iux_left];
                            next_x = arrVec[ez][ey][ex+1][iux_left];
                            inter_x = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey+1][ex][iux_left];
                            third_y = arrVec[ez][ey+2][ex][iux_left];
                            fourth_y = arrVec[ez][ey+3][ex][iux_left];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iux_left];
                            prev_z = arrVec[ez-1][ey][ex][iux_left];
                            inter_z = arrVec[ez][ey][ex][iux_left];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 
                                                 
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {

                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iux_left];
                            next_x = arrVec[ez][ey][ex+1][iux_left];
                            inter_x = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey-1][ex][iux_left];
                            third_y = arrVec[ez][ey-2][ex][iux_left];
                            fourth_y = arrVec[ez][ey-3][ex][iux_left];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_left];    
                            second_z = arrVec[ez+1][ey][ex][iux_left];
                            third_z = arrVec[ez+2][ey][ex][iux_left];
                            fourth_z = arrVec[ez+3][ey][ex][iux_left];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian;                            

                                                      
                        } else if (ez == N[2] - 1) {

                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iux_left];
                            next_x = arrVec[ez][ey][ex+1][iux_left];
                            inter_x = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey-1][ex][iux_left];
                            third_y = arrVec[ez][ey-2][ex][iux_left];
                            fourth_y = arrVec[ez][ey-3][ex][iux_left];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iux_left];    
                            second_z = arrVec[ez-1][ey][ex][iux_left];
                            third_z = arrVec[ez-2][ey][ex][iux_left];
                            fourth_z = arrVec[ez-3][ey][ex][iux_left];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 

                            
                        } else {

                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iux_left];
                            next_x = arrVec[ez][ey][ex+1][iux_left];
                            inter_x = arrVec[ez][ey][ex][iux_left];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iux_left];
                            second_y = arrVec[ez][ey-1][ex][iux_left];
                            third_y = arrVec[ez][ey-2][ex][iux_left];
                            fourth_y = arrVec[ez][ey-3][ex][iux_left];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iux_left];
                            prev_z = arrVec[ez-1][ey][ex][iux_left];
                            inter_z = arrVec[ez][ey][ex][iux_left];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iux_left] = laplacian; 
                                                      
                        }
                    } else if (ez == 0) {

                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iux_left];
                        next_x = arrVec[ez][ey][ex+1][iux_left];
                        inter_x = arrVec[ez][ey][ex][iux_left];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iux_left];
                        next_y = arrVec[ez][ey+1][ex][iux_left];
                        inter_y = arrVec[ez][ey][ex][iux_left];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iux_left];
                        second_z = arrVec[ez+1][ey][ex][iux_left];
                        third_z = arrVec[ez+2][ey][ex][iux_left];
                        fourth_z = arrVec[ez+3][ey][ex][iux_left];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (next_y - 2*inter_y + prev_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iux_left] = laplacian;
                        
                               
                    } else if (ez == N[2] - 1) {

                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iux_left];
                        next_x = arrVec[ez][ey][ex+1][iux_left];
                        inter_x = arrVec[ez][ey][ex][iux_left];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iux_left];
                        next_y = arrVec[ez][ey+1][ex][iux_left];
                        inter_y = arrVec[ez][ey][ex][iux_left];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iux_left];
                        second_z = arrVec[ez-1][ey][ex][iux_left];
                        third_z = arrVec[ez-2][ey][ex][iux_left];
                        fourth_z = arrVec[ez-3][ey][ex][iux_left];

                    

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (next_y - 2*inter_y + prev_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iux_left] = laplacian;
                        
                                           
                    } else {

                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iux_left];
                        next_x = arrVec[ez][ey][ex+1][iux_left];
                        inter_x = arrVec[ez][ey][ex][iux_left];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iux_left];
                        next_y = arrVec[ez][ey+1][ex][iux_left];
                        inter_y = arrVec[ez][ey][ex][iux_left];

                        PetscReal next_z, prev_z, inter_z;
                        next_z = arrVec[ez+1][ey][ex][iux_left];
                        prev_z = arrVec[ez-1][ey][ex][iux_left];
                        inter_z = arrVec[ez][ey][ex][iux_left];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (next_y - 2*inter_y + prev_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                        arrOut[ez][ey][ex][iux_left] = laplacian;
                        
                    }
                }                
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, laplacian);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);


    PetscFunctionReturn(0); 
}

PetscErrorCode laplacian_y_explicit(DM const & dmGrid, Vec & laplacian, Vec const & vec) 
{
    PetscInt iuy_down, iuy_up;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ey == N[1] - 1) {   
                    if (ex == 0) {
                        if (ez == 0) {      
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_up];
                            second_x = arrVec[ez][ey][ex+1][iuy_up];
                            third_x = arrVec[ez][ey][ex+2][iuy_up];
                            fourth_x = arrVec[ez][ey][ex+3][iuy_up];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_up];                        
                            second_y = arrVec[ez][ey-1][ex][iuy_up];
                            third_y = arrVec[ez][ey-2][ex][iuy_up];
                            fourth_y = arrVec[ez][ey-3][ex][iuy_up];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_up];    
                            second_z = arrVec[ez+1][ey][ex][iuy_up];    
                            third_z = arrVec[ez+2][ey][ex][iuy_up];
                            fourth_z = arrVec[ez+3][ey][ex][iuy_up];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_up] = laplacian; 

                            PetscReal first_x_down, second_x_down, third_x_down, fourth_x_down;
                            first_x_down = arrVec[ez][ey][ex][iuy_down];
                            second_x_down = arrVec[ez][ey][ex+1][iuy_down];
                            third_x_down = arrVec[ez][ey][ex+2][iuy_down];
                            fourth_x_down = arrVec[ez][ey][ex+3][iuy_down];

                            PetscReal prev_y_down, next_y_down, inter_y_down;
                            prev_y_down = arrVec[ez][ey-1][ex][iuy_down];
                            next_y_down = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y_down = arrVec[ez][ey][ex][iuy_down];

                            PetscReal first_z_down, second_z_down, third_z_down, fourth_z_down;
                            first_z_down = arrVec[ez][ey][ex][iuy_down];
                            second_z_down = arrVec[ez+1][ey][ex][iuy_down];
                            third_z_down = arrVec[ez+2][ey][ex][iuy_down];
                            fourth_z_down = arrVec[ez+3][ey][ex][iuy_down];

                            PetscReal laplacian_down = (2*first_x_down - 5*second_x_down + 4*third_x_down - fourth_x_down) / (hx*hx) + (next_y_down - 2*inter_y_down + prev_y_down) / (hy*hy) + (2*first_z_down - 5*second_z_down + 4*third_z_down - fourth_z_down) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian_down;          
                                                  
                        } else if (ez == N[2] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_up];
                            second_x = arrVec[ez][ey][ex+1][iuy_up];
                            third_x = arrVec[ez][ey][ex+2][iuy_up];
                            fourth_x = arrVec[ez][ey][ex+3][iuy_up];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_up];                        
                            second_y = arrVec[ez][ey-1][ex][iuy_up];
                            third_y = arrVec[ez][ey-2][ex][iuy_up];
                            fourth_y = arrVec[ez][ey-3][ex][iuy_up];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_up];    
                            second_z = arrVec[ez-1][ey][ex][iuy_up];    
                            third_z = arrVec[ez-2][ey][ex][iuy_up];
                            fourth_z = arrVec[ez-3][ey][ex][iuy_up];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_up] = laplacian;   

                            PetscReal first_x_down, second_x_down, third_x_down, fourth_x_down;
                            first_x_down = arrVec[ez][ey][ex][iuy_down];
                            second_x_down = arrVec[ez][ey][ex+1][iuy_down];
                            third_x_down = arrVec[ez][ey][ex+2][iuy_down];
                            fourth_x_down = arrVec[ez][ey][ex+3][iuy_down];

                            PetscReal prev_y_down, next_y_down, inter_y_down;
                            prev_y_down = arrVec[ez][ey-1][ex][iuy_down];
                            next_y_down = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y_down = arrVec[ez][ey][ex][iuy_down];

                            PetscReal first_z_down, second_z_down, third_z_down, fourth_z_down;
                            first_z_down = arrVec[ez][ey][ex][iuy_down];
                            second_z_down = arrVec[ez-1][ey][ex][iuy_down];
                            third_z_down = arrVec[ez-2][ey][ex][iuy_down];
                            fourth_z_down = arrVec[ez-3][ey][ex][iuy_down];

                            PetscReal laplacian_down = (2*first_x_down - 5*second_x_down + 4*third_x_down - fourth_x_down) / (hx*hx) + (next_y_down - 2*inter_y_down + prev_y_down) / (hy*hy) + (2*first_z_down - 5*second_z_down + 4*third_z_down - fourth_z_down) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian_down;  

                        } else {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_up];
                            second_x = arrVec[ez][ey][ex+1][iuy_up];
                            third_x = arrVec[ez][ey][ex+2][iuy_up];
                            fourth_x = arrVec[ez][ey][ex+3][iuy_up];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_up];                        
                            second_y = arrVec[ez][ey-1][ex][iuy_up];
                            third_y = arrVec[ez][ey-2][ex][iuy_up];
                            fourth_y = arrVec[ez][ey-3][ex][iuy_up];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iuy_up];
                            prev_z = arrVec[ez-1][ey][ex][iuy_up];
                            inter_z = arrVec[ez][ey][ex][iuy_up];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_up] = laplacian;   

                            PetscReal first_x_down, second_x_down, third_x_down, fourth_x_down;
                            first_x_down = arrVec[ez][ey][ex][iuy_down];
                            second_x_down = arrVec[ez][ey][ex+1][iuy_down];
                            third_x_down = arrVec[ez][ey][ex+2][iuy_down];
                            fourth_x_down = arrVec[ez][ey][ex+3][iuy_down];

                            PetscReal prev_y_down, next_y_down, inter_y_down;
                            prev_y_down = arrVec[ez][ey-1][ex][iuy_down];
                            next_y_down = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y_down = arrVec[ez][ey][ex][iuy_down];

                            PetscReal prev_z_down, next_z_down, inter_z_down;
                            prev_z_down = arrVec[ez-1][ey][ex][iuy_down];
                            next_z_down = arrVec[ez+1][ey][ex][iuy_down];
                            inter_z_down = arrVec[ez][ey][ex][iuy_down];

                            PetscReal laplacian_down = (2*first_x_down - 5*second_x_down + 4*third_x_down - fourth_x_down) / (hx*hx) + (next_y_down - 2*inter_y_down + prev_y_down) / (hy*hy) + (next_z_down - 2*inter_z_down + prev_z_down) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian_down;                      
                                                 
                        }
                    } else if (ex == N[0] - 1) {
                        if (ez == 0) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_up];
                            second_x = arrVec[ez][ey][ex-1][iuy_up];
                            third_x = arrVec[ez][ey][ex-2][iuy_up];
                            fourth_x = arrVec[ez][ey][ex-3][iuy_up];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_up];                        
                            second_y = arrVec[ez][ey-1][ex][iuy_up];
                            third_y = arrVec[ez][ey-2][ex][iuy_up];
                            fourth_y = arrVec[ez][ey-3][ex][iuy_up];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_up];    
                            second_z = arrVec[ez+1][ey][ex][iuy_up];    
                            third_z = arrVec[ez+2][ey][ex][iuy_up];
                            fourth_z = arrVec[ez+3][ey][ex][iuy_up];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_up] = laplacian;

                            PetscReal first_x_down, second_x_down, third_x_down, fourth_x_down;
                            first_x_down = arrVec[ez][ey][ex][iuy_down];
                            second_x_down = arrVec[ez][ey][ex-1][iuy_down];
                            third_x_down = arrVec[ez][ey][ex-2][iuy_down];
                            fourth_x_down = arrVec[ez][ey][ex-3][iuy_down];

                            PetscReal prev_y_down, next_y_down, inter_y_down;
                            prev_y_down = arrVec[ez][ey-1][ex][iuy_down];
                            next_y_down = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y_down = arrVec[ez][ey][ex][iuy_down];

                            PetscReal first_z_down, second_z_down, third_z_down, fourth_z_down;
                            first_z_down = arrVec[ez][ey][ex][iuy_down];
                            second_z_down = arrVec[ez+1][ey][ex][iuy_down];
                            third_z_down = arrVec[ez+2][ey][ex][iuy_down];
                            fourth_z_down = arrVec[ez+3][ey][ex][iuy_down];

                            PetscReal laplacian_down = (2*first_x_down - 5*second_x_down + 4*third_x_down - fourth_x_down) / (hx*hx) + (next_y_down - 2*inter_y_down + prev_y_down) / (hy*hy) + (2*first_z_down - 5*second_z_down + 4*third_z_down - fourth_z_down) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian_down; 

                        } else if (ez == N[2] - 1) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_up];
                            second_x = arrVec[ez][ey][ex-1][iuy_up];
                            third_x = arrVec[ez][ey][ex-2][iuy_up];
                            fourth_x = arrVec[ez][ey][ex-3][iuy_up];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_up];                        
                            second_y = arrVec[ez][ey-1][ex][iuy_up];
                            third_y = arrVec[ez][ey-2][ex][iuy_up];
                            fourth_y = arrVec[ez][ey-3][ex][iuy_up];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_up];    
                            second_z = arrVec[ez-1][ey][ex][iuy_up];    
                            third_z = arrVec[ez-2][ey][ex][iuy_up];
                            fourth_z = arrVec[ez-3][ey][ex][iuy_up];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_up] = laplacian;      

                            PetscReal first_x_down, second_x_down, third_x_down, fourth_x_down;
                            first_x_down = arrVec[ez][ey][ex][iuy_down];
                            second_x_down = arrVec[ez][ey][ex-1][iuy_down];
                            third_x_down = arrVec[ez][ey][ex-2][iuy_down];
                            fourth_x_down = arrVec[ez][ey][ex-3][iuy_down];

                            PetscReal prev_y_down, next_y_down, inter_y_down;
                            prev_y_down = arrVec[ez][ey-1][ex][iuy_down];
                            next_y_down = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y_down = arrVec[ez][ey][ex][iuy_down];

                            PetscReal first_z_down, second_z_down, third_z_down, fourth_z_down;
                            first_z_down = arrVec[ez][ey][ex][iuy_down];
                            second_z_down = arrVec[ez-1][ey][ex][iuy_down];
                            third_z_down = arrVec[ez-2][ey][ex][iuy_down];
                            fourth_z_down = arrVec[ez-3][ey][ex][iuy_down];

                            PetscReal laplacian_down = (2*first_x_down - 5*second_x_down + 4*third_x_down - fourth_x_down) / (hx*hx) + (next_y_down - 2*inter_y_down + prev_y_down) / (hy*hy) + (2*first_z_down - 5*second_z_down + 4*third_z_down - fourth_z_down) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian_down;

                        } else {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_up];
                            second_x = arrVec[ez][ey][ex-1][iuy_up];
                            third_x = arrVec[ez][ey][ex-2][iuy_up];
                            fourth_x = arrVec[ez][ey][ex-3][iuy_up];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_up];                        
                            second_y = arrVec[ez][ey-1][ex][iuy_up];
                            third_y = arrVec[ez][ey-2][ex][iuy_up];
                            fourth_y = arrVec[ez][ey-3][ex][iuy_up];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iuy_up];
                            prev_z = arrVec[ez-1][ey][ex][iuy_up];
                            inter_z = arrVec[ez][ey][ex][iuy_up];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_up] = laplacian;  

                            PetscReal first_x_down, second_x_down, third_x_down, fourth_x_down;
                            first_x_down = arrVec[ez][ey][ex][iuy_down];
                            second_x_down = arrVec[ez][ey][ex-1][iuy_down];
                            third_x_down = arrVec[ez][ey][ex-2][iuy_down];
                            fourth_x_down = arrVec[ez][ey][ex-3][iuy_down];

                            PetscReal prev_y_down, next_y_down, inter_y_down;
                            prev_y_down = arrVec[ez][ey-1][ex][iuy_down];
                            next_y_down = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y_down = arrVec[ez][ey][ex][iuy_down];

                            PetscReal prev_z_down, next_z_down, inter_z_down;
                            prev_z_down = arrVec[ez-1][ey][ex][iuy_down];
                            next_z_down = arrVec[ez+1][ey][ex][iuy_down];
                            inter_z_down = arrVec[ez][ey][ex][iuy_down];

                            PetscReal laplacian_down = (2*first_x_down - 5*second_x_down + 4*third_x_down - fourth_x_down) / (hx*hx) + (next_y_down - 2*inter_y_down + prev_y_down) / (hy*hy) + (next_z_down - 2*inter_z_down + prev_z_down) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian_down;

                        }
                    } else if (ez == 0) {

                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuy_up];
                        next_x = arrVec[ez][ey][ex+1][iuy_up];
                        inter_x = arrVec[ez][ey][ex][iuy_up];

                        PetscReal first_y, second_y, third_y, fourth_y;
                        first_y = arrVec[ez][ey][ex][iuy_up];
                        second_y = arrVec[ez][ey-1][ex][iuy_up];
                        third_y = arrVec[ez][ey-2][ex][iuy_up];
                        fourth_y = arrVec[ez][ey-3][ex][iuy_up];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuy_up];
                        second_z = arrVec[ez+1][ey][ex][iuy_up];
                        third_z = arrVec[ez+2][ey][ex][iuy_up];
                        fourth_z = arrVec[ez+3][ey][ex][iuy_up];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_up] = laplacian;

                        PetscReal prev_x_down, next_x_down, inter_x_down;
                        prev_x_down = arrVec[ez][ey][ex-1][iuy_down];
                        next_x_down = arrVec[ez][ey][ex+1][iuy_down];
                        inter_x_down = arrVec[ez][ey][ex][iuy_down];

                        PetscReal prev_y_down, next_y_down, inter_y_down;
                        prev_y_down = arrVec[ez][ey-1][ex][iuy_down];
                        next_y_down = arrVec[ez][ey+1][ex][iuy_down];
                        inter_y_down = arrVec[ez][ey][ex][iuy_down];

                        PetscReal first_z_down, second_z_down, third_z_down, fourth_z_down;
                        first_z_down = arrVec[ez][ey][ex][iuy_down];
                        second_z_down = arrVec[ez+1][ey][ex][iuy_down];
                        third_z_down = arrVec[ez+2][ey][ex][iuy_down];
                        fourth_z_down = arrVec[ez+3][ey][ex][iuy_down];

                        PetscReal laplacian_down = (prev_x_down - 2*inter_x_down + next_x_down) / (hx*hx) + (prev_y_down - 2*inter_y_down + next_y_down) / (hy*hy) + (2*first_z_down - 5*second_z_down + 4*third_z_down - fourth_z_down) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_down] = laplacian_down;
                                             
                    } else if (ez == N[2] - 1) {

                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuy_up];
                        next_x = arrVec[ez][ey][ex+1][iuy_up];
                        inter_x = arrVec[ez][ey][ex][iuy_up];

                        PetscReal first_y, second_y, third_y, fourth_y;
                        first_y = arrVec[ez][ey][ex][iuy_up];
                        second_y = arrVec[ez][ey-1][ex][iuy_up];
                        third_y = arrVec[ez][ey-2][ex][iuy_up];
                        fourth_y = arrVec[ez][ey-3][ex][iuy_up];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuy_up];
                        second_z = arrVec[ez-1][ey][ex][iuy_up];
                        third_z = arrVec[ez-2][ey][ex][iuy_up];
                        fourth_z = arrVec[ez-3][ey][ex][iuy_up];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_up] = laplacian;    

                        PetscReal prev_x_down, next_x_down, inter_x_down;
                        prev_x_down = arrVec[ez][ey][ex-1][iuy_down];
                        next_x_down = arrVec[ez][ey][ex+1][iuy_down];
                        inter_x_down = arrVec[ez][ey][ex][iuy_down];

                        PetscReal prev_y_down, next_y_down, inter_y_down;
                        prev_y_down = arrVec[ez][ey-1][ex][iuy_down];
                        next_y_down = arrVec[ez][ey+1][ex][iuy_down];
                        inter_y_down = arrVec[ez][ey][ex][iuy_down];

                        PetscReal first_z_down, second_z_down, third_z_down, fourth_z_down;
                        first_z_down = arrVec[ez][ey][ex][iuy_down];
                        second_z_down = arrVec[ez-1][ey][ex][iuy_down];
                        third_z_down = arrVec[ez-2][ey][ex][iuy_down];
                        fourth_z_down = arrVec[ez-3][ey][ex][iuy_down];

                        PetscReal laplacian_down = (prev_x_down - 2*inter_x_down + next_x_down) / (hx*hx) + (prev_y_down - 2*inter_y_down + next_y_down) / (hy*hy) + (2*first_z_down - 5*second_z_down + 4*third_z_down - fourth_z_down) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_down] = laplacian_down;
                                                                  
                    } else {

                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuy_up];
                        next_x = arrVec[ez][ey][ex+1][iuy_up];
                        inter_x = arrVec[ez][ey][ex][iuy_up];

                        PetscReal first_y, second_y, third_y, fourth_y;
                        first_y = arrVec[ez][ey][ex][iuy_up];
                        second_y = arrVec[ez][ey-1][ex][iuy_up];
                        third_y = arrVec[ez][ey-2][ex][iuy_up];
                        fourth_y = arrVec[ez][ey-3][ex][iuy_up];

                        PetscReal next_z, prev_z, inter_z;
                        next_z = arrVec[ez+1][ey][ex][iuy_up];
                        prev_z = arrVec[ez-1][ey][ex][iuy_up];
                        inter_z = arrVec[ez][ey][ex][iuy_up];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_up] = laplacian; 

                        PetscReal prev_x_down, next_x_down, inter_x_down;
                        prev_x_down = arrVec[ez][ey][ex-1][iuy_down];
                        next_x_down = arrVec[ez][ey][ex+1][iuy_down];
                        inter_x_down = arrVec[ez][ey][ex][iuy_down];

                        PetscReal prev_y_down, next_y_down, inter_y_down;
                        prev_y_down = arrVec[ez][ey-1][ex][iuy_down];
                        next_y_down = arrVec[ez][ey+1][ex][iuy_down];
                        inter_y_down = arrVec[ez][ey][ex][iuy_down];

                        PetscReal next_z_down, prev_z_down, inter_z_down;
                        next_z_down = arrVec[ez+1][ey][ex][iuy_down];
                        prev_z_down = arrVec[ez-1][ey][ex][iuy_down];
                        inter_z_down = arrVec[ez][ey][ex][iuy_down];

                        PetscReal laplacian_down = (prev_x_down - 2*inter_x_down + next_x_down) / (hx*hx) + (prev_y_down - 2*inter_y_down + next_y_down) / (hy*hy) + (next_z_down - 2*inter_z_down + prev_z_down) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_down] = laplacian_down; 

                    }
                } 
                else if (ey == 0) {

                    if (ex == 0) {
                        if (ez == 0) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex+1][iuy_down];
                            third_x = arrVec[ez][ey][ex+2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex+3][iuy_down];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_down];                        
                            second_y = arrVec[ez][ey+1][ex][iuy_down];
                            third_y = arrVec[ez][ey+2][ex][iuy_down];
                            fourth_y = arrVec[ez][ey+3][ex][iuy_down];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_down];    
                            second_z = arrVec[ez+1][ey][ex][iuy_down];    
                            third_z = arrVec[ez+2][ey][ex][iuy_down];
                            fourth_z = arrVec[ez+3][ey][ex][iuy_down];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian; 

                                                 
                        } else if (ez == N[2] - 1) {
                                
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex+1][iuy_down];
                            third_x = arrVec[ez][ey][ex+2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex+3][iuy_down];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_down];                        
                            second_y = arrVec[ez][ey+1][ex][iuy_down];
                            third_y = arrVec[ez][ey+2][ex][iuy_down];
                            fourth_y = arrVec[ez][ey+3][ex][iuy_down];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_down];    
                            second_z = arrVec[ez-1][ey][ex][iuy_down];    
                            third_z = arrVec[ez-2][ey][ex][iuy_down];
                            fourth_z = arrVec[ez-3][ey][ex][iuy_down];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian; 
    

                        } else {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex+1][iuy_down];
                            third_x = arrVec[ez][ey][ex+2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex+3][iuy_down];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_down];                        
                            second_y = arrVec[ez][ey+1][ex][iuy_down];
                            third_y = arrVec[ez][ey+2][ex][iuy_down];
                            fourth_y = arrVec[ez][ey+3][ex][iuy_down];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iuy_down];
                            prev_z = arrVec[ez-1][ey][ex][iuy_down];
                            inter_z = arrVec[ez][ey][ex][iuy_down];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian; 
                           
                        }
                    } else if (ex == N[0] - 1) {
                        if (ez == 0) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex-1][iuy_down];
                            third_x = arrVec[ez][ey][ex-2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex-3][iuy_down];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_down];                        
                            second_y = arrVec[ez][ey+1][ex][iuy_down];
                            third_y = arrVec[ez][ey+2][ex][iuy_down];
                            fourth_y = arrVec[ez][ey+3][ex][iuy_down];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_down];    
                            second_z = arrVec[ez+1][ey][ex][iuy_down];    
                            third_z = arrVec[ez+2][ey][ex][iuy_down];
                            fourth_z = arrVec[ez+3][ey][ex][iuy_down];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian;              
                          
                        } else if (ez == N[2] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex-1][iuy_down];
                            third_x = arrVec[ez][ey][ex-2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex-3][iuy_down];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_down];                        
                            second_y = arrVec[ez][ey+1][ex][iuy_down];
                            third_y = arrVec[ez][ey+2][ex][iuy_down];
                            fourth_y = arrVec[ez][ey+3][ex][iuy_down];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_down];    
                            second_z = arrVec[ez-1][ey][ex][iuy_down];    
                            third_z = arrVec[ez-2][ey][ex][iuy_down];
                            fourth_z = arrVec[ez-3][ey][ex][iuy_down];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian; 

                        } else {
                                
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex-1][iuy_down];
                            third_x = arrVec[ez][ey][ex-2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex-3][iuy_down];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuy_down];                        
                            second_y = arrVec[ez][ey+1][ex][iuy_down];
                            third_y = arrVec[ez][ey+2][ex][iuy_down];
                            fourth_y = arrVec[ez][ey+3][ex][iuy_down];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iuy_down];
                            prev_z = arrVec[ez-1][ey][ex][iuy_down];
                            inter_z = arrVec[ez][ey][ex][iuy_down];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian; 
                           
                        }
                    } else if (ez == 0) {
                            
                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuy_down];
                        next_x = arrVec[ez][ey][ex+1][iuy_down];
                        inter_x = arrVec[ez][ey][ex][iuy_down];

                        PetscReal first_y, second_y, third_y, fourth_y;
                        first_y = arrVec[ez][ey][ex][iuy_down];
                        second_y = arrVec[ez][ey+1][ex][iuy_down];
                        third_y = arrVec[ez][ey+2][ex][iuy_down];
                        fourth_y = arrVec[ez][ey+3][ex][iuy_down];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuy_down];
                        second_z = arrVec[ez+1][ey][ex][iuy_down];
                        third_z = arrVec[ez+2][ey][ex][iuy_down];
                        fourth_z = arrVec[ez+3][ey][ex][iuy_down];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_down] = laplacian;
                      
                    } else if (ez == N[2] - 1) {

                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuy_down];
                        next_x = arrVec[ez][ey][ex+1][iuy_down];
                        inter_x = arrVec[ez][ey][ex][iuy_down];

                        PetscReal first_y, second_y, third_y, fourth_y;
                        first_y = arrVec[ez][ey][ex][iuy_down];
                        second_y = arrVec[ez][ey+1][ex][iuy_down];
                        third_y = arrVec[ez][ey+2][ex][iuy_down];
                        fourth_y = arrVec[ez][ey+3][ex][iuy_down];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuy_down];
                        second_z = arrVec[ez-1][ey][ex][iuy_down];
                        third_z = arrVec[ez-2][ey][ex][iuy_down];
                        fourth_z = arrVec[ez-3][ey][ex][iuy_down];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_down] = laplacian; 

                     
                    } else {
                            
                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuy_down];
                        next_x = arrVec[ez][ey][ex+1][iuy_down];
                        inter_x = arrVec[ez][ey][ex][iuy_down];

                        PetscReal first_y, second_y, third_y, fourth_y;
                        first_y = arrVec[ez][ey][ex][iuy_down];
                        second_y = arrVec[ez][ey+1][ex][iuy_down];
                        third_y = arrVec[ez][ey+2][ex][iuy_down];
                        fourth_y = arrVec[ez][ey+3][ex][iuy_down];

                        PetscReal next_z, prev_z, inter_z;
                        next_z = arrVec[ez+1][ey][ex][iuy_down];
                        prev_z = arrVec[ez-1][ey][ex][iuy_down];
                        inter_z = arrVec[ez][ey][ex][iuy_down];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_down] = laplacian; 

                    }
                } else {
                   
                    if (ex == 0) {
                        if (ez == 0) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex+1][iuy_down];
                            third_x = arrVec[ez][ey][ex+2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex+3][iuy_down];

                            PetscReal prev_y, next_y, inter_y;
                            prev_y = arrVec[ez][ey-1][ex][iuy_down];
                            next_y = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y = arrVec[ez][ey][ex][iuy_down];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_down];  
                            second_z = arrVec[ez+1][ey][ex][iuy_down];
                            third_z = arrVec[ez+2][ey][ex][iuy_down];
                            fourth_z = arrVec[ez+3][ey][ex][iuy_down];                     

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian;    
                                                     
                        } else if (ez == N[2] - 1) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex+1][iuy_down];
                            third_x = arrVec[ez][ey][ex+2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex+3][iuy_down];

                            PetscReal prev_y, next_y, inter_y;
                            prev_y = arrVec[ez][ey-1][ex][iuy_down];
                            next_y = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y = arrVec[ez][ey][ex][iuy_down];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_down];  
                            second_z = arrVec[ez-1][ey][ex][iuy_down];
                            third_z = arrVec[ez-2][ey][ex][iuy_down];
                            fourth_z = arrVec[ez-3][ey][ex][iuy_down];                     

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian;    
                            
                        } else {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex+1][iuy_down];
                            third_x = arrVec[ez][ey][ex+2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex+3][iuy_down];

                            PetscReal prev_y, next_y, inter_y;
                            prev_y = arrVec[ez][ey-1][ex][iuy_down];
                            next_y = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y = arrVec[ez][ey][ex][iuy_down];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iuy_down];
                            prev_z = arrVec[ez-1][ey][ex][iuy_down];
                            inter_z = arrVec[ez][ey][ex][iuy_down];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian;    
                                                 
                        }
                    } else if (ex == N[0] - 1) {
                        if (ez == 0) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex-1][iuy_down];
                            third_x = arrVec[ez][ey][ex-2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex-3][iuy_down];

                            PetscReal prev_y, next_y, inter_y;
                            prev_y = arrVec[ez][ey-1][ex][iuy_down];
                            next_y = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y = arrVec[ez][ey][ex][iuy_down];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_down];  
                            second_z = arrVec[ez+1][ey][ex][iuy_down];
                            third_z = arrVec[ez+2][ey][ex][iuy_down];
                            fourth_z = arrVec[ez+3][ey][ex][iuy_down];                     

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian;                    
                                                      
                        } else if (ez == N[2] - 1) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex-1][iuy_down];
                            third_x = arrVec[ez][ey][ex-2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex-3][iuy_down];

                            PetscReal prev_y, next_y, inter_y;
                            prev_y = arrVec[ez][ey-1][ex][iuy_down];
                            next_y = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y = arrVec[ez][ey][ex][iuy_down];

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuy_down];  
                            second_z = arrVec[ez-1][ey][ex][iuy_down];
                            third_z = arrVec[ez-2][ey][ex][iuy_down];
                            fourth_z = arrVec[ez-3][ey][ex][iuy_down];                     

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian; 

                            
                        } else {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuy_down];
                            second_x = arrVec[ez][ey][ex-1][iuy_down];
                            third_x = arrVec[ez][ey][ex-2][iuy_down];
                            fourth_x = arrVec[ez][ey][ex-3][iuy_down];

                            PetscReal prev_y, next_y, inter_y;
                            prev_y = arrVec[ez][ey-1][ex][iuy_down];
                            next_y = arrVec[ez][ey+1][ex][iuy_down];
                            inter_y = arrVec[ez][ey][ex][iuy_down];

                            PetscReal next_z, prev_z, inter_z;
                            next_z = arrVec[ez+1][ey][ex][iuy_down];
                            prev_z = arrVec[ez-1][ey][ex][iuy_down];
                            inter_z = arrVec[ez][ey][ex][iuy_down];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuy_down] = laplacian; 
                                                      
                        }
                    } else if (ez == 0) {
                                
                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuy_down];
                        next_x = arrVec[ez][ey][ex+1][iuy_down];
                        inter_x = arrVec[ez][ey][ex][iuy_down];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuy_down];
                        next_y = arrVec[ez][ey+1][ex][iuy_down];
                        inter_y = arrVec[ez][ey][ex][iuy_down];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuy_down];  
                        second_z = arrVec[ez+1][ey][ex][iuy_down];
                        third_z = arrVec[ez+2][ey][ex][iuy_down];
                        fourth_z = arrVec[ez+3][ey][ex][iuy_down];                     

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_down] = laplacian;

                               
                    } else if (ez == N[2] - 1) {
                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuy_down];
                        next_x = arrVec[ez][ey][ex+1][iuy_down];
                        inter_x = arrVec[ez][ey][ex][iuy_down];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuy_down];
                        next_y = arrVec[ez][ey+1][ex][iuy_down];
                        inter_y = arrVec[ez][ey][ex][iuy_down];

                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuy_down];  
                        second_z = arrVec[ez-1][ey][ex][iuy_down];
                        third_z = arrVec[ez-2][ey][ex][iuy_down];
                        fourth_z = arrVec[ez-3][ey][ex][iuy_down];                     

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_down] = laplacian;                                           
                    } else {
                                
                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuy_down];
                        next_x = arrVec[ez][ey][ex+1][iuy_down];
                        inter_x = arrVec[ez][ey][ex][iuy_down];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuy_down];
                        next_y = arrVec[ez][ey+1][ex][iuy_down];
                        inter_y = arrVec[ez][ey][ex][iuy_down];

                        PetscReal next_z, prev_z, inter_z;
                        next_z = arrVec[ez+1][ey][ex][iuy_down];
                        prev_z = arrVec[ez-1][ey][ex][iuy_down];
                        inter_z = arrVec[ez][ey][ex][iuy_down];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (next_z - 2*inter_z + prev_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuy_down] = laplacian; 
                        
                    }
                }                
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, laplacian);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);


    PetscFunctionReturn(0);
}

PetscErrorCode laplacian_z_explicit(DM const & dmGrid, Vec & laplacian, Vec const & vec) 
{
    PetscInt iuz_back, iuz_front;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ez == N[2] - 1) {

                    if (ey == 0) {
                        if (ex == 0) {  
   
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_front];   
                            second_x = arrVec[ez][ey][ex+1][iuz_front];
                            third_x = arrVec[ez][ey][ex+2][iuz_front];
                            fourth_x = arrVec[ez][ey][ex+3][iuz_front];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_front];    
                            second_y = arrVec[ez][ey+1][ex][iuz_front];
                            third_y = arrVec[ez][ey+2][ex][iuz_front];
                            fourth_y = arrVec[ez][ey+3][ex][iuz_front];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_front];    
                            second_z = arrVec[ez-1][ey][ex][iuz_front];
                            third_z = arrVec[ez-2][ey][ex][iuz_front];
                            fourth_z = arrVec[ez-3][ey][ex][iuz_front];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_front] = laplacian;

                            PetscReal first_x_back, second_x_back, third_x_back, fourth_x_back;
                            first_x_back = arrVec[ez][ey][ex][iuz_back];
                            second_x_back = arrVec[ez][ey][ex+1][iuz_back];
                            third_x_back = arrVec[ez][ey][ex+2][iuz_back];
                            fourth_x_back = arrVec[ez][ey][ex+3][iuz_back];

                            PetscReal first_y_back, second_y_back, third_y_back, fourth_y_back;
                            first_y_back = arrVec[ez][ey][ex][iuz_back];
                            second_y_back = arrVec[ez][ey+1][ex][iuz_back];
                            third_y_back = arrVec[ez][ey+2][ex][iuz_back];
                            fourth_y_back = arrVec[ez][ey+3][ex][iuz_back];

                            PetscReal prev_z_back, next_z_back, inter_z_back;
                            prev_z_back = arrVec[ez-1][ey][ex][iuz_back];
                            next_z_back = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z_back = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian_back = (2*first_x_back - 5*second_x_back + 4*third_x_back - fourth_x_back) / (hx*hx) + (2*first_y_back - 5*second_y_back + 4*third_y_back - fourth_y_back) / (hy*hy) + (next_z_back - 2*inter_z_back + prev_z_back) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian_back;                
                                                  
                        } else if (ex == N[0] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_front];   
                            second_x = arrVec[ez][ey][ex-1][iuz_front];
                            third_x = arrVec[ez][ey][ex-2][iuz_front];
                            fourth_x = arrVec[ez][ey][ex-3][iuz_front];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_front];    
                            second_y = arrVec[ez][ey+1][ex][iuz_front];
                            third_y = arrVec[ez][ey+2][ex][iuz_front];
                            fourth_y = arrVec[ez][ey+3][ex][iuz_front];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_front];    
                            second_z = arrVec[ez-1][ey][ex][iuz_front];
                            third_z = arrVec[ez-2][ey][ex][iuz_front];
                            fourth_z = arrVec[ez-3][ey][ex][iuz_front];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_front] = laplacian;  
                            
                            PetscReal first_x_back, second_x_back, third_x_back, fourth_x_back;
                            first_x_back = arrVec[ez][ey][ex][iuz_back];
                            second_x_back = arrVec[ez][ey][ex-1][iuz_back];
                            third_x_back = arrVec[ez][ey][ex-2][iuz_back];
                            fourth_x_back = arrVec[ez][ey][ex-3][iuz_back];

                            PetscReal first_y_back, second_y_back, third_y_back, fourth_y_back;
                            first_y_back = arrVec[ez][ey][ex][iuz_back];
                            second_y_back = arrVec[ez][ey+1][ex][iuz_back];
                            third_y_back = arrVec[ez][ey+2][ex][iuz_back];
                            fourth_y_back = arrVec[ez][ey+3][ex][iuz_back];

                            PetscReal prev_z_back, next_z_back, inter_z_back;
                            prev_z_back = arrVec[ez-1][ey][ex][iuz_back];
                            next_z_back = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z_back = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian_back = (2*first_x_back - 5*second_x_back + 4*third_x_back - fourth_x_back) / (hx*hx) + (2*first_y_back - 5*second_y_back + 4*third_y_back - fourth_y_back) / (hy*hy) + (next_z_back - 2*inter_z_back + prev_z_back) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian_back;                                                    
                        } else {

                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iuz_front];
                            next_x = arrVec[ez][ey][ex+1][iuz_front];
                            inter_x = arrVec[ez][ey][ex][iuz_front];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_front];    
                            second_y = arrVec[ez][ey+1][ex][iuz_front];
                            third_y = arrVec[ez][ey+2][ex][iuz_front];
                            fourth_y = arrVec[ez][ey+3][ex][iuz_front];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_front];    
                            second_z = arrVec[ez-1][ey][ex][iuz_front];
                            third_z = arrVec[ez-2][ey][ex][iuz_front];
                            fourth_z = arrVec[ez-3][ey][ex][iuz_front];         
             
                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_front] = laplacian; 

                            PetscReal prev_x_back, next_x_back, inter_x_back;
                            prev_x_back = arrVec[ez][ey][ex-1][iuz_back];
                            next_x_back = arrVec[ez][ey][ex+1][iuz_back];
                            inter_x_back = arrVec[ez][ey][ex][iuz_back];

                            PetscReal first_y_back, second_y_back, third_y_back, fourth_y_back;
                            first_y_back = arrVec[ez][ey][ex][iuz_back];
                            second_y_back = arrVec[ez][ey+1][ex][iuz_back];
                            third_y_back = arrVec[ez][ey+2][ex][iuz_back];
                            fourth_y_back = arrVec[ez][ey+3][ex][iuz_back];

                            PetscReal prev_z_back, next_z_back, inter_z_back;
                            prev_z_back = arrVec[ez-1][ey][ex][iuz_back];
                            next_z_back = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z_back = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian_back = (prev_x_back - 2*inter_x_back + next_x_back) / (hx*hx) + (2*first_y_back - 5*second_y_back + 4*third_y_back - fourth_y_back) / (hy*hy) + (next_z_back - 2*inter_z_back + prev_z_back) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian_back; 
                                                 
                        }
                    } else if (ey == N[1] - 1) {
                        if (ex == 0) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_front];   
                            second_x = arrVec[ez][ey][ex+1][iuz_front];
                            third_x = arrVec[ez][ey][ex+2][iuz_front];
                            fourth_x = arrVec[ez][ey][ex+3][iuz_front];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_front];    
                            second_y = arrVec[ez][ey-1][ex][iuz_front];
                            third_y = arrVec[ez][ey-2][ex][iuz_front];
                            fourth_y = arrVec[ez][ey-3][ex][iuz_front];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_front];    
                            second_z = arrVec[ez-1][ey][ex][iuz_front];
                            third_z = arrVec[ez-2][ey][ex][iuz_front];
                            fourth_z = arrVec[ez-3][ey][ex][iuz_front];  

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);      
                            arrOut[ez][ey][ex][iuz_front] = laplacian;

                            PetscReal first_x_back, second_x_back, third_x_back, fourth_x_back;
                            first_x_back = arrVec[ez][ey][ex][iuz_back];
                            second_x_back = arrVec[ez][ey][ex+1][iuz_back];
                            third_x_back = arrVec[ez][ey][ex+2][iuz_back];
                            fourth_x_back = arrVec[ez][ey][ex+3][iuz_back];

                            PetscReal first_y_back, second_y_back, third_y_back, fourth_y_back;
                            first_y_back = arrVec[ez][ey][ex][iuz_back];
                            second_y_back = arrVec[ez][ey-1][ex][iuz_back];
                            third_y_back = arrVec[ez][ey-2][ex][iuz_back];
                            fourth_y_back = arrVec[ez][ey-3][ex][iuz_back];

                            PetscReal prev_z_back, next_z_back, inter_z_back;
                            prev_z_back = arrVec[ez-1][ey][ex][iuz_back];
                            next_z_back = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z_back = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian_back = (2*first_x_back - 5*second_x_back + 4*third_x_back - fourth_x_back) / (hx*hx) + (2*first_y_back - 5*second_y_back + 4*third_y_back - fourth_y_back) / (hy*hy) + (next_z_back - 2*inter_z_back + prev_z_back) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian_back; 
                        } else if (ex == N[0] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_front];   
                            second_x = arrVec[ez][ey][ex-1][iuz_front];
                            third_x = arrVec[ez][ey][ex-2][iuz_front];
                            fourth_x = arrVec[ez][ey][ex-3][iuz_front];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_front];    
                            second_y = arrVec[ez][ey-1][ex][iuz_front];
                            third_y = arrVec[ez][ey-2][ex][iuz_front];
                            fourth_y = arrVec[ez][ey-3][ex][iuz_front];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_front];    
                            second_z = arrVec[ez-1][ey][ex][iuz_front];
                            third_z = arrVec[ez-2][ey][ex][iuz_front];
                            fourth_z = arrVec[ez-3][ey][ex][iuz_front];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_front] = laplacian;  

                            PetscReal first_x_back, second_x_back, third_x_back, fourth_x_back;
                            first_x_back = arrVec[ez][ey][ex][iuz_back];
                            second_x_back = arrVec[ez][ey][ex-1][iuz_back];
                            third_x_back = arrVec[ez][ey][ex-2][iuz_back];
                            fourth_x_back = arrVec[ez][ey][ex-3][iuz_back];

                            PetscReal first_y_back, second_y_back, third_y_back, fourth_y_back;
                            first_y_back = arrVec[ez][ey][ex][iuz_back];
                            second_y_back = arrVec[ez][ey-1][ex][iuz_back];
                            third_y_back = arrVec[ez][ey-2][ex][iuz_back];
                            fourth_y_back = arrVec[ez][ey-3][ex][iuz_back];

                            PetscReal prev_z_back, next_z_back, inter_z_back;
                            prev_z_back = arrVec[ez-1][ey][ex][iuz_back];
                            next_z_back = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z_back = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian_back = (2*first_x_back - 5*second_x_back + 4*third_x_back - fourth_x_back) / (hx*hx) + (2*first_y_back - 5*second_y_back + 4*third_y_back - fourth_y_back) / (hy*hy) + (next_z_back - 2*inter_z_back + prev_z_back) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian_back;   

                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iuz_front];
                            next_x = arrVec[ez][ey][ex+1][iuz_front];
                            inter_x = arrVec[ez][ey][ex][iuz_front];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_front];    
                            second_y = arrVec[ez][ey-1][ex][iuz_front];
                            third_y = arrVec[ez][ey-2][ex][iuz_front];
                            fourth_y = arrVec[ez][ey-3][ex][iuz_front];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_front];    
                            second_z = arrVec[ez-1][ey][ex][iuz_front];
                            third_z = arrVec[ez-2][ey][ex][iuz_front];
                            fourth_z = arrVec[ez-3][ey][ex][iuz_front];         
             
                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_front] = laplacian;     

                            PetscReal prev_x_back, next_x_back, inter_x_back;
                            prev_x_back = arrVec[ez][ey][ex-1][iuz_back];
                            next_x_back = arrVec[ez][ey][ex+1][iuz_back];
                            inter_x_back = arrVec[ez][ey][ex][iuz_back];

                            PetscReal first_y_back, second_y_back, third_y_back, fourth_y_back;
                            first_y_back = arrVec[ez][ey][ex][iuz_back];
                            second_y_back = arrVec[ez][ey-1][ex][iuz_back];
                            third_y_back = arrVec[ez][ey-2][ex][iuz_back];
                            fourth_y_back = arrVec[ez][ey-3][ex][iuz_back];

                            PetscReal prev_z_back, next_z_back, inter_z_back;
                            prev_z_back = arrVec[ez-1][ey][ex][iuz_back];
                            next_z_back = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z_back = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian_back = (prev_x_back - 2*inter_x_back + next_x_back) / (hx*hx) + (2*first_y_back - 5*second_y_back + 4*third_y_back - fourth_y_back) / (hy*hy) + (next_z_back - 2*inter_z_back + prev_z_back) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian_back;                             
                                                    
                        }
                    } else if (ex == 0) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iuz_front];   
                        second_x = arrVec[ez][ey][ex+1][iuz_front];
                        third_x = arrVec[ez][ey][ex+2][iuz_front];
                        fourth_x = arrVec[ez][ey][ex+3][iuz_front];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuz_front];
                        next_y = arrVec[ez][ey+1][ex][iuz_front];
                        inter_y = arrVec[ez][ey][ex][iuz_front];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuz_front];    
                        second_z = arrVec[ez-1][ey][ex][iuz_front];
                        third_z = arrVec[ez-2][ey][ex][iuz_front];
                        fourth_z = arrVec[ez-3][ey][ex][iuz_front];         
            
                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_front] = laplacian;   

                        PetscReal first_x_back, second_x_back, third_x_back, fourth_x_back;
                        first_x_back = arrVec[ez][ey][ex][iuz_back];
                        second_x_back = arrVec[ez][ey][ex+1][iuz_back];
                        third_x_back = arrVec[ez][ey][ex+2][iuz_back];
                        fourth_x_back = arrVec[ez][ey][ex+3][iuz_back];

                        PetscReal prev_y_back, next_y_back, inter_y_back;
                        prev_y_back = arrVec[ez][ey-1][ex][iuz_back];
                        next_y_back = arrVec[ez][ey+1][ex][iuz_back];
                        inter_y_back = arrVec[ez][ey][ex][iuz_back];

                        PetscReal prev_z_back, next_z_back, inter_z_back;
                        prev_z_back = arrVec[ez-1][ey][ex][iuz_back];
                        next_z_back = arrVec[ez+1][ey][ex][iuz_back];
                        inter_z_back = arrVec[ez][ey][ex][iuz_back];

                        PetscReal laplacian_back = (2*first_x_back - 5*second_x_back + 4*third_x_back - fourth_x_back) / (hx*hx) + (prev_y_back - 2*inter_y_back + next_y_back) / (hy*hy) + (next_z_back - 2*inter_z_back + prev_z_back) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_back] = laplacian_back;


                                             
                    } else if (ex == N[0] - 1) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iuz_front];   
                        second_x = arrVec[ez][ey][ex-1][iuz_front];
                        third_x = arrVec[ez][ey][ex-2][iuz_front];
                        fourth_x = arrVec[ez][ey][ex-3][iuz_front];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuz_front];
                        next_y = arrVec[ez][ey+1][ex][iuz_front];
                        inter_y = arrVec[ez][ey][ex][iuz_front];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuz_front];    
                        second_z = arrVec[ez-1][ey][ex][iuz_front];
                        third_z = arrVec[ez-2][ey][ex][iuz_front];
                        fourth_z = arrVec[ez-3][ey][ex][iuz_front];         
            
                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_front] = laplacian;    

                        PetscReal first_x_back, second_x_back, third_x_back, fourth_x_back;
                        first_x_back = arrVec[ez][ey][ex][iuz_back];
                        second_x_back = arrVec[ez][ey][ex-1][iuz_back];
                        third_x_back = arrVec[ez][ey][ex-2][iuz_back];
                        fourth_x_back = arrVec[ez][ey][ex-3][iuz_back];

                        PetscReal prev_y_back, next_y_back, inter_y_back;
                        prev_y_back = arrVec[ez][ey-1][ex][iuz_back];
                        next_y_back = arrVec[ez][ey+1][ex][iuz_back];
                        inter_y_back = arrVec[ez][ey][ex][iuz_back];

                        PetscReal prev_z_back, next_z_back, inter_z_back;
                        prev_z_back = arrVec[ez-1][ey][ex][iuz_back];
                        next_z_back = arrVec[ez+1][ey][ex][iuz_back];
                        inter_z_back = arrVec[ez][ey][ex][iuz_back];

                        PetscReal laplacian_back = (2*first_x_back - 5*second_x_back + 4*third_x_back - fourth_x_back) / (hx*hx) + (prev_y_back - 2*inter_y_back + next_y_back) / (hy*hy) + (next_z_back - 2*inter_z_back + prev_z_back) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_back] = laplacian_back;
                                              
                    } else {
                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuz_front];
                        next_x = arrVec[ez][ey][ex+1][iuz_front];
                        inter_x = arrVec[ez][ey][ex][iuz_front];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuz_front];
                        next_y = arrVec[ez][ey+1][ex][iuz_front];
                        inter_y = arrVec[ez][ey][ex][iuz_front];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuz_front];    
                        second_z = arrVec[ez-1][ey][ex][iuz_front];
                        third_z = arrVec[ez-2][ey][ex][iuz_front];
                        fourth_z = arrVec[ez-3][ey][ex][iuz_front];         
            
                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_front] = laplacian; 

                        PetscReal prev_x_back, next_x_back, inter_x_back;
                        prev_x_back = arrVec[ez][ey][ex-1][iuz_back];
                        next_x_back = arrVec[ez][ey][ex+1][iuz_back];
                        inter_x_back = arrVec[ez][ey][ex][iuz_back];

                        PetscReal prev_y_back, next_y_back, inter_y_back;
                        prev_y_back = arrVec[ez][ey-1][ex][iuz_back];
                        next_y_back = arrVec[ez][ey+1][ex][iuz_back];
                        inter_y_back = arrVec[ez][ey][ex][iuz_back];

                        PetscReal prev_z_back, next_z_back, inter_z_back;
                        prev_z_back = arrVec[ez-1][ey][ex][iuz_back];
                        next_z_back = arrVec[ez+1][ey][ex][iuz_back];
                        inter_z_back = arrVec[ez][ey][ex][iuz_back];

                        PetscReal laplacian_back = (prev_x_back - 2*inter_x_back + next_x_back) / (hx*hx) + (prev_y_back - 2*inter_y_back + next_y_back) / (hy*hy) + (next_z_back - 2*inter_z_back + prev_z_back) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_back] = laplacian_back;
                    }
                } 
                else if (ez == 0) {

                    if (ey == 0) {
                        if (ex == 0) { 
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_back];   
                            second_x = arrVec[ez][ey][ex+1][iuz_back];
                            third_x = arrVec[ez][ey][ex+2][iuz_back];
                            fourth_x = arrVec[ez][ey][ex+3][iuz_back];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];    
                            second_y = arrVec[ez][ey+1][ex][iuz_back];
                            third_y = arrVec[ez][ey+2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey+3][ex][iuz_back];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_back];    
                            second_z = arrVec[ez+1][ey][ex][iuz_back];
                            third_z = arrVec[ez+2][ey][ex][iuz_back];
                            fourth_z = arrVec[ez+3][ey][ex][iuz_back];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                           
                       
                                                  
                        } else if (ex == N[0] - 1) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_back];   
                            second_x = arrVec[ez][ey][ex-1][iuz_back];
                            third_x = arrVec[ez][ey][ex-2][iuz_back];
                            fourth_x = arrVec[ez][ey][ex-3][iuz_back];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];    
                            second_y = arrVec[ez][ey+1][ex][iuz_back];
                            third_y = arrVec[ez][ey+2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey+3][ex][iuz_back];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_back];    
                            second_z = arrVec[ez+1][ey][ex][iuz_back];
                            third_z = arrVec[ez+2][ey][ex][iuz_back];
                            fourth_z = arrVec[ez+3][ey][ex][iuz_back];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                                                   
                 
                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iuz_back];
                            next_x = arrVec[ez][ey][ex+1][iuz_back];
                            inter_x = arrVec[ez][ey][ex][iuz_back];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];    
                            second_y = arrVec[ez][ey+1][ex][iuz_back];
                            third_y = arrVec[ez][ey+2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey+3][ex][iuz_back];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_back];    
                            second_z = arrVec[ez+1][ey][ex][iuz_back];
                            third_z = arrVec[ez+2][ey][ex][iuz_back];
                            fourth_z = arrVec[ez+3][ey][ex][iuz_back];         
             
                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                                                  
                        }
                    } else if (ey == N[1] - 1) {
                        if (ex == 0) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_back];   
                            second_x = arrVec[ez][ey][ex+1][iuz_back];
                            third_x = arrVec[ez][ey][ex+2][iuz_back];
                            fourth_x = arrVec[ez][ey][ex+3][iuz_back];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];    
                            second_y = arrVec[ez][ey-1][ex][iuz_back];
                            third_y = arrVec[ez][ey-2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey-3][ex][iuz_back];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_back];    
                            second_z = arrVec[ez+1][ey][ex][iuz_back];
                            third_z = arrVec[ez+2][ey][ex][iuz_back];
                            fourth_z = arrVec[ez+3][ey][ex][iuz_back];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                                                
                                                
                        } else if (ex == N[0] - 1) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_back];   
                            second_x = arrVec[ez][ey][ex-1][iuz_back];
                            third_x = arrVec[ez][ey][ex-2][iuz_back];
                            fourth_x = arrVec[ez][ey][ex-3][iuz_back];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];    
                            second_y = arrVec[ez][ey-1][ex][iuz_back];
                            third_y = arrVec[ez][ey-2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey-3][ex][iuz_back];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_back];    
                            second_z = arrVec[ez+1][ey][ex][iuz_back];
                            third_z = arrVec[ez+2][ey][ex][iuz_back];
                            fourth_z = arrVec[ez+3][ey][ex][iuz_back];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                                                 
                        

                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iuz_back];
                            next_x = arrVec[ez][ey][ex+1][iuz_back];
                            inter_x = arrVec[ez][ey][ex][iuz_back];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];    
                            second_y = arrVec[ez][ey-1][ex][iuz_back];
                            third_y = arrVec[ez][ey-2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey-3][ex][iuz_back];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iuz_back];    
                            second_z = arrVec[ez+1][ey][ex][iuz_back];
                            third_z = arrVec[ez+2][ey][ex][iuz_back];
                            fourth_z = arrVec[ez+3][ey][ex][iuz_back];         
             
                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                                                 
                            
                        }
                    } else if (ex == 0) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iuz_back];   
                        second_x = arrVec[ez][ey][ex+1][iuz_back];
                        third_x = arrVec[ez][ey][ex+2][iuz_back];
                        fourth_x = arrVec[ez][ey][ex+3][iuz_back];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuz_back];
                        next_y = arrVec[ez][ey+1][ex][iuz_back];
                        inter_y = arrVec[ez][ey][ex][iuz_back];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuz_back];    
                        second_z = arrVec[ez+1][ey][ex][iuz_back];
                        third_z = arrVec[ez+2][ey][ex][iuz_back];
                        fourth_z = arrVec[ez+3][ey][ex][iuz_back];         
            
                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_back] = laplacian; 
                                             
                    } else if (ex == N[0] - 1) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iuz_back];   
                        second_x = arrVec[ez][ey][ex-1][iuz_back];
                        third_x = arrVec[ez][ey][ex-2][iuz_back];
                        fourth_x = arrVec[ez][ey][ex-3][iuz_back];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuz_back];
                        next_y = arrVec[ez][ey+1][ex][iuz_back];
                        inter_y = arrVec[ez][ey][ex][iuz_back];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuz_back];    
                        second_z = arrVec[ez+1][ey][ex][iuz_back];
                        third_z = arrVec[ez+2][ey][ex][iuz_back];
                        fourth_z = arrVec[ez+3][ey][ex][iuz_back];         
            
                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_back] = laplacian; 

                                              
                    } else {

                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuz_back];
                        next_x = arrVec[ez][ey][ex+1][iuz_back];
                        inter_x = arrVec[ez][ey][ex][iuz_back];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuz_back];
                        next_y = arrVec[ez][ey+1][ex][iuz_back];
                        inter_y = arrVec[ez][ey][ex][iuz_back];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iuz_back];    
                        second_z = arrVec[ez+1][ey][ex][iuz_back];
                        third_z = arrVec[ez+2][ey][ex][iuz_back];
                        fourth_z = arrVec[ez+3][ey][ex][iuz_back];         
            
                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_back] = laplacian; 

                    }
                } else {
                   
                    if (ey == 0) {
                        if (ex == 0) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_back];
                            second_x = arrVec[ez][ey][ex+1][iuz_back];
                            third_x = arrVec[ez][ey][ex+2][iuz_back];
                            fourth_x = arrVec[ez][ey][ex+3][iuz_back];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];
                            second_y = arrVec[ez][ey+1][ex][iuz_back];
                            third_y = arrVec[ez][ey+2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey+3][ex][iuz_back];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iuz_back];
                            next_z = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                            
                                                      
                        } else if (ex == N[0] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_back];
                            second_x = arrVec[ez][ey][ex-1][iuz_back];
                            third_x = arrVec[ez][ey][ex-2][iuz_back];
                            fourth_x = arrVec[ez][ey][ex-3][iuz_back];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];
                            second_y = arrVec[ez][ey+1][ex][iuz_back];
                            third_y = arrVec[ez][ey+2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey+3][ex][iuz_back];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iuz_back];
                            next_z = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                        
   
                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iuz_back];
                            next_x = arrVec[ez][ey][ex+1][iuz_back];
                            inter_x = arrVec[ez][ey][ex][iuz_back];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];
                            second_y = arrVec[ez][ey+1][ex][iuz_back];
                            third_y = arrVec[ez][ey+2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey+3][ex][iuz_back];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iuz_back];
                            next_z = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                            
                                                 
                        }
                    } else if (ey == N[1] - 1) {
                        if (ex == 0) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_back];
                            second_x = arrVec[ez][ey][ex+1][iuz_back];
                            third_x = arrVec[ez][ey][ex+2][iuz_back];
                            fourth_x = arrVec[ez][ey][ex+3][iuz_back];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];
                            second_y = arrVec[ez][ey-1][ex][iuz_back];
                            third_y = arrVec[ez][ey-2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey-3][ex][iuz_back];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iuz_back];
                            next_z = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                          
                                                                                                      
                        } else if (ex == N[0] - 1) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iuz_back];
                            second_x = arrVec[ez][ey][ex-1][iuz_back];
                            third_x = arrVec[ez][ey][ex-2][iuz_back];
                            fourth_x = arrVec[ez][ey][ex-3][iuz_back];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];
                            second_y = arrVec[ez][ey-1][ex][iuz_back];
                            third_y = arrVec[ez][ey-2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey-3][ex][iuz_back];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iuz_back];
                            next_z = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                         
  
                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iuz_back];
                            next_x = arrVec[ez][ey][ex+1][iuz_back];
                            inter_x = arrVec[ez][ey][ex][iuz_back];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iuz_back];
                            second_y = arrVec[ez][ey-1][ex][iuz_back];
                            third_y = arrVec[ez][ey-2][ex][iuz_back];
                            fourth_y = arrVec[ez][ey-3][ex][iuz_back];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iuz_back];
                            next_z = arrVec[ez+1][ey][ex][iuz_back];
                            inter_z = arrVec[ez][ey][ex][iuz_back];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iuz_back] = laplacian;                         
                                                      
                        }
                    } else if (ex == 0) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iuz_back];
                        second_x = arrVec[ez][ey][ex+1][iuz_back];
                        third_x = arrVec[ez][ey][ex+2][iuz_back];
                        fourth_x = arrVec[ez][ey][ex+3][iuz_back];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuz_back];
                        next_y = arrVec[ez][ey+1][ex][iuz_back];
                        inter_y = arrVec[ez][ey][ex][iuz_back];

                        PetscReal prev_z, next_z, inter_z;
                        prev_z = arrVec[ez-1][ey][ex][iuz_back];
                        next_z = arrVec[ez+1][ey][ex][iuz_back];
                        inter_z = arrVec[ez][ey][ex][iuz_back];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_back] = laplacian;                      
                               
                    } else if (ex == N[0] - 1) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iuz_back];
                        second_x = arrVec[ez][ey][ex-1][iuz_back];
                        third_x = arrVec[ez][ey][ex-2][iuz_back];
                        fourth_x = arrVec[ez][ey][ex-3][iuz_back];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuz_back];
                        next_y = arrVec[ez][ey+1][ex][iuz_back];
                        inter_y = arrVec[ez][ey][ex][iuz_back];

                        PetscReal prev_z, next_z, inter_z;
                        prev_z = arrVec[ez-1][ey][ex][iuz_back];
                        next_z = arrVec[ez+1][ey][ex][iuz_back];
                        inter_z = arrVec[ez][ey][ex][iuz_back];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_back] = laplacian;                        
                    
                    } else {
                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iuz_back];
                        next_x = arrVec[ez][ey][ex+1][iuz_back];
                        inter_x = arrVec[ez][ey][ex][iuz_back];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iuz_back];
                        next_y = arrVec[ez][ey+1][ex][iuz_back];
                        inter_y = arrVec[ez][ey][ex][iuz_back];

                        PetscReal prev_z, next_z, inter_z;
                        prev_z = arrVec[ez-1][ey][ex][iuz_back];
                        next_z = arrVec[ez+1][ey][ex][iuz_back];
                        inter_z = arrVec[ez][ey][ex][iuz_back];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                        arrOut[ez][ey][ex][iuz_back] = laplacian; 
                        
                    }
                }                
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, laplacian);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);


    PetscFunctionReturn(0);
}

PetscErrorCode laplacian_element_explicit(DM const & dmGrid, Vec & laplacian, Vec const & vec) 
{
    PetscInt iue;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iue);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecOutLocal);
    DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ez == N[2] - 1) {

                    if (ey == 0) {
                        if (ex == 0) {  
   
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];   
                            second_x = arrVec[ez][ey][ex+1][iue];
                            third_x = arrVec[ez][ey][ex+2][iue];
                            fourth_x = arrVec[ez][ey][ex+3][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey+1][ex][iue];
                            third_y = arrVec[ez][ey+2][ex][iue];
                            fourth_y = arrVec[ez][ey+3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez-1][ey][ex][iue];
                            third_z = arrVec[ez-2][ey][ex][iue];
                            fourth_z = arrVec[ez-3][ey][ex][iue];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;           
                                                  
                        } else if (ex == N[0] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];   
                            second_x = arrVec[ez][ey][ex-1][iue];
                            third_x = arrVec[ez][ey][ex-2][iue];
                            fourth_x = arrVec[ez][ey][ex-3][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey+1][ex][iue];
                            third_y = arrVec[ez][ey+2][ex][iue];
                            fourth_y = arrVec[ez][ey+3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez-1][ey][ex][iue];
                            third_z = arrVec[ez-2][ey][ex][iue];
                            fourth_z = arrVec[ez-3][ey][ex][iue];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                              
                                                   
                        } else {

                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iue];
                            next_x = arrVec[ez][ey][ex+1][iue];
                            inter_x = arrVec[ez][ey][ex][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey+1][ex][iue];
                            third_y = arrVec[ez][ey+2][ex][iue];
                            fourth_y = arrVec[ez][ey+3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez-1][ey][ex][iue];
                            third_z = arrVec[ez-2][ey][ex][iue];
                            fourth_z = arrVec[ez-3][ey][ex][iue];         
             
                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian; 
                                                 
                        }
                    } else if (ey == N[1] - 1) {
                        if (ex == 0) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];   
                            second_x = arrVec[ez][ey][ex+1][iue];
                            third_x = arrVec[ez][ey][ex+2][iue];
                            fourth_x = arrVec[ez][ey][ex+3][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey-1][ex][iue];
                            third_y = arrVec[ez][ey-2][ex][iue];
                            fourth_y = arrVec[ez][ey-3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez-1][ey][ex][iue];
                            third_z = arrVec[ez-2][ey][ex][iue];
                            fourth_z = arrVec[ez-3][ey][ex][iue];  

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);      
                            arrOut[ez][ey][ex][iue] = laplacian;



                        } else if (ex == N[0] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];   
                            second_x = arrVec[ez][ey][ex-1][iue];
                            third_x = arrVec[ez][ey][ex-2][iue];
                            fourth_x = arrVec[ez][ey][ex-3][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey-1][ex][iue];
                            third_y = arrVec[ez][ey-2][ex][iue];
                            fourth_y = arrVec[ez][ey-3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez-1][ey][ex][iue];
                            third_z = arrVec[ez-2][ey][ex][iue];
                            fourth_z = arrVec[ez-3][ey][ex][iue];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;  
  

                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iue];
                            next_x = arrVec[ez][ey][ex+1][iue];
                            inter_x = arrVec[ez][ey][ex][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey-1][ex][iue];
                            third_y = arrVec[ez][ey-2][ex][iue];
                            fourth_y = arrVec[ez][ey-3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez-1][ey][ex][iue];
                            third_z = arrVec[ez-2][ey][ex][iue];
                            fourth_z = arrVec[ez-3][ey][ex][iue];         
             
                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;     
                        }
                    } else if (ex == 0) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iue];   
                        second_x = arrVec[ez][ey][ex+1][iue];
                        third_x = arrVec[ez][ey][ex+2][iue];
                        fourth_x = arrVec[ez][ey][ex+3][iue];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iue];
                        next_y = arrVec[ez][ey+1][ex][iue];
                        inter_y = arrVec[ez][ey][ex][iue];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iue];    
                        second_z = arrVec[ez-1][ey][ex][iue];
                        third_z = arrVec[ez-2][ey][ex][iue];
                        fourth_z = arrVec[ez-3][ey][ex][iue];         
            
                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iue] = laplacian;   
                                             
                    } else if (ex == N[0] - 1) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iue];   
                        second_x = arrVec[ez][ey][ex-1][iue];
                        third_x = arrVec[ez][ey][ex-2][iue];
                        fourth_x = arrVec[ez][ey][ex-3][iue];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iue];
                        next_y = arrVec[ez][ey+1][ex][iue];
                        inter_y = arrVec[ez][ey][ex][iue];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iue];    
                        second_z = arrVec[ez-1][ey][ex][iue];
                        third_z = arrVec[ez-2][ey][ex][iue];
                        fourth_z = arrVec[ez-3][ey][ex][iue];         
            
                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iue] = laplacian;    
                                              
                    } else {
                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iue];
                        next_x = arrVec[ez][ey][ex+1][iue];
                        inter_x = arrVec[ez][ey][ex][iue];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iue];
                        next_y = arrVec[ez][ey+1][ex][iue];
                        inter_y = arrVec[ez][ey][ex][iue];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iue];    
                        second_z = arrVec[ez-1][ey][ex][iue];
                        third_z = arrVec[ez-2][ey][ex][iue];
                        fourth_z = arrVec[ez-3][ey][ex][iue];         
            
                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iue] = laplacian; 
                    }
                } 
                else if (ez == 0) {

                    if (ey == 0) {
                        if (ex == 0) { 
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];   
                            second_x = arrVec[ez][ey][ex+1][iue];
                            third_x = arrVec[ez][ey][ex+2][iue];
                            fourth_x = arrVec[ez][ey][ex+3][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey+1][ex][iue];
                            third_y = arrVec[ez][ey+2][ex][iue];
                            fourth_y = arrVec[ez][ey+3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez+1][ey][ex][iue];
                            third_z = arrVec[ez+2][ey][ex][iue];
                            fourth_z = arrVec[ez+3][ey][ex][iue];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                           
                       
                                                  
                        } else if (ex == N[0] - 1) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];   
                            second_x = arrVec[ez][ey][ex-1][iue];
                            third_x = arrVec[ez][ey][ex-2][iue];
                            fourth_x = arrVec[ez][ey][ex-3][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey+1][ex][iue];
                            third_y = arrVec[ez][ey+2][ex][iue];
                            fourth_y = arrVec[ez][ey+3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez+1][ey][ex][iue];
                            third_z = arrVec[ez+2][ey][ex][iue];
                            fourth_z = arrVec[ez+3][ey][ex][iue];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                                                   
                 
                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iue];
                            next_x = arrVec[ez][ey][ex+1][iue];
                            inter_x = arrVec[ez][ey][ex][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey+1][ex][iue];
                            third_y = arrVec[ez][ey+2][ex][iue];
                            fourth_y = arrVec[ez][ey+3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez+1][ey][ex][iue];
                            third_z = arrVec[ez+2][ey][ex][iue];
                            fourth_z = arrVec[ez+3][ey][ex][iue];         
             
                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                                                  
                        }
                    } else if (ey == N[1] - 1) {
                        if (ex == 0) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];   
                            second_x = arrVec[ez][ey][ex+1][iue];
                            third_x = arrVec[ez][ey][ex+2][iue];
                            fourth_x = arrVec[ez][ey][ex+3][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey-1][ex][iue];
                            third_y = arrVec[ez][ey-2][ex][iue];
                            fourth_y = arrVec[ez][ey-3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez+1][ey][ex][iue];
                            third_z = arrVec[ez+2][ey][ex][iue];
                            fourth_z = arrVec[ez+3][ey][ex][iue];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                                                
                                                
                        } else if (ex == N[0] - 1) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];   
                            second_x = arrVec[ez][ey][ex-1][iue];
                            third_x = arrVec[ez][ey][ex-2][iue];
                            fourth_x = arrVec[ez][ey][ex-3][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey-1][ex][iue];
                            third_y = arrVec[ez][ey-2][ex][iue];
                            fourth_y = arrVec[ez][ey-3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez+1][ey][ex][iue];
                            third_z = arrVec[ez+2][ey][ex][iue];
                            fourth_z = arrVec[ez+3][ey][ex][iue];         
             
                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                                                 
                        

                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iue];
                            next_x = arrVec[ez][ey][ex+1][iue];
                            inter_x = arrVec[ez][ey][ex][iue];
                
                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];    
                            second_y = arrVec[ez][ey-1][ex][iue];
                            third_y = arrVec[ez][ey-2][ex][iue];
                            fourth_y = arrVec[ez][ey-3][ex][iue];   

                            PetscReal first_z, second_z, third_z, fourth_z;
                            first_z = arrVec[ez][ey][ex][iue];    
                            second_z = arrVec[ez+1][ey][ex][iue];
                            third_z = arrVec[ez+2][ey][ex][iue];
                            fourth_z = arrVec[ez+3][ey][ex][iue];         
             
                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                                                 
                            
                        }
                    } else if (ex == 0) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iue];   
                        second_x = arrVec[ez][ey][ex+1][iue];
                        third_x = arrVec[ez][ey][ex+2][iue];
                        fourth_x = arrVec[ez][ey][ex+3][iue];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iue];
                        next_y = arrVec[ez][ey+1][ex][iue];
                        inter_y = arrVec[ez][ey][ex][iue];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iue];    
                        second_z = arrVec[ez+1][ey][ex][iue];
                        third_z = arrVec[ez+2][ey][ex][iue];
                        fourth_z = arrVec[ez+3][ey][ex][iue];         
            
                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iue] = laplacian; 
                                             
                    } else if (ex == N[0] - 1) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iue];   
                        second_x = arrVec[ez][ey][ex-1][iue];
                        third_x = arrVec[ez][ey][ex-2][iue];
                        fourth_x = arrVec[ez][ey][ex-3][iue];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iue];
                        next_y = arrVec[ez][ey+1][ex][iue];
                        inter_y = arrVec[ez][ey][ex][iue];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iue];    
                        second_z = arrVec[ez+1][ey][ex][iue];
                        third_z = arrVec[ez+2][ey][ex][iue];
                        fourth_z = arrVec[ez+3][ey][ex][iue];         
            
                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iue] = laplacian; 

                                              
                    } else {

                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iue];
                        next_x = arrVec[ez][ey][ex+1][iue];
                        inter_x = arrVec[ez][ey][ex][iue];
            
                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iue];
                        next_y = arrVec[ez][ey+1][ex][iue];
                        inter_y = arrVec[ez][ey][ex][iue];
            
                        PetscReal first_z, second_z, third_z, fourth_z;
                        first_z = arrVec[ez][ey][ex][iue];    
                        second_z = arrVec[ez+1][ey][ex][iue];
                        third_z = arrVec[ez+2][ey][ex][iue];
                        fourth_z = arrVec[ez+3][ey][ex][iue];         
            
                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (2*first_z - 5*second_z + 4*third_z - fourth_z) / (hz*hz);
                        arrOut[ez][ey][ex][iue] = laplacian; 

                    }
                } else {
                   
                    if (ey == 0) {
                        if (ex == 0) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];
                            second_x = arrVec[ez][ey][ex+1][iue];
                            third_x = arrVec[ez][ey][ex+2][iue];
                            fourth_x = arrVec[ez][ey][ex+3][iue];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];
                            second_y = arrVec[ez][ey+1][ex][iue];
                            third_y = arrVec[ez][ey+2][ex][iue];
                            fourth_y = arrVec[ez][ey+3][ex][iue];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iue];
                            next_z = arrVec[ez+1][ey][ex][iue];
                            inter_z = arrVec[ez][ey][ex][iue];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                            
                                                      
                        } else if (ex == N[0] - 1) {

                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];
                            second_x = arrVec[ez][ey][ex-1][iue];
                            third_x = arrVec[ez][ey][ex-2][iue];
                            fourth_x = arrVec[ez][ey][ex-3][iue];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];
                            second_y = arrVec[ez][ey+1][ex][iue];
                            third_y = arrVec[ez][ey+2][ex][iue];
                            fourth_y = arrVec[ez][ey+3][ex][iue];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iue];
                            next_z = arrVec[ez+1][ey][ex][iue];
                            inter_z = arrVec[ez][ey][ex][iue];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                        
   
                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iue];
                            next_x = arrVec[ez][ey][ex+1][iue];
                            inter_x = arrVec[ez][ey][ex][iue];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];
                            second_y = arrVec[ez][ey+1][ex][iue];
                            third_y = arrVec[ez][ey+2][ex][iue];
                            fourth_y = arrVec[ez][ey+3][ex][iue];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iue];
                            next_z = arrVec[ez+1][ey][ex][iue];
                            inter_z = arrVec[ez][ey][ex][iue];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                            
                                                 
                        }
                    } else if (ey == N[1] - 1) {
                        if (ex == 0) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];
                            second_x = arrVec[ez][ey][ex+1][iue];
                            third_x = arrVec[ez][ey][ex+2][iue];
                            fourth_x = arrVec[ez][ey][ex+3][iue];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];
                            second_y = arrVec[ez][ey-1][ex][iue];
                            third_y = arrVec[ez][ey-2][ex][iue];
                            fourth_y = arrVec[ez][ey-3][ex][iue];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iue];
                            next_z = arrVec[ez+1][ey][ex][iue];
                            inter_z = arrVec[ez][ey][ex][iue];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                          
                                                                                                      
                        } else if (ex == N[0] - 1) {
                            PetscReal first_x, second_x, third_x, fourth_x;
                            first_x = arrVec[ez][ey][ex][iue];
                            second_x = arrVec[ez][ey][ex-1][iue];
                            third_x = arrVec[ez][ey][ex-2][iue];
                            fourth_x = arrVec[ez][ey][ex-3][iue];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];
                            second_y = arrVec[ez][ey-1][ex][iue];
                            third_y = arrVec[ez][ey-2][ex][iue];
                            fourth_y = arrVec[ez][ey-3][ex][iue];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iue];
                            next_z = arrVec[ez+1][ey][ex][iue];
                            inter_z = arrVec[ez][ey][ex][iue];

                            PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                         
  
                        } else {
                            PetscReal prev_x, next_x, inter_x;
                            prev_x = arrVec[ez][ey][ex-1][iue];
                            next_x = arrVec[ez][ey][ex+1][iue];
                            inter_x = arrVec[ez][ey][ex][iue];

                            PetscReal first_y, second_y, third_y, fourth_y;
                            first_y = arrVec[ez][ey][ex][iue];
                            second_y = arrVec[ez][ey-1][ex][iue];
                            third_y = arrVec[ez][ey-2][ex][iue];
                            fourth_y = arrVec[ez][ey-3][ex][iue];

                            PetscReal prev_z, next_z, inter_z;
                            prev_z = arrVec[ez-1][ey][ex][iue];
                            next_z = arrVec[ez+1][ey][ex][iue];
                            inter_z = arrVec[ez][ey][ex][iue];

                            PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (2*first_y - 5*second_y + 4*third_y - fourth_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                            arrOut[ez][ey][ex][iue] = laplacian;                         
                                                      
                        }
                    } else if (ex == 0) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iue];
                        second_x = arrVec[ez][ey][ex+1][iue];
                        third_x = arrVec[ez][ey][ex+2][iue];
                        fourth_x = arrVec[ez][ey][ex+3][iue];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iue];
                        next_y = arrVec[ez][ey+1][ex][iue];
                        inter_y = arrVec[ez][ey][ex][iue];

                        PetscReal prev_z, next_z, inter_z;
                        prev_z = arrVec[ez-1][ey][ex][iue];
                        next_z = arrVec[ez+1][ey][ex][iue];
                        inter_z = arrVec[ez][ey][ex][iue];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                        arrOut[ez][ey][ex][iue] = laplacian;                      
                               
                    } else if (ex == N[0] - 1) {
                        PetscReal first_x, second_x, third_x, fourth_x;
                        first_x = arrVec[ez][ey][ex][iue];
                        second_x = arrVec[ez][ey][ex-1][iue];
                        third_x = arrVec[ez][ey][ex-2][iue];
                        fourth_x = arrVec[ez][ey][ex-3][iue];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iue];
                        next_y = arrVec[ez][ey+1][ex][iue];
                        inter_y = arrVec[ez][ey][ex][iue];

                        PetscReal prev_z, next_z, inter_z;
                        prev_z = arrVec[ez-1][ey][ex][iue];
                        next_z = arrVec[ez+1][ey][ex][iue];
                        inter_z = arrVec[ez][ey][ex][iue];

                        PetscReal laplacian = (2*first_x - 5*second_x + 4*third_x - fourth_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                        arrOut[ez][ey][ex][iue] = laplacian;                        
                    
                    } else {
                        PetscReal prev_x, next_x, inter_x;
                        prev_x = arrVec[ez][ey][ex-1][iue];
                        next_x = arrVec[ez][ey][ex+1][iue];
                        inter_x = arrVec[ez][ey][ex][iue];

                        PetscReal prev_y, next_y, inter_y;
                        prev_y = arrVec[ez][ey-1][ex][iue];
                        next_y = arrVec[ez][ey+1][ex][iue];
                        inter_y = arrVec[ez][ey][ex][iue];

                        PetscReal prev_z, next_z, inter_z;
                        prev_z = arrVec[ez-1][ey][ex][iue];
                        next_z = arrVec[ez+1][ey][ex][iue];
                        inter_z = arrVec[ez][ey][ex][iue];

                        PetscReal laplacian = (prev_x - 2*inter_x + next_x) / (hx*hx) + (prev_y - 2*inter_y + next_y) / (hy*hy) + (prev_z - 2*inter_z + next_z) / (hz*hz);
                        arrOut[ez][ey][ex][iue] = laplacian; 
                        
                    }
                }                
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, laplacian);
    DMRestoreLocalVector(dmGrid, &vecOutLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);


    PetscFunctionReturn(0);
}



int main(int argc, char **argv)
{   

    auto start = std::chrono::high_resolution_clock::now();

    PetscInitialize(&argc, &argv, (char*)0, (char*)0);

   

    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Staggered; //need to declare to due to solving linear system laplacian messing things up
    {
        CreateGrid(&dmGrid_Shifted, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered_x, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMClone(dmGrid_Staggered_x, &dmGrid_Staggered_y);
        DMClone(dmGrid_Staggered_x, &dmGrid_Staggered_z);
    }

    Vec U_0, V_0, W_0;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_0);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_0);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_0);
    CreateAnalyticalU(dmGrid_Staggered_x, U_0, 0);
    CreateAnalyticalV(dmGrid_Staggered_y, V_0, 0);
    CreateAnalyticalW(dmGrid_Staggered_z, W_0, 0);

    Vec U_int, V_int, W_int;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_int);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_int);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_int);
    Vec U_pre, V_pre, W_pre;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_pre);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_pre);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_pre);
    Vec P, P_x, P_y, P_z;
    DMCreateGlobalVector(dmGrid_Centered, &P);
    DMCreateGlobalVector(dmGrid_Staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_Staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_Staggered_z, &P_z);
    Vec U_up, V_up, W_up;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_up);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_up);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_up);
    Vec Magnitude;
    DMCreateGlobalVector(dmGrid_Centered, &Magnitude);

    ManageAdvection_z(dt, W_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);
    /*Vec laplacian;
    DMCreateGlobalVector(dmGrid_Staggered_x, &laplacian);
    laplacian_x_explicit(dmGrid_Staggered_x, laplacian, U_0);
    Vec bench;
    DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
    CreateReferenceSolutionTry(dmGrid_Staggered_x, bench, 0);
    CheckSolution(laplacian, bench); */






    /*parabolic_problem_x parabolic_x(dmGrid_Staggered_x, U_0);

    parabolic_x.solve();

    parabolic_x.~parabolic_problem_x();

    parabolic_problem_y parabolic_y(dmGrid_Staggered_y, V_0);

    parabolic_y.solve();

    parabolic_y.~parabolic_problem_y();

    parabolic_problem_z parabolic_z(dmGrid_Staggered_z, W_0);

    parabolic_z.solve();

    parabolic_z.~parabolic_problem_z();*/

    stokes_problem stokes(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Staggered, dmGrid_Centered, dmGrid_Shifted, U_0, V_0, W_0);

    //stokes.solve();

    stokes.~stokes_problem();

    /*euler_problem euler(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Staggered, dmGrid_Centered, dmGrid_Shifted, U_0, V_0, W_0);

    euler.solve();

    euler.~euler_problem();*/


    /*for(size_t i = 0; i < iter; ++i){
        if (i == 0){

            Vec U_0_second, V_0_second, W_0_second;
            Vec U_0_third, V_0_third, W_0_third;
            DMCreateGlobalVector(dmGrid_Staggered_x, &U_0_second);
            DMCreateGlobalVector(dmGrid_Staggered_y, &V_0_second);
            DMCreateGlobalVector(dmGrid_Staggered_z, &W_0_second);
            DMCreateGlobalVector(dmGrid_Staggered_x, &U_0_third);
            DMCreateGlobalVector(dmGrid_Staggered_y, &V_0_third);
            DMCreateGlobalVector(dmGrid_Staggered_z, &W_0_third);
            VecCopy(U_0, U_0_second);
            VecCopy(V_0, V_0_second);
            VecCopy(W_0, W_0_second);
            VecCopy(U_0, U_0_third);
            VecCopy(V_0, V_0_third);
            VecCopy(W_0, W_0_third);
            
            //theta = d*d*(i)*dt;
            theta = (vRef/Re)*k*k*dt*(i);

            ManageAdvection_x(dt, U_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);
            ManageAdvection_y(dt, V_int, U_0_second, V_0_second, W_0_second, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);
            ManageAdvection_z(dt, W_int, U_0_third, V_0_third, W_0_third, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);

            //theta = d*d*(i+1)*dt;

            ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_pre, V_pre, W_pre, U_int, V_int, W_int, theta);

            //std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
            ManagePressure(dmGrid_Centered, dmGrid_Shifted, dmGrid_Staggered, dt, P, U_pre, V_pre, W_pre);
            ManagePressure_x(dmGrid_Staggered_x, dmGrid_Centered, dmGrid_Shifted, P_x, P);
            ManagePressure_y(dmGrid_Staggered_y, dmGrid_Centered, dmGrid_Shifted, P_y, P);
            ManagePressure_z(dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, P_z, P);
            
            //std::cout << "Spirit of Oklahoma completed:   pressure done." << std::endl;
            theta = (vRef/Re)*k*k*dt*(i+1);
            //theta = d*d*(i+1)*dt;

            UpdateVelocity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, U_up, V_up, W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, theta);
            std::cout << "Spirit of California comdmGrid_centpleted: iteration " << i << " done" << std::endl;

            Vec bench;
            DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
            CreateReferenceSolutionTry(dmGrid_Staggered_x, bench, theta);
            CheckSolution(U_up, bench);
            VecDestroy(&bench);

            ComputeMagnitude(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, Magnitude, U_up, V_up, W_up);

            PetscViewer viewer_magnitude;
            DM DM_magnitude;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec magnitude;
            DMStagVecSplitToDMDA(dmGrid_Centered, Magnitude, ELEMENT, 0, &DM_magnitude, &magnitude);
            PetscObjectSetName((PetscObject)magnitude, "magnitude");
            char filename_magnitude[50]; 
            sprintf(filename_magnitude, "results/magnitude%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_magnitude, FILE_MODE_WRITE, &viewer_magnitude);
            VecView(magnitude, viewer_magnitude);
            VecDestroy(&magnitude);
            DMDestroy(&DM_magnitude);
            PetscViewerDestroy(&viewer_magnitude);

            PetscViewer viewer_u;
            DM DM_u;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmGrid_Staggered_x, U_up, LEFT, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "x_component");
            char filename_u[50]; 
            sprintf(filename_u, "results/x_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u); 

            PetscViewer viewer_v;
            DM DM_v;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_y, 0, 0, 1, 0, &DM_v);
            Vec v;
            DMStagVecSplitToDMDA(dmGrid_Staggered_y, V_up, DOWN, 0, &DM_v, &v);
            PetscObjectSetName((PetscObject)v, "y_component");
            char filename_v[50];
            sprintf(filename_v, "results/y_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_y), filename_v, FILE_MODE_WRITE, &viewer_v);
            VecView(v, viewer_v);
            VecDestroy(&v);
            DMDestroy(&DM_v);
            PetscViewerDestroy(&viewer_v);

            PetscViewer viewer_w;
            DM DM_w;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_z, 0, 0, 1, 0, &DM_w);
            Vec w;
            DMStagVecSplitToDMDA(dmGrid_Staggered_z, W_up, BACK, 0, &DM_w, &w);
            PetscObjectSetName((PetscObject)w, "z_component");
            char filename_w[50];
            sprintf(filename_w, "results/z_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_z), filename_w, FILE_MODE_WRITE, &viewer_w);
            VecView(w, viewer_w);
            VecDestroy(&w);
            DMDestroy(&DM_w);
            PetscViewerDestroy(&viewer_w);

            PetscViewer viewer_p;
            DM DM_p;
            //DMStagCreateCompatibleDMStag(dmGrid_Centered, 0, 0, 0, 1, &DM_p);
            Vec p;
            DMStagVecSplitToDMDA(dmGrid_Centered, P, ELEMENT, 0, &DM_p, &p);
            PetscObjectSetName((PetscObject)p, "p");
            char filename_p[50];
            sprintf(filename_p, "results/p%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_p, FILE_MODE_WRITE, &viewer_p);
            VecView(p, viewer_p);
            VecDestroy(&p);
            DMDestroy(&DM_p);
            PetscViewerDestroy(&viewer_p);
            std::cout << "---------------------------------------------------------" << std::endl;

            VecDestroy(&U_0_second);
            VecDestroy(&V_0_second);
            VecDestroy(&W_0_second);
            VecDestroy(&U_0_third);
            VecDestroy(&V_0_third);
            VecDestroy(&W_0_third);

        } else {

            Vec U_up_second, V_up_second, W_up_second;
            Vec U_up_third, V_up_third, W_up_third;
            DMCreateGlobalVector(dmGrid_Staggered_x, &U_up_second);
            DMCreateGlobalVector(dmGrid_Staggered_y, &V_up_second);
            DMCreateGlobalVector(dmGrid_Staggered_z, &W_up_second);
            DMCreateGlobalVector(dmGrid_Staggered_x, &U_up_third);
            DMCreateGlobalVector(dmGrid_Staggered_y, &V_up_third);
            DMCreateGlobalVector(dmGrid_Staggered_z, &W_up_third);
            VecCopy(U_up, U_up_second);
            VecCopy(V_up, V_up_second);
            VecCopy(W_up, W_up_second);
            VecCopy(U_up, U_up_third);
            VecCopy(V_up, V_up_third);
            VecCopy(W_up, W_up_third);                    

            //theta = d*d*(i)*dt;
            theta = (vRef/Re)*k*k*dt*(i);

            ManageAdvection_x(dt, U_int, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);
            ManageAdvection_y(dt, V_int, U_up_second, V_up_second, W_up_second, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);
            ManageAdvection_z(dt, W_int, U_up_third, V_up_third, W_up_third, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);
            //std::cout << "Spirit of Kitty Hawk completed: advection done." << std::endl;

            //theta = d*d*(i+1)*dt;
            ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_pre, V_pre, W_pre, U_int, V_int, W_int, theta);

            //std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
            ManagePressure(dmGrid_Centered, dmGrid_Shifted, dmGrid_Staggered, dt, P, U_pre, V_pre, W_pre);
            ManagePressure_x(dmGrid_Staggered_x, dmGrid_Centered, dmGrid_Shifted, P_x, P);
            ManagePressure_y(dmGrid_Staggered_y, dmGrid_Centered, dmGrid_Shifted, P_y, P);
            ManagePressure_z(dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, P_z, P);
            //std::cout << "Spirit of Oklahoma completed:   pressure done." << std::endl;
            theta = (vRef/Re)*k*k*dt*(i+1);
            //theta = d*d*(i+1)*dt;
            UpdateVelocity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, U_up, V_up, W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, theta);
            //std::cout << "Spirit of California completed: iteration " << i << " done" << std::endl;
            Vec bench;
            DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
            CreateReferenceSolutionTry(dmGrid_Staggered_x, bench, theta);
            CheckSolution(U_up, bench);
            VecDestroy(&bench);


            ComputeMagnitude(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, Magnitude, U_up, V_up, W_up);

            PetscViewer viewer_magnitude;
            DM DM_magnitude;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec magnitude;
            DMStagVecSplitToDMDA(dmGrid_Centered, Magnitude, ELEMENT, 0, &DM_magnitude, &magnitude);
            PetscObjectSetName((PetscObject)magnitude, "magnitude");
            char filename_magnitude[50]; 
            sprintf(filename_magnitude, "results/magnitude%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_magnitude, FILE_MODE_WRITE, &viewer_magnitude);
            VecView(magnitude, viewer_magnitude);
            VecDestroy(&magnitude);
            DMDestroy(&DM_magnitude);
            PetscViewerDestroy(&viewer_magnitude);

            PetscViewer viewer_u;
            DM DM_u;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmGrid_Staggered_x, U_up, LEFT, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "x_component");
            char filename_u[50]; 
            sprintf(filename_u, "results/x_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u); 

            PetscViewer viewer_v;
            DM DM_v;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_y, 0, 0, 1, 0, &DM_v);
            Vec v;
            DMStagVecSplitToDMDA(dmGrid_Staggered_y, V_up, DOWN, 0, &DM_v, &v);
            PetscObjectSetName((PetscObject)v, "y_component");
            char filename_v[50];
            sprintf(filename_v, "results/y_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_y), filename_v, FILE_MODE_WRITE, &viewer_v);
            VecView(v, viewer_v);
            VecDestroy(&v);
            DMDestroy(&DM_v);
            PetscViewerDestroy(&viewer_v);

            PetscViewer viewer_w;
            DM DM_w;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_z, 0, 0, 1, 0, &DM_w);
            Vec w;
            DMStagVecSplitToDMDA(dmGrid_Staggered_z, W_up, BACK, 0, &DM_w, &w);
            PetscObjectSetName((PetscObject)w, "z_component");
            char filename_w[50];
            sprintf(filename_w, "results/z_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_z), filename_w, FILE_MODE_WRITE, &viewer_w);
            VecView(w, viewer_w);
            VecDestroy(&w);
            DMDestroy(&DM_w);
            PetscViewerDestroy(&viewer_w);

            PetscViewer viewer_p;
            DM DM_p;
            //DMStagCreateCompatibleDMStag(dmGrid_Centered, 0, 0, 0, 1, &DM_p);
            Vec p;
            DMStagVecSplitToDMDA(dmGrid_Centered, P, ELEMENT, 0, &DM_p, &p);
            PetscObjectSetName((PetscObject)p, "p");
            char filename_p[50];
            sprintf(filename_p, "results/p%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_p, FILE_MODE_WRITE, &viewer_p);
            VecView(p, viewer_p);
            VecDestroy(&p);
            DMDestroy(&DM_p);
            PetscViewerDestroy(&viewer_p);
            std::cout << "Iteration " << i << " completed." << std::endl;     
            std::cout << "---------------------------------------------------------" << std::endl;

            VecDestroy(&U_up_second);
            VecDestroy(&V_up_second);
            VecDestroy(&W_up_second);
            VecDestroy(&U_up_third);
            VecDestroy(&V_up_third);
            VecDestroy(&W_up_third);                
   
        }
    }*/


    VecDestroy(&U_0);
    VecDestroy(&V_0);
    VecDestroy(&W_0);
    VecDestroy(&U_int);
    VecDestroy(&V_int);
    VecDestroy(&W_int);
    VecDestroy(&U_pre);
    VecDestroy(&V_pre);
    VecDestroy(&W_pre);
    VecDestroy(&P);
    VecDestroy(&P_x);
    VecDestroy(&P_y);
    VecDestroy(&P_z);
    VecDestroy(&U_up);
    VecDestroy(&V_up);
    VecDestroy(&W_up);
    VecDestroy(&Magnitude); 
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_x);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_y);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_z);
    PetscObjectDestroy((PetscObject*)&dmGrid_Centered);
    PetscObjectDestroy((PetscObject*)&dmGrid_Shifted);  
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered);  

    PetscFinalize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Trinity test successfully completed. Ad maiora!" << std::endl;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    PetscFunctionReturn(0); 
}




