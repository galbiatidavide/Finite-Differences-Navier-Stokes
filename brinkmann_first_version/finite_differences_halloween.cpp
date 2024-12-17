//
// Created by dave on 26/10/23.
//

#include "parameters.hpp"
#include "parabolic.hpp"
#include "stokes.hpp"
#include "euler.hpp"
#include "petscvec.h" 

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


PetscErrorCode FirstShiftU_y(DM const & dmGrid, Vec & UShifted, Vec const & vec, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, up_right, down_left, up_left, left, down_right;
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

    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &right);
    DMStagGetLocationSlot(dmGrid, UP_RIGHT, 0, &up_right);
    DMStagGetLocationSlot(dmGrid, DOWN_LEFT, 0, &down_left);
    DMStagGetLocationSlot(dmGrid, UP_LEFT, 0, &up_left);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &left);
    DMStagGetLocationSlot(dmGrid, DOWN_RIGHT, 0, &down_right);

    Vec vec_local;
    PetscReal ****arrvec_local;
    DMCreateLocalVector(dmGrid,&vec_local);
    DMGlobalToLocalBegin(dmGrid,vec,INSERT_VALUES,vec_local);
    DMGlobalToLocalEnd(dmGrid,vec,INSERT_VALUES,vec_local);
    DMStagVecGetArrayRead(dmGrid, vec_local, &arrvec_local);

    Vec UShifted_local;
    VecSet(UShifted, 0.0);
    PetscReal ****arrUShifted_local;
    DMGetLocalVector(dmGrid, &UShifted_local);
    VecSet(UShifted_local, 0.0);
    DMStagVecGetArray(dmGrid, UShifted_local, &arrUShifted_local);    

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey != N[1] - 1 and ex != N[0] - 1) {

                    PetscReal current = arrvec_local[ez][ey][ex][right];
                    PetscReal next = arrvec_local[ez][ey + 1][ex][right];
                    PetscReal inter = (next + current) / 2.0;
                    arrUShifted_local[ez][ey][ex][up_right] = inter;
                }

                if (ey == 0) {
                    PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_down_left[0]], arrCoord[ez][ey][ex][icux_down_left[1]],
                                  arrCoord[ez][ey][ex][icux_down_left[2]], theta);
                    arrUShifted_local[ez][ey][ex][down_left] = inter;                 
                }

                if (ex == 0) {                  
                    PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_up_left[0]], arrCoord[ez][ey][ex][icux_up_left[1]],
                                  arrCoord[ez][ey][ex][icux_up_left[2]], theta);
                    arrUShifted_local[ez][ey][ex][up_left] = inter;                  
                }

                if (ey == N[1] - 1) {
                    PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_up_right[0]], arrCoord[ez][ey][ex][icux_up_right[1]],
                                  arrCoord[ez][ey][ex][icux_up_right[2]], theta);
                    arrUShifted_local[ez][ey][ex][up_right] = inter;
                }

                if (ex == N[0] - 1) {
                    PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_down_right[0]], arrCoord[ez][ey][ex][icux_down_right[1]],
                                  arrCoord[ez][ey][ex][icux_down_right[2]], theta);
                    arrUShifted_local[ez][ey][ex][down_right] = inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagVecRestoreArrayRead(dmGrid, vec_local, &arrvec_local);
    DMStagVecRestoreArray(dmGrid, UShifted_local, &arrUShifted_local);
    DMLocalToGlobal(dmGrid, UShifted_local, ADD_VALUES, UShifted);
    DMRestoreLocalVector(dmGrid, &vec_local);
    DMRestoreLocalVector(dmGrid, &UShifted_local);

    return 0;
}

PetscErrorCode FirstShiftU_z(DM const & dmGrid, Vec & UShifted, Vec const & vec, PetscScalar const & theta) //ok 
{
    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, front_right, back_left, front_left, back_right;
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

    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &right);
    DMStagGetLocationSlot(dmGrid, FRONT_RIGHT, 0, &front_right);
    DMStagGetLocationSlot(dmGrid, BACK_LEFT, 0, &back_left);
    DMStagGetLocationSlot(dmGrid, FRONT_LEFT, 0, &front_left);
    DMStagGetLocationSlot(dmGrid, BACK_RIGHT, 0, &back_right);

    Vec vec_local;
    PetscReal ****arrvec_local;
    DMCreateLocalVector(dmGrid,&vec_local);
    DMGlobalToLocalBegin(dmGrid,vec,INSERT_VALUES,vec_local);
    DMGlobalToLocalEnd(dmGrid,vec,INSERT_VALUES,vec_local);
    DMStagVecGetArrayRead(dmGrid, vec_local, &arrvec_local);

    Vec UShifted_local;
    VecSet(UShifted, 0.0);
    PetscReal ****arrUShifted_local;
    DMGetLocalVector(dmGrid, &UShifted_local);
    VecSet(UShifted_local, 0.0);
    DMStagVecGetArray(dmGrid, UShifted_local, &arrUShifted_local);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez != N[2] - 1 and ex != N[0] - 1) {
                    PetscReal current = arrvec_local[ez][ey][ex][right];
                    PetscReal next = arrvec_local[ez + 1][ey][ex][right];
                    PetscReal inter = (next + current) / 2.0;
                    arrUShifted_local[ez][ey][ex][front_right] = inter;
                }

                if (ez == 0) {
                    PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_back_left[0]], arrCoord[ez][ey][ex][icux_back_left[1]],
                                  arrCoord[ez][ey][ex][icux_back_left[2]], theta);
                    arrUShifted_local[ez][ey][ex][back_left] = inter;
                }

                if (ex == 0) {
                    PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_front_left[0]], arrCoord[ez][ey][ex][icux_front_left[1]],
                                  arrCoord[ez][ey][ex][icux_front_left[2]], theta);
                    arrUShifted_local[ez][ey][ex][front_left] = inter;
                }

                if (ez == N[2] - 1) {
                    PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_front_right[0]], arrCoord[ez][ey][ex][icux_front_right[1]],
                                  arrCoord[ez][ey][ex][icux_front_right[2]], theta);
                    arrUShifted_local[ez][ey][ex][front_right] = inter;
                }

                if (ex == N[0] - 1) {
                    PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_back_right[0]], arrCoord[ez][ey][ex][icux_back_right[1]],
                                  arrCoord[ez][ey][ex][icux_back_right[2]], theta);
                    arrUShifted_local[ez][ey][ex][back_right] = inter;
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagVecRestoreArrayRead(dmGrid, vec_local, &arrvec_local);
    DMStagVecRestoreArray(dmGrid, UShifted_local, &arrUShifted_local);
    DMLocalToGlobal(dmGrid, UShifted_local, ADD_VALUES, UShifted);
    DMRestoreLocalVector(dmGrid, &vec_local);
    DMRestoreLocalVector(dmGrid, &UShifted_local);

    return 0;
}

PetscErrorCode FirstShiftV_y(DM const & dmGrid, Vec & VShifted, Vec const & vec, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, up, up_right, down_left, up_left, down_right;
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

    DMStagGetLocationSlot(dmGrid, UP, 0, &up);
    DMStagGetLocationSlot(dmGrid, UP_RIGHT, 0, &up_right);
    DMStagGetLocationSlot(dmGrid, DOWN_LEFT, 0, &down_left);
    DMStagGetLocationSlot(dmGrid, UP_LEFT, 0, &up_left);
    DMStagGetLocationSlot(dmGrid, DOWN_RIGHT, 0, &down_right);

    Vec vec_local;
    PetscReal ****arrvec_local;
    DMCreateLocalVector(dmGrid,&vec_local);
    DMGlobalToLocalBegin(dmGrid,vec,INSERT_VALUES,vec_local);
    DMGlobalToLocalEnd(dmGrid,vec,INSERT_VALUES,vec_local);
    DMStagVecGetArrayRead(dmGrid, vec_local, &arrvec_local);

    Vec VShifted_local;
    VecSet(VShifted, 0.0);
    PetscReal ****arrVShifted_local;
    DMGetLocalVector(dmGrid, &VShifted_local);
    VecSet(VShifted_local, 0.0);
    DMStagVecGetArray(dmGrid, VShifted_local, &arrVShifted_local);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ex != N[0] - 1 and ey != N[1] - 1) {
                    PetscReal current = arrvec_local[ez][ey][ex][up];
                    PetscReal next = arrvec_local[ez][ey][ex + 1][up];
                    PetscReal inter = (next + current) / 2.0;
                    arrVShifted_local[ez][ey][ex][up_right] = inter;
                }

                if (ey == 0) {
                    PetscReal inter = uyRef(arrCoord[ez][ey][ex][icuy_down_left[0]], arrCoord[ez][ey][ex][icuy_down_left[1]],
                                  arrCoord[ez][ey][ex][icuy_down_left[2]], theta);
                    arrVShifted_local[ez][ey][ex][down_left] = inter;
                }

                if (ex == 0) {
                    PetscReal inter = uyRef(arrCoord[ez][ey][ex][icuy_up_left[0]], arrCoord[ez][ey][ex][icuy_up_left[1]],
                                  arrCoord[ez][ey][ex][icuy_up_left[2]], theta);
                    arrVShifted_local[ez][ey][ex][up_left] = inter;
                }

                if (ey == N[1] - 1) {
                    PetscReal inter = uyRef(arrCoord[ez][ey][ex][icuy_up_right[0]], arrCoord[ez][ey][ex][icuy_up_right[1]],
                                  arrCoord[ez][ey][ex][icuy_up_right[2]], theta);
                    arrVShifted_local[ez][ey][ex][up_right] = inter;
                }

                if (ex == N[0] - 1) {
                    PetscReal inter = uyRef(arrCoord[ez][ey][ex][icuy_down_right[0]], arrCoord[ez][ey][ex][icuy_down_right[1]],
                                  arrCoord[ez][ey][ex][icuy_down_right[2]], theta);
                    arrVShifted_local[ez][ey][ex][down_right] = inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagVecRestoreArrayRead(dmGrid, vec_local, &arrvec_local);
    DMStagVecRestoreArray(dmGrid, VShifted_local, &arrVShifted_local);
    DMLocalToGlobal(dmGrid, VShifted_local, ADD_VALUES, VShifted);
    DMRestoreLocalVector(dmGrid, &vec_local);
    DMRestoreLocalVector(dmGrid, &VShifted_local);

    return 0;
}

PetscErrorCode FirstShiftW_z(DM const & dmGrid, Vec & WShifted, Vec const & vec, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, front, front_right, back_left, front_left, back_right;
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

    DMStagGetLocationSlot(dmGrid, FRONT, 0, &front);
    DMStagGetLocationSlot(dmGrid, FRONT_RIGHT, 0, &front_right);
    DMStagGetLocationSlot(dmGrid, BACK_LEFT, 0, &back_left);
    DMStagGetLocationSlot(dmGrid, FRONT_LEFT, 0, &front_left);
    DMStagGetLocationSlot(dmGrid, BACK_RIGHT, 0, &back_right);


    Vec vec_local;
    PetscReal ****arrvec_local;
    DMCreateLocalVector(dmGrid,&vec_local);
    DMGlobalToLocalBegin(dmGrid,vec,INSERT_VALUES,vec_local);
    DMGlobalToLocalEnd(dmGrid,vec,INSERT_VALUES,vec_local);
    DMStagVecGetArrayRead(dmGrid, vec_local, &arrvec_local);

    Vec WShifted_local;
    VecSet(WShifted, 0.0);
    PetscReal ****arrWShifted_local;
    DMGetLocalVector(dmGrid, &WShifted_local);
    VecSet(WShifted_local, 0.0);
    DMStagVecGetArray(dmGrid, WShifted_local, &arrWShifted_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ez != N[2] - 1 and ex != N[0] - 1) {
                    PetscReal current = arrvec_local[ez][ey][ex][front];
                    PetscReal next = arrvec_local[ez][ey][ex + 1][front];
                    PetscReal inter = (next + current) / 2.0;
                    arrWShifted_local[ez][ey][ex][front_right] = inter;
                }

                if (ez == 0) {
                    PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_back_left[0]], arrCoord[ez][ey][ex][icux_back_left[1]],
                                  arrCoord[ez][ey][ex][icux_back_left[2]], theta);
                    arrWShifted_local[ez][ey][ex][back_left] = inter;   
                }

                if (ex == 0) {
                    PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_front_left[0]], arrCoord[ez][ey][ex][icux_front_left[1]],
                                  arrCoord[ez][ey][ex][icux_front_left[2]], theta);
                    arrWShifted_local[ez][ey][ex][front_left] = inter;
                }

                if (ez == N[2] - 1) {
                    PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_front_right[0]], arrCoord[ez][ey][ex][icux_front_right[1]],
                                  arrCoord[ez][ey][ex][icux_front_right[2]], theta);
                    arrWShifted_local[ez][ey][ex][front_right] = inter;
                }

                if (ex == N[0] - 1) {
                    PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_back_right[0]], arrCoord[ez][ey][ex][icux_back_right[1]],
                                  arrCoord[ez][ey][ex][icux_back_right[2]], theta);
                    arrWShifted_local[ez][ey][ex][back_right] = inter;
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagVecRestoreArrayRead(dmGrid, vec_local, &arrvec_local);
    DMStagVecRestoreArray(dmGrid, WShifted_local, &arrWShifted_local);
    DMLocalToGlobal(dmGrid, WShifted_local, ADD_VALUES, WShifted);
    DMRestoreLocalVector(dmGrid, &vec_local);
    DMRestoreLocalVector(dmGrid, &WShifted_local);

    return 0;
}


PetscErrorCode FirstDerive_y(DM const & dmGrid, Vec & AB_y, Vec const & AB) //ok
{
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, down_left, up_left, left, down_right, up_right, right;
    PetscReal       hy;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pAB_y);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    hy = 1.0 / N[1];
    
    Vec AB_local;
    PetscScalar ****arrAB_local;
    DMCreateLocalVector(dmGrid, &AB_local);
    DMGlobalToLocalBegin(dmGrid, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmGrid, AB, INSERT_VALUES, AB_local);
    DMStagVecGetArrayRead(dmGrid, AB_local, &arrAB_local);

    Vec AB_y_local;
    VecSet(AB_y, 0.0);
    PetscScalar ****arrAB_y_local;
    DMGetLocalVector(dmGrid, &AB_y_local);
    VecSet(AB_y_local, 0.0);
    DMStagVecGetArray(dmGrid, AB_y_local, &arrAB_y_local);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &right);
    DMStagGetLocationSlot(dmGrid, DOWN_LEFT, 0, &down_left);
    DMStagGetLocationSlot(dmGrid, DOWN_RIGHT, 0, &down_right);
    DMStagGetLocationSlot(dmGrid, UP_LEFT, 0, &up_left);
    DMStagGetLocationSlot(dmGrid, UP_RIGHT, 0, &up_right);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal val_down = arrAB_local[ez][ey][ex][down_left];
                PetscReal val_up = arrAB_local[ez][ey][ex][up_left];
                PetscReal der = (val_up - val_down) / hy;
                arrAB_y_local[ez][ey][ex][left] = der;

                if (ex == N[0] - 1) {

                    val_down = arrAB_local[ez][ey][ex][down_right];
                    val_up = arrAB_local[ez][ey][ex][up_right];
                    der = (val_up - val_down) / hy;
                    arrAB_y_local[ez][ey][ex][right] = der;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, AB_local, &arrAB_local);
    DMStagVecRestoreArray(dmGrid, AB_y_local, &arrAB_y_local);
    DMLocalToGlobal(dmGrid, AB_y_local, ADD_VALUES, AB_y);
    DMRestoreLocalVector(dmGrid, &AB_y_local);
    DMRestoreLocalVector(dmGrid, &AB_local);

    return 0;
}

PetscErrorCode FirstDerive_z(DM const & dmGrid, Vec & AB_z, Vec const & AB) //ok
{
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, back_left, front_left, left, back_right, front_right, right;
    PetscReal       hz;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pAB_z);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    hz = 1.0 / N[2];

    Vec AB_local;
    PetscScalar ****arrAB_local;
    DMCreateLocalVector(dmGrid, &AB_local);
    DMGlobalToLocalBegin(dmGrid, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmGrid, AB, INSERT_VALUES, AB_local);
    DMStagVecGetArrayRead(dmGrid, AB_local, &arrAB_local);

    Vec AB_z_local;
    VecSet(AB_z, 0.0);
    PetscScalar ****arrAB_z_local;
    DMGetLocalVector(dmGrid, &AB_z_local);
    VecSet(AB_z_local, 0.0);
    DMStagVecGetArray(dmGrid, AB_z_local, &arrAB_z_local);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &right);
    DMStagGetLocationSlot(dmGrid, BACK_LEFT, 0, &back_left);
    DMStagGetLocationSlot(dmGrid, BACK_RIGHT, 0, &back_right);
    DMStagGetLocationSlot(dmGrid, FRONT_LEFT, 0, &front_left);
    DMStagGetLocationSlot(dmGrid, FRONT_RIGHT, 0, &front_right);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal val_back = arrAB_local[ez][ey][ex][back_left];
                PetscReal val_front = arrAB_local[ez][ey][ex][front_left];
                PetscReal der = (val_front - val_back) / hz;
                arrAB_z_local[ez][ey][ex][left] = der;

                if (ex == N[0] - 1) {

                    val_back = arrAB_local[ez][ey][ex][back_right];
                    val_front = arrAB_local[ez][ey][ex][front_right];
                    der = (val_front - val_back) / hz;
                    arrAB_z_local[ez][ey][ex][right] = der;
                }
                
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, AB_local, &arrAB_local);
    DMStagVecRestoreArray(dmGrid, AB_z_local, &arrAB_z_local);
    DMLocalToGlobal(dmGrid, AB_z_local, ADD_VALUES, AB_z);
    DMRestoreLocalVector(dmGrid, &AB_z_local);
    DMRestoreLocalVector(dmGrid, &AB_local);

    return 0;
}


// Second non-linear members: mixed
PetscErrorCode SecondShiftV_z(DM const & dmGrid, Vec & VShifted, Vec const & vec, PetscScalar const & theta) //ok
{
    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, up, front_up, back_down, front_down, back_up;
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

    DMStagGetLocationSlot(dmGrid, UP, 0, &up);
    DMStagGetLocationSlot(dmGrid, FRONT_UP, 0, &front_up);
    DMStagGetLocationSlot(dmGrid, BACK_DOWN, 0, &back_down);
    DMStagGetLocationSlot(dmGrid, FRONT_DOWN, 0, &front_down);  
    DMStagGetLocationSlot(dmGrid, BACK_UP, 0, &back_up);

    Vec vec_local;
    PetscReal ****arrvec_local;
    DMCreateLocalVector(dmGrid,&vec_local);
    DMGlobalToLocalBegin(dmGrid,vec,INSERT_VALUES,vec_local);
    DMGlobalToLocalEnd(dmGrid,vec,INSERT_VALUES,vec_local);
    DMStagVecGetArrayRead(dmGrid, vec_local, &arrvec_local);

    Vec VShifted_local;
    VecSet(VShifted, 0.0);
    PetscReal ****arrVShifted_local;
    DMGetLocalVector(dmGrid, &VShifted_local);
    VecSet(VShifted_local, 0.0);
    DMStagVecGetArray(dmGrid, VShifted_local, &arrVShifted_local);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ey != N[1] - 1 and ez != N[2] - 1) {
                    PetscReal current = arrvec_local[ez][ey][ex][up];
                    PetscReal next = arrvec_local[ez + 1][ey][ex][up];
                    PetscReal inter = (next + current) / 2.0;
                    arrVShifted_local[ez][ey][ex][front_up] = inter;
                }

                if (ey == 0) {
                    PetscScalar inter = uyRef(arrCoord[ez][ey][ex][icuy_front_down[0]], arrCoord[ez][ey][ex][icuy_front_down[1]],
                                  arrCoord[ez][ey][ex][icuy_front_down[2]], theta);
                    arrVShifted_local[ez][ey][ex][front_down] = inter;
                }

                if (ez == N[2] - 1) {
                    PetscScalar inter = uyRef(arrCoord[ez][ey][ex][icuy_front_up[0]], arrCoord[ez][ey][ex][icuy_front_up[1]],
                                  arrCoord[ez][ey][ex][icuy_front_up[2]], theta);
                    arrVShifted_local[ez][ey][ex][front_up] = inter;
                }

                if (ey == N[1] - 1) {
                    PetscScalar inter = uyRef(arrCoord[ez][ey][ex][icuy_back_up[0]], arrCoord[ez][ey][ex][icuy_back_up[1]],
                                  arrCoord[ez][ey][ex][icuy_back_up[2]], theta);
                    arrVShifted_local[ez][ey][ex][back_up] = inter;
                }

                if (ez == 0) {
                    PetscScalar inter = uyRef(arrCoord[ez][ey][ex][icuy_back_down[0]], arrCoord[ez][ey][ex][icuy_back_down[1]],
                                  arrCoord[ez][ey][ex][icuy_back_down[2]], theta);
                    arrVShifted_local[ez][ey][ex][back_down] = inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagVecRestoreArrayRead(dmGrid, vec_local, &arrvec_local);
    DMStagVecRestoreArray(dmGrid, VShifted_local, &arrVShifted_local);
    DMLocalToGlobal(dmGrid, VShifted_local, ADD_VALUES, VShifted);
    DMRestoreLocalVector(dmGrid, &vec_local);
    DMRestoreLocalVector(dmGrid, &VShifted_local);
    
    return 0;
}

PetscErrorCode SecondShiftW_z(DM const & dmGrid, Vec & WShifted, Vec const & vec, PetscScalar const & theta) //ok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, front, front_up, back_down, front_down, back_up;
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

    DMStagGetLocationSlot(dmGrid, FRONT, 0, &front);
    DMStagGetLocationSlot(dmGrid, FRONT_UP, 0, &front_up);
    DMStagGetLocationSlot(dmGrid, BACK_DOWN, 0, &back_down);
    DMStagGetLocationSlot(dmGrid, FRONT_DOWN, 0, &front_down);
    DMStagGetLocationSlot(dmGrid, BACK_UP, 0, &back_up);


    Vec vec_local;
    PetscReal ****arrvec_local;
    DMCreateLocalVector(dmGrid,&vec_local);
    DMGlobalToLocalBegin(dmGrid,vec,INSERT_VALUES,vec_local);
    DMGlobalToLocalEnd(dmGrid,vec,INSERT_VALUES,vec_local);
    DMStagVecGetArrayRead(dmGrid, vec_local, &arrvec_local);

    Vec WShifted_local;
    VecSet(WShifted, 0.0);
    PetscReal ****arrWShifted_local;
    DMGetLocalVector(dmGrid, &WShifted_local);
    VecSet(WShifted_local, 0.0);
    DMStagVecGetArray(dmGrid, WShifted_local, &arrWShifted_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez != N[2] - 1 and ey != N[1] - 1) {
                    PetscReal current = arrvec_local[ez][ey][ex][front];
                    PetscReal next = arrvec_local[ez][ey + 1][ex][front];
                    PetscReal inter = (next + current) / 2.0;
                    arrWShifted_local[ez][ey][ex][front_up] = inter;
                }

                if (ez == 0) {
                    PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_back_down[0]], arrCoord[ez][ey][ex][icux_back_down[1]],
                                  arrCoord[ez][ey][ex][icux_back_down[2]], theta);
                    arrWShifted_local[ez][ey][ex][back_down] = inter;
                }

                if (ey == 0) {
                    PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_front_down[0]], arrCoord[ez][ey][ex][icux_front_down[1]],
                                  arrCoord[ez][ey][ex][icux_front_down[2]], theta);
                    arrWShifted_local[ez][ey][ex][front_down] = inter;
                }

                if (ez == N[2] - 1) {
                    PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_front_up[0]], arrCoord[ez][ey][ex][icux_front_up[1]],
                                  arrCoord[ez][ey][ex][icux_front_up[2]], theta);
                    arrWShifted_local[ez][ey][ex][front_up] = inter;
                }

                if (ey == N[1] - 1) {
                    PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_back_up[0]], arrCoord[ez][ey][ex][icux_back_up[1]],
                                  arrCoord[ez][ey][ex][icux_back_up[2]], theta);
                    arrWShifted_local[ez][ey][ex][back_up] = inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagVecRestoreArrayRead(dmGrid, vec_local, &arrvec_local);
    DMStagVecRestoreArray(dmGrid, WShifted_local, &arrWShifted_local);
    DMLocalToGlobal(dmGrid, WShifted_local, ADD_VALUES, WShifted);
    DMRestoreLocalVector(dmGrid, &vec_local);
    DMRestoreLocalVector(dmGrid, &WShifted_local);

    return 0;
}

PetscErrorCode SecondDerive_x(DM const & dmGrid, Vec & AB_x, Vec const & AB) //ok
{
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, down, down_left, down_right, up, up_left, up_right;
    PetscReal       hx;
    PetscFunctionBeginUser;
   // DMCreateGlobalVector(dmGrid, pAB_x);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    hx = 1.0 / N[0];

    Vec AB_local;
    PetscScalar ****arrAB_local;
    DMCreateLocalVector(dmGrid, &AB_local);
    DMGlobalToLocalBegin(dmGrid, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmGrid, AB, INSERT_VALUES, AB_local);
    DMStagVecGetArrayRead(dmGrid, AB_local, &arrAB_local);

    Vec AB_x_local;
    VecSet(AB_x, 0.0);
    PetscScalar ****arrAB_x_local;
    DMGetLocalVector(dmGrid, &AB_x_local);
    VecSet(AB_x_local, 0.0);
    DMStagVecGetArray(dmGrid, AB_x_local, &arrAB_x_local);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &down);
    DMStagGetLocationSlot(dmGrid, DOWN_LEFT, 0, &down_left);
    DMStagGetLocationSlot(dmGrid, DOWN_RIGHT, 0, &down_right);
    DMStagGetLocationSlot(dmGrid, UP, 0, &up);
    DMStagGetLocationSlot(dmGrid, UP_LEFT, 0, &up_left);
    DMStagGetLocationSlot(dmGrid, UP_RIGHT, 0, &up_right);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal val_left = arrAB_local[ez][ey][ex][down_left];
                PetscReal val_right = arrAB_local[ez][ey][ex][down_right];
                PetscReal der = (val_right - val_left) / hx;
                arrAB_x_local[ez][ey][ex][down] = der;

                if (ey == N[1] - 1) {
                    val_left = arrAB_local[ez][ey][ex][up_left];
                    val_right = arrAB_local[ez][ey][ex][up_right];
                    der = (val_right - val_left) / hx;
                    arrAB_x_local[ez][ey][ex][up] = der;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, AB_local, &arrAB_local);
    DMStagVecRestoreArray(dmGrid, AB_x_local, &arrAB_x_local);
    DMLocalToGlobal(dmGrid, AB_x_local, ADD_VALUES, AB_x);
    DMRestoreLocalVector(dmGrid, &AB_x_local);
    DMRestoreLocalVector(dmGrid, &AB_local);
    
    return 0;
}

PetscErrorCode SecondDerive_z(DM const & dmGrid, Vec & AB_z, Vec const & AB) //ok
{
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, down, front_down, back_down, up, front_up, back_up;
    PetscReal       hz;
    PetscFunctionBeginUser;
   // DMCreateGlobalVector(dmGrid, pAB_z);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    hz = 1.0 / N[2];

    Vec AB_local;
    PetscScalar ****arrAB_local;
    DMCreateLocalVector(dmGrid, &AB_local);
    DMGlobalToLocalBegin(dmGrid, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmGrid, AB, INSERT_VALUES, AB_local);
    DMStagVecGetArrayRead(dmGrid, AB_local, &arrAB_local);

    Vec AB_z_local;
    VecSet(AB_z, 0.0);
    PetscScalar ****arrAB_z_local;
    DMGetLocalVector(dmGrid, &AB_z_local);
    VecSet(AB_z_local, 0.0);
    DMStagVecGetArray(dmGrid, AB_z_local, &arrAB_z_local);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &down);
    DMStagGetLocationSlot(dmGrid, FRONT_DOWN, 0, &front_down);
    DMStagGetLocationSlot(dmGrid, BACK_DOWN, 0, &back_down);
    DMStagGetLocationSlot(dmGrid, UP, 0, &up);
    DMStagGetLocationSlot(dmGrid, FRONT_UP, 0, &front_up);
    DMStagGetLocationSlot(dmGrid, BACK_UP, 0, &back_up);



    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal val_back = arrAB_local[ez][ey][ex][back_down];
                PetscReal val_front = arrAB_local[ez][ey][ex][front_down];
                PetscReal der = (val_front - val_back) / hz;
                arrAB_z_local[ez][ey][ex][down] = der;

                if (ey == N[1] - 1) {
                    val_back = arrAB_local[ez][ey][ex][back_up];
                    val_front = arrAB_local[ez][ey][ex][front_up];
                    der = (val_front - val_back) / hz;
                    arrAB_z_local[ez][ey][ex][up] = der;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, AB_local, &arrAB_local);
    DMStagVecRestoreArray(dmGrid, AB_z_local, &arrAB_z_local);
    DMLocalToGlobal(dmGrid, AB_z_local, ADD_VALUES, AB_z);
    DMRestoreLocalVector(dmGrid, &AB_z_local);
    DMRestoreLocalVector(dmGrid, &AB_local);


    return 0;
}


PetscErrorCode ThirdDerive_x(DM const & dmGrid, Vec & AB_x, Vec const & AB) //ok
{
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, back_left, back_right, back, front_left, front_right, front;
    PetscReal       hx;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pAB_x);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    hx = 1.0 / N[0];

    Vec AB_local;
    PetscScalar ****arrAB_local;
    DMCreateLocalVector(dmGrid, &AB_local);
    DMGlobalToLocalBegin(dmGrid, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmGrid, AB, INSERT_VALUES, AB_local);
    DMStagVecGetArrayRead(dmGrid, AB_local, &arrAB_local);

    Vec AB_x_local;
    VecSet(AB_x, 0.0);
    PetscScalar ****arrAB_x_local;
    DMGetLocalVector(dmGrid, &AB_x_local);
    VecSet(AB_x_local, 0.0);
    DMStagVecGetArray(dmGrid, AB_x_local, &arrAB_x_local);

    DMStagGetLocationSlot(dmGrid, BACK_LEFT, 0, &back_left);
    DMStagGetLocationSlot(dmGrid, BACK_RIGHT, 0, &back_right);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &back);
    DMStagGetLocationSlot(dmGrid, FRONT_LEFT, 0, &front_left);
    DMStagGetLocationSlot(dmGrid, FRONT_RIGHT, 0, &front_right);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &front);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal val_left = arrAB_local[ez][ey][ex][back_left];
                PetscReal val_right = arrAB_local[ez][ey][ex][back_right];
                PetscReal der = (val_right - val_left) / hx;
                arrAB_x_local[ez][ey][ex][back] = der;

                if (ez == N[2] - 1) {
                    val_left = arrAB_local[ez][ey][ex][front_left];
                    val_right = arrAB_local[ez][ey][ex][front_right];
                    der = (val_right - val_left) / hx;
                    arrAB_x_local[ez][ey][ex][front] = der;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, AB_local, &arrAB_local);
    DMStagVecRestoreArray(dmGrid, AB_x_local, &arrAB_x_local);
    DMLocalToGlobal(dmGrid, AB_x_local, ADD_VALUES, AB_x);
    DMRestoreLocalVector(dmGrid, &AB_x_local);
    DMRestoreLocalVector(dmGrid, &AB_local);


    return 0;
}

PetscErrorCode ThirdDerive_y(DM const & dmGrid, Vec & AB_y, Vec const & AB) //ok
{
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, back_down, back_up, back, front_down, front_up, front;
    PetscReal       hy;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pAB_y);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    hy = 1.0 / N[1];

    Vec AB_local;
    PetscScalar ****arrAB_local;
    DMCreateLocalVector(dmGrid, &AB_local);
    DMGlobalToLocalBegin(dmGrid, AB, INSERT_VALUES, AB_local);
    DMGlobalToLocalEnd(dmGrid, AB, INSERT_VALUES, AB_local);
    DMStagVecGetArrayRead(dmGrid, AB_local, &arrAB_local);

    Vec AB_y_local;
    VecSet(AB_y, 0.0);
    PetscScalar ****arrAB_y_local;
    DMGetLocalVector(dmGrid, &AB_y_local);
    VecSet(AB_y_local, 0.0);
    DMStagVecGetArray(dmGrid, AB_y_local, &arrAB_y_local);

    DMStagGetLocationSlot(dmGrid, BACK_DOWN, 0, &back_down);
    DMStagGetLocationSlot(dmGrid, BACK_UP, 0, &back_up);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &back);
    DMStagGetLocationSlot(dmGrid, FRONT_DOWN, 0, &front_down);
    DMStagGetLocationSlot(dmGrid, FRONT_UP, 0, &front_up);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &front);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {     

                PetscReal val_down = arrAB_local[ez][ey][ex][back_down];
                PetscReal val_up = arrAB_local[ez][ey][ex][back_up];
                PetscReal der = (val_up - val_down) / hy;
                arrAB_y_local[ez][ey][ex][back] = der;

                if (ez == N[2] - 1) {
                    val_down = arrAB_local[ez][ey][ex][front_down];
                    val_up = arrAB_local[ez][ey][ex][front_up];
                    der = (val_up - val_down) / hy;
                    arrAB_y_local[ez][ey][ex][front] = der;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmGrid, AB_local, &arrAB_local);
    DMStagVecRestoreArray(dmGrid, AB_y_local, &arrAB_y_local);
    DMLocalToGlobal(dmGrid, AB_y_local, ADD_VALUES, AB_y);
    DMRestoreLocalVector(dmGrid, &AB_y_local);
    DMRestoreLocalVector(dmGrid, &AB_local);

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
                    arrU[ez][ey][ex][iux_right] = -der;

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
                    arrV[ez][ey][ex][iuy_up] = -der;

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
                    arrW[ez][ey][ex][iuz_front] = -der;
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
PetscErrorCode ManageAdvection_x(PetscReal const & dt, Vec & U_int, Vec const & U_0, Vec const & V_0, Vec const & W_0, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz, PetscReal const & theta, Vec const & U_int_prev)
{
    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered;
    PetscFunctionBegin;
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    }

    Vec U_n, V_n, W_n;
    DMCreateGlobalVector(dmGrid_Staggered, &U_n);
    DMCreateGlobalVector(dmGrid_Staggered, &V_n);
    DMCreateGlobalVector(dmGrid_Staggered, &W_n);
    VecCopy(U_0, U_n);
    VecCopy(V_0, V_n);
    VecCopy(W_0, W_n);

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

    Vec U2_x;
    DMCreateGlobalVector(dmGrid_Centered, &U2_x);
    Derive_x(dmGrid_Centered, U2_x, U_c, theta);
    Vec homoFirst;
    DMCreateGlobalVector(dmGrid_Staggered, &homoFirst);
    DMStagMigrateVec(dmGrid_Centered, U2_x, dmGrid_Staggered, homoFirst);

    VecAXPBY(homoFirst,1,1, mixedFirst);
    VecAXPBY(homoFirst,-0.5, 1.5, U_int_prev); 
    VecCopy(homoFirst, U_int_prev);
    VecAXPBY(U_n,-dt, 1.0, homoFirst);
    
    //VecAXPBYPCZ(U_n, -dt, -dt, 1.0, homoFirst, mixedFirst);
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
    VecDestroy(&U_n);
    VecDestroy(&V_n);
    VecDestroy(&W_n);

    PetscFunctionReturn(0);
}

PetscErrorCode ManageAdvection_y(PetscReal const & dt, Vec & V_int, Vec const & U_0, Vec const & V_0, Vec const & W_0, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz, PetscReal const & theta, Vec const & V_int_prev)
{
    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered;
    PetscFunctionBegin;
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    }

    Vec U_n, V_n, W_n;
    DMCreateGlobalVector(dmGrid_Staggered, &U_n);
    DMCreateGlobalVector(dmGrid_Staggered, &V_n);
    DMCreateGlobalVector(dmGrid_Staggered, &W_n);
    VecCopy(U_0, U_n);
    VecCopy(V_0, V_n);
    VecCopy(W_0, W_n);

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

    VecAXPBY(homoSecond,1,1, mixedSecond);
    VecAXPBY(homoSecond,-0.5, 1.5, V_int_prev);
    VecCopy(homoSecond, V_int_prev);
    VecAXPBY(V_n,-dt, 1.0, homoSecond);


    //VecAXPBYPCZ(V_n, -dt, -dt, 1.0, homoSecond, mixedSecond);
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
    VecDestroy(&U_n);
    VecDestroy(&V_n);
    VecDestroy(&W_n);
    
    PetscFunctionReturn(0);
}

PetscErrorCode ManageAdvection_z(PetscReal const & dt, Vec & W_int, Vec const & U_0, Vec const & V_0, Vec const & W_0, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz, PetscReal const & theta, Vec const & W_int_prev)
{
    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered;
    PetscFunctionBegin;
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    }

    Vec U_n, V_n, W_n;
    DMCreateGlobalVector(dmGrid_Staggered, &U_n);
    DMCreateGlobalVector(dmGrid_Staggered, &V_n);
    DMCreateGlobalVector(dmGrid_Staggered, &W_n);
    VecCopy(U_0, U_n);
    VecCopy(V_0, V_n);
    VecCopy(W_0, W_n);

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

    VecAXPBY(homoThird,1,1, mixedThird);
    VecAXPBY(homoThird,-0.5, 1.5, W_int_prev);
    VecCopy(homoThird, W_int_prev);
    VecAXPBY(W_n,-dt, 1.0, homoThird);


    //VecAXPBYPCZ(W_n, -dt, -dt, 1.0, homoThird, mixedThird);
    VecCopy(W_n, W_int);

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
    VecDestroy(&U_n);
    VecDestroy(&V_n);
    VecDestroy(&W_n);

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

    for (ez = startz; ez < startz + nz; ++ez) { 
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
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iux_element];
                    second = arrVec[ez][ey][ex + 1][iux_element];
                    third = arrVec[ez][ey][ex + 2][iux_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hx);
                    arrOut[ez][ey][ex][iux_left] = inter;
                }        
                if(ex == N[0] - 1){
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iux_element];
                    second = arrVec[ez][ey][ex - 1][iux_element];
                    third = arrVec[ez][ey][ex - 2][iux_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hx);
                    arrOut[ez][ey][ex][iux_right] = -inter;
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
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iuy_element];
                    second = arrVec[ez][ey + 1][ex][iuy_element];
                    third = arrVec[ez][ey + 2][ex][iuy_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hy);
                    arrOut[ez][ey][ex][iuy_down] = inter;
                }
                if(ey == N[1] - 1){
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iuy_element];
                    second = arrVec[ez][ey - 1][ex][iuy_element];
                    third = arrVec[ez][ey - 2][ex][iuy_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hy);
                    arrOut[ez][ey][ex][iuy_up] = -inter;
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
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iuz_element];
                    second = arrVec[ez + 1][ey][ex][iuz_element];
                    third = arrVec[ez + 2][ey][ex][iuz_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hz);
                    arrOut[ez][ey][ex][iuz_back] = inter;
                }
                if(ez == N[2] - 1){
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iuz_element];
                    second = arrVec[ez - 1][ey][ex][iuz_element];
                    third = arrVec[ez - 2][ey][ex][iuz_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hz);
                    arrOut[ez][ey][ex][iuz_front] = -inter;
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


PetscErrorCode ManagePressure(DM const & dmGrid_centered, DM const & dmGrid_shifted, DM const & dmGrid_staggered, PetscReal const & dt, Vec & P, Vec const & U_pre, Vec const & V_pre, Vec const & W_pre)
{
    Mat A;
    Vec rhs;
    KSP ksp;
    PC  pc;

    PetscFunctionBegin;

    Vec U_n, V_n, W_n;
    DMCreateGlobalVector(dmGrid_staggered, &U_n);
    DMCreateGlobalVector(dmGrid_staggered, &V_n);
    DMCreateGlobalVector(dmGrid_staggered, &W_n);
    VecCopy(U_pre, U_n);
    VecCopy(V_pre, V_n);
    VecCopy(W_pre, W_n);

    DMCreateGlobalVector(dmGrid_centered, &rhs);
    DMCreateMatrix(dmGrid_centered, &A);

    Vec div;
    DMCreateGlobalVector(dmGrid_centered, &div);
    ComputeDivergence(dmGrid_centered, dmGrid_shifted, dmGrid_staggered, div, U_n, V_n, W_n, dt); 

    

    /*Vec force;
    DMCreateGlobalVector(dmGrid_centered, &force);
    CreateReferenceSolutionTryForce(dmGrid_centered, force, 0);

    CheckSolution(div, force);*/


    Assemble_P(dmGrid_centered, A, rhs, div);
    PetscReal mean;
    PetscInt size;
    VecSum(rhs, &mean);
    VecGetSize(rhs, &size);
    mean = mean / size;
    VecShift(rhs, -mean);
    VecDestroy(&div);

    //AttachNullspace(dmGrid_centered, A);

    //questo sistema non e' precondizionato: il preconditioner e' nullo. Bisogna farlo ed e' importante: punto piu' lento del codice
    /*KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    //KSPGetPC(ksp, &pc);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, P);*/

    /*KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCHYPRE);
    PCHYPRESetType(pc, "pilut");
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, P);*/

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);

    // Set the preconditioner type to 'none' to disable it
    PCSetType(pc, PCNONE);

    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, P);

    MatDestroy(&A);
    VecDestroy(&rhs);
    KSPDestroy(&ksp);
    VecDestroy(&U_n);
    VecDestroy(&V_n);
    VecDestroy(&W_n);
       
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

    Vec U_int_prev, V_int_prev, W_int_prev;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_int_prev);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_int_prev);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_int_prev);
    VecZeroEntries(U_int_prev);
    VecZeroEntries(V_int_prev);
    VecZeroEntries(W_int_prev);

    /*ManagePressure(dmGrid_Centered, dmGrid_Shifted, dmGrid_Staggered, dt, P, U_0, V_0, W_0);

    ManagePressure_y(dmGrid_Staggered_y, dmGrid_Centered, dmGrid_Shifted, P_y, P);

    Vec prova;
    DMCreateGlobalVector(dmGrid_Staggered_y, &prova);
    CreateReferenceSolutionTry(dmGrid_Staggered_y, prova, 0);
    CheckSolution(P_y, prova);*/

    /*Vec force;
    DMCreateGlobalVector(dmGrid_Centered, &force);
    CreateReferenceSolutionTryForce(dmGrid_Centered, force, 0);*/




    /*PetscViewer viewer_p;
    DM DM_p;
    //DMStagCreateCompatibleDMStag(dmGrid_Centered, 0, 0, 0, 1, &DM_p);
    Vec p;
    DMStagVecSplitToDMDA(dmGrid_centered, P, ELEMENT, 0, &DM_p, &p);
    PetscObjectSetName((PetscObject)p, "shalom");
    char filename_p[50];
    sprintf(filename_p, "results/shalom.vtr");
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_p, FILE_MODE_WRITE, &viewer_p);
    VecView(p, viewer_p);
    VecDestroy(&p);
    DMDestroy(&DM_p);
    PetscViewerDestroy(&viewer_p);*/





    /*ManageAdvection_x(dt, U_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);
    ManageAdvection_y(dt, V_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);
    ManageAdvection_z(dt, W_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta);
    Vec laplacian;
    DMCreateGlobalVector(dmGrid_Staggered_x, &laplacian);
    laplacian_x_explicit(dmGrid_Staggered_x, laplacian, U_0);
    Vec bench;
    DMCreateGlobalVector(dmGrid_Staggered_z, &bench);
    CreateReferenceSolutionTry(dmGrid_Staggered_z, bench, 0);
    CheckSolution(W_int, bench);*/


    /*for(size_t i = 0; i < iter; i++){
            if (i == 0){
            
                theta = (1/Re)*(i+1)*dt*pi*pi/3;

                ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_up, V_up, W_up, U_0, V_0, W_0, theta);




                Vec bench;
                DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
                CreateReferenceSolutionTry(dmGrid_Staggered_x, bench, theta);

                CheckSolution(U_up, bench);
                VecDestroy(&bench);

                            PetscViewer viewer_u;
            DM DM_u;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmGrid_Staggered_x, U_up, LEFT, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "solution");
            char filename_u[50]; 
            sprintf(filename_u, "results/solution%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u); 



            } else {

                std::cout<<i<<std::endl;
                
                theta = (1/Re)*(i+1)*dt*pi*pi/3;


                //theta = d*d*(i+1)*dt;
                ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_up, V_up, W_up, U_up, V_up, W_up, theta);

                //std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;

                
            //std::cout << "Spirit of Oklahoma completed:   pressure done." << std::endl;
                //theta = (1/Re)*k*k*dt*(i+1);
                //std::cout << "Spirit of California completed: iteration " << i << " done" << std::endl;
                Vec bench;
                DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
                CreateReferenceSolutionTry(dmGrid_Staggered_x, bench, theta);

                CheckSolution(U_up, bench);
                VecDestroy(&bench);

                                            PetscViewer viewer_u;
            DM DM_u;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmGrid_Staggered_x, U_up, LEFT, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "solution");
            char filename_u[50]; 
            sprintf(filename_u, "results/solution%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u); 



                  
    
            }
        }*/



    /*parabolic_problem_x parabolic_x(dmGrid_Staggered_x, U_0);

    parabolic_x.solve();

    parabolic_x.~parabolic_problem_x();

    parabolic_problem_y parabolic_y(dmGrid_Staggered_y, V_0);

    parabolic_y.solve();

    parabolic_y.~parabolic_problem_y();

    parabolic_problem_z parabolic_z(dmGrid_Staggered_z, W_0);

    parabolic_z.solve();

    parabolic_z.~parabolic_problem_z();

    stokes_problem stokes(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Staggered, dmGrid_Centered, dmGrid_Shifted, U_0, V_0, W_0);

    stokes.solve();

    stokes.~stokes_problem();*/

    /*euler_problem euler(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Staggered, dmGrid_Centered, dmGrid_Shifted, U_0, V_0, W_0);

    euler.solve();

    euler.~euler_problem();*/
    


    for(size_t i = 0; i < iter; ++i){
        if (i == 0){
            
            theta = d*d*(i)*dt;
            //theta = 3*(vRef*A/Re)*k*k*dt*(i);

            ManageAdvection_x(dt, U_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, U_int_prev);
            ManageAdvection_y(dt, V_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, V_int_prev);
            ManageAdvection_z(dt, W_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, W_int_prev);

            //theta = d*d*(i+1)*dt;

            ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_pre, V_pre, W_pre, U_int, V_int, W_int, theta);

            //std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
            ManagePressure(dmGrid_Centered, dmGrid_Shifted, dmGrid_Staggered, dt, P, U_pre, V_pre, W_pre);
            ManagePressure_x(dmGrid_Staggered_x, dmGrid_Centered, dmGrid_Shifted, P_x, P);
            ManagePressure_y(dmGrid_Staggered_y, dmGrid_Centered, dmGrid_Shifted, P_y, P);
            ManagePressure_z(dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, P_z, P);
            
            //std::cout << "Spirit of Oklahoma completed:   pressure done." << std::endl;
            //theta = 3*(vRef*A/Re)*k*k*dt*(i+1);
            theta = d*d*(i+1)*dt;

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


        } else {
                

            theta = d*d*(i)*dt;
            //theta = 3*(vRef*A/Re)*k*k*dt*(i);

            ManageAdvection_x(dt, U_int, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, U_int_prev);
            ManageAdvection_y(dt, V_int, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, V_int_prev);
            ManageAdvection_z(dt, W_int, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, W_int_prev);    

            //std::cout << "Spirit of Kitty Hawk completed: advection done." << std::endl;

            //theta = d*d*(i+1)*dt;
            ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_pre, V_pre, W_pre, U_int, V_int, W_int, theta);

            //std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
            ManagePressure(dmGrid_Centered, dmGrid_Shifted, dmGrid_Staggered, dt, P, U_pre, V_pre, W_pre);
            Vec pressure_bench;
            DMCreateGlobalVector(dmGrid_Centered, &pressure_bench);
            CreateReferencePressure(dmGrid_Centered, pressure_bench, theta);
            //CheckSolution(P, pressure_bench);
            VecDestroy(&pressure_bench);

            ManagePressure_x(dmGrid_Staggered_x, dmGrid_Centered, dmGrid_Shifted, P_x, P);
            ManagePressure_y(dmGrid_Staggered_y, dmGrid_Centered, dmGrid_Shifted, P_y, P);
            ManagePressure_z(dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, P_z, P);

            //std::cout << "Spirit of Oklahoma completed:   pressure done." << std::endl;
            //theta = 3*(vRef*A/Re)*k*k*dt*(i+1);
            theta = d*d*(i+1)*dt;
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
              
   
        }
    }


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
    VecDestroy(&U_int_prev);
    VecDestroy(&V_int_prev);
    VecDestroy(&W_int_prev);

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



