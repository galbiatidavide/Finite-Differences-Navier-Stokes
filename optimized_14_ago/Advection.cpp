#include "Setting.cpp"

#define pi 3.1415926535

#define RIGHT            DMSTAG_RIGHT
#define UP_LEFT          DMSTAG_UP_LEFT
#define UP_RIGHT         DMSTAG_UP_RIGHT
#define DOWN_LEFT        DMSTAG_DOWN_LEFT
#define DOWN_RIGHT       DMSTAG_DOWN_RIGHT
#define BACK_LEFT        DMSTAG_BACK_LEFT
#define BACK_RIGHT       DMSTAG_BACK_RIGHT
#define FRONT_LEFT       DMSTAG_FRONT_LEFT
#define FRONT_RIGHT      DMSTAG_FRONT_RIGHT
#define FRONT_DOWN       DMSTAG_FRONT_DOWN
#define LEFT             DMSTAG_LEFT
#define BACK_DOWN        DMSTAG_BACK_DOWN
#define BACK             DMSTAG_BACK
#define BACK_UP          DMSTAG_BACK_UP
#define DOWN             DMSTAG_DOWN
#define ELEMENT          DMSTAG_ELEMENT
#define UP               DMSTAG_UP
#define FRONT            DMSTAG_FRONT
#define FRONT_UP         DMSTAG_FRONT_UP

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

PetscErrorCode FirstDerive_y(DM const & dmGrid, Vec & AB_y, Vec const & vec) 
{
    PetscInt iux_up_right, iux_down_right, iux_up_left, iux_down_left, iux_left, iux_right;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecULocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrU;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hy = 1.0 / N[1];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid, UP_LEFT, 0, &iux_up_left);
    DMStagGetLocationSlot(dmGrid, UP_RIGHT, 0, &iux_up_right);
    DMStagGetLocationSlot(dmGrid, DOWN_LEFT, 0, &iux_down_left);
    DMStagGetLocationSlot(dmGrid, DOWN_RIGHT, 0, &iux_down_right);

    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecULocal);
    DMStagVecGetArray(dmGrid, vecULocal, &arrU);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal inter, next, prev;
                next = arrVec[ez][ey][ex][iux_up_left];
                prev = arrVec[ez][ey][ex][iux_down_left];
                inter = (next - prev) / hy;
                arrU[ez][ey][ex][iux_left] = inter;

                if (ex == N[0] - 1) {
                    PetscReal inter, next, prev;
                    next = arrVec[ez][ey][ex][iux_up_right];
                    prev = arrVec[ez][ey][ex][iux_down_right];
                    inter = (next - prev) / hy; 
                    arrU[ez][ey][ex][iux_right] = inter;
                }

            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecULocal, &arrU);
    DMLocalToGlobal(dmGrid, vecULocal, INSERT_VALUES, AB_y);
    DMRestoreLocalVector(dmGrid, &vecULocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode FirstDerive_z(DM const & dmGrid, Vec & AB_z, Vec const & vec) 
{
    PetscInt iux_left, iux_right, iux_back_left, iux_back_right, iux_front_left, iux_front_right;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecULocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrU;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid, BACK_LEFT, 0, &iux_back_left);
    DMStagGetLocationSlot(dmGrid, BACK_RIGHT, 0, &iux_back_right);
    DMStagGetLocationSlot(dmGrid, FRONT_LEFT, 0, &iux_front_left);
    DMStagGetLocationSlot(dmGrid, FRONT_RIGHT, 0, &iux_front_right);
    
    DMCreateLocalVector(dmGrid, &vecLocal);
    DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid, &vecULocal);
    DMStagVecGetArray(dmGrid, vecULocal, &arrU);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal inter, next, prev;
                next = arrVec[ez][ey][ex][iux_front_left];
                prev = arrVec[ez][ey][ex][iux_back_left];
                inter = (next - prev) / hz;
                arrU[ez][ey][ex][iux_left] = inter;

                if (ex == N[0] - 1) {
                    PetscReal inter, next, prev;
                    next = arrVec[ez][ey][ex][iux_front_right];
                    prev = arrVec[ez][ey][ex][iux_back_right];
                    inter = (next - prev) / hz;
                    arrU[ez][ey][ex][iux_right] = inter;
                }
                
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid, vecULocal, &arrU);
    DMLocalToGlobal(dmGrid, vecULocal, INSERT_VALUES, AB_z);
    DMRestoreLocalVector(dmGrid, &vecULocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
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

// Third non-linear members: mixed
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
                    next = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]], theta);
                    inter = ((prev + next) / 2.0) * ((prev + next) / 2.0);
                    arrU[ez][ey][ex][iux_element] = inter;
                } else if(ex == 0) {
                    PetscReal inter, prev, next;
                    prev = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]], arrCoord[ez][ey][ex][icux_left[2]], theta);
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
                    next = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]], theta);
                    inter = ((prev + next) / 2.0) * ((prev + next) / 2.0);
                    arrV[ez][ey][ex][iuy_element] = inter;
                 } else if(ey == 0) {
                    PetscReal inter, prev, next;
                    prev = uyRef(arrCoord[ez][ey][ex][icuy_down[0]], arrCoord[ez][ey][ex][icuy_down[1]], arrCoord[ez][ey][ex][icuy_down[2]], theta);
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
                    next = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]], theta);
                    inter = ((prev + next) / 2.0) * ((prev + next) / 2.0);
                    arrW[ez][ey][ex][iuz_element] = inter;
                } else if(ez == 0) {
                    PetscReal inter, prev, next;
                    prev = uzRef(arrCoord[ez][ey][ex][icuz_back[0]], arrCoord[ez][ey][ex][icuz_back[1]], arrCoord[ez][ey][ex][icuz_back[2]], theta);
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
                    PetscReal inter, prev, next;
                    next = arrVec[ez][ey][ex][iux_element];
                    prev = uxRef(arrCoord[ez][ey][ex][icux_element[0]] - hx, arrCoord[ez][ey][ex][icux_element[1]], arrCoord[ez][ey][ex][icux_element[2]], theta);
                    inter = (next - prev) / hx;
                    arrU[ez][ey][ex][iux_left] = inter;
                }        
                if(ex == N[0] - 1){
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex][iux_element];
                    next = uxRef(arrCoord[ez][ey][ex][icux_element[0]] + hx, arrCoord[ez][ey][ex][icux_element[1]], arrCoord[ez][ey][ex][icux_element[2]], theta);
                    inter = (next - prev) / hx;
                    arrU[ez][ey][ex][iux_right] = inter;
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
                    PetscReal inter, prev, next;
                    next = arrVec[ez][ey][ex][iuy_element];
                    prev = uyRef(arrCoord[ez][ey][ex][icuy_element[0]], arrCoord[ez][ey][ex][icuy_element[1]] - hy, arrCoord[ez][ey][ex][icuy_element[2]], theta);
                    inter = (next - prev) / hy;
                    arrV[ez][ey][ex][iuy_down] = inter;
                }

                if(ey == N[1] - 1){
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex][iuy_element];
                    next = uyRef(arrCoord[ez][ey][ex][icuy_element[0]], arrCoord[ez][ey][ex][icuy_element[1]] + hy, arrCoord[ez][ey][ex][icuy_element[2]], theta);
                    inter = (next - prev) / hy;
                    arrV[ez][ey][ex][iuy_up] = inter;

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
                    PetscReal inter, prev, next;
                    next = arrVec[ez][ey][ex][iuz_element];
                    prev = uzRef(arrCoord[ez][ey][ex][icuz_element[0]], arrCoord[ez][ey][ex][icuz_element[1]], arrCoord[ez][ey][ex][icuz_element[2]] - hz, theta);
                    inter = (next - prev) / hz;
                    arrW[ez][ey][ex][iuz_back] = inter;
                }
                if(ez == N[2] - 1){
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex][iuz_element];
                    next = uzRef(arrCoord[ez][ey][ex][icuz_element[0]], arrCoord[ez][ey][ex][icuz_element[1]], arrCoord[ez][ey][ex][icuz_element[2]] + hz, theta);
                    inter = (next - prev) / hz;
                    arrW[ez][ey][ex][iuz_front] = inter;
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

    VecAXPBYPCZ(W_n, -dt, -dt, 1.0, homoThird, mixedThird);
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

    PetscFunctionReturn(0);
}

