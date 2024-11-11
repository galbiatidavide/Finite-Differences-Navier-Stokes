#include <petscdmstag.h>
#include <petscksp.h>
#include <petscdm.h>

#include "Setting.hpp"

#ifndef TRANSPORT_PROBLEM_X_HPP
#define TRANSPORT_PROBLEM_X_HPP

class transport_problem_x {
protected:

    DM const & dmGrid_Shifted;
    DM const & dmGrid_Centered;
    DM const & dmGrid_Staggered;

    PetscReal const & dt;
    PetscReal const & iter;

    Vec U_0, V_0, W_0;
    Vec U_int;

    PetscErrorCode const FirstShiftU_y(Vec & UShifted, Vec const & vec, PetscScalar const & theta) //okok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, up_right, down_left, up_left, left, down_right;
        PetscInt icux_right[3], icux_up_left[3], icux_up_right[3], icux_down_left[3], icux_down_right[3], icux_left[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hy = 1.0 / N[1];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icux_up_left[d]);
            DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icux_up_right[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icux_down_left[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icux_down_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &right);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_RIGHT, 0, &up_right);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_LEFT, 0, &down_left);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_LEFT, 0, &up_left);
        DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &left);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_RIGHT, 0, &down_right);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec UShifted_local;
        VecSet(UShifted, 0.0);
        PetscReal ****arrUShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &UShifted_local);
        VecSet(UShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);    

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
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]] - hy,
                                    arrCoord[ez][ey][ex][icux_left[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][left];
                        PetscReal inter = (out + inner) / 2.0;
                        arrUShifted_local[ez][ey][ex][down_left] = inter;                 
                    }

                    if (ex == 0) {                  
                        PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_up_left[0]], arrCoord[ez][ey][ex][icux_up_left[1]],
                                    arrCoord[ez][ey][ex][icux_up_left[2]], theta);
                        arrUShifted_local[ez][ey][ex][up_left] = inter;                  
                    }

                    if (ey == N[1] - 1) {
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]] + hy,
                                    arrCoord[ez][ey][ex][icux_right[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][right];
                        PetscReal inter = (out + inner) / 2.0;
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
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, UShifted_local, ADD_VALUES, UShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &UShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstShiftU_z(Vec & UShifted, Vec const & vec, PetscScalar const & theta) //okok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, front_right, back_left, front_left, back_right, left;
        PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3], icux_left[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hz = 1.0 / N[2];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
            DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &left);
        DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &right);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_RIGHT, 0, &front_right);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_LEFT, 0, &back_left);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_LEFT, 0, &front_left);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_RIGHT, 0, &back_right);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec UShifted_local;
        VecSet(UShifted, 0.0);
        PetscReal ****arrUShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &UShifted_local);
        VecSet(UShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);


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
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]],
                                    arrCoord[ez][ey][ex][icux_left[2]] - hz, theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][left];
                        PetscReal inter = (out + inner) / 2.0;
                        arrUShifted_local[ez][ey][ex][back_left] = inter;
                    }

                    if (ex == 0) {
                        PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_front_left[0]], arrCoord[ez][ey][ex][icux_front_left[1]],
                                    arrCoord[ez][ey][ex][icux_front_left[2]], theta);
                        arrUShifted_local[ez][ey][ex][front_left] = inter;
                    }

                    if (ez == N[2] - 1) {
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]],
                                    arrCoord[ez][ey][ex][icux_right[2]] + hz, theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][right];
                        PetscReal inter = (out + inner) / 2.0;                  
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
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, UShifted_local, ADD_VALUES, UShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &UShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstShiftV_y(Vec & VShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBeginUser;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, up, up_right, down_left, up_left, down_right, down;
        PetscInt icuy[3], icuy_up_left[3], icuy_up_right[3], icuy_down_left[3], icuy_down_right[3], icuy_up[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hx = 1.0 / N[0];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
            DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
            DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icuy_up_left[d]);
            DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icuy_up_right[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icuy_down_left[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icuy_down_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &down);
        DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &up);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_RIGHT, 0, &up_right);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_LEFT, 0, &down_left);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_LEFT, 0, &up_left);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_RIGHT, 0, &down_right);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec VShifted_local;
        VecSet(VShifted, 0.0);
        PetscReal ****arrVShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &VShifted_local);
        VecSet(VShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);


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
                        PetscReal out = uyRef(arrCoord[ez][ey][ex][icuy_up[0]] - hx, arrCoord[ez][ey][ex][icuy_up[1]],
                                    arrCoord[ez][ey][ex][icuy_up[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][up];
                        PetscReal inter = (out + inner) / 2.0;                    
                        arrVShifted_local[ez][ey][ex][up_left] = inter;
                    }

                    if (ey == N[1] - 1) {
                        PetscReal inter = uyRef(arrCoord[ez][ey][ex][icuy_up_right[0]], arrCoord[ez][ey][ex][icuy_up_right[1]],
                                    arrCoord[ez][ey][ex][icuy_up_right[2]], theta);
                        arrVShifted_local[ez][ey][ex][up_right] = inter;
                    }

                    if (ex == N[0] - 1) {
                        PetscReal out = uyRef(arrCoord[ez][ey][ex][icuy[0]] + hx, arrCoord[ez][ey][ex][icuy[1]],
                                    arrCoord[ez][ey][ex][icuy[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][down];
                        PetscReal inter = (out + inner) / 2.0;
                        arrVShifted_local[ez][ey][ex][down_right] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, VShifted_local, ADD_VALUES, VShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &VShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, front, front_right, back_left, front_left, back_right, back;
        PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3], icux_front[3], icux_back[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hx = 1.0 / N[0];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, BACK, d, &icux_back[d]);
            DMStagGetLocationSlot(dmCoord, FRONT, d, &icux_front[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
            DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &back);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &front);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_RIGHT, 0, &front_right);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_LEFT, 0, &back_left);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_LEFT, 0, &front_left);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_RIGHT, 0, &back_right);


        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec WShifted_local;
        VecSet(WShifted, 0.0);
        PetscReal ****arrWShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &WShifted_local);
        VecSet(WShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);

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
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_front[0]] - hx, arrCoord[ez][ey][ex][icux_front[1]],
                                    arrCoord[ez][ey][ex][icux_front[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][front];
                        PetscReal inter = (out + inner) / 2.0;                    
                        arrWShifted_local[ez][ey][ex][front_left] = inter;
                    }

                    if (ez == N[2] - 1) {
                        PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_front_right[0]], arrCoord[ez][ey][ex][icux_front_right[1]],
                                    arrCoord[ez][ey][ex][icux_front_right[2]], theta);
                        arrWShifted_local[ez][ey][ex][front_right] = inter;
                    }

                    if (ex == N[0] - 1) {
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_back[0]] + hx, arrCoord[ez][ey][ex][icux_back[1]],
                                    arrCoord[ez][ey][ex][icux_back[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][back];
                        PetscReal inter = (out + inner) / 2.0;
                        arrWShifted_local[ez][ey][ex][back_right] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, WShifted_local, ADD_VALUES, WShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &WShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstDerive_y(Vec & AB_y, Vec const & AB) //ok
    {
        PetscFunctionBegin;
        PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, down_left, up_left, left, down_right, up_right, right;
        PetscReal       hy;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        hy = 1.0 / N[1];
        Vec AB_local;
        PetscScalar ****arrAB_local;
        DMCreateLocalVector(dmGrid_Shifted, &AB_local);
        DMGlobalToLocalBegin(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMGlobalToLocalEnd(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);
        Vec AB_y_local;
        VecSet(AB_y, 0.0);
        PetscScalar ****arrAB_y_local;
        DMGetLocalVector(dmGrid_Shifted, &AB_y_local);
        VecSet(AB_y_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, AB_y_local, &arrAB_y_local);
        DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &left);
        DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &right);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_LEFT, 0, &down_left);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_RIGHT, 0, &down_right);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_LEFT, 0, &up_left);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_RIGHT, 0, &up_right);


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

        DMStagVecRestoreArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);
        DMStagVecRestoreArray(dmGrid_Shifted, AB_y_local, &arrAB_y_local);
        DMLocalToGlobal(dmGrid_Shifted, AB_y_local, ADD_VALUES, AB_y);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_y_local);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstDerive_z(Vec & AB_z, Vec const & AB) //ok
    {
        PetscFunctionBegin;
        PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, back_left, front_left, left, back_right, front_right, right;
        PetscReal       hz;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        hz = 1.0 / N[2];

        Vec AB_local;
        PetscScalar ****arrAB_local;
        DMCreateLocalVector(dmGrid_Shifted, &AB_local);
        DMGlobalToLocalBegin(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMGlobalToLocalEnd(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);

        Vec AB_z_local;
        VecSet(AB_z, 0.0);
        PetscScalar ****arrAB_z_local;
        DMGetLocalVector(dmGrid_Shifted, &AB_z_local);
        VecSet(AB_z_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, AB_z_local, &arrAB_z_local);

        DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &left);
        DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &right);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_LEFT, 0, &back_left);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_RIGHT, 0, &back_right);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_LEFT, 0, &front_left);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_RIGHT, 0, &front_right);

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

        DMStagVecRestoreArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);
        DMStagVecRestoreArray(dmGrid_Shifted, AB_z_local, &arrAB_z_local);
        DMLocalToGlobal(dmGrid_Shifted, AB_z_local, ADD_VALUES, AB_z);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_z_local);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const CenterU(Vec & UCenter, Vec const & vec, PetscReal const & theta)
    {
        PetscFunctionBegin;
        PetscInt icux_left[3], icux_right[3], iux_left, iux_right, iux_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
        DM dmCoord;
        Vec vecLocal, vecULocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec, ****arrU;    
        DMStagGetCorners(dmGrid_Centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Centered, &N[0], &N[1], &N[2]);
        DMGetCoordinateDM(dmGrid_Centered, &dmCoord);   
        DMGetCoordinates(dmGrid_Centered, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        }  
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid_Centered, LEFT, 0, &iux_left);
        DMStagGetLocationSlot(dmGrid_Centered, RIGHT, 0, &iux_right);
        DMStagGetLocationSlot(dmGrid_Centered, ELEMENT, 0, &iux_element);

        DMCreateLocalVector(dmGrid_Centered, &vecLocal);
        DMGlobalToLocalBegin(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid_Centered, vecLocal, &arrVec);

        DMGetLocalVector(dmGrid_Centered, &vecULocal);
        DMStagVecGetArray(dmGrid_Centered, vecULocal, &arrU);

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
        DMStagVecRestoreArrayRead(dmGrid_Centered, vecLocal, &arrVec);
        DMStagVecRestoreArray(dmGrid_Centered, vecULocal, &arrU);
        DMLocalToGlobal(dmGrid_Centered, vecULocal, INSERT_VALUES, UCenter);
        DMRestoreLocalVector(dmGrid_Centered, &vecULocal);
        DMRestoreLocalVector(dmGrid_Centered, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const Derive_x(Vec & U2_x, Vec const & vec, PetscReal const & theta)
    {
        PetscFunctionBegin;
        PetscInt icux_element[3], iux_left, iux_right, iux_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
        DM dmCoord;
        Vec vecLocal, vecULocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec, ****arrU;    
        DMStagGetCorners(dmGrid_Centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Centered, &N[0], &N[1], &N[2]);
        PetscReal const hx = 1.0 / N[0];
        DMGetCoordinateDM(dmGrid_Centered, &dmCoord);
        DMGetCoordinates(dmGrid_Centered, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icux_element[d]);
        }  
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid_Centered, LEFT, 0, &iux_left);
        DMStagGetLocationSlot(dmGrid_Centered, RIGHT, 0, &iux_right);
        DMStagGetLocationSlot(dmGrid_Centered, ELEMENT, 0, &iux_element);

        DMCreateLocalVector(dmGrid_Centered, &vecLocal);
        DMGlobalToLocalBegin(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid_Centered, vecLocal, &arrVec);

        DMGetLocalVector(dmGrid_Centered, &vecULocal);
        DMStagVecGetArray(dmGrid_Centered, vecULocal, &arrU);
    
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
        DMStagVecRestoreArrayRead(dmGrid_Centered, vecLocal, &arrVec);
        DMStagVecRestoreArray(dmGrid_Centered, vecULocal, &arrU);
        DMLocalToGlobal(dmGrid_Centered, vecULocal, INSERT_VALUES, U2_x);
        DMRestoreLocalVector(dmGrid_Centered, &vecULocal);
        DMRestoreLocalVector(dmGrid_Centered, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

public:

    transport_problem_x(ProblemSetting<Transport> const & setting, Vec const & U_0_, Vec const & V_0_, Vec const & W_0_) : 
    dmGrid_Shifted(setting.dmGrid_Shifted), dmGrid_Staggered(setting.dmGrid_Staggered), dmGrid_Centered(setting.dmGrid_Centered), dt(setting.dt),
    iter(setting.iter)
    {
        DMCreateGlobalVector(dmGrid_Staggered, &U_0);
        DMCreateGlobalVector(dmGrid_Staggered, &V_0);
        DMCreateGlobalVector(dmGrid_Staggered, &W_0);
        DMCreateGlobalVector(dmGrid_Staggered, &U_int);
        VecCopy(U_0_, U_0);
        VecCopy(V_0_, V_0);
        VecCopy(W_0_, W_0);
    }

    PetscErrorCode const solve_step_x(PetscScalar const & theta)
    {
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

        Vec UV_y, UW_z;
        DMCreateGlobalVector(dmGrid_Shifted, &UV_y);
        DMCreateGlobalVector(dmGrid_Shifted, &UW_z);
        Vec mixedFirst;
        DMCreateGlobalVector(dmGrid_Staggered, &mixedFirst);

        {        
            Vec U_y, U_z, V_y, W_z, UV, UW;
            DMCreateGlobalVector(dmGrid_Shifted, &U_y);
            DMCreateGlobalVector(dmGrid_Shifted, &U_z);
            DMCreateGlobalVector(dmGrid_Shifted, &V_y);
            DMCreateGlobalVector(dmGrid_Shifted, &W_z);
            DMCreateGlobalVector(dmGrid_Shifted, &UV);
            DMCreateGlobalVector(dmGrid_Shifted, &UW);
            FirstShiftU_y(U_y, U_shift, theta);
            FirstShiftU_z(U_z, U_shift, theta);
            FirstShiftV_y(V_y, V_shift, theta);
            FirstShiftW_z(W_z, W_shift, theta); 
            VecPointwiseMult(UV, U_y, V_y);
            VecPointwiseMult(UW, U_z, W_z);
            FirstDerive_y(UV_y, UV);
            FirstDerive_z(UW_z, UW);
            VecAXPY(UV_y, 1.0, UW_z);
            PetscObjectDestroy((PetscObject*)&U_y);
            PetscObjectDestroy((PetscObject*)&U_z);
            PetscObjectDestroy((PetscObject*)&V_y);
            PetscObjectDestroy((PetscObject*)&W_z);
            PetscObjectDestroy((PetscObject*)&UV);
            PetscObjectDestroy((PetscObject*)&UW);
        }           

        DMStagMigrateVec(dmGrid_Shifted, UV_y, dmGrid_Staggered, mixedFirst);

        Vec U_center;
        DMCreateGlobalVector(dmGrid_Centered, &U_center);
        DMStagMigrateVec(dmGrid_Staggered, U_n, dmGrid_Centered, U_center);

        Vec U_c;
        DMCreateGlobalVector(dmGrid_Centered, &U_c);
        CenterU(U_c, U_center, theta);

        Vec U2_x;
        DMCreateGlobalVector(dmGrid_Centered, &U2_x);
        Derive_x(U2_x, U_c, theta);

        Vec homoFirst;
        DMCreateGlobalVector(dmGrid_Staggered, &homoFirst);
        DMStagMigrateVec(dmGrid_Centered, U2_x, dmGrid_Staggered, homoFirst);
        VecAXPBYPCZ(U_n, -dt, -dt, 1.0, homoFirst, mixedFirst);
        VecCopy(U_n, U_int);
        //VecCopy(U_n, U_0);

        PetscObjectDestroy((PetscObject*)&U_shift);
        PetscObjectDestroy((PetscObject*)&V_shift);
        PetscObjectDestroy((PetscObject*)&W_shift);
        PetscObjectDestroy((PetscObject*)&UV_y);
        PetscObjectDestroy((PetscObject*)&UW_z);
        PetscObjectDestroy((PetscObject*)&mixedFirst);
        PetscObjectDestroy((PetscObject*)&U_c);
        PetscObjectDestroy((PetscObject*)&U2_x);
        PetscObjectDestroy((PetscObject*)&U_center);
        PetscObjectDestroy((PetscObject*)&homoFirst);

        VecDestroy(&U_n);
        VecDestroy(&V_n);
        VecDestroy(&W_n);

        PetscFunctionReturn(0);        
    }

    const Vec & get_U() const { return U_int; }


    ~transport_problem_x()
    {
        VecDestroy(&U_0);
        VecDestroy(&V_0);
        VecDestroy(&W_0);
        VecDestroy(&U_int);
   }

};

#endif // TRANSPORT_PROBLEM_X_HPP


#ifndef TRANSPORT_PROBLEM_Y_HPP
#define TRANSPORT_PROBLEM_Y_HPP

class transport_problem_y {
protected:

    DM const & dmGrid_Shifted;
    DM const & dmGrid_Centered;
    DM const & dmGrid_Staggered;

    PetscReal const & dt;
    PetscReal const & iter;

    Vec U_0, V_0, W_0;
    Vec V_int;

    PetscErrorCode const FirstShiftU_y(Vec & UShifted, Vec const & vec, PetscScalar const & theta) //okok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, up_right, down_left, up_left, left, down_right;
        PetscInt icux_right[3], icux_up_left[3], icux_up_right[3], icux_down_left[3], icux_down_right[3], icux_left[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hy = 1.0 / N[1];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icux_up_left[d]);
            DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icux_up_right[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icux_down_left[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icux_down_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &right);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_RIGHT, 0, &up_right);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_LEFT, 0, &down_left);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_LEFT, 0, &up_left);
        DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &left);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_RIGHT, 0, &down_right);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec UShifted_local;
        VecSet(UShifted, 0.0);
        PetscReal ****arrUShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &UShifted_local);
        VecSet(UShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);    

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
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]] - hy,
                                    arrCoord[ez][ey][ex][icux_left[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][left];
                        PetscReal inter = (out + inner) / 2.0;
                        arrUShifted_local[ez][ey][ex][down_left] = inter;                 
                    }

                    if (ex == 0) {                  
                        PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_up_left[0]], arrCoord[ez][ey][ex][icux_up_left[1]],
                                    arrCoord[ez][ey][ex][icux_up_left[2]], theta);
                        arrUShifted_local[ez][ey][ex][up_left] = inter;                  
                    }

                    if (ey == N[1] - 1) {
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]] + hy,
                                    arrCoord[ez][ey][ex][icux_right[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][right];
                        PetscReal inter = (out + inner) / 2.0;
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
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, UShifted_local, ADD_VALUES, UShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &UShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstShiftU_z(Vec & UShifted, Vec const & vec, PetscScalar const & theta) //okok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, front_right, back_left, front_left, back_right, left;
        PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3], icux_left[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hz = 1.0 / N[2];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
            DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &left);
        DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &right);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_RIGHT, 0, &front_right);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_LEFT, 0, &back_left);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_LEFT, 0, &front_left);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_RIGHT, 0, &back_right);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec UShifted_local;
        VecSet(UShifted, 0.0);
        PetscReal ****arrUShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &UShifted_local);
        VecSet(UShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);


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
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]],
                                    arrCoord[ez][ey][ex][icux_left[2]] - hz, theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][left];
                        PetscReal inter = (out + inner) / 2.0;
                        arrUShifted_local[ez][ey][ex][back_left] = inter;
                    }

                    if (ex == 0) {
                        PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_front_left[0]], arrCoord[ez][ey][ex][icux_front_left[1]],
                                    arrCoord[ez][ey][ex][icux_front_left[2]], theta);
                        arrUShifted_local[ez][ey][ex][front_left] = inter;
                    }

                    if (ez == N[2] - 1) {
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]],
                                    arrCoord[ez][ey][ex][icux_right[2]] + hz, theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][right];
                        PetscReal inter = (out + inner) / 2.0;                  
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
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, UShifted_local, ADD_VALUES, UShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &UShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstShiftV_y(Vec & VShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBeginUser;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, up, up_right, down_left, up_left, down_right, down;
        PetscInt icuy[3], icuy_up_left[3], icuy_up_right[3], icuy_down_left[3], icuy_down_right[3], icuy_up[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hx = 1.0 / N[0];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
            DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
            DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icuy_up_left[d]);
            DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icuy_up_right[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icuy_down_left[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icuy_down_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &down);
        DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &up);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_RIGHT, 0, &up_right);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_LEFT, 0, &down_left);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_LEFT, 0, &up_left);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_RIGHT, 0, &down_right);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec VShifted_local;
        VecSet(VShifted, 0.0);
        PetscReal ****arrVShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &VShifted_local);
        VecSet(VShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);


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
                        PetscReal out = uyRef(arrCoord[ez][ey][ex][icuy_up[0]] - hx, arrCoord[ez][ey][ex][icuy_up[1]],
                                    arrCoord[ez][ey][ex][icuy_up[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][up];
                        PetscReal inter = (out + inner) / 2.0;                    
                        arrVShifted_local[ez][ey][ex][up_left] = inter;
                    }

                    if (ey == N[1] - 1) {
                        PetscReal inter = uyRef(arrCoord[ez][ey][ex][icuy_up_right[0]], arrCoord[ez][ey][ex][icuy_up_right[1]],
                                    arrCoord[ez][ey][ex][icuy_up_right[2]], theta);
                        arrVShifted_local[ez][ey][ex][up_right] = inter;
                    }

                    if (ex == N[0] - 1) {
                        PetscReal out = uyRef(arrCoord[ez][ey][ex][icuy[0]] + hx, arrCoord[ez][ey][ex][icuy[1]],
                                    arrCoord[ez][ey][ex][icuy[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][down];
                        PetscReal inter = (out + inner) / 2.0;
                        arrVShifted_local[ez][ey][ex][down_right] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, VShifted_local, ADD_VALUES, VShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &VShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, front, front_right, back_left, front_left, back_right, back;
        PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3], icux_front[3], icux_back[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hx = 1.0 / N[0];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, BACK, d, &icux_back[d]);
            DMStagGetLocationSlot(dmCoord, FRONT, d, &icux_front[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
            DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &back);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &front);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_RIGHT, 0, &front_right);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_LEFT, 0, &back_left);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_LEFT, 0, &front_left);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_RIGHT, 0, &back_right);


        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec WShifted_local;
        VecSet(WShifted, 0.0);
        PetscReal ****arrWShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &WShifted_local);
        VecSet(WShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);

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
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_front[0]] - hx, arrCoord[ez][ey][ex][icux_front[1]],
                                    arrCoord[ez][ey][ex][icux_front[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][front];
                        PetscReal inter = (out + inner) / 2.0;                    
                        arrWShifted_local[ez][ey][ex][front_left] = inter;
                    }

                    if (ez == N[2] - 1) {
                        PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_front_right[0]], arrCoord[ez][ey][ex][icux_front_right[1]],
                                    arrCoord[ez][ey][ex][icux_front_right[2]], theta);
                        arrWShifted_local[ez][ey][ex][front_right] = inter;
                    }

                    if (ex == N[0] - 1) {
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_back[0]] + hx, arrCoord[ez][ey][ex][icux_back[1]],
                                    arrCoord[ez][ey][ex][icux_back[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][back];
                        PetscReal inter = (out + inner) / 2.0;
                        arrWShifted_local[ez][ey][ex][back_right] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, WShifted_local, ADD_VALUES, WShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &WShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const SecondShiftV_z(Vec & VShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, up, front_up, back_down, front_down, back_up, down;
        PetscInt icuy[3], icuy_front_up[3], icuy_front_down[3], icuy_back_up[3], icuy_back_down[3], icuy_up[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hz = 1.0 / N[2];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
            DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icuy_front_up[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icuy_front_down[d]);
            DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icuy_back_up[d]);
            DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icuy_back_down[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &down);
        DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &up);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_UP, 0, &front_up);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_DOWN, 0, &back_down);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_DOWN, 0, &front_down);  
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_UP, 0, &back_up);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec VShifted_local;
        VecSet(VShifted, 0.0);
        PetscReal ****arrVShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &VShifted_local);
        VecSet(VShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);


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
                        PetscScalar out = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]],
                                    arrCoord[ez][ey][ex][icuy_up[2]] + hz, theta);
                        PetscScalar inner = arrvec_local[ez][ey][ex][up];
                        PetscScalar inter = (out + inner) / 2.0;
                        arrVShifted_local[ez][ey][ex][front_up] = inter;
                    }

                    if (ey == N[1] - 1) {
                        PetscScalar inter = uyRef(arrCoord[ez][ey][ex][icuy_back_up[0]], arrCoord[ez][ey][ex][icuy_back_up[1]],
                                    arrCoord[ez][ey][ex][icuy_back_up[2]], theta);
                        arrVShifted_local[ez][ey][ex][back_up] = inter;
                    }

                    if (ez == 0) {
                        PetscScalar out = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]],
                                    arrCoord[ez][ey][ex][icuy[2]] - hz, theta);
                        PetscScalar inner = arrvec_local[ez][ey][ex][down];
                        PetscScalar inter = (out + inner) / 2.0;
                        arrVShifted_local[ez][ey][ex][back_down] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, VShifted_local, ADD_VALUES, VShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &VShifted_local);
        
        PetscFunctionReturn(0);
    }

    PetscErrorCode const SecondShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBeginUser;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, front, front_up, back_down, front_down, back_up, back;
        PetscInt icux_right[3], icux_back_up[3], icux_back_down[3], icux_front_up[3], icux_front_down[3], icux_front[3], icux_back[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hy = 1.0 / N[1];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, BACK, d, &icux_back[d]);
            DMStagGetLocationSlot(dmCoord, FRONT, d, &icux_front[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icux_front_up[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icux_front_down[d]);
            DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icux_back_up[d]);
            DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icux_back_down[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &back);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &front);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_UP, 0, &front_up);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_DOWN, 0, &back_down);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_DOWN, 0, &front_down);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_UP, 0, &back_up);


        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec WShifted_local;
        VecSet(WShifted, 0.0);
        PetscReal ****arrWShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &WShifted_local);
        VecSet(WShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);

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
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_front[0]], arrCoord[ez][ey][ex][icux_front[1]] - hy,
                                    arrCoord[ez][ey][ex][icux_front[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][front];
                        PetscReal inter = (out + inner) / 2.0;
                        arrWShifted_local[ez][ey][ex][front_down] = inter;
                    }

                    if (ez == N[2] - 1) {
                        PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_front_up[0]], arrCoord[ez][ey][ex][icux_front_up[1]],
                                    arrCoord[ez][ey][ex][icux_front_up[2]], theta);
                        arrWShifted_local[ez][ey][ex][front_up] = inter;
                    }

                    if (ey == N[1] - 1) {
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_back[0]], arrCoord[ez][ey][ex][icux_back[1]] + hy,
                                    arrCoord[ez][ey][ex][icux_back[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][back];
                        PetscReal inter = (out + inner) / 2.0;
                        arrWShifted_local[ez][ey][ex][back_up] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, WShifted_local, ADD_VALUES, WShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &WShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const SecondDerive_x(Vec & AB_x, Vec const & AB) //ok
    {
        PetscFunctionBeginUser;
        PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, down, down_left, down_right, up, up_left, up_right;
        PetscReal       hx;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        hx = 1.0 / N[0];
        Vec AB_local;
        PetscScalar ****arrAB_local;
        DMCreateLocalVector(dmGrid_Shifted, &AB_local);
        DMGlobalToLocalBegin(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMGlobalToLocalEnd(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);

        Vec AB_x_local;
        VecSet(AB_x, 0.0);
        PetscScalar ****arrAB_x_local;
        DMGetLocalVector(dmGrid_Shifted, &AB_x_local);
        VecSet(AB_x_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, AB_x_local, &arrAB_x_local);

        DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &down);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_LEFT, 0, &down_left);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_RIGHT, 0, &down_right);
        DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &up);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_LEFT, 0, &up_left);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_RIGHT, 0, &up_right);


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

        DMStagVecRestoreArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);
        DMStagVecRestoreArray(dmGrid_Shifted, AB_x_local, &arrAB_x_local);
        DMLocalToGlobal(dmGrid_Shifted, AB_x_local, ADD_VALUES, AB_x);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_x_local);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_local);
        
        PetscFunctionReturn(0);
    }

    PetscErrorCode const SecondDerive_z(Vec & AB_z, Vec const & AB) //ok
    {
        PetscFunctionBegin;
        PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, down, front_down, back_down, up, front_up, back_up;
        PetscReal       hz;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        hz = 1.0 / N[2];

        Vec AB_local;
        PetscScalar ****arrAB_local;
        DMCreateLocalVector(dmGrid_Shifted, &AB_local);
        DMGlobalToLocalBegin(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMGlobalToLocalEnd(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);

        Vec AB_z_local;
        VecSet(AB_z, 0.0);
        PetscScalar ****arrAB_z_local;
        DMGetLocalVector(dmGrid_Shifted, &AB_z_local);
        VecSet(AB_z_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, AB_z_local, &arrAB_z_local);

        DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &down);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_DOWN, 0, &front_down);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_DOWN, 0, &back_down);
        DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &up);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_UP, 0, &front_up);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_UP, 0, &back_up);



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

        DMStagVecRestoreArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);
        DMStagVecRestoreArray(dmGrid_Shifted, AB_z_local, &arrAB_z_local);
        DMLocalToGlobal(dmGrid_Shifted, AB_z_local, ADD_VALUES, AB_z);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_z_local);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const CenterV(Vec & VCenter, Vec const & vec, PetscReal const & theta) 
    {
        PetscFunctionBegin;
        PetscInt icuy_down[3], icuy_up[3], iuy_up, iuy_down, iuy_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
        DM dmCoord;
        Vec vecLocal, vecVLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec, ****arrV;    
        DMStagGetCorners(dmGrid_Centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Centered, &N[0], &N[1], &N[2]);
        DMGetCoordinateDM(dmGrid_Centered, &dmCoord);
        DMGetCoordinates(dmGrid_Centered, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy_down[d]);
            DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        } 
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid_Centered, DOWN, 0, &iuy_down);
        DMStagGetLocationSlot(dmGrid_Centered, UP, 0, &iuy_up);
        DMStagGetLocationSlot(dmGrid_Centered, ELEMENT, 0, &iuy_element);

        DMCreateLocalVector(dmGrid_Centered, &vecLocal);
        DMGlobalToLocalBegin(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid_Centered, vecLocal, &arrVec);

        DMGetLocalVector(dmGrid_Centered, &vecVLocal);
        DMStagVecGetArray(dmGrid_Centered, vecVLocal, &arrV); 

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
        DMStagVecRestoreArrayRead(dmGrid_Centered, vecLocal, &arrVec);
        DMStagVecRestoreArray(dmGrid_Centered, vecVLocal, &arrV);
        DMLocalToGlobal(dmGrid_Centered, vecVLocal, INSERT_VALUES, VCenter);
        DMRestoreLocalVector(dmGrid_Centered, &vecVLocal);
        DMRestoreLocalVector(dmGrid_Centered, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const Derive_y(Vec & V2_y, Vec const & vec, PetscReal const & theta)
    {
        PetscFunctionBegin;
        PetscInt icuy_element[3], iuy_up, iuy_down, iuy_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
        DM dmCoord;
        Vec vecLocal, vecVLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec, ****arrV;    
        DMStagGetCorners(dmGrid_Centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Centered, &N[0], &N[1], &N[2]);
        PetscReal const hy = 1.0 / N[1];
        DMGetCoordinateDM(dmGrid_Centered, &dmCoord);
        DMGetCoordinates(dmGrid_Centered, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icuy_element[d]);
        } 
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid_Centered, DOWN, 0, &iuy_down);
        DMStagGetLocationSlot(dmGrid_Centered, UP, 0, &iuy_up);
        DMStagGetLocationSlot(dmGrid_Centered, ELEMENT, 0, &iuy_element);

        DMCreateLocalVector(dmGrid_Centered, &vecLocal);
        DMGlobalToLocalBegin(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid_Centered, vecLocal, &arrVec);

        DMGetLocalVector(dmGrid_Centered, &vecVLocal);
        DMStagVecGetArray(dmGrid_Centered, vecVLocal, &arrV); 

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
        DMStagVecRestoreArrayRead(dmGrid_Centered, vecLocal, &arrVec);
        DMStagVecRestoreArray(dmGrid_Centered, vecVLocal, &arrV);
        DMLocalToGlobal(dmGrid_Centered, vecVLocal, INSERT_VALUES, V2_y);
        DMRestoreLocalVector(dmGrid_Centered, &vecVLocal);
        DMRestoreLocalVector(dmGrid_Centered, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

public:

    transport_problem_y(ProblemSetting<Transport> const & setting, Vec const & U_0_, Vec const & V_0_, Vec const & W_0_) : 
    dmGrid_Shifted(setting.dmGrid_Shifted), dmGrid_Staggered(setting.dmGrid_Staggered), dmGrid_Centered(setting.dmGrid_Centered), dt(setting.dt),
    iter(setting.iter)
    {
        DMCreateGlobalVector(dmGrid_Staggered, &U_0);
        DMCreateGlobalVector(dmGrid_Staggered, &V_0);
        DMCreateGlobalVector(dmGrid_Staggered, &W_0);
        DMCreateGlobalVector(dmGrid_Staggered, &V_int);
        VecCopy(U_0_, U_0);
        VecCopy(V_0_, V_0);
        VecCopy(W_0_, W_0);
    }

    PetscErrorCode const solve_step_y(PetscScalar const & theta)
    {
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

        Vec VU_x, VW_z;
        DMCreateGlobalVector(dmGrid_Shifted, &VU_x);
        DMCreateGlobalVector(dmGrid_Shifted, &VW_z);
        Vec mixedSecond;
        DMCreateGlobalVector(dmGrid_Staggered, &mixedSecond);

        {
            Vec V_x, V_z, U_x, W_z, VU, VW;
            DMCreateGlobalVector(dmGrid_Shifted, &V_x);
            DMCreateGlobalVector(dmGrid_Shifted, &V_z);
            DMCreateGlobalVector(dmGrid_Shifted, &U_x);
            DMCreateGlobalVector(dmGrid_Shifted, &W_z);
            DMCreateGlobalVector(dmGrid_Shifted, &VU);
            DMCreateGlobalVector(dmGrid_Shifted, &VW);
            FirstShiftV_y(V_x, V_shift, theta);
            SecondShiftV_z(V_z, V_shift, theta);
            FirstShiftU_y(U_x, U_shift, theta);
            SecondShiftW_z(W_z, W_shift, theta);
            VecPointwiseMult(VU, V_x, U_x);
            VecPointwiseMult(VW, V_z, W_z);
            SecondDerive_x(VU_x, VU);
            SecondDerive_z(VW_z, VW);
            VecAXPY(VU_x, 1.0, VW_z);
            PetscObjectDestroy((PetscObject*)&V_x);
            PetscObjectDestroy((PetscObject*)&V_z);
            PetscObjectDestroy((PetscObject*)&U_x);
            PetscObjectDestroy((PetscObject*)&W_z);
            PetscObjectDestroy((PetscObject*)&VU);
            PetscObjectDestroy((PetscObject*)&VW);
        }

        DMStagMigrateVec(dmGrid_Shifted, VU_x, dmGrid_Staggered, mixedSecond);

        Vec V_center;
        DMCreateGlobalVector(dmGrid_Centered, &V_center);
        DMStagMigrateVec(dmGrid_Staggered, V_n, dmGrid_Centered, V_center);

        Vec V_c;
        DMCreateGlobalVector(dmGrid_Centered, &V_c);
        CenterV(V_c, V_center, theta);

        Vec V2_y;
        DMCreateGlobalVector(dmGrid_Centered, &V2_y);
        Derive_y(V2_y, V_c, theta);

        Vec homoSecond;
        DMCreateGlobalVector(dmGrid_Staggered, &homoSecond);
        DMStagMigrateVec(dmGrid_Centered, V2_y, dmGrid_Staggered, homoSecond);

        VecAXPBYPCZ(V_n, -dt, -dt, 1.0, homoSecond, mixedSecond);
        VecCopy(V_n, V_int);
        //VecCopy(V_n, V_0);

        PetscObjectDestroy((PetscObject*)&U_shift);
        PetscObjectDestroy((PetscObject*)&V_shift);
        PetscObjectDestroy((PetscObject*)&W_shift);
        
        PetscObjectDestroy((PetscObject*)&VU_x);
        PetscObjectDestroy((PetscObject*)&VW_z);
        PetscObjectDestroy((PetscObject*)&mixedSecond);

        PetscObjectDestroy((PetscObject*)&V_c);
        PetscObjectDestroy((PetscObject*)&V2_y);
        PetscObjectDestroy((PetscObject*)&V_center);
        PetscObjectDestroy((PetscObject*)&homoSecond);

        VecDestroy(&U_n);
        VecDestroy(&V_n);
        VecDestroy(&W_n);

        PetscFunctionReturn(0);
    }

    const Vec & get_V() const { return V_int; }

    ~transport_problem_y()
    {
        VecDestroy(&U_0);
        VecDestroy(&V_0);
        VecDestroy(&W_0);
        VecDestroy(&V_int);
    }

};

#endif // TRANSPORT_PROBLEM__Y_HPP


#ifndef TRANSPORT_PROBLEM_Z_HPP
#define TRANSPORT_PROBLEM_Z_HPP

class transport_problem_z {
protected:

    DM const & dmGrid_Shifted;
    DM const & dmGrid_Centered;
    DM const & dmGrid_Staggered;

    PetscReal const & dt;
    PetscReal const & iter;

    Vec U_0, V_0, W_0;
    Vec W_int;

    PetscErrorCode const FirstShiftU_y(Vec & UShifted, Vec const & vec, PetscScalar const & theta) //okok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, up_right, down_left, up_left, left, down_right;
        PetscInt icux_right[3], icux_up_left[3], icux_up_right[3], icux_down_left[3], icux_down_right[3], icux_left[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hy = 1.0 / N[1];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icux_up_left[d]);
            DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icux_up_right[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icux_down_left[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icux_down_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &right);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_RIGHT, 0, &up_right);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_LEFT, 0, &down_left);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_LEFT, 0, &up_left);
        DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &left);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_RIGHT, 0, &down_right);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec UShifted_local;
        VecSet(UShifted, 0.0);
        PetscReal ****arrUShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &UShifted_local);
        VecSet(UShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);    

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
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]] - hy,
                                    arrCoord[ez][ey][ex][icux_left[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][left];
                        PetscReal inter = (out + inner) / 2.0;
                        arrUShifted_local[ez][ey][ex][down_left] = inter;                 
                    }

                    if (ex == 0) {                  
                        PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_up_left[0]], arrCoord[ez][ey][ex][icux_up_left[1]],
                                    arrCoord[ez][ey][ex][icux_up_left[2]], theta);
                        arrUShifted_local[ez][ey][ex][up_left] = inter;                  
                    }

                    if (ey == N[1] - 1) {
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]] + hy,
                                    arrCoord[ez][ey][ex][icux_right[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][right];
                        PetscReal inter = (out + inner) / 2.0;
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
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, UShifted_local, ADD_VALUES, UShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &UShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstShiftU_z(Vec & UShifted, Vec const & vec, PetscScalar const & theta) //okok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, front_right, back_left, front_left, back_right, left;
        PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3], icux_left[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hz = 1.0 / N[2];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
            DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &left);
        DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &right);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_RIGHT, 0, &front_right);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_LEFT, 0, &back_left);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_LEFT, 0, &front_left);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_RIGHT, 0, &back_right);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec UShifted_local;
        VecSet(UShifted, 0.0);
        PetscReal ****arrUShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &UShifted_local);
        VecSet(UShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);


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
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]],
                                    arrCoord[ez][ey][ex][icux_left[2]] - hz, theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][left];
                        PetscReal inter = (out + inner) / 2.0;
                        arrUShifted_local[ez][ey][ex][back_left] = inter;
                    }

                    if (ex == 0) {
                        PetscReal inter = uxRef(arrCoord[ez][ey][ex][icux_front_left[0]], arrCoord[ez][ey][ex][icux_front_left[1]],
                                    arrCoord[ez][ey][ex][icux_front_left[2]], theta);
                        arrUShifted_local[ez][ey][ex][front_left] = inter;
                    }

                    if (ez == N[2] - 1) {
                        PetscReal out = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]],
                                    arrCoord[ez][ey][ex][icux_right[2]] + hz, theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][right];
                        PetscReal inter = (out + inner) / 2.0;                  
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
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, UShifted_local, &arrUShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, UShifted_local, ADD_VALUES, UShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &UShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstShiftV_y(Vec & VShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBeginUser;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, up, up_right, down_left, up_left, down_right, down;
        PetscInt icuy[3], icuy_up_left[3], icuy_up_right[3], icuy_down_left[3], icuy_down_right[3], icuy_up[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hx = 1.0 / N[0];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
            DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
            DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icuy_up_left[d]);
            DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icuy_up_right[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icuy_down_left[d]);
            DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icuy_down_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &down);
        DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &up);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_RIGHT, 0, &up_right);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_LEFT, 0, &down_left);
        DMStagGetLocationSlot(dmGrid_Shifted, UP_LEFT, 0, &up_left);
        DMStagGetLocationSlot(dmGrid_Shifted, DOWN_RIGHT, 0, &down_right);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec VShifted_local;
        VecSet(VShifted, 0.0);
        PetscReal ****arrVShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &VShifted_local);
        VecSet(VShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);


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
                        PetscReal out = uyRef(arrCoord[ez][ey][ex][icuy_up[0]] - hx, arrCoord[ez][ey][ex][icuy_up[1]],
                                    arrCoord[ez][ey][ex][icuy_up[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][up];
                        PetscReal inter = (out + inner) / 2.0;                    
                        arrVShifted_local[ez][ey][ex][up_left] = inter;
                    }

                    if (ey == N[1] - 1) {
                        PetscReal inter = uyRef(arrCoord[ez][ey][ex][icuy_up_right[0]], arrCoord[ez][ey][ex][icuy_up_right[1]],
                                    arrCoord[ez][ey][ex][icuy_up_right[2]], theta);
                        arrVShifted_local[ez][ey][ex][up_right] = inter;
                    }

                    if (ex == N[0] - 1) {
                        PetscReal out = uyRef(arrCoord[ez][ey][ex][icuy[0]] + hx, arrCoord[ez][ey][ex][icuy[1]],
                                    arrCoord[ez][ey][ex][icuy[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][down];
                        PetscReal inter = (out + inner) / 2.0;
                        arrVShifted_local[ez][ey][ex][down_right] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, VShifted_local, ADD_VALUES, VShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &VShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const FirstShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, front, front_right, back_left, front_left, back_right, back;
        PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3], icux_front[3], icux_back[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hx = 1.0 / N[0];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, BACK, d, &icux_back[d]);
            DMStagGetLocationSlot(dmCoord, FRONT, d, &icux_front[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
            DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &back);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &front);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_RIGHT, 0, &front_right);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_LEFT, 0, &back_left);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_LEFT, 0, &front_left);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_RIGHT, 0, &back_right);


        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec WShifted_local;
        VecSet(WShifted, 0.0);
        PetscReal ****arrWShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &WShifted_local);
        VecSet(WShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);

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
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_front[0]] - hx, arrCoord[ez][ey][ex][icux_front[1]],
                                    arrCoord[ez][ey][ex][icux_front[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][front];
                        PetscReal inter = (out + inner) / 2.0;                    
                        arrWShifted_local[ez][ey][ex][front_left] = inter;
                    }

                    if (ez == N[2] - 1) {
                        PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_front_right[0]], arrCoord[ez][ey][ex][icux_front_right[1]],
                                    arrCoord[ez][ey][ex][icux_front_right[2]], theta);
                        arrWShifted_local[ez][ey][ex][front_right] = inter;
                    }

                    if (ex == N[0] - 1) {
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_back[0]] + hx, arrCoord[ez][ey][ex][icux_back[1]],
                                    arrCoord[ez][ey][ex][icux_back[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][back];
                        PetscReal inter = (out + inner) / 2.0;
                        arrWShifted_local[ez][ey][ex][back_right] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, WShifted_local, ADD_VALUES, WShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &WShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const SecondShiftV_z(Vec & VShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBegin;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, up, front_up, back_down, front_down, back_up, down;
        PetscInt icuy[3], icuy_front_up[3], icuy_front_down[3], icuy_back_up[3], icuy_back_down[3], icuy_up[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hz = 1.0 / N[2];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
            DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icuy_front_up[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icuy_front_down[d]);
            DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icuy_back_up[d]);
            DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icuy_back_down[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &down);
        DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &up);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_UP, 0, &front_up);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_DOWN, 0, &back_down);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_DOWN, 0, &front_down);  
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_UP, 0, &back_up);

        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec VShifted_local;
        VecSet(VShifted, 0.0);
        PetscReal ****arrVShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &VShifted_local);
        VecSet(VShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);


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
                        PetscScalar out = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]],
                                    arrCoord[ez][ey][ex][icuy_up[2]] + hz, theta);
                        PetscScalar inner = arrvec_local[ez][ey][ex][up];
                        PetscScalar inter = (out + inner) / 2.0;
                        arrVShifted_local[ez][ey][ex][front_up] = inter;
                    }

                    if (ey == N[1] - 1) {
                        PetscScalar inter = uyRef(arrCoord[ez][ey][ex][icuy_back_up[0]], arrCoord[ez][ey][ex][icuy_back_up[1]],
                                    arrCoord[ez][ey][ex][icuy_back_up[2]], theta);
                        arrVShifted_local[ez][ey][ex][back_up] = inter;
                    }

                    if (ez == 0) {
                        PetscScalar out = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]],
                                    arrCoord[ez][ey][ex][icuy[2]] - hz, theta);
                        PetscScalar inner = arrvec_local[ez][ey][ex][down];
                        PetscScalar inter = (out + inner) / 2.0;
                        arrVShifted_local[ez][ey][ex][back_down] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, VShifted_local, &arrVShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, VShifted_local, ADD_VALUES, VShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &VShifted_local);
        
        PetscFunctionReturn(0);
    }

    PetscErrorCode const SecondShiftW_z(Vec & WShifted, Vec const & vec, PetscScalar const & theta) //ok
    {
        PetscFunctionBeginUser;
        Vec coordLocal;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, front, front_up, back_down, front_down, back_up, back;
        PetscInt icux_right[3], icux_back_up[3], icux_back_down[3], icux_front_up[3], icux_front_down[3], icux_front[3], icux_back[3];
        DM dmCoord;
        PetscScalar ****arrCoord;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        PetscReal hy = 1.0 / N[1];
        DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);
        DMGetCoordinatesLocal(dmGrid_Shifted, &coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, BACK, d, &icux_back[d]);
            DMStagGetLocationSlot(dmCoord, FRONT, d, &icux_front[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icux_front_up[d]);
            DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icux_front_down[d]);
            DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icux_back_up[d]);
            DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icux_back_down[d]);
        }

        DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &back);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &front);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_UP, 0, &front_up);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_DOWN, 0, &back_down);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_DOWN, 0, &front_down);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_UP, 0, &back_up);


        Vec vec_local;
        PetscReal ****arrvec_local;
        DMCreateLocalVector(dmGrid_Shifted,&vec_local);
        DMGlobalToLocalBegin(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMGlobalToLocalEnd(dmGrid_Shifted,vec,INSERT_VALUES,vec_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);

        Vec WShifted_local;
        VecSet(WShifted, 0.0);
        PetscReal ****arrWShifted_local;
        DMGetLocalVector(dmGrid_Shifted, &WShifted_local);
        VecSet(WShifted_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);

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
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_front[0]], arrCoord[ez][ey][ex][icux_front[1]] - hy,
                                    arrCoord[ez][ey][ex][icux_front[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][front];
                        PetscReal inter = (out + inner) / 2.0;
                        arrWShifted_local[ez][ey][ex][front_down] = inter;
                    }

                    if (ez == N[2] - 1) {
                        PetscReal inter = uzRef(arrCoord[ez][ey][ex][icux_front_up[0]], arrCoord[ez][ey][ex][icux_front_up[1]],
                                    arrCoord[ez][ey][ex][icux_front_up[2]], theta);
                        arrWShifted_local[ez][ey][ex][front_up] = inter;
                    }

                    if (ey == N[1] - 1) {
                        PetscReal out = uzRef(arrCoord[ez][ey][ex][icux_back[0]], arrCoord[ez][ey][ex][icux_back[1]] + hy,
                                    arrCoord[ez][ey][ex][icux_back[2]], theta);
                        PetscReal inner = arrvec_local[ez][ey][ex][back];
                        PetscReal inter = (out + inner) / 2.0;
                        arrWShifted_local[ez][ey][ex][back_up] = inter;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid_Shifted, vec_local, &arrvec_local);
        DMStagVecRestoreArray(dmGrid_Shifted, WShifted_local, &arrWShifted_local);
        DMLocalToGlobal(dmGrid_Shifted, WShifted_local, ADD_VALUES, WShifted);
        DMRestoreLocalVector(dmGrid_Shifted, &vec_local);
        DMRestoreLocalVector(dmGrid_Shifted, &WShifted_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const ThirdDerive_x(Vec & AB_x, Vec const & AB) //ok
    {
        PetscFunctionBegin;
        PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, back_left, back_right, back, front_left, front_right, front;
        PetscReal       hx;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        hx = 1.0 / N[0];
        Vec AB_local;
        PetscScalar ****arrAB_local;
        DMCreateLocalVector(dmGrid_Shifted, &AB_local);
        DMGlobalToLocalBegin(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMGlobalToLocalEnd(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);

        Vec AB_x_local;
        VecSet(AB_x, 0.0);
        PetscScalar ****arrAB_x_local;
        DMGetLocalVector(dmGrid_Shifted, &AB_x_local);
        VecSet(AB_x_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, AB_x_local, &arrAB_x_local);

        DMStagGetLocationSlot(dmGrid_Shifted, BACK_LEFT, 0, &back_left);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_RIGHT, 0, &back_right);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &back);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_LEFT, 0, &front_left);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_RIGHT, 0, &front_right);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &front);

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

        DMStagVecRestoreArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);
        DMStagVecRestoreArray(dmGrid_Shifted, AB_x_local, &arrAB_x_local);
        DMLocalToGlobal(dmGrid_Shifted, AB_x_local, ADD_VALUES, AB_x);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_x_local);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const ThirdDerive_y(Vec & AB_y, Vec const & AB) //ok
    {
        PetscFunctionBegin;
        PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, back_down, back_up, back, front_down, front_up, front;
        PetscReal       hy;
        DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
        hy = 1.0 / N[1];
        Vec AB_local;
        PetscScalar ****arrAB_local;
        DMCreateLocalVector(dmGrid_Shifted, &AB_local);
        DMGlobalToLocalBegin(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMGlobalToLocalEnd(dmGrid_Shifted, AB, INSERT_VALUES, AB_local);
        DMStagVecGetArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);

        Vec AB_y_local;
        VecSet(AB_y, 0.0);
        PetscScalar ****arrAB_y_local;
        DMGetLocalVector(dmGrid_Shifted, &AB_y_local);
        VecSet(AB_y_local, 0.0);
        DMStagVecGetArray(dmGrid_Shifted, AB_y_local, &arrAB_y_local);

        DMStagGetLocationSlot(dmGrid_Shifted, BACK_DOWN, 0, &back_down);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK_UP, 0, &back_up);
        DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &back);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_DOWN, 0, &front_down);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT_UP, 0, &front_up);
        DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &front);

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

        DMStagVecRestoreArrayRead(dmGrid_Shifted, AB_local, &arrAB_local);
        DMStagVecRestoreArray(dmGrid_Shifted, AB_y_local, &arrAB_y_local);
        DMLocalToGlobal(dmGrid_Shifted, AB_y_local, ADD_VALUES, AB_y);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_y_local);
        DMRestoreLocalVector(dmGrid_Shifted, &AB_local);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const CenterW(Vec & WCenter, Vec const & vec, PetscReal const & theta) 
    {
        PetscFunctionBegin;
        PetscInt icuz_back[3], icuz_front[3], iuz_back, iuz_front, iuz_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
        DM dmCoord;
        Vec vecLocal, vecWLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec, ****arrW;    
        DMStagGetCorners(dmGrid_Centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Centered, &N[0], &N[1], &N[2]);
        DMGetCoordinateDM(dmGrid_Centered, &dmCoord);
        DMGetCoordinates(dmGrid_Centered, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, BACK, d, &icuz_back[d]);
            DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
        }  
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid_Centered, BACK, 0, &iuz_back);
        DMStagGetLocationSlot(dmGrid_Centered, FRONT, 0, &iuz_front); 
        DMStagGetLocationSlot(dmGrid_Centered, ELEMENT, 0, &iuz_element);
        
        DMCreateLocalVector(dmGrid_Centered, &vecLocal);
        DMGlobalToLocalBegin(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid_Centered, vecLocal, &arrVec);

        DMGetLocalVector(dmGrid_Centered, &vecWLocal);
        DMStagVecGetArray(dmGrid_Centered, vecWLocal, &arrW); 

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
        DMStagVecRestoreArrayRead(dmGrid_Centered, vecLocal, &arrVec);
        DMStagVecRestoreArray(dmGrid_Centered, vecWLocal, &arrW);
        DMLocalToGlobal(dmGrid_Centered, vecWLocal, INSERT_VALUES, WCenter);
        DMRestoreLocalVector(dmGrid_Centered, &vecWLocal);
        DMRestoreLocalVector(dmGrid_Centered, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

    PetscErrorCode const Derive_z(Vec & W2_z, Vec const & vec, PetscReal const & theta)
    {
        PetscFunctionBegin;
        PetscInt icuz_element[3], iuz_back, iuz_front, iuz_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
        DM dmCoord;
        Vec vecLocal, vecWLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec, ****arrW;    
        DMStagGetCorners(dmGrid_Centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid_Centered, &N[0], &N[1], &N[2]);
        PetscReal const hz = 1.0 / N[2];
        DMGetCoordinateDM(dmGrid_Centered, &dmCoord);
        DMGetCoordinates(dmGrid_Centered, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icuz_element[d]);
        }  
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid_Centered, BACK, 0, &iuz_back);
        DMStagGetLocationSlot(dmGrid_Centered, FRONT, 0, &iuz_front); 
        DMStagGetLocationSlot(dmGrid_Centered, ELEMENT, 0, &iuz_element);
        
        DMCreateLocalVector(dmGrid_Centered, &vecLocal);
        DMGlobalToLocalBegin(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid_Centered, vec, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid_Centered, vecLocal, &arrVec);

        DMGetLocalVector(dmGrid_Centered, &vecWLocal);
        DMStagVecGetArray(dmGrid_Centered, vecWLocal, &arrW); 

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
        DMStagVecRestoreArrayRead(dmGrid_Centered, vecLocal, &arrVec);
        DMStagVecRestoreArray(dmGrid_Centered, vecWLocal, &arrW);
        DMLocalToGlobal(dmGrid_Centered, vecWLocal, INSERT_VALUES, W2_z);
        DMRestoreLocalVector(dmGrid_Centered, &vecWLocal);
        DMRestoreLocalVector(dmGrid_Centered, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

public:

    transport_problem_z(ProblemSetting<Transport> const & setting, Vec const & U_0_, Vec const & V_0_, Vec const & W_0_) : 
    dmGrid_Shifted(setting.dmGrid_Shifted), dmGrid_Staggered(setting.dmGrid_Staggered), dmGrid_Centered(setting.dmGrid_Centered), dt(setting.dt),
    iter(setting.iter)
    {
        DMCreateGlobalVector(dmGrid_Staggered, &U_0);
        DMCreateGlobalVector(dmGrid_Staggered, &V_0);
        DMCreateGlobalVector(dmGrid_Staggered, &W_0);
        DMCreateGlobalVector(dmGrid_Staggered, &W_int);
        VecCopy(U_0_, U_0);
        VecCopy(V_0_, V_0);
        VecCopy(W_0_, W_0);
    }

    PetscErrorCode const solve_step_z(PetscScalar const & theta)
    {

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

        Vec WU_x, WV_y;
        DMCreateGlobalVector(dmGrid_Shifted, &WU_x);
        DMCreateGlobalVector(dmGrid_Shifted, &WV_y);
        Vec mixedThird;
        DMCreateGlobalVector(dmGrid_Staggered, &mixedThird);
        {
            Vec W_x, W_y, U_x, V_y, WU, WV;
            DMCreateGlobalVector(dmGrid_Shifted, &W_x);
            DMCreateGlobalVector(dmGrid_Shifted, &W_y);
            DMCreateGlobalVector(dmGrid_Shifted, &U_x);
            DMCreateGlobalVector(dmGrid_Shifted, &V_y);
            DMCreateGlobalVector(dmGrid_Shifted, &WU);
            DMCreateGlobalVector(dmGrid_Shifted, &WV);
            FirstShiftW_z(W_x, W_shift, theta);// ==FirsShiftW_z
            SecondShiftW_z(W_y, W_shift, theta);// ==SecondShiftW_z
            FirstShiftU_z(U_x, U_shift, theta);// ==FirsShiftU_z
            SecondShiftV_z(V_y, V_shift, theta);// ==SecondShiftV_z
            VecPointwiseMult(WU, W_x, U_x);
            VecPointwiseMult(WV, W_y, V_y);
            ThirdDerive_x(WU_x, WU);
            ThirdDerive_y(WV_y, WV);
            VecAXPY(WU_x, 1.0, WV_y);
            PetscObjectDestroy((PetscObject*)&W_x);
            PetscObjectDestroy((PetscObject*)&W_y);
            PetscObjectDestroy((PetscObject*)&U_x);
            PetscObjectDestroy((PetscObject*)&V_y);
            PetscObjectDestroy((PetscObject*)&WU);
            PetscObjectDestroy((PetscObject*)&WV);
        }

        DMStagMigrateVec(dmGrid_Shifted, WU_x, dmGrid_Staggered, mixedThird);

        Vec W_center;
        DMCreateGlobalVector(dmGrid_Centered, &W_center);
        DMStagMigrateVec(dmGrid_Staggered, W_n, dmGrid_Centered, W_center);

        Vec W_c;
        DMCreateGlobalVector(dmGrid_Centered, &W_c);
        CenterW(W_c, W_center, theta);

        Vec W2_z;
        DMCreateGlobalVector(dmGrid_Centered, &W2_z);
        Derive_z(W2_z, W_c, theta);

        Vec homoThird;
        DMCreateGlobalVector(dmGrid_Staggered, &homoThird);
        DMStagMigrateVec(dmGrid_Centered, W2_z, dmGrid_Staggered, homoThird);

        VecAXPBYPCZ(W_n, -dt, -dt, 1.0, homoThird, mixedThird);
        VecCopy(W_n, W_int);
        //VecCopy(W_n, W_0);

        PetscObjectDestroy((PetscObject*)&U_shift);
        PetscObjectDestroy((PetscObject*)&V_shift);
        PetscObjectDestroy((PetscObject*)&W_shift);

        PetscObjectDestroy((PetscObject*)&WU_x);
        PetscObjectDestroy((PetscObject*)&WV_y);
        PetscObjectDestroy((PetscObject*)&mixedThird);

        PetscObjectDestroy((PetscObject*)&W_c);
        PetscObjectDestroy((PetscObject*)&W2_z);
        PetscObjectDestroy((PetscObject*)&W_center);
        PetscObjectDestroy((PetscObject*)&homoThird);

        VecDestroy(&U_n);
        VecDestroy(&V_n);
        VecDestroy(&W_n);

        PetscFunctionReturn(0);
    }

    const Vec & get_W() const { return W_int; }

    ~transport_problem_z()
    {
        VecDestroy(&U_0);
        VecDestroy(&V_0);
        VecDestroy(&W_0);
        VecDestroy(&W_int);
    }

};

#endif // TRANSPORT_PROBLEM_Z_HPP
