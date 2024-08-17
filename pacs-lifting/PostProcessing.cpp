
#include <petscdm.h>

#define ELEMENT          DMSTAG_ELEMENT
#define RIGHT            DMSTAG_RIGHT
#define LEFT             DMSTAG_LEFT
#define UP               DMSTAG_UP
#define DOWN             DMSTAG_DOWN
#define FRONT            DMSTAG_FRONT
#define BACK             DMSTAG_BACK

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
