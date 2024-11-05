//
// Created by dave on 26/10/23.
//

/*WHAT TO BE DONE
1. Indentare tutto correttamente
2. Eliminare tutti i warning e le variabili inutili indicate dai warning
3. Ragionare se alcune copie di vettori siano necessarie (potebbe essere molto dispendioso)
4. Trovare un benchmark che vada bene perche' Taylor Green non rispetta ne Neumman ne Dirichlet. A sto punto forse l'esponenziale per il seno e' la cosa migliore
5. Capire come integrare questa roba qua in Lifex (ad esempio il Parsing dei parametri e i sistemi lineari vengono gestiti da Lifex di default)
6. Capire come integrare la vena (lol)
7. Trovare uno stencil migliore per le condizioni di Neumman

8. non sono sicuro sia giusto il non lineare di V

NOTA BENE: Il programma prende come input U0, V0, W0 da uxRef, uyRef, UzRef e lo stesso per le condizioni al contorno (di fatto assegnate nel termine convettivo)

*/
static char help[] = "Solve a toy 3D problem on a staggered grid\n\n";

#include <petscdm.h>
#include <petscksp.h>
#include <petscdmstag.h>
#include <iostream>
#include <chrono>
#include "petscvec.h" 




#define pi 3.1415926535



#define BACK_DOWN_LEFT   DMSTAG_BACK_DOWN_LEFT
#define BACK_DOWN        DMSTAG_BACK_DOWN
#define BACK_DOWN_RIGHT  DMSTAG_BACK_DOWN_RIGHT
#define BACK_LEFT        DMSTAG_BACK_LEFT
#define BACK             DMSTAG_BACK
#define BACK_RIGHT       DMSTAG_BACK_RIGHT
#define BACK_UP_LEFT     DMSTAG_BACK_UP_LEFT
#define BACK_UP          DMSTAG_BACK_UP
#define BACK_UP_RIGHT    DMSTAG_BACK_UP_RIGHT
#define DOWN_LEFT        DMSTAG_DOWN_LEFT
#define DOWN             DMSTAG_DOWN
#define DOWN_RIGHT       DMSTAG_DOWN_RIGHT
#define LEFT             DMSTAG_LEFT
#define ELEMENT          DMSTAG_ELEMENT
#define RIGHT            DMSTAG_RIGHT
#define UP_LEFT          DMSTAG_UP_LEFT
#define UP               DMSTAG_UP
#define UP_RIGHT         DMSTAG_UP_RIGHT
#define FRONT_DOWN_LEFT  DMSTAG_FRONT_DOWN_LEFT
#define FRONT_DOWN       DMSTAG_FRONT_DOWN
#define FRONT_DOWN_RIGHT DMSTAG_FRONT_DOWN_RIGHT
#define FRONT_LEFT       DMSTAG_FRONT_LEFT
#define FRONT            DMSTAG_FRONT
#define FRONT_RIGHT      DMSTAG_FRONT_RIGHT
#define FRONT_UP_LEFT    DMSTAG_FRONT_UP_LEFT
#define FRONT_UP         DMSTAG_FRONT_UP
#define FRONT_UP_RIGHT   DMSTAG_FRONT_UP_RIGHT

static PetscScalar a = pi/4;
static PetscScalar d = pi/2;
static PetscScalar par = 1;
//static PetscScalar A = 4*sqrt(2)/(3*sqrt(3));
/*static PetscScalar k = 2*pi;
static PetscScalar c = (5/6)*pi;
static PetscScalar d = (1/6)*pi;
static PetscScalar A = 1;
static PetscScalar B = 1;
static PetscScalar C = 1;*/


static PetscScalar uxRef(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar theta)

{
    //return  cos(pi*x)*sin(2*pi*y);
    //return  -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y))*exp(-theta)*par;
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    //return (sin(k*x - c)*cos(k*y - d)*sin(k*z) - cos(k*z - c)*sin(k*x - d)*sin(k*y))*exp(-theta)*1e-4;
    //return (A*sin(k*z) + C*cos(k*y))*exp(-theta);
    return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z;

}

static PetscScalar uyRef(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar theta)
{
    //return  -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z))*exp(-theta)*par;
    //return (sin(k*y - c)*cos(k*z - d)*sin(k*x) - cos(k*x - c)*sin(k*y - d)*sin(k*z))*exp(-theta)*1e-4;
    return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return -cos(2*pi*x)*sin(2*pi*y)*exp(-theta);
    //return (B*sin(k*x) + A*cos(k*z))*exp(-theta);
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);


}

static PetscScalar uzRef(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar theta)
{
    return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return (sin(k*z - c)*cos(k*x - d)*sin(k*y) - cos(k*y - c)*sin(k*z - d)*sin(k*x))*exp(-theta)*1e-4;
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    //return 0.0*x + 0.0*y + 0.0*z;
    //return (C*sin(k*y) + B*cos(k*x))*exp(-theta);
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);

}


static PetscScalar solution(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar theta)
{
    //return (sin(k*x - c)*cos(k*y - d)*sin(k*z) - cos(k*z - c)*sin(k*x - d)*sin(k*y))*exp(-theta)*1e-4;
    //return -84*pi*pi*sin(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    //return 0.0*x + 0.0*y + 0.0*z;
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z;
    //return (A*sin(k*z) + C*cos(k*y))*exp(-theta);
    return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) + 0.01*4*pi*cos(2*pi*x)*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)*cos(2*pi*z)*sin(2*pi*y) + 0.01*4*pi*cos(2*pi*x)*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*y)*cos(2*pi*z)*sin(2*pi*z) + 0.01*4*pi*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*y)*cos(2*pi*z)*cos(2*pi*z)*sin(2*pi*x);

}



static PetscErrorCode CreateReferenceSolutionFirst(DM dmGrid, Vec *pSolRef, PetscScalar theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmGrid, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmGrid, &solRefLocal);
    DMStagVecGetArray(dmGrid, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[1] + n[1] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iux] = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, solRefLocal, &arrSol);
    DMLocalToGlobal(dmGrid, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &solRefLocal);

    return 0;
}

static PetscErrorCode CreateReferenceSolutionSecond(DM dmGrid, Vec *pSolRef, PetscScalar theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmGrid, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmGrid, &solRefLocal);
    DMStagVecGetArray(dmGrid, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[0] + n[0] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iuy] = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, solRefLocal, &arrSol);
    DMLocalToGlobal(dmGrid, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &solRefLocal);

    return 0;
}

static PetscErrorCode CreateReferenceSolutionThird(DM dmGrid, Vec *pSolRef, PetscScalar theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmGrid, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmGrid, &solRefLocal);
    DMStagVecGetArray(dmGrid, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[0] + n[0] && ey < start[1] + n[1]) 
                arrSol[ez][ey][ex][iuz] = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, solRefLocal, &arrSol);
    DMLocalToGlobal(dmGrid, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &solRefLocal);

    return 0;
}

static PetscErrorCode CreateReferenceSolutionTry(DM dmGrid, Vec *pSolRef, PetscScalar theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3], iue, icue[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmGrid, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iue);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icue[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmGrid, &solRefLocal);
    DMStagVecGetArray(dmGrid, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[1] + n[1] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iuy] = solution(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, solRefLocal, &arrSol);
    DMLocalToGlobal(dmGrid, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &solRefLocal);

    return 0;
}

// Create Domain routines
static PetscErrorCode CreateGrid(DM* pdmGrid, PetscInt dof1, PetscInt dof2, PetscInt dof3, PetscScalar nx, PetscScalar ny, PetscScalar nz, PetscReal Lx_0, PetscReal Lx, PetscReal Ly_0, PetscReal Ly, PetscReal Lz_0, PetscReal Lz)
{
    DM dmGrid;
    const PetscInt dof0 = 0;
    const PetscInt stencilWidth = 1;
    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, nx, ny, nz, PETSC_DECIDE,
                   PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL,
                   NULL, pdmGrid);
    dmGrid = *pdmGrid;
    DMSetFromOptions(dmGrid);
    DMSetUp(dmGrid);
    DMStagSetUniformCoordinatesExplicit(dmGrid, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    //PetscObjectDestroy((PetscObject*)&dmGrid);

    return 0;
}

PetscErrorCode FirstShiftU_y(DM const & dmGrid, Vec & UShifted, Vec const & vec, PetscScalar const & theta) //okok
{

    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, up_right, down_left, up_left, left, down_right;
    PetscInt icux_right[3], icux_up_left[3], icux_up_right[3], icux_down_left[3], icux_down_right[3], icux_left[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pUShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal hy = 1.0 / N[1];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
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

    DMStagVecRestoreArrayRead(dmGrid, vec_local, &arrvec_local);
    DMStagVecRestoreArray(dmGrid, UShifted_local, &arrUShifted_local);
    DMLocalToGlobal(dmGrid, UShifted_local, ADD_VALUES, UShifted);
    DMRestoreLocalVector(dmGrid, &vec_local);
    DMRestoreLocalVector(dmGrid, &UShifted_local);

    return 0;
}

PetscErrorCode FirstShiftU_z(DM const & dmGrid, Vec & UShifted, Vec const & vec, PetscScalar const & theta) //okok
{
    Vec coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, right, front_right, back_left, front_left, back_right, left;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3], icux_left[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pUShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal hz = 1.0 / N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
        DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
    }

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &left);
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

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, up, up_right, down_left, up_left, down_right, down;
    PetscInt icuy[3], icuy_up_left[3], icuy_up_right[3], icuy_down_left[3], icuy_down_right[3], icuy_up[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pVShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal hx = 1.0 / N[0];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icuy_up_left[d]);
        DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icuy_up_right[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icuy_down_left[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icuy_down_right[d]);
    }

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &down);
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

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, front, front_right, back_left, front_left, back_right, back;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3], icux_front[3], icux_back[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pWShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal hx = 1.0 / N[0];
    
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
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

    DMStagGetLocationSlot(dmGrid, BACK, 0, &back);
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

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, up, front_up, back_down, front_down, back_up, down;
    PetscInt icuy[3], icuy_front_up[3], icuy_front_down[3], icuy_back_up[3], icuy_back_down[3], icuy_up[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pVShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal hz = 1.0 / N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icuy_front_up[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icuy_front_down[d]);
        DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icuy_back_up[d]);
        DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icuy_back_down[d]);
    }

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &down);
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

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d, front, front_up, back_down, front_down, back_up, back;
    PetscInt icux_right[3], icux_back_up[3], icux_back_down[3], icux_front_up[3], icux_front_down[3], icux_front[3], icux_back[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmGrid, pWShifted);

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal hy = 1.0 / N[1];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
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

    DMStagGetLocationSlot(dmGrid, BACK, 0, &back);
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
static PetscErrorCode ManageAdvection_x(PetscScalar dt, Vec & U_int, Vec const & U_n, Vec const & V_n, Vec const & W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz, PetscScalar theta)
{

    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered;
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

        FirstShiftU_y(dmGrid_Shifted, U_y, U_shift, theta);
        FirstShiftU_z(dmGrid_Shifted, U_z, U_shift, theta);
        FirstShiftV_y(dmGrid_Shifted, V_y, V_shift, theta);
        FirstShiftW_z(dmGrid_Shifted, W_z, W_shift, theta); 

        VecPointwiseMult(UV, U_y, V_y);
        VecPointwiseMult(UW, U_z, W_z);

        FirstDerive_y(dmGrid_Shifted, UV_y, UV);
        FirstDerive_z(dmGrid_Shifted, UW_z, UW);

        VecAXPY(UV_y, 1.0, UW_z);

        PetscObjectDestroy((PetscObject*)&U_y);
        PetscObjectDestroy((PetscObject*)&U_z);
        PetscObjectDestroy((PetscObject*)&V_y);
        PetscObjectDestroy((PetscObject*)&W_z);
        PetscObjectDestroy((PetscObject*)&UV);
        PetscObjectDestroy((PetscObject*)&UW);

    }

    DMStagMigrateVec(dmGrid_Shifted, UV_y, dmGrid_Staggered, mixedFirst);


    // Managing homogenous components of non-linear term
    //////////////////////////////////////////////////////////////
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

 
  

    // Copy to output

    VecCopy(U_n, U_int);


    PetscObjectDestroy((PetscObject*)&dmGrid_Shifted);
    PetscObjectDestroy((PetscObject*)&dmGrid_Centered);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered);
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
    
    return 0;
}

static PetscErrorCode ManageAdvection_y(PetscScalar dt, Vec & V_int, Vec const & U_n, Vec const & V_n, Vec const & W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz, PetscScalar theta)
{

    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered;
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


        FirstShiftV_y(dmGrid_Shifted, V_x, V_shift, theta);
        SecondShiftV_z(dmGrid_Shifted, V_z, V_shift, theta);

        FirstShiftU_y(dmGrid_Shifted, U_x, U_shift, theta);
        SecondShiftW_z(dmGrid_Shifted, W_z, W_shift, theta);


        VecPointwiseMult(VU, V_x, U_x);
        VecPointwiseMult(VW, V_z, W_z);



        SecondDerive_x(dmGrid_Shifted, VU_x, VU);
        SecondDerive_z(dmGrid_Shifted, VW_z, VW);

        VecAXPY(VU_x, 1.0, VW_z);

        PetscObjectDestroy((PetscObject*)&V_x);
        PetscObjectDestroy((PetscObject*)&V_z);
        PetscObjectDestroy((PetscObject*)&U_x);
        PetscObjectDestroy((PetscObject*)&W_z);
        PetscObjectDestroy((PetscObject*)&VU);
        PetscObjectDestroy((PetscObject*)&VW);
    }
    DMStagMigrateVec(dmGrid_Shifted, VU_x, dmGrid_Staggered, mixedSecond);

    ///////////////////////////////

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



    /*VecAXPBY(homoFirst, -1.0, -1.0, mixedFirst);

    VecScale(homoFirst, dt);

    VecAXPY(homoFirst, 1.0, U_n);*/    

    // Copy to output

    VecCopy(V_n, V_int);


    PetscObjectDestroy((PetscObject*)&dmGrid_Shifted);
    PetscObjectDestroy((PetscObject*)&dmGrid_Centered);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered);
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

    
    return 0;
}

static PetscErrorCode ManageAdvection_z(PetscScalar dt, Vec & W_int, Vec const & U_n, Vec const & V_n, Vec const & W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz, PetscScalar theta)
{

    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered;
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

        FirstShiftW_z(dmGrid_Shifted, W_x, W_shift, theta);// ==FirsShiftW_z
        SecondShiftW_z(dmGrid_Shifted, W_y, W_shift, theta);// ==SecondShiftW_z
        FirstShiftU_z(dmGrid_Shifted, U_x, U_shift, theta);// ==FirsShiftU_z
        SecondShiftV_z(dmGrid_Shifted, V_y, V_shift, theta);// ==SecondShiftV_z

        VecPointwiseMult(WU, W_x, U_x);
        VecPointwiseMult(WV, W_y, V_y);

        ThirdDerive_x(dmGrid_Shifted, WU_x, WU);
        ThirdDerive_y(dmGrid_Shifted, WV_y, WV);

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
    CenterW(dmGrid_Centered, W_c, W_center, theta);

    Vec W2_z;
    DMCreateGlobalVector(dmGrid_Centered, &W2_z);
    Derive_z(dmGrid_Centered, W2_z, W_c, theta);

    Vec homoThird;
    DMCreateGlobalVector(dmGrid_Staggered, &homoThird);
    DMStagMigrateVec(dmGrid_Centered, W2_z, dmGrid_Staggered, homoThird);



    // Assembling whole non-linear
    VecAXPBYPCZ(W_n, -dt, -dt, 1.0, homoThird, mixedThird);

    VecCopy(W_n, W_int);

    PetscObjectDestroy((PetscObject*)&dmGrid_Shifted);
    PetscObjectDestroy((PetscObject*)&dmGrid_Centered);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered);
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
    
    return 0;
}

PetscErrorCode CheckSolution(Vec const & sol, Vec const & solRef)
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
    PetscPrintf(PETSC_COMM_WORLD, "Error (abs): %g\nError (rel): %g\n", (double)errAbs, (double)errRel);
    VecDestroy(&diff);
    PetscFunctionReturn(0);
}




int main(int argc, char **argv)
{
    auto start = std::chrono::high_resolution_clock::now();

    PetscInitialize(&argc, &argv, (char*)0, "Help text for this program");
    // Get the rank of the current process

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    //std::cout<<rank<<std::endl;

    // Get the size of the communicator (total number of processes)
    int size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    //std::cout<<size<<std::endl;

    static PetscScalar T   = 1.0;
    static PetscInt nt  = 320;
    static PetscInt nx  = 40;
    static PetscInt ny  = 40;
    static PetscInt nz  = 40;

    static PetscReal Lx_0  = 0;
    static PetscReal Ly_0  = 0;
    static PetscReal Lz_0  = 0;
    static PetscReal Lx = 1;
    static PetscReal Ly = 1;
    static PetscReal Lz = 1;
    static PetscReal Re = 10;


    PetscScalar dt = 0.01; //T / nt;
    std::cout << "dt: " << dt << std::endl;

    PetscFunctionBeginUser;
    PetscInitialize(&argc, &argv, (char *)0, help);

    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z; //need to declare to due to solving linear system laplacian messing things up
    {
        CreateGrid(&dmGrid_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered_x, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMClone(dmGrid_Staggered_x, &dmGrid_Staggered_y);
        DMClone(dmGrid_Staggered_x, &dmGrid_Staggered_z);
    }



    //PetscObjectDestroy((PetscObject*)&bench);

    Vec U_0, V_0, W_0;
    CreateReferenceSolutionFirst(dmGrid_Staggered_x, &U_0, 0);
    CreateReferenceSolutionSecond(dmGrid_Staggered_y, &V_0, 0);
    CreateReferenceSolutionThird(dmGrid_Staggered_z, &W_0, 0);
    Vec U_int, V_int, W_int;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_int);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_int);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_int);

    //ManageAdvection_x(dt, U_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, 0);
    ManageAdvection_y(dt, V_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, 0);
    //ManageAdvection_z(dt, W_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, 0);

    Vec bench_x;
    DMCreateGlobalVector(dmGrid_Staggered_y, &bench_x);
    CreateReferenceSolutionTry(dmGrid_Staggered_y, &bench_x, 0);
    CheckSolution(V_int, bench_x);
            
    VecDestroy(&U_0);
    VecDestroy(&V_0);
    VecDestroy(&W_0);
    VecDestroy(&U_int);
    VecDestroy(&V_int);
    VecDestroy(&W_int);
    

    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_x);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_y);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_z);
    PetscObjectDestroy((PetscObject*)&dmGrid_Centered);
    PetscObjectDestroy((PetscObject*)&dmGrid_Shifted);    

    PetscFinalize();
        // Get the ending timestamp
    // Get the ending timestamp
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration in seconds
    std::chrono::duration<double> duration = end - start;

    std::cout << "Trinity test successfully completed. Ad maiora!" << std::endl;
    // Output the duration
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;


    return 0;
}




