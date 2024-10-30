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


static PetscScalar uxRef(PetscScalar x, PetscScalar y, PetscScalar z)

{
    //return  cos(pi*x)*sin(2*pi*y);
    //return  -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y));
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    return sin(2*pi*x)*cos(2*pi*y);
}

static PetscScalar uyRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return  -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z));
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    return -cos(2*pi*x)*sin(2*pi*y);

}

static PetscScalar uzRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return  -a*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x));
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    return 0.0*x + 0.0*y + 0.0*z;
}

static PetscScalar solution(PetscScalar x, PetscScalar y, PetscScalar z)
{
    return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z) + 0.003*4*pi*cos(2*pi*x)*cos(4*pi*y)*cos(4*pi*y)*cos(8*pi*z)*cos(8*pi*z)*sin(2*pi*x) + 0.003*8*pi*cos(2*pi*x)*cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z)*cos(8*pi*z)*sin(4*pi*y) + 0.003*16*pi*cos(2*pi*x)*cos(2*pi*x)*cos(4*pi*y)*cos(4*pi*y)*cos(8*pi*z)*sin(8*pi*z);
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
}

static PetscErrorCode CreateReferenceSolutionFirst(DM dmSol, Vec *pSolRef)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmSol, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinates(dmSol, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmSol, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmSol, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmSol, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmSol, &solRefLocal);
    DMStagVecGetArray(dmSol, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[1] + n[1] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iux] = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmSol, solRefLocal, &arrSol);
    DMLocalToGlobal(dmSol, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmSol, &solRefLocal);

    return 0;
}

static PetscErrorCode CreateReferenceSolutionSecond(DM dmSol, Vec *pSolRef)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmSol, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinates(dmSol, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmSol, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmSol, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmSol, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmSol, &solRefLocal);
    DMStagVecGetArray(dmSol, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[0] + n[0] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iuy] = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmSol, solRefLocal, &arrSol);
    DMLocalToGlobal(dmSol, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmSol, &solRefLocal);

    return 0;
}

static PetscErrorCode CreateReferenceSolutionThird(DM dmSol, Vec *pSolRef)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmSol, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinates(dmSol, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmSol, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmSol, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmSol, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmSol, &solRefLocal);
    DMStagVecGetArray(dmSol, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[0] + n[0] && ey < start[1] + n[1]) 
                arrSol[ez][ey][ex][iuz] = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmSol, solRefLocal, &arrSol);
    DMLocalToGlobal(dmSol, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmSol, &solRefLocal);

    return 0;
}

static PetscErrorCode CreateReferenceSolutionTry(DM dmSol, Vec *pSolRef)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmSol, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinates(dmSol, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmSol, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmSol, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmSol, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmSol, &solRefLocal);
    DMStagVecGetArray(dmSol, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[1] + n[1] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iuz] = solution(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmSol, solRefLocal, &arrSol);
    DMLocalToGlobal(dmSol, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmSol, &solRefLocal);

    return 0;
}


// Create Domain routines
static PetscErrorCode CreateGrid(DM* pdmSol, PetscInt dof1, PetscInt dof2, PetscInt dof3, PetscScalar nx, PetscScalar ny, PetscScalar nz, PetscReal Lx_0, PetscReal Lx, PetscReal Ly_0, PetscReal Ly, PetscReal Lz_0, PetscReal Lz)
{
    DM dmSol;
    const PetscInt dof0 = 0;
    const PetscInt stencilWidth = 1;
    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, nx, ny, nz, PETSC_DECIDE,
                   PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL,
                   NULL, pdmSol);
    dmSol = *pdmSol;
    DMSetFromOptions(dmSol);
    DMSetUp(dmSol);
    DMStagSetUniformCoordinatesExplicit(dmSol, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    //PetscObjectDestroy((PetscObject*)&dmSol);

    return 0;
}

// First non-linear members: mixed
static PetscErrorCode FirstShiftU_y(DM dmSol, Vec *pUShifted, Vec solRef) //ok
{

    Vec UShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_up_left[3], icux_up_right[3], icux_down_left[3], icux_down_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pUShifted);
    UShifted = *pUShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icux_up_left[d]);
        DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icux_up_right[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icux_down_left[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icux_down_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

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

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);


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
                                  arrCoord[ez][ey][ex][icux_down_left[2]]);
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_up_left[2]]);
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_up_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_down_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);


                }
            }
        }
    }



    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    VecAssemblyBegin(UShifted);
    VecAssemblyEnd(UShifted);


    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);

    return 0;
}

static PetscErrorCode FirstShiftU_z(DM dmSol, Vec *pUShifted, Vec solRef) //ok
{

    Vec UShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pUShifted);
    UShifted = *pUShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
        DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

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

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_back_left[2]]);
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_front_left[2]]);
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_front_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_back_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, UShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(UShifted);
    VecAssemblyEnd(UShifted);
    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);

    PetscObjectDestroy((PetscObject*)&l);

    return 0;
}

static PetscErrorCode FirstShiftV_y(DM dmSol, Vec *pVShifted, Vec solRef) //ok
{

    Vec VShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up_left[3], icuy_up_right[3], icuy_down_left[3], icuy_down_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pVShifted);
    VShifted = *pVShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icuy_up_left[d]);
        DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icuy_up_right[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icuy_down_left[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icuy_down_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol, &l);
    DMGlobalToLocalBegin(dmSol, solRef, INSERT_VALUES, l);


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

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icuy_down_left[2]]);
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icuy_up_left[2]]);
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icuy_up_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icuy_down_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VShifted);
    VecAssemblyEnd(VShifted);
    DMGlobalToLocalEnd(dmSol, solRef, INSERT_VALUES, l);

    PetscObjectDestroy((PetscObject*)&l);



    return 0;
}

static PetscErrorCode FirstShiftW_z(DM dmSol, Vec *pWShifted, Vec solRef) //ok
{

    Vec WShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pWShifted);
    WShifted = *pWShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
        DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, &icux_back_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, &icux_front_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

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

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_back_left[2]]);
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_front_left[2]]);
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_front_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_back_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WShifted);
    VecAssemblyEnd(WShifted);
    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);

    PetscObjectDestroy((PetscObject*)&l);


    return 0;
}

static PetscErrorCode FirstDerive_y(DM dmSol, Vec *pAB_y, Vec AB) //ok
{
    Vec             AB_y, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hy;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pAB_y);
    AB_y = *pAB_y;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hy = 1.0 / N[1];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

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
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);


    PetscObjectDestroy((PetscObject*)&AB_local);


    return 0;
}

static PetscErrorCode FirstDerive_z(DM dmSol, Vec *pAB_z, Vec AB) //ok
{
    Vec             AB_z, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hz;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pAB_z);
    AB_z = *pAB_z;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hz = 1.0 / N[2];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

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

    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);

    PetscObjectDestroy((PetscObject*)&AB_local);


    return 0;
}

// Second non-linear members: mixed
static PetscErrorCode SecondShiftV_z(DM dmSol, Vec *pVShifted, Vec solRef) //ok
{

    Vec VShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_front_up[3], icuy_front_down[3], icuy_back_up[3], icuy_back_down[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pVShifted);
    VShifted = *pVShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icuy_front_up[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icuy_front_down[d]);
        DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icuy_back_up[d]);
        DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icuy_back_down[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol, &l);
    DMGlobalToLocalBegin(dmSol, solRef, INSERT_VALUES, l);


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

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icuy_front_down[2]]);
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icuy_front_up[2]]);
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icuy_back_up[2]]);
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icuy_back_down[2]]);
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VShifted);
    VecAssemblyEnd(VShifted);
    DMGlobalToLocalEnd(dmSol, solRef, INSERT_VALUES, l);
    PetscObjectDestroy((PetscObject*)&l);
    


    return 0;
}

static PetscErrorCode SecondShiftW_z(DM dmSol, Vec *pWShifted, Vec solRef) //ok
{

    Vec WShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_up[3], icux_back_down[3], icux_front_up[3], icux_front_down[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pWShifted);
    WShifted = *pWShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_UP, d, &icux_front_up[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, &icux_front_down[d]);
        DMStagGetLocationSlot(dmCoord, BACK_UP, d, &icux_back_up[d]);
        DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, &icux_back_down[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

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

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_current, &current);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_next, &next);
                    inter = (next + current) / 2.0;
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_back_down[2]]);
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_front_down[2]]);
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_front_up[2]]);
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
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
                                  arrCoord[ez][ey][ex][icux_back_up[2]]);
                    DMStagVecSetValuesStencil(dmSol, WShifted, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WShifted);
    VecAssemblyEnd(WShifted);
    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);

    PetscObjectDestroy((PetscObject*)&l);


    return 0;
}

static PetscErrorCode SecondDerive_x(DM dmSol, Vec *pAB_x, Vec AB) //ok
{
    Vec             AB_x, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hx;
    PetscFunctionBeginUser;
   // DMCreateGlobalVector(dmSol, pAB_x);
    AB_x = *pAB_x;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hx = 1.0 / N[0];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

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
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);
    PetscObjectDestroy((PetscObject*)&AB_local);

    return 0;
}

static PetscErrorCode SecondDerive_z(DM dmSol, Vec *pAB_z, Vec AB) //ok
{
    Vec             AB_z, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hz;
    PetscFunctionBeginUser;
   // DMCreateGlobalVector(dmSol, pAB_z);
    AB_z = *pAB_z;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hz = 1.0 / N[2];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

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
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);
    PetscObjectDestroy((PetscObject*)&AB_local);


    return 0;
}

// Third non-linear members: mixed
static PetscErrorCode ThirdDerive_x(DM dmSol, Vec *pAB_x, Vec AB) //ok
{
    Vec             AB_x, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hx;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pAB_x);
    AB_x = *pAB_x;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hx = 1.0 / N[0];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

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
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);
    PetscObjectDestroy((PetscObject*)&AB_local);

    return 0;
}

static PetscErrorCode ThirdDerive_y(DM dmSol, Vec *pAB_y, Vec AB) //ok
{
    Vec             AB_y, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    PetscReal       hy;
    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pAB_y);
    AB_y = *pAB_y;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hy = 1.0 / N[1];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

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
    DMGlobalToLocalEnd(dmSol, AB, INSERT_VALUES, AB_local);
    PetscObjectDestroy((PetscObject*)&AB_local);

    return 0;
}

// non-linear members: homogeneous
// non-linear members: homogeneous
static PetscErrorCode CenterU(DM dmSol, Vec *pUCenter, Vec solRef)
{

    Vec UCenter, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux[3], icux_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pUCenter);
    UCenter = *pUCenter;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ex == N[0] - 1) {
                    PetscScalar inter, prev;
                    DMStagStencil row_prev;
                    row_prev.i = ex;
                    row_prev.j = ey;
                    row_prev.k = ez;
                    row_prev.loc = LEFT;
                    row_prev.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_prev, &prev);
                    inter = ((prev + uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]],
                                           arrCoord[ez][ey][ex][icux_right[2]])) / 2.0) * ((prev + uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]],
                                                                                                         arrCoord[ez][ey][ex][icux_right[2]])) / 2.0);
                    DMStagVecSetValuesStencil(dmSol, UCenter, 1, &row, &inter, INSERT_VALUES);
                } else if(ex == 0) {
                    PetscScalar inter, prev;
                    DMStagStencil row_prev;
                    row_prev.i = ex;
                    row_prev.j = ey;
                    row_prev.k = ez;
                    row_prev.loc = RIGHT;
                    row_prev.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_prev, &prev);
                    inter = ((prev + uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]],
                                           arrCoord[ez][ey][ex][icux[2]])) / 2.0) * ((prev + uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]],
                                                                                                   arrCoord[ez][ey][ex][icux[2]])) / 2.0);
                    DMStagVecSetValuesStencil(dmSol, UCenter, 1, &row, &inter, INSERT_VALUES);

                } else {
                    PetscScalar inter, next, prev;
                    DMStagStencil row_prev;
                    row_prev.i = ex;
                    row_prev.j = ey;
                    row_prev.k = ez;
                    row_prev.loc = LEFT;
                    row_prev.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey;
                    row_next.k = ez;
                    row_next.loc = RIGHT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;


                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_prev, &prev);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_next, &next);
                    inter = ((next + prev) / 2.0) * ((next + prev) / 2.0);
                    DMStagVecSetValuesStencil(dmSol, UCenter, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(UCenter);
    VecAssemblyEnd(UCenter);
    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);

    

    return 0;
}

static PetscErrorCode CenterV(DM dmSol, Vec *pVCenter, Vec solRef) 
{

    Vec VCenter, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pVCenter);
    VCenter = *pVCenter;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
   
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey == N[1] - 1) {
                    PetscScalar inter, prev;
                    DMStagStencil row_prev;
                    row_prev.i = ex;
                    row_prev.j = ey;
                    row_prev.k = ez;
                    row_prev.loc = DOWN;
                    row_prev.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_prev, &prev);
                    inter = ((prev + uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]],
                                           arrCoord[ez][ey][ex][icuy_up[2]])) / 2.0) * ((prev + uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]],
                                                                                                      arrCoord[ez][ey][ex][icuy_up[2]])) / 2.0);
                    DMStagVecSetValuesStencil(dmSol, VCenter, 1, &row, &inter, INSERT_VALUES);
                } else if(ey == 0) {
                    PetscScalar inter, prev;
                    DMStagStencil row_prev;
                    row_prev.i = ex;
                    row_prev.j = ey;
                    row_prev.k = ez;
                    row_prev.loc = UP;
                    row_prev.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_prev, &prev);
                    inter = ((prev + uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]],
                                           arrCoord[ez][ey][ex][icuy[2]])) / 2.0) * ((prev + uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]],
                                                                                                   arrCoord[ez][ey][ex][icuy[2]])) / 2.0);
                    DMStagVecSetValuesStencil(dmSol, VCenter, 1, &row, &inter, INSERT_VALUES);

                } else {
                    PetscScalar inter, next, prev;
                    DMStagStencil row_prev;
                    row_prev.i = ex;
                    row_prev.j = ey;
                    row_prev.k = ez;
                    row_prev.loc = DOWN;
                    row_prev.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey;
                    row_next.k = ez;
                    row_next.loc = UP;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;


                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_prev, &prev);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_next, &next);
                    inter = ((next + prev) / 2.0) * ((next + prev) / 2.0);
                    DMStagVecSetValuesStencil(dmSol, VCenter, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VCenter);
    VecAssemblyEnd(VCenter);
    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);


    return 0;
}

static PetscErrorCode CenterW(DM dmSol, Vec *pWCenter, Vec solRef) 
{

    Vec WCenter, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuz[3], icuz_front[3];
    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pWCenter);
    WCenter = *pWCenter;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
  
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez == N[2] - 1) {
                    PetscScalar inter, prev;
                    DMStagStencil row_prev;
                    row_prev.i = ex;
                    row_prev.j = ey;
                    row_prev.k = ez;
                    row_prev.loc = BACK;
                    row_prev.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_prev, &prev);
                    inter = ((prev + uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]],
                                           arrCoord[ez][ey][ex][icuz_front[2]])) / 2.0) * ((prev + uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]],
                                                                                                         arrCoord[ez][ey][ex][icuz_front[2]])) / 2.0);
                    DMStagVecSetValuesStencil(dmSol, WCenter, 1, &row, &inter, INSERT_VALUES);
                } else if(ez == 0) {
                    PetscScalar inter, prev;
                    DMStagStencil row_prev;
                    row_prev.i = ex;
                    row_prev.j = ey;
                    row_prev.k = ez;
                    row_prev.loc = FRONT;
                    row_prev.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_prev, &prev);
                    inter = ((prev + uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]],
                                           arrCoord[ez][ey][ex][icuz[2]])) / 2.0) * ((prev + uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]],
                                                                                                   arrCoord[ez][ey][ex][icuz[2]])) / 2.0);
                    DMStagVecSetValuesStencil(dmSol, WCenter, 1, &row, &inter, INSERT_VALUES);

                } else {
                    PetscScalar inter, next, prev;
                    DMStagStencil row_prev;
                    row_prev.i = ex;
                    row_prev.j = ey;
                    row_prev.k = ez;
                    row_prev.loc = BACK;
                    row_prev.c = 0;
                    DMStagStencil row_next;
                    row_next.i = ex;
                    row_next.j = ey;
                    row_next.k = ez;
                    row_next.loc = FRONT;
                    row_next.c = 0;
                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;


                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_prev, &prev);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_next, &next);
                    inter = ((next + prev) / 2.0) * ((next + prev) / 2.0);
                    DMStagVecSetValuesStencil(dmSol, WCenter, 1, &row, &inter, INSERT_VALUES);
                }
            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WCenter);
    VecAssemblyEnd(WCenter);
    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);


    return 0;
}

static PetscErrorCode Derive_x(DM dmSol, Vec *pU2_x, Vec U2)
{
    Vec             U2_x, coordLocal;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hx;
    DM              dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pU2_x);
    U2_x = *pU2_x;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hx = 1.0 / N[0];

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    PetscInt icux[3], icux_right[3], icux_e[3];
    
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icux_e[d]);

    }


    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,U2,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                /*if (ex != N[0] - 1) {


                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex + 1;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = ELEMENT;
                    row_right.c = 0;

                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = ELEMENT;
                    row_left.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    der = (val_right - val_left) / hx;

                    DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row, &der, INSERT_VALUES);
                    
                }

                if(ex == 0) {
                        DMStagStencil row_left;
                        PetscScalar val_left;
                        row_left.i = ex;
                        row_left.j = ey;
                        row_left.k = ez;
                        row_left.loc = ELEMENT;
                        row_left.c = 0;

                        DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                        DMStagStencil row_l;
                        PetscScalar val_l;
                        row_l.i = ex;
                        row_l.j = ey;
                        row_l.k = ez;
                        row_l.loc = LEFT;
                        row_l.c = 0;

                        val_l = (val_left - uxRef(arrCoord[ez][ey][ex][icux[0]] - hx, arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]])*uxRef(arrCoord[ez][ey][ex][icux[0]] - hx, arrCoord[ez][ey][ex][icux[1]],
                                            arrCoord[ez][ey][ex][icux[2]]))/hx;

                        DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row_l, &val_l, INSERT_VALUES);
                    }

                

                if(ex == N[0] - 1){
                        DMStagStencil row_right;
                        PetscScalar val_right;
                        row_right.i = ex;
                        row_right.j = ey;
                        row_right.k = ez;
                        row_right.loc = ELEMENT;
                        row_right.c = 0;

                        DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);

                        DMStagStencil row_r;
                        PetscScalar val_r;
                        row_r.i = ex;
                        row_r.j = ey;
                        row_r.k = ez;
                        row_r.loc = RIGHT;
                        row_r.c = 0;

                        val_r = (uxRef(arrCoord[ez][ey][ex][icux_right[0]] + hx, arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]])*uxRef(arrCoord[ez][ey][ex][icux_right[0]] + hx, arrCoord[ez][ey][ex][icux_right[1]],
                                            arrCoord[ez][ey][ex][icux_right[2]]) - val_right)/hx;
                        DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row_r, &val_r, INSERT_VALUES);
                }*/

                if (ex != 0) {


                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex - 1;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = ELEMENT;
                    row_right.c = 0;

                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = ELEMENT;
                    row_left.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = LEFT;
                    row.c = 0;
                    der = (val_left - val_right) / hx;

                    DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row, &der, INSERT_VALUES);                    
                }

                if(ex == 0) {
                        DMStagStencil row_left;
                        PetscScalar val_left;
                        row_left.i = ex;
                        row_left.j = ey;
                        row_left.k = ez;
                        row_left.loc = ELEMENT;
                        row_left.c = 0;

                        DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                        DMStagStencil row_l;
                        PetscScalar val_l;
                        row_l.i = ex;
                        row_l.j = ey;
                        row_l.k = ez;
                        row_l.loc = LEFT;
                        row_l.c = 0;

                        val_l = (val_left - uxRef(arrCoord[ez][ey][ex][icux_e[0]] - hx, arrCoord[ez][ey][ex][icux_e[1]], arrCoord[ez][ey][ex][icux_e[2]])*uxRef(arrCoord[ez][ey][ex][icux_e[0]] - hx, arrCoord[ez][ey][ex][icux_e[1]],
                                            arrCoord[ez][ey][ex][icux_e[2]]))/hx;
                        DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row_l, &val_l, INSERT_VALUES);

                    }

                

                if(ex == N[0] - 1){
                        DMStagStencil row_right;
                        PetscScalar val_right;
                        row_right.i = ex;
                        row_right.j = ey;
                        row_right.k = ez;
                        row_right.loc = ELEMENT;
                        row_right.c = 0;

                        DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);

                        DMStagStencil row_r;
                        PetscScalar val_r;
                        row_r.i = ex;
                        row_r.j = ey;
                        row_r.k = ez;
                        row_r.loc = RIGHT;
                        row_r.c = 0;



                        val_r = (uxRef(arrCoord[ez][ey][ex][icux_e[0]] + hx, arrCoord[ez][ey][ex][icux_e[1]], arrCoord[ez][ey][ex][icux_e[2]])*uxRef(arrCoord[ez][ey][ex][icux_e[0]] + hx, arrCoord[ez][ey][ex][icux_e[1]],
                                            arrCoord[ez][ey][ex][icux_e[2]]) - val_right)/hx;


                        DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row_r, &val_r, INSERT_VALUES);

                }


            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(U2_x);
    VecAssemblyEnd(U2_x);
    DMGlobalToLocalEnd(dmSol,U2,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);




    return 0;
}

static PetscErrorCode Derive_y(DM dmSol, Vec *pV2_y, Vec V2)
{
    Vec             V2_y, coordLocal;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hy;
    DM              dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pV2_y);
    V2_y = *pV2_y;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hy = 1.0 / N[1];

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,V2,INSERT_VALUES,l);

    PetscInt icuy[3], icuy_up[3], icuy_e[3];
    
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icuy_e[d]);
    }

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if(ey != 0){
                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey - 1;
                    row_right.k = ez;
                    row_right.loc = ELEMENT;
                    row_right.c = 0;

                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = ELEMENT;
                    row_left.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = DOWN;
                    row.c = 0;
                    der = (val_left - val_right) / hy;

                    DMStagVecSetValuesStencil(dmSol, V2_y, 1, &row, &der, INSERT_VALUES);
                }

                if(ey == 0) {
                        DMStagStencil row_left;
                        PetscScalar val_left;
                        row_left.i = ex;
                        row_left.j = ey;
                        row_left.k = ez;
                        row_left.loc = ELEMENT;
                        row_left.c = 0;

                        DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                        DMStagStencil row_l;
                        PetscScalar val_l;
                        row_l.i = ex;
                        row_l.j = ey;
                        row_l.k = ez;
                        row_l.loc = DOWN;
                        row_l.c = 0;
                        val_l = (val_left - uyRef(arrCoord[ez][ey][ex][icuy_e[0]], arrCoord[ez][ey][ex][icuy_e[1]] - hy, arrCoord[ez][ey][ex][icuy_e[2]])*uyRef(arrCoord[ez][ey][ex][icuy_e[0]], arrCoord[ez][ey][ex][icuy_e[1]] - hy, arrCoord[ez][ey][ex][icuy_e[2]]))/hy;

                        DMStagVecSetValuesStencil(dmSol, V2_y, 1, &row_l, &val_l, INSERT_VALUES);
                }

                if(ey == N[1] - 1){
                        DMStagStencil row_right;
                        PetscScalar val_right;
                        row_right.i = ex;
                        row_right.j = ey;
                        row_right.k = ez;
                        row_right.loc = ELEMENT;
                        row_right.c = 0;

                        DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);

                        DMStagStencil row_r;
                        PetscScalar val_r;
                        row_r.i = ex;
                        row_r.j = ey;
                        row_r.k = ez;
                        row_r.loc = UP;
                        row_r.c = 0;
                        val_r = (uyRef(arrCoord[ez][ey][ex][icuy_e[0]], arrCoord[ez][ey][ex][icuy_e[1]] + hy, arrCoord[ez][ey][ex][icuy_e[2]])*uyRef(arrCoord[ez][ey][ex][icuy_e[0]], arrCoord[ez][ey][ex][icuy_e[1]] + hy, arrCoord[ez][ey][ex][icuy_e[2]]) - val_right)/hy;

                        DMStagVecSetValuesStencil(dmSol, V2_y, 1, &row_r, &val_r, INSERT_VALUES);
                }


            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(V2_y);
    VecAssemblyEnd(V2_y);
    DMGlobalToLocalEnd(dmSol,V2,INSERT_VALUES,l);

    PetscObjectDestroy((PetscObject*)&l);

    return 0;
}

static PetscErrorCode Derive_z(DM dmSol, Vec *pW2_z, Vec W2)
{
    Vec             W2_z, coordLocal;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hz;
    DM              dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pW2_z);
    W2_z = *pW2_z;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hz = 1.0 / N[2];

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,W2,INSERT_VALUES,l);


    PetscInt icuz[3], icuz_front[3], icuz_e[3];
    
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icuz_e[d]);
    }

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ez != 0) {
                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey;
                    row_right.k = ez - 1;
                    row_right.loc = ELEMENT;
                    row_right.c = 0;

                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = ELEMENT;
                    row_left.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = BACK;
                    row.c = 0;
                    der = (val_left - val_right) / hz;

                    DMStagVecSetValuesStencil(dmSol, W2_z, 1, &row, &der, INSERT_VALUES);
                }

                if(ez == 0) {
                        DMStagStencil row_left;
                        PetscScalar val_left;
                        row_left.i = ex;
                        row_left.j = ey;
                        row_left.k = ez;
                        row_left.loc = ELEMENT;
                        row_left.c = 0;

                        DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                        DMStagStencil row_l;
                        PetscScalar val_l;
                        row_l.i = ex;
                        row_l.j = ey;
                        row_l.k = ez;
                        row_l.loc = BACK;
                        row_l.c = 0;
                        val_l = (val_left - uzRef(arrCoord[ez][ey][ex][icuz_e[0]], arrCoord[ez][ey][ex][icuz_e[1]], arrCoord[ez][ey][ex][icuz_e[2]] - hz)*uzRef(arrCoord[ez][ey][ex][icuz_e[0]], arrCoord[ez][ey][ex][icuz_e[1]], arrCoord[ez][ey][ex][icuz_e[2]] - hz))/hz;

                        DMStagVecSetValuesStencil(dmSol, W2_z, 1, &row_l, &val_l, INSERT_VALUES);
                }

                if(ez == N[2] - 1){
                        DMStagStencil row_right;
                        PetscScalar val_right;
                        row_right.i = ex;
                        row_right.j = ey;
                        row_right.k = ez;
                        row_right.loc = ELEMENT;
                        row_right.c = 0;

                        DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);

                        DMStagStencil row_r;
                        PetscScalar val_r;
                        row_r.i = ex;
                        row_r.j = ey;
                        row_r.k = ez;
                        row_r.loc = FRONT;
                        row_r.c = 0;
                        val_r = (uzRef(arrCoord[ez][ey][ex][icuz_e[0]], arrCoord[ez][ey][ex][icuz_e[1]], arrCoord[ez][ey][ex][icuz_e[2]] + hz)*uzRef(arrCoord[ez][ey][ex][icuz_e[0]], arrCoord[ez][ey][ex][icuz_e[1]],
                                            arrCoord[ez][ey][ex][icuz_e[2]] + hz) - val_right)/hz;
                        DMStagVecSetValuesStencil(dmSol, W2_z, 1, &row_r, &val_r, INSERT_VALUES);
                }

            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(W2_z);
    VecAssemblyEnd(W2_z);
    DMGlobalToLocalEnd(dmSol,W2,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);



    return 0;
}

// Assembling advection term
static PetscErrorCode ManageAdvection_x(PetscScalar dt, Vec* U_int, Vec U_n, Vec V_n, Vec W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz)
{

    // Create necessary grids
    DM dmSol_Shifted, dmSol_Centered, dmSol_Staggered;
    {
        CreateGrid(&dmSol_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    }

    Vec U_shift;
    DMCreateGlobalVector(dmSol_Shifted, &U_shift);
    DMStagMigrateVec(dmSol_Staggered, U_n, dmSol_Shifted, U_shift);

    Vec V_shift;
    DMCreateGlobalVector(dmSol_Shifted, &V_shift);
    DMStagMigrateVec(dmSol_Staggered, V_n, dmSol_Shifted, V_shift);

    Vec W_shift;
    DMCreateGlobalVector(dmSol_Shifted, &W_shift);
    DMStagMigrateVec(dmSol_Staggered, W_n, dmSol_Shifted, W_shift);

    // Managing first mixed non-linear term, please DO NOT remove any comment
    Vec UV_y, UW_z;
    DMCreateGlobalVector(dmSol_Shifted, &UV_y);
    DMCreateGlobalVector(dmSol_Shifted, &UW_z);
    Vec mixedFirst;
    DMCreateGlobalVector(dmSol_Staggered, &mixedFirst);
    {
        Vec U_y, U_z, V_y, W_z, UV, UW;
        DMCreateGlobalVector(dmSol_Shifted, &U_y);
        DMCreateGlobalVector(dmSol_Shifted, &U_z);
        DMCreateGlobalVector(dmSol_Shifted, &V_y);
        DMCreateGlobalVector(dmSol_Shifted, &W_z);
        DMCreateGlobalVector(dmSol_Shifted, &UV);
        DMCreateGlobalVector(dmSol_Shifted, &UW);

        FirstShiftU_y(dmSol_Shifted, &U_y, U_shift);
        FirstShiftU_z(dmSol_Shifted, &U_z, U_shift);
        FirstShiftV_y(dmSol_Shifted, &V_y, V_shift);
        FirstShiftW_z(dmSol_Shifted, &W_z, W_shift); 

        VecPointwiseMult(UV, U_y, V_y);
        VecPointwiseMult(UW, U_z, W_z);

        FirstDerive_y(dmSol_Shifted, &UV_y, UV);
        FirstDerive_z(dmSol_Shifted, &UW_z, UW);

        VecAXPY(UV_y, 1.0, UW_z);

        PetscObjectDestroy((PetscObject*)&U_y);
        PetscObjectDestroy((PetscObject*)&U_z);
        PetscObjectDestroy((PetscObject*)&V_y);
        PetscObjectDestroy((PetscObject*)&W_z);
        PetscObjectDestroy((PetscObject*)&UV);
        PetscObjectDestroy((PetscObject*)&UW);

    }

    DMStagMigrateVec(dmSol_Shifted, UV_y, dmSol_Staggered, mixedFirst);

    // Managing homogenous components of non-linear term
    //////////////////////////////////////////////////////////////
    Vec U_center;
    DMCreateGlobalVector(dmSol_Centered, &U_center);
    DMStagMigrateVec(dmSol_Staggered, U_n, dmSol_Centered, U_center);

    Vec U_c;
    DMCreateGlobalVector(dmSol_Centered, &U_c);
    CenterU(dmSol_Centered, &U_c, U_center);

    Vec U2_x;
    DMCreateGlobalVector(dmSol_Centered, &U2_x);
    Derive_x(dmSol_Centered, &U2_x, U_c);

    Vec homoFirst;
    DMCreateGlobalVector(dmSol_Staggered, &homoFirst);
    DMStagMigrateVec(dmSol_Centered, U2_x, dmSol_Staggered, homoFirst);

    VecAXPBYPCZ(U_n, -dt, -dt, 1.0, homoFirst, mixedFirst);
  

    // Copy to output

    VecCopy(U_n, *U_int);

    PetscObjectDestroy((PetscObject*)&dmSol_Shifted);
    PetscObjectDestroy((PetscObject*)&dmSol_Centered);
    PetscObjectDestroy((PetscObject*)&dmSol_Staggered);
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

static PetscErrorCode ManageAdvection_y(PetscScalar dt, Vec* V_int, Vec U_n, Vec V_n, Vec W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz)
{

    // Create necessary grids
    DM dmSol_Shifted, dmSol_Centered, dmSol_Staggered;
    {
        CreateGrid(&dmSol_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    }


    Vec U_shift;
    DMCreateGlobalVector(dmSol_Shifted, &U_shift);
    DMStagMigrateVec(dmSol_Staggered, U_n, dmSol_Shifted, U_shift);

    Vec V_shift;
    DMCreateGlobalVector(dmSol_Shifted, &V_shift);
    DMStagMigrateVec(dmSol_Staggered, V_n, dmSol_Shifted, V_shift);

    Vec W_shift;
    DMCreateGlobalVector(dmSol_Shifted, &W_shift);
    DMStagMigrateVec(dmSol_Staggered, W_n, dmSol_Shifted, W_shift);

    
    // Managing second mixed non-linear term, please DO NOT remove any comment
    Vec VU_x, VW_z;
    DMCreateGlobalVector(dmSol_Shifted, &VU_x);
    DMCreateGlobalVector(dmSol_Shifted, &VW_z);
    Vec mixedSecond;
    DMCreateGlobalVector(dmSol_Staggered, &mixedSecond);

    {
        Vec V_x, V_z, U_x, W_z, VU, VW;
        DMCreateGlobalVector(dmSol_Shifted, &V_x);
        DMCreateGlobalVector(dmSol_Shifted, &V_z);
        DMCreateGlobalVector(dmSol_Shifted, &U_x);
        DMCreateGlobalVector(dmSol_Shifted, &W_z);
        DMCreateGlobalVector(dmSol_Shifted, &VU);
        DMCreateGlobalVector(dmSol_Shifted, &VW);


        FirstShiftV_y(dmSol_Shifted, &V_x, V_shift);
        SecondShiftV_z(dmSol_Shifted, &V_z, V_shift);

        FirstShiftU_y(dmSol_Shifted, &U_x, U_shift);
        SecondShiftW_z(dmSol_Shifted, &W_z, W_shift);


        VecPointwiseMult(VU, V_x, U_x);
        VecPointwiseMult(VW, V_z, W_z);



        SecondDerive_x(dmSol_Shifted, &VU_x, VU);
        SecondDerive_z(dmSol_Shifted, &VW_z, VW);

        VecAXPY(VU_x, 1.0, VW_z);

        PetscObjectDestroy((PetscObject*)&V_x);
        PetscObjectDestroy((PetscObject*)&V_z);
        PetscObjectDestroy((PetscObject*)&U_x);
        PetscObjectDestroy((PetscObject*)&W_z);
        PetscObjectDestroy((PetscObject*)&VU);
        PetscObjectDestroy((PetscObject*)&VW);
    }
    DMStagMigrateVec(dmSol_Shifted, VU_x, dmSol_Staggered, mixedSecond);

    ///////////////////////////////

    Vec V_center;
    DMCreateGlobalVector(dmSol_Centered, &V_center);
    DMStagMigrateVec(dmSol_Staggered, V_n, dmSol_Centered, V_center);

    Vec V_c;
    DMCreateGlobalVector(dmSol_Centered, &V_c);
    CenterV(dmSol_Centered, &V_c, V_center);

    Vec V2_y;
    DMCreateGlobalVector(dmSol_Centered, &V2_y);
    Derive_y(dmSol_Centered, &V2_y, V_c);

    Vec homoSecond;
    DMCreateGlobalVector(dmSol_Staggered, &homoSecond);
    DMStagMigrateVec(dmSol_Centered, V2_y, dmSol_Staggered, homoSecond);

    VecAXPBYPCZ(V_n, -dt, -dt, 1.0, homoSecond, mixedSecond);



    /*VecAXPBY(homoFirst, -1.0, -1.0, mixedFirst);

    VecScale(homoFirst, dt);

    VecAXPY(homoFirst, 1.0, U_n);*/    

    // Copy to output

    VecCopy(V_n, *V_int);


    PetscObjectDestroy((PetscObject*)&dmSol_Shifted);
    PetscObjectDestroy((PetscObject*)&dmSol_Centered);
    PetscObjectDestroy((PetscObject*)&dmSol_Staggered);
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

static PetscErrorCode ManageAdvection_z(PetscScalar dt, Vec* W_int, Vec U_n, Vec V_n, Vec W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz)
{

    // Create necessary grids
    DM dmSol_Shifted, dmSol_Centered, dmSol_Staggered;
    {
        CreateGrid(&dmSol_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    }


    Vec U_shift;
    DMCreateGlobalVector(dmSol_Shifted, &U_shift);
    DMStagMigrateVec(dmSol_Staggered, U_n, dmSol_Shifted, U_shift);

    Vec V_shift;
    DMCreateGlobalVector(dmSol_Shifted, &V_shift);
    DMStagMigrateVec(dmSol_Staggered, V_n, dmSol_Shifted, V_shift);

    Vec W_shift;
    DMCreateGlobalVector(dmSol_Shifted, &W_shift);
    DMStagMigrateVec(dmSol_Staggered, W_n, dmSol_Shifted, W_shift);

    
    // Managing third mixed non-linear term, please DO NOT remove any comment, ESPECCIALLY the "==[] ones"
    Vec WU_x, WV_y;
    DMCreateGlobalVector(dmSol_Shifted, &WU_x);
    DMCreateGlobalVector(dmSol_Shifted, &WV_y);
    Vec mixedThird;
    DMCreateGlobalVector(dmSol_Staggered, &mixedThird);
    {
        Vec W_x, W_y, U_x, V_y, WU, WV;
        DMCreateGlobalVector(dmSol_Shifted, &W_x);
        DMCreateGlobalVector(dmSol_Shifted, &W_y);
        DMCreateGlobalVector(dmSol_Shifted, &U_x);
        DMCreateGlobalVector(dmSol_Shifted, &V_y);
        DMCreateGlobalVector(dmSol_Shifted, &WU);
        DMCreateGlobalVector(dmSol_Shifted, &WV);

        FirstShiftW_z(dmSol_Shifted, &W_x, W_shift);// ==FirsShiftW_z
        SecondShiftW_z(dmSol_Shifted, &W_y, W_shift);// ==SecondShiftW_z
        FirstShiftU_z(dmSol_Shifted, &U_x, U_shift);// ==FirsShiftU_z
        SecondShiftV_z(dmSol_Shifted, &V_y, V_shift);// ==SecondShiftV_z

        VecPointwiseMult(WU, W_x, U_x);
        VecPointwiseMult(WV, W_y, V_y);

        ThirdDerive_x(dmSol_Shifted, &WU_x, WU);
        ThirdDerive_y(dmSol_Shifted, &WV_y, WV);

        VecAXPY(WU_x, 1.0, WV_y);

        PetscObjectDestroy((PetscObject*)&W_x);
        PetscObjectDestroy((PetscObject*)&W_y);
        PetscObjectDestroy((PetscObject*)&U_x);
        PetscObjectDestroy((PetscObject*)&V_y);
        PetscObjectDestroy((PetscObject*)&WU);
        PetscObjectDestroy((PetscObject*)&WV);
    }
    DMStagMigrateVec(dmSol_Shifted, WU_x, dmSol_Staggered, mixedThird);




    Vec W_center;
    DMCreateGlobalVector(dmSol_Centered, &W_center);
    DMStagMigrateVec(dmSol_Staggered, W_n, dmSol_Centered, W_center);



    Vec W_c;
    DMCreateGlobalVector(dmSol_Centered, &W_c);
    CenterW(dmSol_Centered, &W_c, W_center);

    Vec W2_z;
    DMCreateGlobalVector(dmSol_Centered, &W2_z);
    Derive_z(dmSol_Centered, &W2_z, W_c);

    Vec homoThird;
    DMCreateGlobalVector(dmSol_Staggered, &homoThird);
    DMStagMigrateVec(dmSol_Centered, W2_z, dmSol_Staggered, homoThird);



    // Assembling whole non-linear
    VecAXPBYPCZ(W_n, -dt, -dt, 1.0, homoThird, mixedThird);

    VecCopy(W_n, *W_int);

    PetscObjectDestroy((PetscObject*)&dmSol_Shifted);
    PetscObjectDestroy((PetscObject*)&dmSol_Centered);
    PetscObjectDestroy((PetscObject*)&dmSol_Staggered);
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

// Second step: viscosity members to solve implicit laplacian
static PetscErrorCode Assemble_x(DM dmSol, Mat *pA, Vec *pRhs, Vec solRef, PetscScalar dt, PetscScalar Re) {

    Vec rhs, coordLocal;
    Mat A;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux[3], icuy[3], icuz[3], icux_right[3], icuy_up[3], icuz_front[3], icux_back_left[3], icux_down_left[3], icux_front_left[3], icux_up_left[3];
    PetscReal hx, hy, hz;
    PetscScalar Ret = Re/dt;

    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateMatrix(dmSol, pA);
    A = *pA;
    DMCreateGlobalVector(dmSol, pRhs);
    rhs = *pRhs;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);

    hx = 1.0 / N[0];
    hy = 1.0 / N[1];
    hz = 1.0 / N[2];
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {
                    /* Right Boundary velocity Dirichlet */
                    DMStagStencil row;
                    PetscScalar valRhs;
                    const PetscScalar valA = 1.0;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]],
                                   arrCoord[ez][ey][ex][icux_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }        

                /* Equation on left face of this element */
                if (ex == 0) {
                    /* Left velocity Dirichlet */
                    DMStagStencil row;
                    PetscScalar valRhs;
                    const PetscScalar valA = 1.0;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = LEFT;
                    row.c = 0;

                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]],
                                   arrCoord[ez][ey][ex][icux[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                } else {
                    /* X-momentum interior equation : (u_xx + u_yy + u_zz) - p_x = f^x */

                    DMStagStencil row, col[7];
                    PetscScalar valA[7], valRhs;
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            PetscScalar bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]-hy, arrCoord[ez][ey][ex][icux[2]]);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz);
                            valRhs = -Ret*valRhs - bc_1/(hy*hy) - bc_2/(hz*hz);
                            //ok
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                               

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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            PetscScalar bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] - hy, arrCoord[ez][ey][ex][icux[2]]);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]+hz);
                            valRhs = -Ret*valRhs - bc_1/(hy*hy) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            PetscScalar bc_2;

                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] - hy, arrCoord[ez][ey][ex][icux[2]]);

                            valRhs = -Ret*valRhs - bc_2/(hy*hy);                           
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                            /* Missing up term */
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
                            /* Missing back entry */
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            PetscScalar bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]+hy, arrCoord[ez][ey][ex][icux[2]]);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz);
                            valRhs = -Ret*valRhs -bc_1/(hy*hy) - bc_2/(hz*hz);                         
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            

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
                            /* Missing up term */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            PetscScalar bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]+hy, arrCoord[ez][ey][ex][icux[2]]);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]+hz);
                            
                            valRhs = -Ret*valRhs -bc_1/(hy*hy) - bc_2/(hz*hz);       

                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                     
                            /* Missing front term */
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
                            /* Missing up term */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            PetscScalar bc_2;
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]+hy, arrCoord[ez][ey][ex][icux[2]]);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                        /* Missing back term */
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = LEFT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        PetscScalar bc_1;
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz);
                        valRhs = -Ret*valRhs  - bc_1/(hz*hz);

                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                        PetscScalar bc_1;
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]+hz);
                        valRhs = -Ret*valRhs - bc_1/(hz*hz);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        valRhs = -Ret*valRhs;
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);

                }
                
            }
        }
    }
    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    VecDestroy(&l);

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyEnd(rhs);

    return 0;



}

static PetscErrorCode Assemble_y(DM dmSol, Mat *pA, Vec *pRhs, Vec solRef, PetscScalar dt, PetscScalar Re) {

    Vec rhs, coordLocal;
    Mat A;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up[3];
    PetscReal hx, hy, hz;
    PetscReal Ret = Re/dt;

    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateMatrix(dmSol, pA);
    A = *pA;
    DMCreateGlobalVector(dmSol, pRhs);
    rhs = *pRhs;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);

    hx = 1.0 / N[0];
    hy = 1.0 / N[1];
    hz = 1.0 / N[2];
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        
    }


    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
  
                if (ey == N[1] - 1) {
                    /* Top boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = UP;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }


                
                /* Equation on bottom face of this element */
                if (ey == 0) {
                    /* Bottom boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = DOWN;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);

                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                } else {
                    /* Y-momentum equation, (v_xx + v_yy + v_zz) - p_y = f^y */
                    DMStagStencil row, col[7];
                    PetscScalar   valA[7], valRhs;
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
                            /* Left term missing */
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            /* Back term missing */
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            PetscScalar bc_1, bc_2;
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

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
                            /* Left term missing */
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
                            PetscScalar bc_1, bc_2;
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz);
                            valRhs = -Ret*valRhs -bc_1/(hx*hx) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

                            /* Front term missing */
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
                            /* Left term missing */
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
                            PetscScalar bc_1;

                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);
                            valRhs = -Ret*valRhs -bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

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
                            /* Right term missing */
                            /* Back term missing */
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            PetscScalar bc_1, bc_2;
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                         
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
                            /* Right term missing */
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            /* Front term missing */
                            PetscScalar bc_1, bc_2;
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);


                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx) - bc_2/(hz*hz);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);


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
                            /* Right term missing */
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
                            PetscScalar bc_1;
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
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
                        /* Back term missing */
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = DOWN;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        PetscScalar bc_2;
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                        
                        bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz);
                        valRhs = -Ret*valRhs - bc_2/(hz*hz);


                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                        PetscScalar bc_1;
                        bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz);
                        valRhs = -Ret*valRhs - bc_1/(hz*hz);

                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                    
                        /* Front term missing */
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        valRhs = -Ret*valRhs;
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);

                }
                
            }
        }
    }

    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    VecDestroy(&l);


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);

    return 0;    

}

static PetscErrorCode Assemble_z(DM dmSol, Mat *pA, Vec *pRhs, Vec solRef, PetscScalar dt, PetscScalar Re) {
    Vec rhs, coordLocal;
    Mat A;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuz[3], icuz_front[3];
    PetscReal hx, hy, hz;
    PetscReal Ret = Re/dt;

    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateMatrix(dmSol, pA);
    A = *pA;
    DMCreateGlobalVector(dmSol, pRhs);
    rhs = *pRhs;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);

    hx = 1.0 / N[0];
    hy = 1.0 / N[1];
    hz = 1.0 / N[2];
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                
         
                if (ez == N[2] - 1) {
                    /* Front boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = FRONT;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }
                
                /* Equation on back face of this element */
                if (ez == 0) {
                    /* Back boundary velocity Dirichlet */
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = BACK;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                } else {
                    /* Z-momentum equation, (w_xx + w_yy + w_zz) - p_z = f^z */
                    DMStagStencil row, col[7];
                    PetscScalar   valA[7], valRhs;
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
                            /* Down term missing */
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            /* Left term missing */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            PetscScalar bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]]- hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]]);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);                            
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES); 

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
                            /* Up term missing */
                            /* Left term missing */
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

                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            PetscScalar bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]]-hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]]+hy, arrCoord[ez][ey][ex][icuz[2]]);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                             
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
                            /* Left term missing */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            PetscScalar bc_1;

                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES); 

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
                            /* Down term missing */
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
                            /* Right term missing */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);

                            PetscScalar bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] +hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] -hy, arrCoord[ez][ey][ex][icuz[2]]);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);       
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES); 

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
                            /* Up term missing */
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            /* Right term missing */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            PetscScalar bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]]+hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]]+hy, arrCoord[ez][ey][ex][icuz[2]]);
                            valRhs = -Ret*valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);

                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES); 

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
                            /* Right term missing */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            PetscScalar bc_1;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]]+hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
                            valRhs = -Ret*valRhs - bc_1/(hx*hx);

                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES); 

                        }
                    } else if (ey == 0) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        /* Down term missing */
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        PetscScalar bc_2;
                        bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]]);
                        valRhs = -Ret*valRhs - bc_2/(hy*hy);

                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                         
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
                        /* Up term missing */
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        PetscScalar bc_2;
                        bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]]+hy, arrCoord[ez][ey][ex][icuz[2]]);
                        valRhs = -Ret*valRhs - bc_2/(hy*hy);

                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES); 

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

                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        valRhs = -Ret*valRhs;
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES); 
                        
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);

                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    VecDestroy(&l);

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);
    return 0;
}

static PetscErrorCode CheckSolution(Vec sol, Vec solRef)
{
    Vec       diff;
    PetscReal normsolRef, errAbs, errRel;

    PetscFunctionBeginUser;
    VecDuplicate(sol, &diff);
    VecCopy(sol, diff);
    VecAXPY(diff, -1.0, solRef);
    VecNorm(diff, NORM_2, &errAbs);
    VecNorm(solRef, NORM_2, &normsolRef);
    errRel = errAbs / normsolRef;
    PetscPrintf(PETSC_COMM_WORLD, "Error (abs): %g\nError (rel): %g\n", (double)errAbs, (double)errRel);
    VecDestroy(&diff);
    return 0;
}

// Assembling Viscosity term
static PetscErrorCode ManageViscosity(PetscScalar dt, PetscScalar Re, Vec* pU_pre, Vec* pV_pre, Vec* pW_pre, Vec U_int, Vec V_int, Vec W_int, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz){

    DM dmSol_staggered;
    Mat A_x, A_y, A_z;
    Vec rhs_x, rhs_y, rhs_z;
    KSP       ksp_x, ksp_y, ksp_z;
    PC        pc_x, pc_y, pc_z;

    PetscFunctionBeginUser;
    
    CreateGrid(&dmSol_staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    {
    Assemble_x(dmSol_staggered, &A_x, &rhs_x, U_int, dt, Re);
    Vec U_pre;
    DMCreateGlobalVector(dmSol_staggered, &U_pre);
    KSPCreate(PETSC_COMM_WORLD, &ksp_x);
    KSPSetType(ksp_x, KSPCG);
    KSPSetOperators(ksp_x, A_x, A_x);
    KSPGetPC(ksp_x, &pc_x);
    PCSetType(pc_x, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc_x, PETSC_TRUE);
    KSPSetFromOptions(ksp_x);
    KSPSolve(ksp_x, rhs_x, U_pre);

    VecCopy(U_pre, *pU_pre);
    PetscObjectDestroy((PetscObject*)&ksp_x);
    PetscObjectDestroy((PetscObject*)&U_pre);
    PetscObjectDestroy((PetscObject*)&A_x);
    PetscObjectDestroy((PetscObject*)&rhs_x);
    }

    {
    Assemble_y(dmSol_staggered, &A_y, &rhs_y, V_int, dt, Re);
    Vec V_pre;
    DMCreateGlobalVector(dmSol_staggered, &V_pre);
    KSPCreate(PETSC_COMM_WORLD, &ksp_y);
    KSPSetType(ksp_y, KSPCG);
    KSPSetOperators(ksp_y, A_y, A_y);
    KSPGetPC(ksp_y, &pc_y);
    PCSetType(pc_y, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc_y, PETSC_TRUE);
    KSPSetFromOptions(ksp_y);
    KSPSolve(ksp_y, rhs_y, V_pre);

    VecCopy(V_pre, *pV_pre);
    PetscObjectDestroy((PetscObject*)&ksp_y);
    PetscObjectDestroy((PetscObject*)&V_pre);
    PetscObjectDestroy((PetscObject*)&A_y);
    PetscObjectDestroy((PetscObject*)&rhs_y);
    }
    
    {
    Assemble_z(dmSol_staggered, &A_z, &rhs_z, W_int, dt, Re);
    Vec W_pre;
    DMCreateGlobalVector(dmSol_staggered, &W_pre);
    KSPCreate(PETSC_COMM_WORLD, &ksp_z);
    KSPSetType(ksp_z, KSPCG);
    KSPSetOperators(ksp_z, A_z, A_z);
    KSPGetPC(ksp_z, &pc_z);
    PCSetType(pc_z, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc_z, PETSC_TRUE);
    KSPSetFromOptions(ksp_z);
    KSPSolve(ksp_z, rhs_z, W_pre);

    VecCopy(W_pre, *pW_pre);
    PetscObjectDestroy((PetscObject*)&ksp_z);
    PetscObjectDestroy((PetscObject*)&W_pre);
    PetscObjectDestroy((PetscObject*)&A_z);
    PetscObjectDestroy((PetscObject*)&rhs_z);
    }

      
    /*DMCreateGlobalVector(dmSol_staggered, pU_pre);
    DMCreateGlobalVector(dmSol_staggered, pV_pre);
    DMCreateGlobalVector(dmSol_staggered, pW_pre);*/

    PetscObjectDestroy((PetscObject*)&dmSol_staggered);

    return 0;
}

// Assembling Pressure

static PetscErrorCode AttachNullspace(DM dmSol, Mat A)
{
  DM           dmPressure;
  Vec          constantPressure, basis;
  PetscReal    nrm;
  MatNullSpace matNullSpace;

  PetscFunctionBeginUser;
  DMStagCreateCompatibleDMStag(dmSol, 0, 0, 1, 0, &dmPressure);
  DMGetGlobalVector(dmPressure, &constantPressure);
  VecSet(constantPressure, 1.0);
  VecNorm(constantPressure, NORM_2, &nrm);
  VecScale(constantPressure, 1.0 / nrm);
  DMCreateGlobalVector(dmSol, &basis);
  DMStagMigrateVec(dmPressure, constantPressure, dmSol, basis);
  MatNullSpaceCreate(PetscObjectComm((PetscObject)dmSol), PETSC_FALSE, 1, &basis, &matNullSpace);
  VecDestroy(&basis);
  MatSetNullSpace(A, matNullSpace);
  MatNullSpaceDestroy(&matNullSpace);
  DMRestoreGlobalVector(dmPressure, &constantPressure);
  DMDestroy(&dmPressure);
  return 0;
}

// Assembling Pressure

static PetscErrorCode Assemble_P(DM dmSol, Mat *pA, Vec *pRhs, Vec rhs_input) {
    Vec rhs, coordLocal;
    Mat A;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux[3], icuy[3], icuz[3], icux_right[3], icuy_up[3], icuz_front[3], icux_back_left[3], icux_down_left[3], icux_front_left[3], icux_up_left[3];
    PetscReal hx, hy, hz;

    DM dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateMatrix(dmSol, pA);
    A = *pA;
    DMCreateGlobalVector(dmSol, pRhs);
    rhs = *pRhs;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);

    hx = 1.0 / N[0];
    hy = 1.0 / N[1];
    hz = 1.0 / N[2];
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icux[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,rhs_input,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {

                    DMStagStencil row, col[6];
                    PetscScalar valA[6], valRhs;
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            

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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);                          
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                            /* Missing up term */
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            
                            /* Missing back entry */
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);                      
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            

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
                            /* Missing up term */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                     
                            /* Missing front term */
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
                            /* Missing up term */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                        
                        /* Missing back term */
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              

                } 

                /* Equation on left face of this element */
                else if (ex == 0) {

                    DMStagStencil row, col[6];
                    PetscScalar valA[6], valRhs;
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                               

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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);                          
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                            /* Missing up term */
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            
                            /* Missing back entry */
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);                      
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            

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
                            /* Missing up term */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                     
                            /* Missing front term */
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
                            /* Missing up term */
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                        
                        /* Missing back term */
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              

                } else {
                    /* X-momentum interior equation : (u_xx + u_yy + u_zz) - p_x = f^x */

                    DMStagStencil row, col[7];
                    PetscScalar valA[7], valRhs;
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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                               

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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);                          
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                            /* Missing back entry */
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);                      
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            

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
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                     
                            /* Missing front term */
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
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                        /* Missing back term */
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                        DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);

                }
                
            }
        }
    }
    DMGlobalToLocalEnd(dmSol,rhs_input,INSERT_VALUES,l);
    VecDestroy(&l);

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyEnd(rhs);

    return 0;

}

static PetscErrorCode AssembleDivergence(DM dmSol, Vec *pDiv, Vec U_x, Vec U_y, Vec U_z, PetscScalar dt) 
{

    Vec div, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    PetscScalar ****arrCoord;
    PetscReal hx, hy, hz;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pDiv);
    div = *pDiv;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
   
    hy = 1.0 / N[1];
    hx = 1.0 / N[0];
    hz = 1.0 / N[2];



    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    /*

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);

    }*/

    DM dmSol_y, dmSol_z;

    DMClone(dmSol, &dmSol_y);
    DMClone(dmSol, &dmSol_z);


    Vec l_x;
    DMCreateLocalVector(dmSol,&l_x);
    DMGlobalToLocalBegin(dmSol,U_x,INSERT_VALUES,l_x);

    Vec l_y;
    DMCreateLocalVector(dmSol_y,&l_y);
    DMGlobalToLocalBegin(dmSol_y,U_y,INSERT_VALUES,l_y);
    
    Vec l_z;
    DMCreateLocalVector(dmSol_z,&l_z);
    DMGlobalToLocalBegin(dmSol_z,U_z,INSERT_VALUES,l_z);

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {


                PetscScalar inter, left, right, up, down, front, back;
                DMStagStencil row_up;
                row_up.i = ex;
                row_up.j = ey;
                row_up.k = ez;
                row_up.loc = UP;
                row_up.c = 0;
                DMStagStencil row_down;
                row_down.i = ex;
                row_down.j = ey;
                row_down.k = ez;
                row_down.loc = DOWN;
                row_down.c = 0;
                DMStagStencil row_left;
                row_left.i = ex;
                row_left.j = ey;
                row_left.k = ez;
                row_left.loc = LEFT;
                row_left.c = 0;
                DMStagStencil row_right;
                row_right.i = ex;
                row_right.j = ey;
                row_right.k = ez;
                row_right.loc = RIGHT;
                row_right.c = 0;
                DMStagStencil row_front;
                row_front.i = ex;
                row_front.j = ey;
                row_front.k = ez;
                row_front.loc = FRONT;
                row_front.c = 0;
                DMStagStencil row_back;
                row_back.i = ex;
                row_back.j = ey;
                row_back.k = ez;
                row_back.loc = BACK;
                row_back.c = 0;


                    DMStagStencil row;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;


                DMStagVecGetValuesStencil(dmSol_y, l_y, 1, &row_up, &up);
                DMStagVecGetValuesStencil(dmSol_y, l_y, 1, &row_down, &down);
                DMStagVecGetValuesStencil(dmSol, l_x, 1, &row_left, &left);
                DMStagVecGetValuesStencil(dmSol, l_x, 1, &row_right, &right);
                DMStagVecGetValuesStencil(dmSol_z, l_z, 1, &row_front, &front);
                DMStagVecGetValuesStencil(dmSol_z, l_z, 1, &row_back, &back);

                    inter = ((up - down) / hy + (right - left) / hx + (front - back) / hz)/dt;
                    DMStagVecSetValuesStencil(dmSol, div, 1, &row, &inter, INSERT_VALUES);

            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(div);
    VecAssemblyEnd(div);



    DMGlobalToLocalEnd(dmSol,U_x,INSERT_VALUES,l_x);
    PetscObjectDestroy((PetscObject*)&l_x);
    DMGlobalToLocalEnd(dmSol_y,U_y,INSERT_VALUES,l_y);
    PetscObjectDestroy((PetscObject*)&l_y);
    DMGlobalToLocalEnd(dmSol_z,U_z,INSERT_VALUES,l_z);
    PetscObjectDestroy((PetscObject*)&l_z);
    DMDestroy(&dmSol_y);
    DMDestroy(&dmSol_z);



    return 0;
}

static PetscErrorCode ComputeDivergence(Vec* pDiv, Vec U_n, Vec V_n, Vec W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz, PetscScalar dt) 
{
        DM dmSol_centered, dmSol_shifted, dmSol_staggered;
        CreateGrid(&dmSol_centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_shifted, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

        Vec U_shifted, V_shifted, W_shifted;

        DMCreateGlobalVector(dmSol_shifted, &U_shifted);
        DMCreateGlobalVector(dmSol_shifted, &V_shifted);
        DMCreateGlobalVector(dmSol_shifted, &W_shifted);

        DMStagMigrateVec(dmSol_staggered, U_n, dmSol_shifted, U_shifted);
        DMStagMigrateVec(dmSol_staggered, V_n, dmSol_shifted, V_shifted);
        DMStagMigrateVec(dmSol_staggered, W_n, dmSol_shifted, W_shifted);

        Vec Div;
        AssembleDivergence(dmSol_shifted, &Div, U_shifted, V_shifted, W_shifted, dt);

        DMCreateGlobalVector(dmSol_centered, pDiv);
        DMStagMigrateVec(dmSol_shifted, Div, dmSol_centered, *pDiv);

        PetscObjectDestroy((PetscObject*)&dmSol_staggered);
        PetscObjectDestroy((PetscObject*)&dmSol_shifted);
        PetscObjectDestroy((PetscObject*)&dmSol_centered);
        PetscObjectDestroy((PetscObject*)&U_shifted);
        PetscObjectDestroy((PetscObject*)&V_shifted);
        PetscObjectDestroy((PetscObject*)&W_shifted);
        PetscObjectDestroy((PetscObject*)&Div);



        return 0;
}

static PetscErrorCode Derive_x_P(DM dmSol, Vec *pU2_x, Vec U2)
{
    Vec             U2_x, coordLocal;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hx;
    DM              dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pU2_x);
    U2_x = *pU2_x;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hx = 1.0 / N[0];

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    PetscInt icux[3], icux_right[3];
    
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
    }


    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,U2,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex != N[0] - 1) {
                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex + 1;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = ELEMENT;
                    row_right.c = 0;

                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = ELEMENT;
                    row_left.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    der = (val_right - val_left) / hx;

                    DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row, &der, INSERT_VALUES);
                    
                }

                if(ex == 0) {
                       

                        DMStagStencil row_l;
                        PetscScalar val_l;
                        row_l.i = ex;
                        row_l.j = ey;
                        row_l.k = ez;
                        row_l.loc = LEFT;
                        row_l.c = 0;
                        val_l = 0;

                        DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row_l, &val_l, INSERT_VALUES);
                    }

                

                if(ex == N[0] - 1){
            

                        DMStagStencil row_r;
                        PetscScalar val_r;
                        row_r.i = ex;
                        row_r.j = ey;
                        row_r.k = ez;
                        row_r.loc = RIGHT;
                        row_r.c = 0;
                        val_r = 0;
                        DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row_r, &val_r, INSERT_VALUES);
                }



            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(U2_x);
    VecAssemblyEnd(U2_x);
    DMGlobalToLocalEnd(dmSol,U2,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);




    return 0;
}

static PetscErrorCode Derive_y_P(DM dmSol, Vec *pV2_y, Vec V2)
{
    Vec             V2_y, coordLocal;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hy;
    DM              dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pV2_y);
    V2_y = *pV2_y;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hy = 1.0 / N[1];

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,V2,INSERT_VALUES,l);

    PetscInt icuy[3], icuy_up[3];
    
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
    }

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ey != N[1] - 1) {
                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey + 1;
                    row_right.k = ez;
                    row_right.loc = ELEMENT;
                    row_right.c = 0;

                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = ELEMENT;
                    row_left.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = UP;
                    row.c = 0;
                    der = (val_right - val_left) / hy;

                    DMStagVecSetValuesStencil(dmSol, V2_y, 1, &row, &der, INSERT_VALUES);
                }

                if(ey == 0) {
                        

                        DMStagStencil row_l;
                        PetscScalar val_l;
                        row_l.i = ex;
                        row_l.j = ey;
                        row_l.k = ez;
                        row_l.loc = DOWN;
                        row_l.c = 0;
                        val_l = 0;

                        DMStagVecSetValuesStencil(dmSol, V2_y, 1, &row_l, &val_l, INSERT_VALUES);
                }
                if(ey == N[1] - 1){
              

                        DMStagStencil row_r;
                        PetscScalar val_r;
                        row_r.i = ex;
                        row_r.j = ey;
                        row_r.k = ez;
                        row_r.loc = UP;
                        row_r.c = 0;
                        val_r = 0;
                        DMStagVecSetValuesStencil(dmSol, V2_y, 1, &row_r, &val_r, INSERT_VALUES);
                }


            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(V2_y);
    VecAssemblyEnd(V2_y);

    DMGlobalToLocalEnd(dmSol,V2,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);



    return 0;
}

static PetscErrorCode Derive_z_P(DM dmSol, Vec *pW2_z, Vec W2)
{
    Vec             W2_z, coordLocal;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hz;
    DM              dmCoord;
    PetscScalar ****arrCoord;

    PetscFunctionBeginUser;
    //DMCreateGlobalVector(dmSol, pW2_z);
    W2_z = *pW2_z;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    hz = 1.0 / N[2];

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,W2,INSERT_VALUES,l);


    PetscInt icuz[3], icuz_front[3];
    
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ez != N[2] - 1) {
                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey;
                    row_right.k = ez + 1;
                    row_right.loc = ELEMENT;
                    row_right.c = 0;

                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = ELEMENT;
                    row_left.c = 0;

                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_right, &val_right);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row_left, &val_left);

                    DMStagStencil row;
                    PetscScalar der;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = FRONT;
                    row.c = 0;
                    der = (val_right - val_left) / hz;

                    DMStagVecSetValuesStencil(dmSol, W2_z, 1, &row, &der, INSERT_VALUES);
                }

                if(ez == 0) {
                      

                        DMStagStencil row_l;
                        PetscScalar val_l;
                        row_l.i = ex;
                        row_l.j = ey;
                        row_l.k = ez;
                        row_l.loc = BACK;
                        row_l.c = 0;
                        val_l = 0;

                        DMStagVecSetValuesStencil(dmSol, W2_z, 1, &row_l, &val_l, INSERT_VALUES);
                }

                if(ez == N[2] - 1){
              

                        DMStagStencil row_r;
                        PetscScalar val_r;
                        row_r.i = ex;
                        row_r.j = ey;
                        row_r.k = ez;
                        row_r.loc = FRONT;
                        row_r.c = 0;
                        val_r = 0;
                        DMStagVecSetValuesStencil(dmSol, W2_z, 1, &row_r, &val_r, INSERT_VALUES);
                }

            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(W2_z);
    VecAssemblyEnd(W2_z);

    DMGlobalToLocalEnd(dmSol,W2,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);



    return 0;
}

static PetscErrorCode ManagePressure(PetscScalar dt, Vec* pP, Vec* pP_x, Vec*  pP_y, Vec* pP_z, Vec U_n, Vec V_n, Vec W_n, PetscScalar nx, PetscScalar ny, PetscScalar nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz)
{
    DM dmSol_centered, dmSol_shifted, dmSol_staggered;
    Mat A;
    Vec rhs;
    KSP ksp;
    PC  pc;

    CreateGrid(&dmSol_centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmSol_shifted, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmSol_staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    // Create Divergence vector centered from input (staggered)
    Vec Div;
    ComputeDivergence(&Div, U_n, V_n, W_n, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, dt); 

    /*Vec DivRef;
    CreateReferenceSolutionDivergence(dmSol_centered, &DivRef);
    CheckSolution(Div, DivRef);*/

    //VecScale(Div, 1/dt);

    Assemble_P(dmSol_centered, &A, &rhs, Div);

    PetscObjectDestroy((PetscObject*)&Div);


    AttachNullspace(dmSol_centered, A);

    Vec P;
    DMCreateGlobalVector(dmSol_centered, &P);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, A, A);
    //KSPGetPC(ksp, &pc);
    //PCSetType(pc, PCFIELDSPLIT);
    //PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, P);
    PetscObjectDestroy((PetscObject*)&A);
    PetscObjectDestroy((PetscObject*)&rhs);
    PetscObjectDestroy((PetscObject*)&ksp);

        
    //DMCreateGlobalVector(dmSol_centered, pP);
    VecCopy(P, *pP);
    Vec P_shifted;
    PetscObjectDestroy((PetscObject*)&P);

    DMCreateGlobalVector(dmSol_shifted, &P_shifted);
    DMStagMigrateVec(dmSol_centered, *pP, dmSol_shifted, P_shifted);

    
    Vec P_x, P_y, P_z;
    DMCreateGlobalVector(dmSol_shifted, &P_x);
    DMCreateGlobalVector(dmSol_shifted, &P_y);
    DMCreateGlobalVector(dmSol_shifted, &P_z);

    Derive_x_P(dmSol_shifted, &P_x, P_shifted);
    Derive_y_P(dmSol_shifted, &P_y, P_shifted);
    Derive_z_P(dmSol_shifted, &P_z, P_shifted);

    PetscObjectDestroy((PetscObject*)&P_shifted);

    /*DMCreateGlobalVector(dmSol_staggered, pP_x);
    DMCreateGlobalVector(dmSol_staggered, pP_y);
    DMCreateGlobalVector(dmSol_staggered, pP_z);*/

    DMStagMigrateVec(dmSol_shifted, P_x, dmSol_staggered, *pP_x);
    DMStagMigrateVec(dmSol_shifted, P_y, dmSol_staggered, *pP_y);
    DMStagMigrateVec(dmSol_shifted, P_z, dmSol_staggered, *pP_z);

    PetscObjectDestroy((PetscObject*)&P_x);
    PetscObjectDestroy((PetscObject*)&P_y);
    PetscObjectDestroy((PetscObject*)&P_z);

    PetscObjectDestroy((PetscObject*)&dmSol_centered);
    PetscObjectDestroy((PetscObject*)&dmSol_shifted);
    PetscObjectDestroy((PetscObject*)&dmSol_staggered);
    return 0;
}

// Update velocity

static PetscErrorCode UpdateVelocity(PetscScalar dt, Vec* pU_up, Vec* pV_up, Vec* pW_up, Vec P_x, Vec P_y, Vec P_z, Vec U_pre, Vec V_pre, Vec W_pre, PetscScalar nx, PetscScalar ny, PetscScalar nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz)
{
    VecAXPY(U_pre, -dt, P_x);
    VecAXPY(V_pre, -dt, P_y);
    VecAXPY(W_pre, -dt, P_z);

    DM dmSol_staggered;

    CreateGrid(&dmSol_staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    VecCopy(U_pre, *pU_up);
    VecCopy(V_pre, *pV_up);
    VecCopy(W_pre, *pW_up);

    PetscObjectDestroy((PetscObject*)&dmSol_staggered);


    return 0;
}

int main(int argc, char **argv)
{
    auto start = std::chrono::high_resolution_clock::now();

    PetscInitialize(&argc, &argv, (char*)0, "Help text for this program");
    // Get the rank of the current process

    static PetscScalar T   = 1.0;
    static PetscInt nt  = 10;
    static PetscInt nx  = 60;
    static PetscInt ny  = 60;
    static PetscInt nz  = 60;

    static PetscReal Lx_0   = 0;
    static PetscReal Ly_0  = 0;
    static PetscReal Lz_0  = 0;
    static PetscReal Lx = 1.0;
    static PetscReal Ly = 1.0;
    static PetscReal Lz = 1.0;
    static PetscReal Re   = 80.0;

    PetscScalar dt = 0.01; //T / nt;
    std::cout << "dt: " << dt << std::endl;

    PetscFunctionBeginUser;
    PetscInitialize(&argc, &argv, (char *)0, help);

    // Create necessary grids
    DM dmSol_Shifted, dmSol_Centered, dmSol_Staggered_x, dmSol_Staggered_y, dmSol_Staggered_z; //need to declare to due to solving linear system laplacian messing things up
    {
        CreateGrid(&dmSol_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Staggered_x, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMClone(dmSol_Staggered_x, &dmSol_Staggered_y);
        DMClone(dmSol_Staggered_x, &dmSol_Staggered_z);
    }

    //PetscObjectDestroy((PetscObject*)&bench);

    Vec U_0, V_0, W_0;
    CreateReferenceSolutionFirst(dmSol_Staggered_x, &U_0);
    CreateReferenceSolutionSecond(dmSol_Staggered_y, &V_0);
    CreateReferenceSolutionThird(dmSol_Staggered_z, &W_0);
    Vec U_int, V_int, W_int;
    DMCreateGlobalVector(dmSol_Staggered_x, &U_int);
    DMCreateGlobalVector(dmSol_Staggered_y, &V_int);
    DMCreateGlobalVector(dmSol_Staggered_z, &W_int);
    Vec U_pre, V_pre, W_pre;
    DMCreateGlobalVector(dmSol_Staggered_x, &U_pre);
    DMCreateGlobalVector(dmSol_Staggered_y, &V_pre);
    DMCreateGlobalVector(dmSol_Staggered_z, &W_pre);

    Vec P, P_x, P_y, P_z;
    DMCreateGlobalVector(dmSol_Centered, &P);
    DMCreateGlobalVector(dmSol_Staggered_x, &P_x);
    DMCreateGlobalVector(dmSol_Staggered_y, &P_y);
    DMCreateGlobalVector(dmSol_Staggered_z, &P_z);

    Vec U_up, V_up, W_up;
    DMCreateGlobalVector(dmSol_Staggered_x, &U_up);
    DMCreateGlobalVector(dmSol_Staggered_y, &V_up);
    DMCreateGlobalVector(dmSol_Staggered_z, &W_up);

    for(size_t i = 0; i < 300; i++){
        if (i == 0){
            ManageAdvection_x(dt, &U_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            ManageAdvection_y(dt, &V_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            ManageAdvection_z(dt, &W_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            std::cout << "Spirit of Kitty Hawk completed: advection done." << std::endl;
            /*Vec bench;
            DMCreateGlobalVector(dmSol_Staggered_z, &bench);
            CreateReferenceSolutionTry(dmSol_Staggered_z, &bench);
            CheckSolution(W_int, bench);*/
            ManageViscosity(dt, Re, &U_pre, &V_pre, &W_pre, U_int, V_int, W_int, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
            ManagePressure(dt, &P, &P_x, &P_y, &P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            std::cout << "Spirit of Oklahoma completed:   pressure done." << std::endl;
            /*PetscViewer viewer_Pi;
            DM da_solution_P;
            DMStagCreateCompatibleDMStag(dmSol_Staggered_z, 0, 0, 1, 0, &da_solution_P);
            Vec P_grid;
            DMStagVecSplitToDMDA(dmSol_Staggered_z, P_z, BACK, 0, &da_solution_P, &P_grid);
            PetscObjectSetName((PetscObject)P_grid, "p_simple_k");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "p_simple_k.vtr", FILE_MODE_WRITE, &viewer_Pi);
            VecView(P_grid, viewer_Pi);
            PetscViewerDestroy(&viewer_Pi); */
            UpdateVelocity(dt, &U_up, &V_up, &W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            std::cout << "Spirit of California completed: iteration " << i << " done" << std::endl;
            std::cout << "---------------------------------------------------------" << std::endl;

            PetscViewer viewer_u;
            DM DM_u;
            //DMStagCreateCompatibleDMStag(dmSol_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmSol_Staggered_x, U_up, LEFT, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "x_component");
            char filename_u[50]; 
            sprintf(filename_u, "results/x_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmSol_Staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u); 

            PetscViewer viewer_v;
            DM DM_v;
            //DMStagCreateCompatibleDMStag(dmSol_Staggered_y, 0, 0, 1, 0, &DM_v);
            Vec v;
            DMStagVecSplitToDMDA(dmSol_Staggered_y, V_up, DOWN, 0, &DM_v, &v);
            PetscObjectSetName((PetscObject)v, "y_component");
            char filename_v[50];
            sprintf(filename_v, "results/y_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmSol_Staggered_y), filename_v, FILE_MODE_WRITE, &viewer_v);
            VecView(v, viewer_v);
            VecDestroy(&v);
            DMDestroy(&DM_v);
            PetscViewerDestroy(&viewer_v);

            PetscViewer viewer_w;
            DM DM_w;
            //DMStagCreateCompatibleDMStag(dmSol_Staggered_z, 0, 0, 1, 0, &DM_w);
            Vec w;
            DMStagVecSplitToDMDA(dmSol_Staggered_z, W_up, BACK, 0, &DM_w, &w);
            PetscObjectSetName((PetscObject)w, "z_component");
            char filename_w[50];
            sprintf(filename_w, "results/z_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmSol_Staggered_z), filename_w, FILE_MODE_WRITE, &viewer_w);
            VecView(w, viewer_w);
            VecDestroy(&w);
            DMDestroy(&DM_w);
            PetscViewerDestroy(&viewer_w);

            PetscViewer viewer_p;
            DM DM_p;
            //DMStagCreateCompatibleDMStag(dmSol_Centered, 0, 0, 0, 1, &DM_p);
            Vec p;
            DMStagVecSplitToDMDA(dmSol_Centered, P, ELEMENT, 0, &DM_p, &p);
            PetscObjectSetName((PetscObject)p, "p");
            char filename_p[50];
            sprintf(filename_p, "results/p%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmSol_Centered), filename_p, FILE_MODE_WRITE, &viewer_p);
            VecView(p, viewer_p);
            VecDestroy(&p);
            DMDestroy(&DM_p);
            PetscViewerDestroy(&viewer_p);


        } else {
            ManageAdvection_x(dt, &U_int, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            ManageAdvection_y(dt, &V_int, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            ManageAdvection_z(dt, &W_int, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            std::cout << "Spirit of Kitty Hawk completed: advection done." << std::endl;
            ManageViscosity(dt, Re, &U_pre, &V_pre, &W_pre, U_int, V_int, W_int, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
            ManagePressure(dt, &P, &P_x, &P_y, &P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            std::cout << "Spirit of Oklahoma completed:   pressure done." << std::endl;
            UpdateVelocity(dt, &U_up, &V_up, &W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            std::cout << "Spirit of California completed: iteration " << i << " done" << std::endl;
            std::cout << "---------------------------------------------------------" << std::endl;

            PetscViewer viewer_u;
            DM DM_u;
            //DMStagCreateCompatibleDMStag(dmSol_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmSol_Staggered_x, U_up, LEFT, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "x_component");
            char filename_u[50]; 
            sprintf(filename_u, "results/x_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmSol_Staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u); 

            PetscViewer viewer_v;
            DM DM_v;
            //DMStagCreateCompatibleDMStag(dmSol_Staggered_y, 0, 0, 1, 0, &DM_v);
            Vec v;
            DMStagVecSplitToDMDA(dmSol_Staggered_y, V_up, DOWN, 0, &DM_v, &v);
            PetscObjectSetName((PetscObject)v, "y_component");
            char filename_v[50];
            sprintf(filename_v, "results/y_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmSol_Staggered_y), filename_v, FILE_MODE_WRITE, &viewer_v);
            VecView(v, viewer_v);
            VecDestroy(&v);
            DMDestroy(&DM_v);
            PetscViewerDestroy(&viewer_v);

            PetscViewer viewer_w;
            DM DM_w;
            //DMStagCreateCompatibleDMStag(dmSol_Staggered_z, 0, 0, 1, 0, &DM_w);
            Vec w;
            DMStagVecSplitToDMDA(dmSol_Staggered_z, W_up, BACK, 0, &DM_w, &w);
            PetscObjectSetName((PetscObject)w, "z_component");
            char filename_w[50];
            sprintf(filename_w, "results/z_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmSol_Staggered_z), filename_w, FILE_MODE_WRITE, &viewer_w);
            VecView(w, viewer_w);
            VecDestroy(&w);
            DMDestroy(&DM_w);
            PetscViewerDestroy(&viewer_w);

            PetscViewer viewer_p;
            DM DM_p;
            //DMStagCreateCompatibleDMStag(dmSol_Centered, 0, 0, 0, 1, &DM_p);
            Vec p;
            DMStagVecSplitToDMDA(dmSol_Centered, P, ELEMENT, 0, &DM_p, &p);
            PetscObjectSetName((PetscObject)p, "p");
            char filename_p[50];
            sprintf(filename_p, "results/p%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmSol_Centered), filename_p, FILE_MODE_WRITE, &viewer_p);
            VecView(p, viewer_p);
            VecDestroy(&p);
            DMDestroy(&DM_p);
            PetscViewerDestroy(&viewer_p);
   
        }
    }

    /*for(size_t i = 0; i < 3; i++){
        if (i == 0){
            ManageViscosity(dt, Re, &U_pre, &V_pre, &W_pre, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz); 
        } else {
            ManageViscosity(dt, Re, &U_pre, &V_pre, &W_pre, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);           
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

    PetscObjectDestroy((PetscObject*)&dmSol_Staggered_x);
    PetscObjectDestroy((PetscObject*)&dmSol_Staggered_y);
    PetscObjectDestroy((PetscObject*)&dmSol_Staggered_z);
    PetscObjectDestroy((PetscObject*)&dmSol_Centered);
    PetscObjectDestroy((PetscObject*)&dmSol_Shifted);    

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



