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

PetscScalar a = pi/4;
PetscScalar d = pi/2;


static PetscScalar uxRef(PetscScalar x, PetscScalar y, PetscScalar z)

{
    //return  cos(pi*x)*sin(2*pi*y);
    //return  -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y));
    return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    //return sqrt(0.5)*x;
}

static PetscScalar square(PetscScalar x, PetscScalar y, PetscScalar z)

{
    //return  cos(pi*x)*sin(2*pi*y);
    //return  -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y));
    return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z)*cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    //return 0.5*x*x;
}

static PetscScalar uyRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return  -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z));
    return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);

}

static PetscScalar uzRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return  -a*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x));
    return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
}

static PetscScalar solution(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return  cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z) - 0.5*(-4*pi*cos(2*pi*x)*cos(4*pi*y)*cos(4*pi*y)*cos(8*pi*z)*cos(8*pi*z)*sin(2*pi*x) - 8*pi*cos(2*pi*x)*cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z)*cos(8*pi*z)*sin(4*pi*y) - 16*pi*cos(2*pi*x)*cos(2*pi*x)*cos(4*pi*y)*cos(4*pi*y)*cos(8*pi*z)*sin(8*pi*z));
    //return -16*pi*cos(2*pi*x)*cos(2*pi*x)*cos(4*pi*y)*cos(4*pi*y)*cos(8*pi*z)*sin(8*pi*z);
    //return -4*pi*cos(2*pi*x)*cos(4*pi*y)*cos(4*pi*y)*cos(8*pi*z)*cos(8*pi*z)*sin(2*pi*x) - 16*pi*cos(2*pi*x)*cos(2*pi*x)*cos(4*pi*y)*cos(4*pi*y)*cos(8*pi*z)*sin(8*pi*z);
    return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z) + 0.01*4*pi*cos(2*pi*x)*cos(4*pi*y)*cos(4*pi*y)*cos(8*pi*z)*cos(8*pi*z)*sin(2*pi*x) +  0.01*8*pi*cos(2*pi*x)*cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z)*cos(8*pi*z)*sin(4*pi*y) + 0.01*16*pi*cos(2*pi*x)*cos(2*pi*x)*cos(4*pi*y)*cos(4*pi*y)*cos(8*pi*z)*sin(8*pi*z);
    //return x;
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
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3], iue, icue[3];
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
    DMStagGetLocationSlot(dmSol, ELEMENT, 0, &iue);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icue[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmSol, &solRefLocal);
    DMStagVecGetArray(dmSol, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[1] + n[1] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iux] = solution(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]);

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


static PetscErrorCode CreateReferenceSquare(DM dmSol, Vec *pSolRef)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3], iue, icue[3];
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
    DMStagGetLocationSlot(dmSol, ELEMENT, 0, &iue);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icue[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmSol, &solRefLocal);
    DMStagVecGetArray(dmSol, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[1] + n[1] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iue] = square(arrCoord[ez][ey][ex][icue[0]], arrCoord[ez][ey][ex][icue[1]], arrCoord[ez][ey][ex][icue[2]]);

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
    DMCreateGlobalVector(dmSol, pUShifted);
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
    DMCreateGlobalVector(dmSol, pUShifted);
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
    DMCreateGlobalVector(dmSol, pVShifted);
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
    DMCreateGlobalVector(dmSol, pWShifted);
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
    DMCreateGlobalVector(dmSol, pAB_y);
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
    DMCreateGlobalVector(dmSol, pAB_z);
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
    DMCreateGlobalVector(dmSol, pVShifted);
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
    DMCreateGlobalVector(dmSol, pWShifted);
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
    DMCreateGlobalVector(dmSol, pAB_x);
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
    DMCreateGlobalVector(dmSol, pAB_z);
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
    DMCreateGlobalVector(dmSol, pAB_x);
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
    DMCreateGlobalVector(dmSol, pAB_y);
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
    DMCreateGlobalVector(dmSol, pUCenter);
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
    DMCreateGlobalVector(dmSol, pVCenter);
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
    DMCreateGlobalVector(dmSol, pWCenter);
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
    DMCreateGlobalVector(dmSol, pU2_x);
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
    DMCreateGlobalVector(dmSol, pV2_y);
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
    DMCreateGlobalVector(dmSol, pW2_z);
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
                        val_l = (val_left - uzRef(arrCoord[ez][ey][ex][icuz_e[0]], arrCoord[ez][ey][ex][icuz_e[1]], arrCoord[ez][ey][ex][icuz_e[2]] - hz)*uzRef(arrCoord[ez][ey][ex][icuz_e[0]], arrCoord[ez][ey][ex][icuz_e[1]],
                                            arrCoord[ez][ey][ex][icuz_e[2]]) - hz)/hz;

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
static PetscErrorCode ManageAdvection(PetscScalar dt, Vec* U_int, Vec* V_int, Vec* W_int, Vec U_n, Vec V_n, Vec W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz)
{

    // Create necessary grids
    DM dmSol_Shifted, dmSol_Centered, dmSol_Staggered, dmSol_cent;
    {
        CreateGrid(&dmSol_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Centered, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_cent, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

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
    Vec mixedFirst;
    DMCreateGlobalVector(dmSol_Staggered, &mixedFirst);
    {
        Vec U_y, U_z, V_y, W_z, UV, UW;

        FirstShiftU_y(dmSol_Shifted, &U_y, U_shift);
        FirstShiftU_z(dmSol_Shifted, &U_z, U_shift);
        FirstShiftV_y(dmSol_Shifted, &V_y, V_shift);
        FirstShiftW_z(dmSol_Shifted, &W_z, W_shift);

 

        DMCreateGlobalVector(dmSol_Shifted, &UV);
        DMCreateGlobalVector(dmSol_Shifted, &UW);

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


int main(int argc, char **argv)
{
    //std::cout<<size<<std::endl;

    static PetscScalar T   = 0.5;
    static PetscInt nt  = 3;
    static PetscInt nx  = 80;
    static PetscInt ny  = 80;
    static PetscInt nz  = 80;

    static PetscReal Lx_0  = 0;
    static PetscReal Ly_0  = 0;
    static PetscReal Lz_0  = 0;
    static PetscReal Lx = 1.0;
    static PetscReal Ly = 1.0;
    static PetscReal Lz = 1.0;

    PetscScalar dt = 0.01;

    PetscFunctionBeginUser;
    PetscInitialize(&argc, &argv, (char *)0, help);

    // Create necessary grids
    DM dmSol_Shifted, dmSol_Centered, dmSol_Staggered_x, dmSol_Staggered_y, dmSol_Staggered_z;
    {
        CreateGrid(&dmSol_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Staggered_x, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMClone(dmSol_Staggered_x, &dmSol_Staggered_y);
        DMClone(dmSol_Staggered_x, &dmSol_Staggered_z);
    }

    Vec U_0, V_0, W_0;
    
    CreateReferenceSolutionFirst(dmSol_Staggered_x, &U_0);
    CreateReferenceSolutionSecond(dmSol_Staggered_y, &V_0);
    CreateReferenceSolutionThird(dmSol_Staggered_z, &W_0);

    Vec U_int, V_int, W_int;
    DMCreateGlobalVector(dmSol_Staggered_x, &U_int);
    DMCreateGlobalVector(dmSol_Staggered_y, &V_int);
    DMCreateGlobalVector(dmSol_Staggered_z, &W_int);
    ManageAdvection(dt, &U_int, &V_int, &W_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    Vec my_bench;
    DMCreateGlobalVector(dmSol_Staggered_x, &my_bench);
    CreateReferenceSolutionTry(dmSol_Staggered_x, &my_bench);
    CheckSolution(U_int, my_bench);


    PetscObjectDestroy((PetscObject*)&dmSol_Staggered_x);
    PetscObjectDestroy((PetscObject*)&dmSol_Staggered_y);
    PetscObjectDestroy((PetscObject*)&dmSol_Staggered_z);
    PetscObjectDestroy((PetscObject*)&dmSol_Centered);
    PetscObjectDestroy((PetscObject*)&dmSol_Shifted);

    PetscObjectDestroy((PetscObject*)&U_0);
    PetscObjectDestroy((PetscObject*)&V_0); 
    PetscObjectDestroy((PetscObject*)&W_0);
    PetscObjectDestroy((PetscObject*)&my_bench);
    PetscObjectDestroy((PetscObject*)&U_int);
    PetscObjectDestroy((PetscObject*)&V_int);
    PetscObjectDestroy((PetscObject*)&W_int);
   

    PetscFinalize();
    return 0;
}




