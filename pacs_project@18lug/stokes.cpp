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
    return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
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
    return -0.001*(-1000*cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z) - 84*pi*pi*cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z));
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
                arrSol[ez][ey][ex][iuy] = solution(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);
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
    static PetscInt nx  = 80;
    static PetscInt ny  = 80;
    static PetscInt nz  = 80;

    static PetscReal Lx_0   = 0;
    static PetscReal Ly_0  = 0;
    static PetscReal Lz_0  = 0;
    static PetscReal Lx = 1.0;
    static PetscReal Ly = 1.0;
    static PetscReal Lz = 1.0;
    static PetscReal Re   = 200;

    PetscScalar dt = 0.1; //T / nt;
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


    /*for(size_t i = 0; i < 1; i++){
        if(i == 0){
            Vec U_0, V_0, W_0;
            //Vec bench;
            CreateReferenceSolutionFirst(dmSol_Staggered_x, &U_0);
            CreateReferenceSolutionSecond(dmSol_Staggered_y, &V_0);
            CreateReferenceSolutionThird(dmSol_Staggered_z, &W_0);
            Vec U_pre, V_pre, W_pre;
            DMCreateGlobalVector(dmSol_Staggered_x, &U_pre);
            DMCreateGlobalVector(dmSol_Staggered_y, &V_pre);
            DMCreateGlobalVector(dmSol_Staggered_z, &W_pre);
            //DMCreateGlobalVector(dmSol_Staggered_y, &bench);
            //CreateReferenceSolutionTry(dmSol_Staggered_y, &bench);
            ManageViscosity(dt, Re, &U_pre, &V_pre, &W_pre, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            PetscObjectDestroy((PetscObject*)&U_0);
            PetscObjectDestroy((PetscObject*)&V_0); 
            PetscObjectDestroy((PetscObject*)&W_0);
            //CheckSolution(V_pre, V_0); 
            Vec P, P_x, P_y, P_z;
            DMCreateGlobalVector(dmSol_Centered, &P);
            DMCreateGlobalVector(dmSol_Staggered_x, &P_x);
            DMCreateGlobalVector(dmSol_Staggered_y, &P_y);
            DMCreateGlobalVector(dmSol_Staggered_z, &P_z);
            ManagePressure(dt, &P, &P_x, &P_y, &P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            Vec U_pre_final;
            Vec V_pre_final;
            Vec W_pre_final;
            DMCreateGlobalVector(dmSol_Staggered_x, &U_pre_final);
            DMCreateGlobalVector(dmSol_Staggered_y, &V_pre_final);
            DMCreateGlobalVector(dmSol_Staggered_z, &W_pre_final);
            VecCopy(U_pre, U_pre_final);
            VecCopy(V_pre, V_pre_final);
            VecCopy(W_pre, W_pre_final);
            PetscObjectDestroy((PetscObject*)&U_pre);
            PetscObjectDestroy((PetscObject*)&V_pre);
            PetscObjectDestroy((PetscObject*)&W_pre);
            Vec U_up, V_up, W_up;
            DMCreateGlobalVector(dmSol_Staggered_x, &U_up);
            DMCreateGlobalVector(dmSol_Staggered_y, &V_up);
            DMCreateGlobalVector(dmSol_Staggered_z, &W_up);
            UpdateVelocity(dt, &U_up, &V_up, &W_up, P_x, P_y, P_z, U_pre_final, V_pre_final, W_pre_final, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            VecDestroy(&U_pre_final);
            VecDestroy(&V_pre_final);
            VecDestroy(&W_pre_final);
            VecDestroy(&P_x);
            VecDestroy(&P_y);
            VecDestroy(&P_z);
            VecDestroy(&P);
                                    PetscViewer viewer_P;
                    DM da_solution_P;
                    DMStagCreateCompatibleDMStag(dmSol_Staggered_x, 0, 0, 1, 0, &da_solution_P); 
                    Vec P_grid;
                    DMStagVecSplitToDMDA(dmSol_Staggered_x, W_up, BACK, 0, &da_solution_P, &P_grid);
                    PetscObjectSetName((PetscObject)P_grid, "ips");
                    PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "ips.vtr", FILE_MODE_WRITE, &viewer_P);            
                    VecView(P_grid, viewer_P);
                    PetscViewerDestroy(&viewer_P); 
                    std::cout<<"ciao";   
            VecDestroy(&U_up);
            VecDestroy(&V_up);
            VecDestroy(&W_up);

        } else {

            Vec U_pre, V_pre, W_pre;
            DMCreateGlobalVector(dmSol_Staggered_x, &U_pre);
            DMCreateGlobalVector(dmSol_Staggered_y, &V_pre);
            DMCreateGlobalVector(dmSol_Staggered_z, &W_pre);
            //DMCreateGlobalVector(dmSol_Staggered_y, &bench);
            //CreateReferenceSolutionTry(dmSol_Staggered_y, &bench);
            ManageViscosity(dt, Re, &U_pre, &V_pre, &W_pre, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            PetscObjectDestroy((PetscObject*)&U_up);
            PetscObjectDestroy((PetscObject*)&V_up); 
            PetscObjectDestroy((PetscObject*)&W_up);
            //CheckSolution(V_pre, V_0); 
            Vec P, P_x, P_y, P_z;
            DMCreateGlobalVector(dmSol_Centered, &P);
            DMCreateGlobalVector(dmSol_Staggered_x, &P_x);
            DMCreateGlobalVector(dmSol_Staggered_y, &P_y);
            DMCreateGlobalVector(dmSol_Staggered_z, &P_z);
            ManagePressure(dt, &P, &P_x, &P_y, &P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            Vec U_pre_final;
            Vec V_pre_final;
            Vec W_pre_final;
            DMCreateGlobalVector(dmSol_Staggered_x, &U_pre_final);
            DMCreateGlobalVector(dmSol_Staggered_y, &V_pre_final);
            DMCreateGlobalVector(dmSol_Staggered_z, &W_pre_final);
            VecCopy(U_pre, U_pre_final);
            VecCopy(V_pre, V_pre_final);
            VecCopy(W_pre, W_pre_final);
            PetscObjectDestroy((PetscObject*)&U_pre);
            PetscObjectDestroy((PetscObject*)&V_pre);
            PetscObjectDestroy((PetscObject*)&W_pre);
            Vec U_up, V_up, W_up;
            DMCreateGlobalVector(dmSol_Staggered_x, &U_up);
            DMCreateGlobalVector(dmSol_Staggered_y, &V_up);
            DMCreateGlobalVector(dmSol_Staggered_z, &W_up);
            UpdateVelocity(dt, &U_up, &V_up, &W_up, P_x, P_y, P_z, U_pre_final, V_pre_final, W_pre_final, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            VecDestroy(&U_pre_final);
            VecDestroy(&V_pre_final);
            VecDestroy(&W_pre_final);
            VecDestroy(&P_x);
            VecDestroy(&P_y);
            VecDestroy(&P_z);
            VecDestroy(&P);
 
            VecDestroy(&U_up);
            VecDestroy(&V_up);
            VecDestroy(&W_up);    
        }
    }*/

    //PetscObjectDestroy((PetscObject*)&bench);

    Vec U_0, V_0, W_0;
    CreateReferenceSolutionFirst(dmSol_Staggered_x, &U_0);
    CreateReferenceSolutionSecond(dmSol_Staggered_y, &V_0);
    CreateReferenceSolutionThird(dmSol_Staggered_z, &W_0);
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

    for(size_t i = 0; i < 3; i++){
        if (i == 0){
            ManageViscosity(dt, Re, &U_pre, &V_pre, &W_pre, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            ManagePressure(dt, &P, &P_x, &P_y, &P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            UpdateVelocity(dt, &U_up, &V_up, &W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            PetscViewer viewer_P;
            DM da_solution_P;
            DMStagCreateCompatibleDMStag(dmSol_Staggered_x, 0, 0, 1, 0, &da_solution_P); 
            Vec P_grid;
            DMStagVecSplitToDMDA(dmSol_Staggered_x, U_up, LEFT, 0, &da_solution_P, &P_grid);
            PetscObjectSetName((PetscObject)P_grid, "ips");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "ips.vtr", FILE_MODE_WRITE, &viewer_P);            
            VecView(P_grid, viewer_P);
            PetscViewerDestroy(&viewer_P); 
            std::cout<<"ciao";   
        } else {
            ManageViscosity(dt, Re, &U_pre, &V_pre, &W_pre, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            ManagePressure(dt, &P, &P_x, &P_y, &P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            UpdateVelocity(dt, &U_up, &V_up, &W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            PetscViewer viewer_P;
            DM da_solution_P;
            DMStagCreateCompatibleDMStag(dmSol_Staggered_x, 0, 0, 1, 0, &da_solution_P); 
            Vec P_grid;
            DMStagVecSplitToDMDA(dmSol_Staggered_x, U_up, LEFT, 0, &da_solution_P, &P_grid);
            PetscObjectSetName((PetscObject)P_grid, "ips2");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "ips2.vtr", FILE_MODE_WRITE, &viewer_P);            
            VecView(P_grid, viewer_P);
            PetscViewerDestroy(&viewer_P); 
            std::cout<<"ciao";          
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

    // Output the duration
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;


    return 0;
}




