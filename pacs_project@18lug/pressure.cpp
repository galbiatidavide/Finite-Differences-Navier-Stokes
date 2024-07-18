//
// Created by dave on 26/10/23.
//

/*WHAT TO BE DONE
1. Indentare tutto correttamente
3. Ragionare se alcune copie di vettori siano necessarie (potebbe essere molto dispendioso)
5. Capire come integrare questa roba qua in Lifex (ad esempio il Parsing dei parametri e i sistemi lineari vengono gestiti da Lifex di default)
6. Capire come integrare la vena (lol)

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

/*static PetscErrorCode CreateReferenceSolutionFirst(DM, Vec *);
static PetscErrorCode CreateReferenceSolutionSecond(DM, Vec *);
static PetscErrorCode CreateReferenceSolutionThird(DM, Vec *);

static PetscErrorCode FirstShiftU_y(DM, Vec*, Vec);
static PetscErrorCode FirstShiftU_z(DM, Vec*, Vec);
static PetscErrorCode FirstShiftV_y(DM, Vec*, Vec);
static PetscErrorCode FirstShiftW_z(DM, Vec*, Vec);
static PetscErrorCode FirstDerive_y(DM, Vec*, Vec);
static PetscErrorCode FirstDerive_z(DM, Vec*, Vec);
//first eq. checked

//static PetscErrorCode SecondShiftV_x(DM, Vec*, Vec); // ==FirstShiftV_y
static PetscErrorCode SecondShiftV_z(DM, Vec*, Vec); //ok
//static PetscErrorCode SecondShiftU_x(DM, Vec*, Vec); // ==FirstShiftU_y
static PetscErrorCode SecondShiftW_z(DM, Vec*, Vec);
static PetscErrorCode SecondDerive_x(DM, Vec*, Vec);
static PetscErrorCode SecondDerive_z(DM, Vec*, Vec);
//second eq. checked

//static PetscErrorCode ThirdShiftW_x(DM, Vec*, Vec);// ==FirsShiftW_z
//static PetscErrorCode ThirdShiftW_y(DM, Vec*, Vec);// ==SecondShiftW_z
//static PetscErrorCode ThirdShiftU_x(DM, Vec*, Vec);// ==FirsShiftU_z
//static PetscErrorCode ThirdShiftV_y(DM, Vec*, Vec);// ==SecondShiftV_z
static PetscErrorCode ThirdDerive_x(DM, Vec*, Vec);
static PetscErrorCode ThirdDerive_y(DM, Vec*, Vec);

static PetscErrorCode CenterU(DM, Vec*, Vec);
static PetscErrorCode CenterV(DM, Vec*, Vec);
static PetscErrorCode CenterW(DM, Vec*, Vec);

static PetscErrorCode Derive_x(DM dmSol, Vec*, Vec);
static PetscErrorCode Derive_y(DM dmSol, Vec*, Vec);
static PetscErrorCode Derive_z(DM dmSol, Vec*, Vec);


//static PetscErrorCode CreateReferenceSolution_y(DM dmSol, Vec *);
//static PetscErrorCode CheckSolution(Vec, Vec);
//static PetscErrorCode CreateReferenceSolution_new(DM, Vec*);
*/

static PetscScalar uxRef(PetscScalar x, PetscScalar y, PetscScalar z)

{
    return  cos(pi*x)*sin(2*pi*y);
    //return  0*x + 0*y;

}

static PetscScalar uyRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    return  -sin(pi*x)*cos(2*pi*y);
    //return  0*x + 0*y;
}

static PetscScalar uzRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    return  0;
}

static PetscScalar solution(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return 0.5*pi*sin(2*pi*x)*cos(2*pi*y)*cos(2*pi*y) + cos(pi*x)*sin(2*pi*y);
    //return 0.5*cos(2*pi*y)*(2*pi*sin(pi*x)*sin(pi*x)*sin(2*pi*y)-2*sin(pi*x) + pi*sin(2*pi*y));

    return  -2*pi*sin(2*pi*x)*cos(pi*y);
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


// Create Domain routines
static PetscErrorCode CreateGrid(DM* pdmSol, PetscInt dof1, PetscInt dof2, PetscInt dof3, PetscInt nx, PetscInt ny, PetscInt nz, PetscReal Lx_0, PetscReal Lx, PetscReal Ly_0, PetscReal Ly, PetscReal Lz_0, PetscReal Lz)
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



// Assembling Pressure

static PetscErrorCode Assemble_P(DM dmSol, Mat *pA, Vec *pRhs, Vec solRef) 
{
    Vec rhs, coordLocal;
    Mat A;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icp[3];
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
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icp[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,solRef,INSERT_VALUES,l);

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {


                    DMStagStencil row, col[6];
                    PetscScalar   valA[6], valRhs;
                    PetscInt      nEntries;
                    row.i   = ex;
                    row.j   = ey;
                    row.k   = ez;
                    row.loc = ELEMENT;
                    row.c   = 0;
                    if(ey == 0) {
                        if(ez == 0){
                            nEntries   = 4;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                        }
                        else if(ez == N[2] - 1) {
                            nEntries   = 4;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                        } else {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries   = 4;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                        }
                        else if (ez == N[2] - 1) {
                            nEntries   = 4;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                        } else {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries   = 5;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = ELEMENT;
                        col[0].c   = 0;
                        valA[0]    = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = ELEMENT;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = ELEMENT;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = ELEMENT;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries   = 5;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = ELEMENT;
                        col[0].c   = 0;
                        valA[0]    = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = ELEMENT;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = ELEMENT;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = ELEMENT;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                    } else {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = ELEMENT;
                        col[0].c   = 0;
                        valA[0]    = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = ELEMENT;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = ELEMENT;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = ELEMENT;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

                }

                else if (ex == 0) {
                    DMStagStencil row, col[6];
                    PetscScalar   valA[6], valRhs;
                    PetscInt      nEntries;
                    row.i   = ex;
                    row.j   = ey;
                    row.k   = ez;
                    row.loc = ELEMENT;
                    row.c   = 0;
                    if(ey == 0) {
                        if(ez == 0){
                            nEntries   = 4;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                        }
                        else if(ez == N[2] - 1) {
                            nEntries   = 4;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                        } else {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries   = 4;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                        }
                        else if (ez == N[2] - 1) {
                            nEntries   = 4;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                        } else {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries   = 5;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = ELEMENT;
                        col[0].c   = 0;
                        valA[0]    = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = ELEMENT;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = ELEMENT;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex + 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = ELEMENT;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries   = 5;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = ELEMENT;
                        col[0].c   = 0;
                        valA[0]    = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = ELEMENT;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = ELEMENT;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex + 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = ELEMENT;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                    } else {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = ELEMENT;
                        col[0].c   = 0;
                        valA[0]    = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = ELEMENT;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = ELEMENT;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex + 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = ELEMENT;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                } else {
                    /* X-momentum interior equation : (u_xx + u_yy + u_zz) - p_x = f^x */
                    DMStagStencil row, col[7];
                    PetscScalar   valA[7], valRhs;
                    PetscInt      nEntries;

                    row.i   = ex;
                    row.j   = ey;
                    row.k   = ez;
                    row.loc = ELEMENT;
                    row.c   = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            /* Missing down term */
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            /* Missing back term */
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);

                        } else if (ez == N[2] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            /* Missing down term */
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            /* Missing front term */
                        } else {//ok
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            /* Missing down term */
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {// cambiato: originale era -2/hy*hy in valA[0]
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            /* Missing up term */
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            /* Missing back entry */
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {// cambiato: originale era -2/hy*hy in valA[0]
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            /* Missing up term */
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            /* Missing front term */
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = ELEMENT;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = ELEMENT;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            /* Missing up term */
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = ELEMENT;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = ELEMENT;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = ELEMENT;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = ELEMENT;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = ELEMENT;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = ELEMENT;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = ELEMENT;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        /* Missing back term */
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = ELEMENT;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = ELEMENT;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = ELEMENT;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = ELEMENT;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = ELEMENT;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez - 1;
                        col[5].loc = ELEMENT;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        /* Missing front term */
                    } else {
                        nEntries   = 7;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = ELEMENT;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = ELEMENT;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = ELEMENT;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = ELEMENT;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = ELEMENT;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez - 1;
                        col[5].loc = ELEMENT;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        col[6].i   = ex;
                        col[6].j   = ey;
                        col[6].k   = ez + 1;
                        col[6].loc = ELEMENT;
                        col[6].c   = 0;
                        valA[6]    = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    DMStagVecGetValuesStencil(dmSol, l, 1, &row, &valRhs);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);

    VecAssemblyEnd(rhs);
    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    PetscObjectDestroy((PetscObject*)&l);
    //DMDestroy(&dmCoord);

    return 0;
}

static PetscErrorCode AssembleDivergence(DM dmSol, Vec *pDiv, Vec U_x, Vec U_y, Vec U_z) 
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

                    inter = (up - down) / hy + (right - left) / hx + (front - back) / hz;
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

static PetscErrorCode ComputeDivergence(Vec* pDiv, Vec U_n, Vec V_n, Vec W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz) 
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
        AssembleDivergence(dmSol_shifted, &Div, U_shifted, V_shifted, W_shifted);

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
    DMCreateGlobalVector(dmSol, pU2_x);
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

static PetscErrorCode ManagePressure(PetscScalar dt, Vec* pP, Vec* pP_x, Vec*  pP_y, Vec* pP_z, Vec U_n, Vec V_n, Vec W_n, PetscInt nx, PetscInt ny, PetscInt nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz)
{
    DM dmSol_centered, dmSol_shifted, dmSol_staggered;
    Mat A;
    Vec rhs;
    KSP ksp;
    //PC  pc;

    CreateGrid(&dmSol_centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmSol_shifted, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmSol_staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    // Create Divergence vector centered from input (staggered)
    Vec Div;
    ComputeDivergence(&Div, U_n, V_n, W_n, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz); 
    VecScale(Div, 1/dt);

    DMCreateGlobalVector(dmSol_centered, pP);

    Assemble_P(dmSol_centered, &A, &rhs, Div);

    Vec P;
    DMCreateGlobalVector(dmSol_centered, &P);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, A, A);
    KSPSetFromOptions(ksp);

    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);
    KSPSolve(ksp, rhs, P);        

    VecCopy(P, *pP);

    Vec P_shifted;
    DMCreateGlobalVector(dmSol_shifted, &P_shifted);
    DMStagMigrateVec(dmSol_centered, *pP, dmSol_shifted, P_shifted);    
    Vec P_x, P_y, P_z;

    Derive_x_P(dmSol_shifted, &P_x, P_shifted);
    Derive_y_P(dmSol_shifted, &P_y, P_shifted);
    Derive_z_P(dmSol_shifted, &P_z, P_shifted);

    DMCreateGlobalVector(dmSol_staggered, pP_x);
    DMCreateGlobalVector(dmSol_staggered, pP_y);
    DMCreateGlobalVector(dmSol_staggered, pP_z);

    DMStagMigrateVec(dmSol_shifted, P_x, dmSol_staggered, *pP_x);
    DMStagMigrateVec(dmSol_shifted, P_y, dmSol_staggered, *pP_y);
    DMStagMigrateVec(dmSol_shifted, P_z, dmSol_staggered, *pP_z);

    
    PetscObjectDestroy((PetscObject*)&Div);
    PetscObjectDestroy((PetscObject*)&dmSol_centered);
    PetscObjectDestroy((PetscObject*)&dmSol_shifted);
    PetscObjectDestroy((PetscObject*)&dmSol_staggered);
    PetscObjectDestroy((PetscObject*)&A);
    PetscObjectDestroy((PetscObject*)&rhs);
    PetscObjectDestroy((PetscObject*)&P);
    PetscObjectDestroy((PetscObject*)&ksp);
    PetscObjectDestroy((PetscObject*)&P_shifted);
    PetscObjectDestroy((PetscObject*)&P_x);
    PetscObjectDestroy((PetscObject*)&P_y);
    PetscObjectDestroy((PetscObject*)&P_z);


    return 0;
}






int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, (char*)0, "Help text for this program");
    // Get the rank of the current process
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    //std::cout<<rank<<std::endl;

    // Get the size of the communicator (total number of processes)
    int size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    //std::cout<<size<<std::endl;
    static PetscScalar T   = 0.05;
    static PetscInt nt  = 5;
    static PetscInt nx  = 50;
    static PetscInt ny  = 50;
    static PetscInt nz  = 50;

    static PetscReal Lx_0   = 0.0;
    static PetscReal Ly_0  = 0.0;
    static PetscReal Lz_0  = 0.0;
    static PetscReal Lx = 1.0;
    static PetscReal Ly = 1.0;
    static PetscReal Lz = 1.0;
    static PetscReal Re   = 10;

    PetscScalar dt = 0.01;

    PetscFunctionBeginUser;
    PetscInitialize(&argc, &argv, (char *)0, help);

    // Create necessary grids
    DM dmSol_Staggered;
    {
        CreateGrid(&dmSol_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);


    }


    Vec U_0, V_0, W_0;

  
    Vec P, P_x, P_y, P_z;
    
    CreateReferenceSolutionFirst(dmSol_Staggered, &U_0);
    CreateReferenceSolutionSecond(dmSol_Staggered, &V_0);
    CreateReferenceSolutionThird(dmSol_Staggered, &W_0);



    std::cout<<"Time is 0"<<std::endl;        


    ManagePressure(dt, &P, &P_x, &P_y, &P_z, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    //UpdateVelocity(dt, &U_up, &V_up, &W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    if(rank == 0){
    PetscObjectDestroy((PetscObject*)&U_0);
    PetscObjectDestroy((PetscObject*)&V_0); 
    PetscObjectDestroy((PetscObject*)&P);
    PetscObjectDestroy((PetscObject*)&P_x);
    PetscObjectDestroy((PetscObject*)&P_y);
    PetscObjectDestroy((PetscObject*)&P_z);}
    

    PetscObjectDestroy((PetscObject*)&dmSol_Staggered);
    PetscObjectDestroy((PetscObject*)&W_0);







    

    PetscFinalize();
    return 0;
}



