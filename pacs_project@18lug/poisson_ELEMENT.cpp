//
// Created by dave on 20/10/23.
//
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

static PetscErrorCode CreateReferenceSolution(DM, Vec *);
static PetscErrorCode CreateSystem(DM, Mat *, Vec *);
static PetscErrorCode AttachNullspace(DM, Mat);
static PetscErrorCode CheckSolution(Vec, Vec);

static PetscScalar pRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return exp(-10 * x * x -10 * y * y -10 * z * z);
    return  cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);


} /* zero integral */

static PetscScalar f(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return 20.0 * exp(-10 * x * x -10 * y * y -10 * z * z) * (-3 + 20 * x * x + 20 * y * y + 20 * z * z);
    return  -84*pi*pi*cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);

}


int main(int argc, char **argv)
{
    DM        dmSol;
    Vec       sol, solRef, rhs;
    Mat       A;
    KSP       ksp;
    PC        pc;

    PetscFunctionBeginUser;
    PetscInitialize(&argc, &argv, (char *)0, help);

    {
        const PetscInt dof0 = 0, dof1 = 0, dof2 = 0, dof3 = 1; /* 1 dof on each face and element center */
        const PetscInt stencilWidth = 1;
        DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, 100, 100, 100, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, &dmSol);
        DMSetFromOptions(dmSol);
        DMSetUp(dmSol);
        DMStagSetUniformCoordinatesExplicit(dmSol, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    }

    CreateReferenceSolution(dmSol, &solRef);

    CreateSystem(dmSol, &A, &rhs);

    AttachNullspace(dmSol, A);

    DMCreateGlobalVector(dmSol, &sol);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPFGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, sol);

    CheckSolution(sol, solRef);

    KSPDestroy(&ksp);
    VecDestroy(&sol);
    VecDestroy(&solRef);
    VecDestroy(&rhs);
    MatDestroy(&A);
    DMDestroy(&dmSol);
    PetscFinalize();
    return 0;
}

//compute numerical reference solution
static PetscErrorCode CreateReferenceSolution(DM dmSol, Vec *pSolRef)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        ip, icp[3];
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
    DMStagGetLocationSlot(dmSol, ELEMENT, 0, &ip);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icp[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmSol, &solRefLocal);
    DMStagVecGetArray(dmSol, solRefLocal, &arrSol);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                if (ex < start[0] + n[0] && ey < start[1] + n[1] && ez < start[2] + n[2]) arrSol[ez][ey][ex][ip] = pRef(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]], arrCoord[ez][ey][ex][icp[2]]);
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

static PetscErrorCode CreateSystem(DM dmSol, Mat *pA, Vec *pRhs) {
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
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hx = 1.0 / N[0];
    hy = 1.0 / N[1];
    hz = 1.0 / N[2];
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icp[d]);
    }

    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {


                    DMStagStencil row, col[7];
                    PetscScalar   valA[7], valRhs;
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
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]], arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

                }

                if (ex == 0) {
                    /* Left velocity Dirichlet */

/*
                    DMStagStencil     row;
                    PetscScalar       valRhs;
                    const PetscScalar valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = ELEMENT;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                    valRhs = pRef(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]], arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    */

                    DMStagStencil row, col[7];
                    PetscScalar   valA[7], valRhs;
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
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]], arrCoord[ez][ey][ex][icp[2]]);
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
                        if (ez == 0) {//ok
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

                        } else if (ez == N[2] - 1) {//ok
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
                        } else {// cambiato: originale era -2/hy*hy in valA[0]
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
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]], arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(rhs);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
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

static PetscErrorCode AttachNullspace(DM dmSol, Mat A)
{
    DM           dmPressure;
    Vec          constantPressure, basis;
    PetscReal    nrm;
    MatNullSpace matNullSpace;

    PetscFunctionBeginUser;
    DMStagCreateCompatibleDMStag(dmSol, 0, 0, 0, 1, &dmPressure);
    DMGetGlobalVector(dmPressure, &constantPressure);
    VecSet(constantPressure, 1.0);
    VecNorm(constantPressure, NORM_2, &nrm);
    VecScale(constantPressure, 1.0 / nrm);
    DMCreateGlobalVector(dmSol, &basis);
    DMStagMigrateVec(dmPressure, constantPressure, dmSol, basis);
    MatNullSpaceCreate(PetscObjectComm((PetscObject)dmSol), PETSC_FALSE, 1, &basis, &matNullSpace);
    VecDestroy(&basis);
    VecDestroy(&constantPressure);
    MatSetNullSpace(A, matNullSpace);
    MatNullSpaceDestroy(&matNullSpace);
    return 0;
}












