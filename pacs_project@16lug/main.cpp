//
// Created by dave on 20/10/23.
//
static char help[] = "Solve a toy 3D problem on a staggered grid\n\n";

#include <petscdm.h>
#include <petscksp.h>
#include <petscdmstag.h>
#include <iostream>


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
    return exp(-10 * x * x -10 * y * y -10 * z * z);

} /* zero integral */

static PetscScalar f(PetscScalar x, PetscScalar y, PetscScalar z)
{
    return 20.0 * exp(-10 * x * x -10 * y * y -10 * z * z) * (-3 + 20 * x * x + 20 * y * y + 20 * z * z);
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
        DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, 50, 50, 50, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, &dmSol);
        DMSetFromOptions(dmSol);
        DMSetUp(dmSol);
        DMStagSetUniformCoordinatesExplicit(dmSol, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    }

    CreateReferenceSolution(dmSol, &solRef);

    std::cout<<"cambia"<<std::endl;

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

// Primo elemento dietro in basso a sinistra (Un solo elemento)
// Faccio la derivata solo su di me e sul punto dopo nelle tre direzioni
// Non considero proprio i punti -1 infatti 4 entries

                if(ex == 0 and ey == 0 and ez == 0)
                {
                    DMStagStencil row, col[4];
                    PetscScalar valA[4], valRhs;
                    PetscInt nEntries = 4;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[3].i = ez + 1;
                    col[3].j = ey;
                    col[3].k = ez;
                    col[3].loc = ELEMENT;
                    col[3].c = 0;
                    valA[3] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                // Elemento in alto a sx dietro
                // Vado in +1 su x e z e in -1 su y
                if(ex == 0 and ey == N[1] - 1 and ez == 0)
                {
                    DMStagStencil row, col[4];
                    PetscScalar valA[4], valRhs;
                    PetscInt nEntries = 4;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[3].i = ez + 1;
                    col[3].j = ey;
                    col[3].k = ez;
                    col[3].loc = ELEMENT;
                    col[3].c = 0;
                    valA[3] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //elemento in alto a sx davanti
                // -1 z e y e, +1 x
                if(ex == 0 and ey == N[1] - 1 and ez == N[2] - 1)
                {
                    DMStagStencil row, col[4];
                    PetscScalar valA[4], valRhs;
                    PetscInt nEntries = 4;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[3].i = ez - 1;
                    col[3].j = ey;
                    col[3].k = ez;
                    col[3].loc = ELEMENT;
                    col[3].c = 0;
                    valA[3] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //elemento in basso a sx davanti
                // +1 x e y, -1 z
                if(ex == 0 and ey == 0 and ez == N[2] - 1)
                {
                    DMStagStencil row, col[4];
                    PetscScalar valA[4], valRhs;
                    PetscInt nEntries = 4;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[3].i = ez - 1;
                    col[3].j = ey;
                    col[3].k = ez;
                    col[3].loc = ELEMENT;
                    col[3].c = 0;
                    valA[3] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //elemento dietro basso a dx
                // -1 x, +1 y e z
                if(ex == N[0] - 1 and ey == 0 and ez == 0)
                {
                    DMStagStencil row, col[4];
                    PetscScalar valA[4], valRhs;
                    PetscInt nEntries = 4;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[3].i = ez + 1;
                    col[3].j = ey;
                    col[3].k = ez;
                    col[3].loc = ELEMENT;
                    col[3].c = 0;
                    valA[3] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }


                //elemento in alto a sx dietro
                // -1 x e y, +1 z
                if(ex == N[0] - 1 and ey == N[1] - 1 and ez == 0)
                {
                    DMStagStencil row, col[4];
                    PetscScalar valA[4], valRhs;
                    PetscInt nEntries = 4;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[3].i = ez + 1;
                    col[3].j = ey;
                    col[3].k = ez;
                    col[3].loc = ELEMENT;
                    col[3].c = 0;
                    valA[3] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //elemento in alto davanti a dx
                // tutti meno
                if(ex == N[0] - 1 and ey == N[1] - 1 and ez == N[2] - 1)
                {
                    DMStagStencil row, col[4];
                    PetscScalar valA[4], valRhs;
                    PetscInt nEntries = 4;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[3].i = ez - 1;
                    col[3].j = ey;
                    col[3].k = ez;
                    col[3].loc = ELEMENT;
                    col[3].c = 0;
                    valA[3] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //elemento davanti in basso a sx
                // +1 y, -1 x e z
                if(ex == N[0] - 1 and ey == 0 and ez == N[2] - 1)
                {
                    DMStagStencil row, col[4];
                    PetscScalar valA[4], valRhs;
                    PetscInt nEntries = 4;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[3].i = ez - 1;
                    col[3].j = ey;
                    col[3].k = ez;
                    col[3].loc = ELEMENT;
                    col[3].c = 0;
                    valA[3] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //elementi in basso a sx non vertici
                //La z può fare tutto normale
                //x e y possono solo fare +1
                //quindi 5 ingressi: ijk, i+1, j+1, k+1, k-1
                if(ex == 0 and ey == 0 and ez != 0 and ez != N[2] - 1)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -2.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato verticale dietro a sx
                //La y può fare tutto normale
                //x e z possono solo fare +1
                //quindi 5 ingressi: ijk, i+1, k+1, j+1, j-1
                if(ex == 0 and ey != 0 and ey != N[1] - 1 and ez == 0)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato obliquo a sx i alto
                //La z può fare tutto normale
                //x può fare +1
                //y può fare -1
                //quindi 5 ingressi: ijk, i+1, j-1, k+1, k-1
                if(ex == 0 and ey == N[1] - 1 and ez != 0 and ez != N[2] - 1)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -2.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato verticale davanti a sx
                //La y può fare tutto normale
                //x può fare +1
                //z può fare -1
                //quindi 5 ingressi: ijk, i+1, j-1, j+1, k-1
                if(ex == 0 and ey != 0 and ey != N[1] - 1 and ez == N[2] - 1)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato obliquo a dx in basso
                //La z può fare tutto normale
                //x può fare -1
                //y può fare +1
                //quindi 5 ingressi: ijk, i-1, j+1, k+1, k-1
                if(ex == N[0] - 1 and ey == 0 and ez != 0 and ez != N[2] - 1)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -2.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato verticale a dx dietro
                //La y può fare tutto normale
                //x può fare -1
                //z può fare +1
                //quindi 5 ingressi: ijk, i-1, j-1, j+1, k+1
                if(ex == N[0] - 1 and ey != 0 and ey != N[1] - 1 and ez == 0)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato obliquo a dx in alto
                //La z può fare tutto normale
                //x può fare -1
                //y può fare -1
                //quindi 5 ingressi: ijk, i-1, j-1, k+1, k-1
                if(ex == N[0] - 1 and ey == N[1] - 1 and ez != 0 and ez != N[2] - 1)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) + -2.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato verticale a dx davanti
                //La y può fare tutto normale
                //x può fare -1
                //z può fare -1
                //quindi 5 ingressi: ijk, i-1, j-1, j+1, k-1
                if(ex == N[0] - 1 and ey != 0 and ey != N[1] - 1 and ez == N[2] - 1)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato orizzontale dietro in basso
                //La x può fare tutto normale
                //y può fare +1
                //z può fare +1
                //quindi 5 ingressi: ijk, i+1, i-1, j+1, k+1
                if(ex != 0 and ex != N[0] - 1 and ey == 0 and ez == 0)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato orizzontale in alto dietro
                //La x può fare tutto normale
                //y può fare -1
                //z può fare +1
                //quindi 5 ingressi: ijk, i+1, i-1, j-1, k+1
                if(ex != 0 and ex != N[0] - 1 and ey == N[1] - 1 and ez == 0)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato orizzotale in alto davanti
                //La x può fare tutto normale
                //y può fare -1
                //z può fare -1
                //quindi 5 ingressi: ijk, i+1, i-1, j-1, k-1
                if(ex != 0 and ex != N[0] - 1 and ey == N[1] - 1 and ez == N[2] - 1)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //Elementi lato orizzontale in basso davanti
                //La x può fare tutto normale
                //y può fare +1
                //z può fare -1
                //quindi 5 ingressi: ijk, i+1, i-1, j+1, k-1
                if(ex != 0 and ex != N[0] - 1 and ey == 0 and ez == N[2] - 1)
                {
                    DMStagStencil row, col[5];
                    PetscScalar valA[5], valRhs;
                    PetscInt nEntries = 5;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //faccia centrale a sx
                //y e z quello che vogliono
                //x solo +1
                //numero di ingressi 6: i+1, j+1, j-1, k+1, k-1
                if(ex == 0 and ey != 0 and ey != N[1] - 1 and ez != 0 and ez != N[2] - 1)
                {
                    DMStagStencil row, col[6];
                    PetscScalar valA[6], valRhs;
                    PetscInt nEntries = 6;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) + -2.0 / (hz * hz);
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
                    col[5].i = ex;
                    col[5].j = ey;
                    col[5].k = ez + 1;
                    col[5].loc = ELEMENT;
                    col[5].c = 0;
                    valA[5] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //faccia centrale a dx
                //y e z quello che vogliono
                //x solo -1
                //numero di ingressi 6: i-1, j+1, j-1, k+1, k-1
                if(ex == N[0] - 1 and ey != 0 and ey != N[1] - 1 and ez != 0 and ez != N[2] - 1)
                {
                    DMStagStencil row, col[6];
                    PetscScalar valA[6], valRhs;
                    PetscInt nEntries = 6;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) + -2.0 / (hz * hz);
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
                    col[5].i = ex;
                    col[5].j = ey;
                    col[5].k = ez + 1;
                    col[5].loc = ELEMENT;
                    col[5].c = 0;
                    valA[5] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //faccia centrale sotto
                //x e z quello che vogliono
                //y solo +1
                //numero di ingressi 6: i+1, i-1, j+1, k+1, k-1
                if(ey == 0 and ex != 0 and ex != N[0] - 1 and ez != 0 and ez != N[2] - 1)
                {
                    DMStagStencil row, col[6];
                    PetscScalar valA[6], valRhs;
                    PetscInt nEntries = 6;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) + -2.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //faccia centrale sopra
                //x e z quello che vogliono
                //y solo -1
                //numero di ingressi 6: i+1, i-1, j-1, k+1, k-1
                if(ey == N[1] - 1 and ex != 0 and ex != N[0] - 1 and ez != 0 and ez != N[2] - 1)
                {
                    DMStagStencil row, col[6];
                    PetscScalar valA[6], valRhs;
                    PetscInt nEntries = 6;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) + -2.0 / (hz * hz);
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //faccia centrale dietro
                //y e x quello che vogliono
                //z solo +1
                //numero di ingressi 6: i+1, i-1, j-1, j+1, k+1
                if(ez == 0 and ex != 0 and ex != N[0] - 1 and ey != 0 and ey != N[1] - 1)
                {
                    DMStagStencil row, col[6];
                    PetscScalar valA[6], valRhs;
                    PetscInt nEntries = 6;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[4].i = ex - 1;
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                //faccia centrale davanti
                //x e y quello che vogliono
                //z solo -1
                //numero di ingressi 6: i+1, i-1, j-1, j+1, k-1
                if(ez == N[2] - 1 and ex != 0 and ex != N[0] - 1 and ey != 0 and ey != N[1] - 1)
                {
                    DMStagStencil row, col[6];
                    PetscScalar valA[6], valRhs;
                    PetscInt nEntries = 6;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) + -1.0 / (hz * hz);
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
                    col[4].i = ex - 1;
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
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
                    DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                }

                if(ez != 0 and ez == N[2] - 1 and ex != 0 and ex != N[0] - 1 and ey != 0 and ey != N[1] - 1)
                {
                    //da qui e' quello tutto in mezzo completissimo
                    DMStagStencil row, col[7];
                    PetscScalar valA[7], valRhs;
                    PetscInt nEntries = 7;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    col[0].i = ex;
                    col[0].j = ey;
                    col[0].k = ez;
                    col[0].loc = ELEMENT;
                    col[0].c = 0;
                    valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) + -2.0 / (hz * hz);
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
                    col[6].i = ex;
                    col[6].j = ey;
                    col[6].k = ez + 1;
                    col[6].loc = ELEMENT;
                    col[6].c = 0;
                    valA[6] = 1.0 / (hz * hz);
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    valRhs = f(arrCoord[ez][ey][ex][icp[0]], arrCoord[ez][ey][ex][icp[1]],
                               arrCoord[ez][ey][ex][icp[2]]);
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
