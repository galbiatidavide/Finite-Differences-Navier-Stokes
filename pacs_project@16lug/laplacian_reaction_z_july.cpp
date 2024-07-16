//
// Created by dave on 20/10/23.
//
static char help[] = "Solve a toy 3D problem on a staggered grid\n\n";

#include <petscdm.h>
#include <petscksp.h>
#include <petscdmstag.h>
#include <iostream>
#include <chrono>

#define pi 3.1415926535

//execute with mpirun --bind-to none --map-by slot -n 5 --oversubscribe ./lifex_test_laplacian_reaction_z_july


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
static PetscErrorCode CreateReferenceRHS(DM, Vec *);

static PetscErrorCode CreateSystem(DM, Mat *, Vec *, Vec);
static PetscErrorCode CheckSolution(Vec, Vec);




static PetscScalar uzRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    return sin(2.0*pi*x) * sin(4.0*pi*y) * sin(8.0*pi*z);
}



static PetscScalar fz(PetscScalar x, PetscScalar y, PetscScalar z)
{
    return -84 * pi * pi * sin(2.0*pi*x) * sin(4.0*pi*y) * sin(8.0*pi*z) - 1000 * (sin(2.0*pi*x) * sin(4.0*pi*y) * sin(8.0*pi*z));
    //return 0;
}




int main(int argc, char **argv)
{
    auto start = std::chrono::high_resolution_clock::now();

    DM        dmSol;
    Vec       sol, solRef, rhs, rhs_input;
    Mat       A;
    KSP       ksp;
    PC        pc;
    PetscInt startx, starty, startz, nx, ny, nz;


    PetscFunctionBeginUser;
    PetscInitialize(&argc, &argv, (char *)0, help);

    {
        const PetscInt dof0 = 0, dof1 = 0, dof2 = 1, dof3 = 0; /* 1 dof on each face and element center */
        const PetscInt stencilWidth = 1;
        DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, 140, 140, 140, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, &dmSol);
        DMSetFromOptions(dmSol);
        DMSetUp(dmSol);
        DMStagSetUniformCoordinatesExplicit(dmSol, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    }

    CreateReferenceSolution(dmSol, &solRef);
    CreateReferenceRHS(dmSol, &rhs_input);
    CreateSystem(dmSol, &A, &rhs, rhs_input);

    DMCreateGlobalVector(dmSol, &sol);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, sol);


    PetscViewer viewer_P;
    DM da_solution_P;
    DMStagCreateCompatibleDMStag(dmSol, 0, 0, 1, 0, &da_solution_P);
    Vec P_grid;
    DMStagVecSplitToDMDA(dmSol, solRef, BACK, 0, &da_solution_P, &P_grid);
    PetscObjectSetName((PetscObject)P_grid, "p_simple");
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "p_simple.vtr", FILE_MODE_WRITE, &viewer_P);
    VecView(P_grid, viewer_P);
    PetscViewerDestroy(&viewer_P);

    /* Check Solution */
    CheckSolution(sol, solRef);

    /* Clean up and finalize PETSc */
    KSPDestroy(&ksp);
    VecDestroy(&sol);
    VecDestroy(&rhs);
    MatDestroy(&A);
    VecDestroy(&rhs_input);
    VecDestroy(&solRef);

    DMDestroy(&dmSol);
    PetscFinalize();
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate duration
    std::chrono::duration<double> duration = end - start;

    // Output duration in seconds
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    return 0;
}

//compute numerical reference solution
static PetscErrorCode CreateReferenceSolution(DM dmSol, Vec *pSolRef)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iuz, icuz[3];
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
    DMStagGetLocationSlot(dmSol, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmSol, &solRefLocal);
    DMStagVecGetArray(dmSol, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
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

//compute numerical reference solution
static PetscErrorCode CreateReferenceRHS(DM dmSol, Vec *pSolRef)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iuz, icuz[3];
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
    DMStagGetLocationSlot(dmSol, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmSol, &solRefLocal);
    DMStagVecGetArray(dmSol, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrSol[ez][ey][ex][iuz] = fz(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
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

static PetscErrorCode CreateSystem(DM dmSol, Mat *pA, Vec *pRhs, Vec rhs_input) {
    Vec rhs, coordLocal;
    Mat A;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuz[3], icuz_front[3];
    PetscReal hx, hy, hz;
    PetscReal Ret = 1000;

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
    DMGlobalToLocalBegin(dmSol,rhs_input,INSERT_VALUES,l);

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
                            valRhs = valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);                            
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
                            valRhs = valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);
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
                            valRhs = valRhs - bc_1/(hx*hx);
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
                            valRhs = valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);       
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
                            valRhs = valRhs - bc_2/(hy*hy) - bc_1/(hx*hx);

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
                            valRhs = valRhs - bc_1/(hx*hx);

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
                        valRhs = valRhs - bc_2/(hy*hy);

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
                        valRhs = valRhs - bc_2/(hy*hy);

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
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES); 
                        
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);

                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);

    DMGlobalToLocalEnd(dmSol,rhs_input,INSERT_VALUES,l);
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









