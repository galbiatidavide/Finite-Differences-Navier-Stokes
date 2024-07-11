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
static PetscErrorCode CreateReferenceRHS(DM, Vec *);

static PetscErrorCode CreateSystem(DM, Mat *, Vec *, Vec);
static PetscErrorCode CheckSolution(Vec, Vec);

PetscScalar a = pi/4;
PetscScalar d = pi/2;
static PetscScalar uxRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return cos(2.0*pi*x) * cos(4.0*pi*y) * sin(8.0*pi*z);
    //return -a*(exp(a*x)*sin(a*y+d*z) + exp(a*z)*cos(a*x+d*y));
    return cos(2.0*pi*x)cos(4.0*pi*y);
    
    

}


static PetscScalar fx(PetscScalar x, PetscScalar y, PetscScalar z)
{
    //return -84 * pi * pi * cos(2.0*pi*x) * cos(4.0*pi*y) * sin(8.0*pi*z) - 1000 * (cos(2.0*pi*x) * cos(4.0*pi*y) * sin(8.0*pi*z));
    //return -a*(d*d*exp(a*z)*cos(a*x+d*y) - a*a*exp(a*x)*sin(a*y+d*z)) - a*(-a*a*exp(a*z)*cos(a*x+d*y) - a*a*exp(a*x)*sin(a*y+d*z)) - a*(a*a*exp(a*z)*cos(a*x+d*y) - d*d*exp(a*x)*sin(a*y+d*z)) - 1000*(-a*(exp(a*x)*sin(a*y+d*z) + exp(a*z)*cos(a*x+d*y)));
    return -16.0 * pi * pi * cos(4.0*pi*y) - 1000 * cos(4.0*pi*y);
}





int main(int argc, char **argv)
{
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
        DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, 80, 80, 80, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, &dmSol);
        DMSetFromOptions(dmSol);
        DMSetUp(dmSol);
        DMStagSetUniformCoordinatesExplicit(dmSol, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        /*
        DMStagGetGhostCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz);
        std::cout<<"nx is"<<startx;
         */
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

    /* Check Solution */
    CheckSolution(sol, solRef);

    PetscViewer viewer_P;
    DM da_solution_P;
    DMStagCreateCompatibleDMStag(dmSol, 0, 0, 1, 0, &da_solution_P);
    Vec P_grid;
    DMStagVecSplitToDMDA(dmSol, sol, LEFT, 0, &da_solution_P, &P_grid);
    PetscObjectSetName((PetscObject)P_grid, "p_simple");
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "p_simple.vtr", FILE_MODE_WRITE, &viewer_P);
    VecView(P_grid, viewer_P);
    PetscViewerDestroy(&viewer_P); 

    /* Clean up and finalize PETSc */
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
    PetscInt        ip, iux, iuy, iuz, icp[3], icux[3], icuy[3], icuz[3];
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
//compute numerical reference solution
static PetscErrorCode CreateReferenceRHS(DM dmSol, Vec *pSolRef)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        ip, iux, iuy, iuz, icp[3], icux[3], icuy[3], icuz[3];
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
                arrSol[ez][ey][ex][iux] = fx(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]);
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
    PetscInt icux[3], icuy[3], icuz[3], icux_right[3], icuy_up[3], icuz_front[3], icux_back_left[3], icux_down_left[3], icux_front_left[3], icux_up_left[3];
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
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hx = 1.0 / N[0];
    hy = 1.0 / N[1];
    hz = 1.0 / N[2];
    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
        DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, &icux_back_left[d]);
        DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icux_down_left[d]);
        DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, &icux_front_left[d]);
        DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icux_up_left[d]);
    }

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,rhs_input,INSERT_VALUES,l);

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
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz) - Ret;
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
                            bc_1 = uxRef(arrCoord[ez][ey-1][ex][icux_back_left[0]], arrCoord[ez][ey-1][ex][icux_back_left[1]], arrCoord[ez][ey-1][ex][icux_back_left[2]]);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux_down_left[0]], arrCoord[ez][ey][ex][icux_down_left[1]], arrCoord[ez][ey][ex][icux_down_left[2]]);
                            valRhs = valRhs -2*bc_2/(hy*hy) - 2*bc_1/(hz*hz);
                            //ok
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            

                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz) - Ret;
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
                            bc_1 = uxRef(arrCoord[ez][ey-1][ex][icux_front_left[0]], arrCoord[ez][ey-1][ex][icux_front_left[1]], arrCoord[ez][ey-1][ex][icux_front_left[2]]);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux_down_left[0]], arrCoord[ez][ey][ex][icux_down_left[1]], arrCoord[ez][ey][ex][icux_down_left[2]]);
                            //ok
                            valRhs = valRhs -2*bc_2/(hy*hy) - 2*bc_1/(hz*hz);
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);

                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
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

                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux_down_left[0]], arrCoord[ez][ey][ex][icux_down_left[1]], arrCoord[ez][ey][ex][icux_down_left[2]]);
                            //ok

                            valRhs = valRhs -2*bc_2/(hy*hy);                            
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
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz) - Ret;
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
                            bc_1 = uxRef(arrCoord[ez][ey+1][ex][icux_back_left[0]], arrCoord[ez][ey+1][ex][icux_back_left[1]], arrCoord[ez][ey+1][ex][icux_back_left[2]]);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux_up_left[0]], arrCoord[ez][ey][ex][icux_up_left[1]], arrCoord[ez][ey][ex][icux_up_left[2]]);
                            //ok
                            valRhs = valRhs -2*bc_2/(hy*hy) - 2*bc_1/(hz*hz);                            
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            

                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 3.0 / (hz * hz) - Ret;
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
                            bc_1 = uxRef(arrCoord[ez+1][ey][ex][icux_front_left[0]], arrCoord[ez+1][ey][ex][icux_front_left[1]], arrCoord[ez+1][ey][ex][icux_front_left[2]]);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux_up_left[0]], arrCoord[ez][ey][ex][icux_up_left[1]], arrCoord[ez][ey][ex][icux_up_left[2]]);
                            //ok
                            
                            valRhs = valRhs -2*bc_2/(hy*hy) - 2*bc_1/(hz*hz);       

                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                     
                            /* Missing front term */
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -3.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
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
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux_up_left[0]], arrCoord[ez][ey][ex][icux_up_left[1]], arrCoord[ez][ey][ex][icux_up_left[2]]);
                            valRhs = valRhs -2*bc_2/(hy*hy);
                            //ok
                            DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
                        }
                    } else if (ez == 0) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 3.0 / (hz * hz) - Ret;
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
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux_back_left[0]], arrCoord[ez][ey][ex][icux_back_left[1]], arrCoord[ez][ey][ex][icux_back_left[2]]);
                        valRhs = valRhs  - 2*bc_1/(hz*hz);

                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
                    } else if (ez == N[2] - 1) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 3.0 / (hz * hz) - Ret;
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
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux_front_left[0]], arrCoord[ez][ey][ex][icux_front_left[1]], arrCoord[ez][ey][ex][icux_front_left[2]]);
                        valRhs = valRhs- 2*bc_1/(hz*hz);
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
                        DMStagVecSetValuesStencil(dmSol, rhs, 1, &row, &valRhs, INSERT_VALUES);
                    }
                    DMStagMatSetValuesStencil(dmSol, A, 1, &row, nEntries, col, valA, INSERT_VALUES);

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









