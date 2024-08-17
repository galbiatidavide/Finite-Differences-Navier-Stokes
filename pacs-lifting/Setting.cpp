/*
annotazioni varie

creare template per le tre funzioni della sol analitica
creare unique ptrs per la griglia; CONCETTO: DICHIARO GRIGLIA E LA DEFINISCO PUNTANDOLA CON UNA VARIABILE DI TIPO UNIQUE PTR

*/

#include <iostream>
#include <cmath>
#include <limits>
#include <petscdmstag.h>
#include <petscksp.h>




PetscReal constexpr pi = 3.14159265358979323846;
PetscReal constexpr eps=std::numeric_limits<float>::max();

//Define the domain & grid
PetscInt constexpr nx{5};
PetscInt constexpr ny{6};
PetscInt constexpr nz{5};
PetscReal constexpr Lx_0{0};
PetscReal constexpr Ly_0{0};
PetscReal constexpr Lz_0{0};

PetscReal constexpr Lx{1.0};
PetscReal constexpr Ly{1.0};
PetscReal constexpr Lz{1.0};

//Define the time
PetscReal constexpr dt{0.1};
PetscReal constexpr iter{1};

//Define flow parameters
PetscReal constexpr A{4*sqrt(2)/(3*sqrt(3))};
PetscReal constexpr k{2*pi};
PetscReal constexpr c{(5/6)*pi};
PetscReal constexpr d{(1/6)*pi};
PetscReal constexpr vNorm{1e-2*(1/A)};
PetscReal constexpr vRef{vNorm*A};
//PetscReal constexpr vRef = 1.0;

PetscReal constexpr Re{1};
PetscReal constexpr LRef{1};
PetscReal constexpr nu{vRef*LRef/Re};
PetscReal theta{0};

PetscReal constexpr A_f = 1;
PetscReal constexpr B = 1;
PetscReal constexpr C = 1;



//Define macros (to change)
#define BACK             DMSTAG_BACK
#define DOWN             DMSTAG_DOWN
#define LEFT             DMSTAG_LEFT
#define ELEMENT          DMSTAG_ELEMENT

//Define analystical solution (also for bc's)

//inline effect

constexpr PetscReal uxRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    return vNorm*A*(sin(k*x - c)*cos(k*y - d)*sin(k*z) - cos(k*z - c)*sin(k*x - d)*sin(k*y))*exp(-theta);
    //return (A_f*sin(k*z) + C*cos(k*y))*exp(-theta);
}

constexpr PetscReal uyRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    return vNorm*A*(sin(k*y - c)*cos(k*z - d)*sin(k*x) - cos(k*x - c)*sin(k*y - d)*sin(k*z))*exp(-theta);
    //return (B*sin(k*x) + A_f*cos(k*z))*exp(-theta);

}

constexpr PetscReal uzRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    return vNorm*A*(sin(k*z - c)*cos(k*x - d)*sin(k*y) - cos(k*y - c)*sin(k*z - d)*sin(k*x))*exp(-theta);
    //return (C*sin(k*y) + B*cos(k*x))*exp(-theta);

}

constexpr PetscReal solution(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    return vNorm*A*(sin(k*x - c)*cos(k*y - d)*sin(k*z) - cos(k*z - c)*sin(k*x - d)*sin(k*y))*exp(-theta);
    //return (A_f*sin(k*z) + C*cos(k*y))*exp(-theta);

}

PetscErrorCode CreateAnalyticalU(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iux, icux[3];
    DM              dmCoord;
    Vec             vecLocal, coord, coordLocal;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, LEFT, 0, &icux[0]);
    DMStagGetLocationSlot(dmCoord, LEFT, 1, &icux[1]);
    DMStagGetLocationSlot(dmCoord, LEFT, 2, &icux[2]);  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iux] = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode CreateAnalyticalV(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iuy, icuy[3];
    Vec             vecLocal, coord, coordLocal;
    DM              dmCoord;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, DOWN, 0, &icuy[0]);
    DMStagGetLocationSlot(dmCoord, DOWN, 1, &icuy[1]);
    DMStagGetLocationSlot(dmCoord, DOWN, 2, &icuy[2]);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iuy] = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode CreateAnalyticalW(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iuz, icuz[3];
    Vec             vecLocal, coord, coordLocal;
    DM              dmCoord;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, BACK, 0, &icuz[0]);
    DMStagGetLocationSlot(dmCoord, BACK, 1, &icuz[1]);
    DMStagGetLocationSlot(dmCoord, BACK, 2, &icuz[2]);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iuz] = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}


PetscErrorCode CreateReferenceSolutionTry(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iux, icux[3], iuy, icuy[3], iuz, icuz[3], iue, icue[3];
    DM              dmCoord;
    Vec             vecLocal, coord, coordLocal;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, LEFT, 0, &icux[0]);
    DMStagGetLocationSlot(dmCoord, LEFT, 1, &icux[1]);
    DMStagGetLocationSlot(dmCoord, LEFT, 2, &icux[2]); 
    DMStagGetLocationSlot(dmCoord, DOWN, 0, &icuy[0]);
    DMStagGetLocationSlot(dmCoord, DOWN, 1, &icuy[1]);
    DMStagGetLocationSlot(dmCoord, DOWN, 2, &icuy[2]);
    DMStagGetLocationSlot(dmCoord, BACK, 0, &icuz[0]);
    DMStagGetLocationSlot(dmCoord, BACK, 1, &icuz[1]);
    DMStagGetLocationSlot(dmCoord, BACK, 2, &icuz[2]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 0, &icue[0]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 1, &icue[1]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 2, &icue[2]);     
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iue);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iux] = solution(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}

// Create Domain routines
PetscErrorCode CreateGrid(DM * dmGrid, PetscInt const & dof1, PetscInt const & dof2, PetscInt const & dof3, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz)
{
    const PetscInt dof0 = 0;
    const PetscInt stencilWidth = 1;

    PetscFunctionBegin;

    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, dmGrid);
    DMSetFromOptions(*dmGrid);
    DMSetUp(*dmGrid);
    DMStagSetUniformCoordinatesExplicit(*dmGrid, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    PetscFunctionReturn(0);
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
