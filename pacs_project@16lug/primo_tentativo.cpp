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

static PetscErrorCode CreateReferenceSolutionFirst(DM, Vec *);
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


static PetscScalar uxRef(PetscScalar x, PetscScalar y, PetscScalar z)

{
    return  sin(2*pi*x)*sin(4*pi*y)*sin(8*pi*z);
}

static PetscScalar uyRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    return  sin(5*pi*y)*sin(27*pi*x)*sin(2*pi*z);
}

static PetscScalar solution(PetscScalar x, PetscScalar y, PetscScalar z)
{
    return 4*pi*sin(2*pi*x)*sin(4*pi*y)*sin(8*pi*z)*sin(4*pi*y)*sin(8*pi*z)*cos(2*pi*x);

    //return 2*pi*sin(2*pi*x)*cos(2*pi*x)*sin(2*pi*z)*sin(2*pi*z)*(cos(2*pi*y)*cos(2*pi*y)-sin(2*pi*y)*sin(2*pi*y));
}


static PetscScalar uzRef(PetscScalar x, PetscScalar y, PetscScalar z)
{
    return  sin(2*pi*y)*sin(2*pi*x)*sin(2*pi*z);
}

static PetscErrorCode CreateReferenceSolution(DM dmSol, Vec *pVShifted)
{

    Vec VShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_right[3], icuy_up_left[3], icuy_up_right[3], icuy_down_left[3], icuy_down_right[3], icuy_e[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pVShifted);
    VShifted = *pVShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icuy_right[d]);
        DMStagGetLocationSlot(dmCoord, UP_LEFT, d, &icuy_up_left[d]); //
        DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, &icuy_up_right[d]);//
        DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, &icuy_down_left[d]);//
        DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, &icuy_down_right[d]);//
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icuy_e[d]);//


    }



    for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {


                DMStagStencil row;
                PetscScalar inter;
                row.i = ex;
                row.j = ey;
                row.k = ez;
                row.loc = LEFT;
                row.c = 0;
                inter = uxRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]],
                                 arrCoord[ez][ey][ex][icuy[2]]);
                DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);

                if (ex == N[0]-1){
                    DMStagStencil row;
                    PetscScalar inter;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    inter = uxRef(arrCoord[ez][ey][ex][icuy_right[0]], arrCoord[ez][ey][ex][icuy_right[1]],
                                  arrCoord[ez][ey][ex][icuy_right[2]]);
                    DMStagVecSetValuesStencil(dmSol, VShifted, 1, &row, &inter, INSERT_VALUES);
                }





            }
        }
    }




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



// Compute reference solution: initial input and b.c.
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
                if (ex < start[1] + n[1] && ey < start[2] + n[2]) arrSol[ez][ey][ex][iux] = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]);
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
                if (ex < start[0] + n[0] && ey < start[2] + n[2]) arrSol[ez][ey][ex][iuy] = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]);
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
                if (ex < start[0] + n[0] && ey < start[1] + n[1]) arrSol[ez][ey][ex][iuz] = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]]);
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

static PetscErrorCode DumpToVTK(DM dmSol_Staggered, DM dmSol_Centered, Vec U, Vec V, Vec W, Vec P) 
{
    
    PetscViewer viewer_U, viewer_V, viewer_W, viewer_P;

    DM da_solution_U;
    DM da_solution_V;
    DM da_solution_W;
    DM da_solution_P;

    DMStagCreateCompatibleDMStag(dmSol_Staggered, 0, 0, 1, 0, &da_solution_U); 
    DMStagCreateCompatibleDMStag(dmSol_Staggered, 0, 0, 1, 0, &da_solution_V); 
    DMStagCreateCompatibleDMStag(dmSol_Staggered, 0, 0, 1, 0, &da_solution_W); 
    DMStagCreateCompatibleDMStag(dmSol_Centered, 0, 0, 0, 1, &da_solution_P); 


    Vec U_grid, V_grid, W_grid, P_grid;
    DMStagVecSplitToDMDA(dmSol_Staggered, U, LEFT, 0, &da_solution_U, &U_grid);
    DMStagVecSplitToDMDA(dmSol_Staggered, V, DOWN, 0, &da_solution_V, &V_grid); 
    DMStagVecSplitToDMDA(dmSol_Staggered, W, BACK, 0, &da_solution_W, &W_grid); 
    DMStagVecSplitToDMDA(dmSol_Centered, P, ELEMENT, 0, &da_solution_P, &P_grid); 

 
    PetscObjectSetName((PetscObject)U_grid, "U");
    PetscObjectSetName((PetscObject)V_grid, "V");
    PetscObjectSetName((PetscObject)W_grid, "W");
    PetscObjectSetName((PetscObject)P_grid, "P");


    PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_U), "U.vtr", FILE_MODE_WRITE, &viewer_U);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_V), "V.vtr", FILE_MODE_WRITE, &viewer_V);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_W), "W.vtr", FILE_MODE_WRITE, &viewer_W);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "P.vtr", FILE_MODE_WRITE, &viewer_P);


    VecView(U_grid, viewer_U);
    VecView(V_grid, viewer_V);
    VecView(W_grid, viewer_W);
    VecView(P_grid, viewer_P);

    PetscViewerDestroy(&viewer_U);    
    PetscViewerDestroy(&viewer_V);    
    PetscViewerDestroy(&viewer_W);
    PetscViewerDestroy(&viewer_P);       

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

    return 0;
}

// First non-linear members: mixed
static PetscErrorCode FirstShiftU_y(DM dmSol, Vec *pUShifted, Vec solRef) 
{

    Vec UShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux[3], icux_right[3], icux_up_left[3], icux_up_right[3], icux_down_left[3], icux_down_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pUShifted);
    UShifted = *pUShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

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


    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(UShifted);
    VecAssemblyEnd(UShifted);

    return 0;
}

static PetscErrorCode FirstShiftU_z(DM dmSol, Vec *pUShifted, Vec solRef)
{

    Vec UShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pUShifted);
    UShifted = *pUShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

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


    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(UShifted);
    VecAssemblyEnd(UShifted);

    return 0;
}

static PetscErrorCode FirstShiftV_y(DM dmSol, Vec *pVShifted, Vec solRef) 
{

    Vec VShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_right[3], icuy_up_left[3], icuy_up_right[3], icuy_down_left[3], icuy_down_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pVShifted);
    VShifted = *pVShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

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


    DMGlobalToLocalEnd(dmSol, solRef, INSERT_VALUES, l);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VShifted);
    VecAssemblyEnd(VShifted);

    return 0;
}

static PetscErrorCode FirstShiftW_z(DM dmSol, Vec *pWShifted, Vec solRef) 
{

    Vec WShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_left[3], icux_back_right[3], icux_front_left[3], icux_front_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pWShifted);
    WShifted = *pWShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

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


    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WShifted);
    VecAssemblyEnd(WShifted);

    return 0;
}

static PetscErrorCode FirstDerive_y(DM dmSol, Vec *pAB_y, Vec AB)
{
    Vec             AB_y, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hy;
    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pAB_y);
    AB_y = *pAB_y;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if(N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm((PetscObject)dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hy = 1.0 / N[1];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
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
                }
            }
        }
    }

    VecAssemblyBegin(AB_y);
    VecAssemblyEnd(AB_y);

    return 0;
}

static PetscErrorCode FirstDerive_z(DM dmSol, Vec *pAB_z, Vec AB)
{
    Vec             AB_z, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hz;
    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pAB_z);
    AB_z = *pAB_z;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if(N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm((PetscObject)dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hz = 1.0 / N[2];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {


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

                } else {
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
                }
            }
        }
    }

    VecAssemblyBegin(AB_z);
    VecAssemblyEnd(AB_z);

    return 0;
}

// Second non-linear members: mixed
static PetscErrorCode SecondShiftV_z(DM dmSol, Vec *pVShifted, Vec solRef) 
{

    Vec VShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_right[3], icuy_front_up[3], icuy_front_down[3], icuy_back_up[3], icuy_back_down[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pVShifted);
    VShifted = *pVShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

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


    DMGlobalToLocalEnd(dmSol, solRef, INSERT_VALUES, l);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VShifted);
    VecAssemblyEnd(VShifted);

    return 0;
}

static PetscErrorCode SecondShiftW_z(DM dmSol, Vec *pWShifted, Vec solRef) 
{

    Vec WShifted, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux_right[3], icux_back_up[3], icux_back_down[3], icux_front_up[3], icux_front_down[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pWShifted);
    WShifted = *pWShifted;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

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


    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WShifted);
    VecAssemblyEnd(WShifted);

    return 0;
}

static PetscErrorCode SecondDerive_x(DM dmSol, Vec *pAB_x, Vec AB)
{
    Vec             AB_x, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hx;
    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pAB_x);
    AB_x = *pAB_x;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if(N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm((PetscObject)dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hx = 1.0 / N[0];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
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
                }
            }
        }
    }

    VecAssemblyBegin(AB_x);
    VecAssemblyEnd(AB_x);
    return 0;
}

static PetscErrorCode SecondDerive_z(DM dmSol, Vec *pAB_z, Vec AB)
{
    Vec             AB_z, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hz;
    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pAB_z);
    AB_z = *pAB_z;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if(N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm((PetscObject)dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hz = 1.0 / N[2];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {


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

                } else {
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
                }
            }
        }
    }

    VecAssemblyBegin(AB_z);
    VecAssemblyEnd(AB_z);

    return 0;
}

// Third non-linear members: mixed

static PetscErrorCode ThirdDerive_x(DM dmSol, Vec *pAB_x, Vec AB)
{
    Vec             AB_x, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hx;
    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pAB_x);
    AB_x = *pAB_x;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if(N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm((PetscObject)dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hx = 1.0 / N[0];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
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

                } else {
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
                }
            }
        }
    }

    VecAssemblyBegin(AB_x);
    VecAssemblyEnd(AB_x);
    return 0;
}

static PetscErrorCode ThirdDerive_y(DM dmSol, Vec *pAB_y, Vec AB)
{
    Vec             AB_y, AB_local;
    PetscInt        startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscReal       hy;
    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pAB_y);
    AB_y = *pAB_y;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if(N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm((PetscObject)dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hy = 1.0 / N[1];
    DMCreateLocalVector(dmSol, &AB_local);
    DMGlobalToLocalBegin(dmSol, AB, INSERT_VALUES, AB_local);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
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

                } else {
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
                }
            }
        }
    }

    VecAssemblyBegin(AB_y);
    VecAssemblyEnd(AB_y);

    return 0;
}


// non-linear members: homogeneous
static PetscErrorCode CenterU(DM dmSol, Vec *pUCenter, Vec solRef)
{

    Vec UCenter, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux[3], icux_right[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pUCenter);
    UCenter = *pUCenter;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

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


    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(UCenter);
    VecAssemblyEnd(UCenter);

    return 0;
}

static PetscErrorCode CenterV(DM dmSol, Vec *pVCenter, Vec solRef) 
{

    Vec VCenter, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pVCenter);
    VCenter = *pVCenter;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

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


    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(VCenter);
    VecAssemblyEnd(VCenter);

    return 0;
}

static PetscErrorCode CenterW(DM dmSol, Vec *pWCenter, Vec solRef) 
{

    Vec WCenter, coordLocal;

    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuz[3], icuz_front[3];
    DM dmCoord;
    PetscScalar ****arrCoord, ****array;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmSol, pWCenter);
    WCenter = *pWCenter;

    DMStagGetCorners(dmSol, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmSol, &N[0], &N[1], &N[2]);
    if (N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm(
            (PetscObject) dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");

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


    DMGlobalToLocalEnd(dmSol,solRef,INSERT_VALUES,l);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(WCenter);
    VecAssemblyEnd(WCenter);

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
    if(N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm((PetscObject)dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hx = 1.0 / N[0];

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

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
                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = LEFT;
                    row_left.c = 0;
                    val_left = 0;
                    DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row_left, &val_left, INSERT_VALUES);
                }

                if(ex == N[0] - 1){
                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = RIGHT;
                    row_right.c = 0;
                    val_right = 0;
                    DMStagVecSetValuesStencil(dmSol, U2_x, 1, &row_right, &val_right, INSERT_VALUES);
                }



            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(U2_x);
    VecAssemblyEnd(U2_x);

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
    if(N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm((PetscObject)dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hy = 1.0 / N[1];

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,V2,INSERT_VALUES,l);

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
                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = DOWN;
                    row_left.c = 0;
                    val_left = 0;
                    DMStagVecSetValuesStencil(dmSol, V2_y, 1, &row_left, &val_left, INSERT_VALUES);
                }

                if(ey == N[1] - 1){
                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = UP;
                    row_right.c = 0;
                    val_right = 0;
                    DMStagVecSetValuesStencil(dmSol, V2_y, 1, &row_right, &val_right, INSERT_VALUES);
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(V2_y);
    VecAssemblyEnd(V2_y);

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
    if(N[0] >= 2 && N[1] >= 2 && N[2] >= 2, PetscObjectComm((PetscObject)dmSol), PETSC_ERR_ARG_SIZ, "This example requires at least two elements in each dimensions");
    hz = 1.0 / N[2];

    DMGetCoordinateDM(dmSol, &dmCoord);
    DMGetCoordinatesLocal(dmSol, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    Vec l;
    DMCreateLocalVector(dmSol,&l);
    DMGlobalToLocalBegin(dmSol,W2,INSERT_VALUES,l);

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
                    DMStagStencil row_left;
                    PetscScalar val_left;
                    row_left.i = ex;
                    row_left.j = ey;
                    row_left.k = ez;
                    row_left.loc = BACK;
                    row_left.c = 0;
                    val_left = 0;
                    DMStagVecSetValuesStencil(dmSol, W2_z, 1, &row_left, &val_left, INSERT_VALUES);
                }

                if(ez == N[2] - 1){
                    DMStagStencil row_right;
                    PetscScalar val_right;
                    row_right.i = ex;
                    row_right.j = ey;
                    row_right.k = ez;
                    row_right.loc = FRONT;
                    row_right.c = 0;
                    val_right = 0;
                    DMStagVecSetValuesStencil(dmSol, W2_z, 1, &row_right, &val_right, INSERT_VALUES);
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    VecAssemblyBegin(W2_z);
    VecAssemblyEnd(W2_z);

    return 0;
}

// Assembling advection term
static PetscErrorCode ManageAdvection(PetscScalar dt, Vec* U_int, Vec* V_int, Vec* W_int, Vec U_n, Vec V_n, Vec W_n, PetscScalar nx, PetscScalar ny, PetscScalar nz, PetscScalar Lx_0, PetscScalar Lx, PetscScalar Ly_0, PetscScalar Ly, PetscScalar Lz_0, PetscScalar Lz)
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

    }

    DMStagMigrateVec(dmSol_Shifted, UV_y, dmSol_Staggered, mixedFirst);

    // Managing second mixed non-linear term, please DO NOT remove any comment
    Vec VU_x, VW_z;
    Vec mixedSecond;
    DMCreateGlobalVector(dmSol_Staggered, &mixedSecond);

    {
        Vec V_x, V_z, U_x, W_z, VU, VW;


        FirstShiftV_y(dmSol_Shifted, &V_x, V_shift);
        SecondShiftV_z(dmSol_Shifted, &V_z, V_shift);

        FirstShiftU_y(dmSol_Shifted, &U_x, U_shift);
        SecondShiftW_z(dmSol_Shifted, &W_z, W_shift);


        DMCreateGlobalVector(dmSol_Shifted, &VU);
        DMCreateGlobalVector(dmSol_Shifted, &VW);

        VecPointwiseMult(VU, V_x, U_x);
        VecPointwiseMult(VW, V_z, W_z);

        SecondDerive_x(dmSol_Shifted, &VU_x, VU);
        SecondDerive_z(dmSol_Shifted, &VW_z, VW);

        VecAXPY(VU_x, 1.0, VW_z);
    }
    DMStagMigrateVec(dmSol_Shifted, VU_x, dmSol_Staggered, mixedSecond);

    // Managing third mixed non-linear term, please DO NOT remove any comment, ESPECCIALLY the "==[] ones"
    Vec WU_x, WV_y;
    Vec mixedThird;
    DMCreateGlobalVector(dmSol_Staggered, &mixedThird);
    {
        Vec W_x, W_y, U_x, V_y, WU, WV;
        FirstShiftW_z(dmSol_Shifted, &W_x, W_shift);// ==FirsShiftW_z
        SecondShiftW_z(dmSol_Shifted, &W_y, W_shift);// ==SecondShiftW_z
        FirstShiftU_z(dmSol_Shifted, &U_x, U_shift);// ==FirsShiftU_z
        SecondShiftV_z(dmSol_Shifted, &V_y, V_shift);// ==SecondShiftV_z

        DMCreateGlobalVector(dmSol_Shifted, &WU);
        DMCreateGlobalVector(dmSol_Shifted, &WV);

        VecPointwiseMult(WU, W_x, U_x);
        VecPointwiseMult(WV, W_y, V_y);

        //DMCreateGlobalVector(dmSol_staggered, &WU_x);
        //DMCreateGlobalVector(dmSol_staggered, &WV_y);

        ThirdDerive_x(dmSol_Shifted, &WU_x, WU);
        ThirdDerive_y(dmSol_Shifted, &WV_y, WV);

        VecAXPY(WU_x, 1.0, WV_y);
    }
    DMStagMigrateVec(dmSol_Shifted, WU_x, dmSol_Staggered, mixedThird);

    // Managing homogenous components of non-linear term
    Vec U_c, V_c, W_c, U2_x, V2_y, W2_z;

    Vec U_center;
    DMCreateGlobalVector(dmSol_Centered, &U_center);
    DMStagMigrateVec(dmSol_Staggered, U_n, dmSol_Centered, U_center);

    Vec V_center;
    DMCreateGlobalVector(dmSol_Centered, &V_center);
    DMStagMigrateVec(dmSol_Staggered, V_n, dmSol_Centered, V_center);

    Vec W_center;
    DMCreateGlobalVector(dmSol_Centered, &W_center);
    DMStagMigrateVec(dmSol_Staggered, W_n, dmSol_Centered, W_center);

    CenterU(dmSol_Centered, &U_c, U_center);
    CenterV(dmSol_Centered, &V_c, V_center);
    CenterW(dmSol_Centered, &W_c, W_center);

    Derive_x(dmSol_Centered, &U2_x, U_c);
    Derive_y(dmSol_Centered, &V2_y, V_c);
    Derive_z(dmSol_Centered, &W2_z, W_c);

    Vec homoFirst, homoSecond, homoThird;
    DMCreateGlobalVector(dmSol_Staggered, &homoFirst);
    DMCreateGlobalVector(dmSol_Staggered, &homoSecond);
    DMCreateGlobalVector(dmSol_Staggered, &homoThird);

    DMStagMigrateVec(dmSol_Centered, U2_x, dmSol_Staggered, homoFirst);
    DMStagMigrateVec(dmSol_Centered, V2_y, dmSol_Staggered, homoSecond);
    DMStagMigrateVec(dmSol_Centered, W2_z, dmSol_Staggered, homoThird);

    Vec ref;
    CreateReferenceSolution(dmSol_Staggered, &ref);

    CheckSolution(homoFirst, ref);

    // DO NOT eliminate this comment
    /*
    Vec sol;
    CreateReferenceSolution_y(dmSol_Centered, &sol);
    CheckSolution(U2_x, sol);
    */

    // Assembling whole non-linear
    VecAXPBYPCZ(U_n, -dt, -dt, 1.0, homoFirst, mixedFirst);
    VecAXPBYPCZ(V_n, -dt, -dt, 1.0, homoSecond, mixedSecond);
    VecAXPBYPCZ(W_n, -dt, -dt, 1.0, homoThird, mixedThird);

    // Copy to output
    DMCreateGlobalVector(dmSol_Staggered, U_int);
    DMCreateGlobalVector(dmSol_Staggered, V_int);
    DMCreateGlobalVector(dmSol_Staggered, W_int);
    VecCopy(U_n, *U_int);
    VecCopy(V_n, *V_int);
    VecCopy(W_n, *W_int);

    return 0;
}



int main(int argc, char **argv)
{
    static PetscReal T   = 1.0;
    static PetscInt nt  = 1;
    static PetscInt nx  = 80;
    static PetscInt ny  = 80;
    static PetscInt nz  = 80;

    static PetscReal Lx_0   = 0.0;
    static PetscReal Ly_0  = 0.0;
    static PetscReal Lz_0  = 0.0;
    static PetscReal Lx = 1.0;
    static PetscReal Ly = 1.0;
    static PetscReal Lz = 1.0;
    static PetscReal tol  = 1e-4;
    static PetscReal Re   = 100000;

    PetscReal dt = T / nt;

    PetscFunctionBeginUser;
    PetscInitialize(&argc, &argv, (char *)0, help);

    // Create necessary grids
    DM dmSol_Shifted, dmSol_Centered, dmSol_Staggered;
    {
        CreateGrid(&dmSol_Shifted, 1, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmSol_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    }


    Vec U_0, V_0, W_0;
    Vec U_n, V_n, W_n;

    DMCreateGlobalVector(dmSol_Staggered, &U_n);
    DMCreateGlobalVector(dmSol_Staggered, &V_n);
    DMCreateGlobalVector(dmSol_Staggered, &W_n);

    Vec U_int, V_int, W_int;
    Vec U_pre, V_pre, W_pre;
    Vec P, P_x, P_y, P_z;
    Vec U_up, V_up, W_up;
    
    CreateReferenceSolutionFirst(dmSol_Staggered, &U_0);
    CreateReferenceSolutionSecond(dmSol_Staggered, &V_0);
    CreateReferenceSolutionThird(dmSol_Staggered, &W_0);



/*
    PetscViewer viewer;

    VecAssemblyBegin(U_n);
    VecAssemblyEnd(U_n);
    DM da_vel_avg;
    Vec vec_vel_avg;
    DMStagVecSplitToDMDA(dmSol_Staggered, U_0, LEFT, 0, &da_vel_avg, &vec_vel_avg); 
    PetscObjectSetName((PetscObject)vec_vel_avg, "Velocity new");

  {
    PetscViewer viewer;

    PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_vel_avg), "porco_diaz.vtr", FILE_MODE_WRITE, &viewer);
    VecView(vec_vel_avg, viewer);
    PetscViewerDestroy(&viewer);
    std::cout<<"I come";

  }

*/




/*
    //VecView(U_0, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "writing vector in binary to vector.dat ...\n");
    //PetscViewerBinaryOpen(PETSC_COMM_WORLD, "vector.dat", FILE_MODE_WRITE, &viewer);

    //PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output",&viewer);
    PetscViewerSetType(viewer, PETSCVIEWERVTK);    
    PetscViewerVTKOpen(PETSC_COMM_WORLD, "mat.output", FILE_MODE_WRITE, &viewer);
    DMDAVTKWriteAll_VTS(dmSol_Staggered, viewer);

    //VecView(U_0, viewer);
    PetscViewerDestroy(&viewer);
    //PetscOptionsSetValue(NULL, "-viewer_binary_mpiio", "");
/*
/*

    PetscViewer viewer;


    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERVTK);
  PetscViewerFileSetName(viewer, "sol.vtu");
 PetscViewerVTKOpen(PETSC_COMM_WORLD, "sol.vtu", FILE_MODE_WRITE, &viewer);

  VecView(U_n, viewer);
  PetscViewerDestroy(&viewer);
*/
    //for(size_t i = 1; i <= nt+1; i++)
   // {
        //if(i == 1)
        {
            /* Managing first step: find the intermediate velocity field by approximating the nonlinear term using a grid shifted on the edges and in the centers*/
            ManageAdvection(dt, &U_int, &V_int, &W_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

          //  VecCopy(U_pre, U_n);
           // VecCopy(V_pre, V_n);
           // VecCopy(W_pre, W_n);

            //std::cout<<"Step    "<<i<<"completed"<<std::endl;
        }
       // else 
        /*{
            ManageAdvection(dt, &U_int, &V_int, &W_int, U_n, V_n, W_n, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            ManageViscosity(dt, Re, &U_pre, &V_pre, &W_pre, U_int, V_int, W_int, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            ManagePressure(dt, &P, &P_x, &P_y, &P_z, U_n, V_n, W_n, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
            UpdateVelocity(dt, &U_up, &V_up, &W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

            VecCopy(U_pre, U_n);
            VecCopy(V_pre, V_n);
            VecCopy(W_pre, W_n);
            
            std::cout<<"Step    "<<i<<"completed"<<std::endl;
        }*/
   // }
    
/*
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    
    PetscViewerSetType(viewer, PETSCVIEWERVTK);

    PetscViewerVTKOpen(PETSC_COMM_WORLD,"prova.vtk",FILE_MODE_WRITE,&viewer);
    VecView(U_n,viewer);*/

    
    
    
/*
    PetscViewer viewer;

    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    
    PetscViewerSetType(viewer, PETSCVIEWERVTK);
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);
    PetscViewerFileSetName(viewer, "exodus.vtk");
    
    DMView(dmSol_Staggered, viewer);
    VecView(U_n, viewer);
    
*/
    PetscFinalize();
    return 0;
}











