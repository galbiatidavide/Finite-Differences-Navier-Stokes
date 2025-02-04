#include "parabolic.hpp"


const PetscErrorCode parabolic_problem_x::assemble_rhs(PetscReal const & theta)
{
    PetscFunctionBegin;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux[3], icux_right[3], iux_left, iux_right;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    PetscReal const hy = D_y /N[1];
    PetscReal const hz = D_z/ N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);

    VecScale(U_up, -Re/dt);

    Vec local;
    PetscReal ****arrVec;
    DMCreateLocalVector(dmGrid, &local);
    DMGlobalToLocalBegin(dmGrid, U_up, INSERT_VALUES, local);
    DMGlobalToLocalEnd(dmGrid, U_up, INSERT_VALUES, local);
    DMStagVecGetArrayRead(dmGrid, local, &arrVec);

    Vec local_rhs;
    PetscReal ****arrRhs;
    DMGetLocalVector(dmGrid, &local_rhs);
    VecSet(local_rhs, 0.0);
    DMStagVecGetArray(dmGrid, local_rhs, &arrRhs);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid, mask_U, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid, mask_U, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid, penalizationLocal, &arrPenalization);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {
                    if(arrPenalization[ez][ey][ex][iux_right] == 0) {
                        arrRhs[ez][ey][ex][iux_right] = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iux_right] = 0.0;
                    }
                }      
                if (ex == 0) {
                    if(arrPenalization[ez][ey][ex][iux_left] == 0) {
                        arrRhs[ez][ey][ex][iux_left] = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iux_left] = 0.0;
                    }

                } else {
                    PetscReal valRhs;
                    if (ey == 0) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]-hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz, theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hy*hy) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iux_left] = valRhs;                           
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] - hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]+hz, theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hy*hy) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iux_left] = valRhs;
                        } else {
                            PetscReal bc_2;
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] - hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_2/(hy*hy);       
                            arrRhs[ez][ey][ex][iux_left] = valRhs;                  
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]+hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz, theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hy*hy) - bc_2/(hz*hz); 
                            arrRhs[ez][ey][ex][iux_left] = valRhs;         
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] + hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]] + hz, theta);                            
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hy*hy) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iux_left] = valRhs;
                        } else {
                            PetscReal bc_2;
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]+hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_2/(hy*hy);  
                            arrRhs[ez][ey][ex][iux_left] = valRhs;                   
                        }
                    } else if (ez == 0) {
                        PetscReal bc_1;
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz, theta);
                        valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hz*hz);
                        arrRhs[ez][ey][ex][iux_left] = valRhs;                      
                    } else if (ez == N[2] - 1) {
                        PetscReal bc_1;
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]+hz, theta);
                        valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hz*hz);
                        arrRhs[ez][ey][ex][iux_left] = valRhs;                        
                    } else {
                        valRhs = arrVec[ez][ey][ex][iux_left];
                        arrRhs[ez][ey][ex][iux_left] = valRhs;
                    }
                }
                
            }
        }
    }    
    
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, local, &arrVec);
    DMRestoreLocalVector(dmGrid, &local);
    DMLocalToGlobal(dmGrid, local_rhs, INSERT_VALUES, rhs);
    DMStagVecRestoreArray(dmGrid, local_rhs, &arrRhs);
    DMRestoreLocalVector(dmGrid, &local_rhs);
    DMStagVecRestoreArrayRead(dmGrid, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalizationLocal);
    


    PetscFunctionReturn(0); 
}

const PetscErrorCode parabolic_problem_x::assemble_rhs(PetscReal const & theta, Vec const & U_up)
{
    PetscFunctionBegin;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux[3], icux_right[3], iux_left, iux_right;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    
    PetscReal const hy = D_y /N[1];
    PetscReal const hz = D_z/ N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);

    VecScale(U_up, -Re/dt);

    Vec local;
    PetscReal ****arrVec;
    DMCreateLocalVector(dmGrid, &local);
    DMGlobalToLocalBegin(dmGrid, U_up, INSERT_VALUES, local);
    DMGlobalToLocalEnd(dmGrid, U_up, INSERT_VALUES, local);
    DMStagVecGetArrayRead(dmGrid, local, &arrVec);

    Vec local_rhs;
    PetscReal ****arrRhs;
    DMGetLocalVector(dmGrid, &local_rhs);
    VecSet(local_rhs, 0.0);
    DMStagVecGetArray(dmGrid, local_rhs, &arrRhs);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid, mask_U, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid, mask_U, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid, penalizationLocal, &arrPenalization);


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {
                    if(arrPenalization[ez][ey][ex][iux_right] == 0) {
                        arrRhs[ez][ey][ex][iux_right] = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iux_right] = 0.0;
                    }
                }      
                if (ex == 0) {
                    if(arrPenalization[ez][ey][ex][iux_left] == 0) {
                        arrRhs[ez][ey][ex][iux_left] = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iux_left] = 0.0;
                    }

                } else {
                    PetscReal valRhs;
                    if (ey == 0) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]-hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz, theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hy*hy) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iux_left] = valRhs;                           
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] - hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]+hz, theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hy*hy) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iux_left] = valRhs;
                        } else {
                            PetscReal bc_2;
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] - hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_2/(hy*hy);       
                            arrRhs[ez][ey][ex][iux_left] = valRhs;                  
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]+hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz, theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hy*hy) - bc_2/(hz*hz); 
                            arrRhs[ez][ey][ex][iux_left] = valRhs;         
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]] + hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]] + hz, theta);                            
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hy*hy) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iux_left] = valRhs;
                        } else {
                            PetscReal bc_2;
                            bc_2 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]]+hy, arrCoord[ez][ey][ex][icux[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iux_left] - bc_2/(hy*hy);  
                            arrRhs[ez][ey][ex][iux_left] = valRhs;                   
                        }
                    } else if (ez == 0) {
                        PetscReal bc_1;
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]-hz, theta);
                        valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hz*hz);
                        arrRhs[ez][ey][ex][iux_left] = valRhs;                      
                    } else if (ez == N[2] - 1) {
                        PetscReal bc_1;
                        bc_1 = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]]+hz, theta);
                        valRhs = arrVec[ez][ey][ex][iux_left] - bc_1/(hz*hz);
                        arrRhs[ez][ey][ex][iux_left] = valRhs;                        
                    } else {
                        valRhs = arrVec[ez][ey][ex][iux_left];
                        arrRhs[ez][ey][ex][iux_left] = valRhs;
                    }
                }
                
            }
        }
    }    
    
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, local, &arrVec);
    DMRestoreLocalVector(dmGrid, &local);
    DMLocalToGlobal(dmGrid, local_rhs, INSERT_VALUES, rhs);
    DMStagVecRestoreArray(dmGrid, local_rhs, &arrRhs);
    DMRestoreLocalVector(dmGrid, &local_rhs);
    DMStagVecRestoreArrayRead(dmGrid, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalizationLocal);
    


    PetscFunctionReturn(0); 
}

PetscErrorCode parabolic_problem_x::exodus(size_t const & i)
{
    PetscFunctionBegin;
    PetscViewer viewer_u;
    DM DM_u;
    //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
    Vec u;
    DMStagVecSplitToDMDA(dmGrid, U_up, LEFT, 0, &DM_u, &u);
    PetscObjectSetName((PetscObject)u, "x_component");
    char filename_u[50]; 
    sprintf(filename_u, "%sx_component%03zu.vtr", base_path, i);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid), filename_u, FILE_MODE_WRITE, &viewer_u);
    VecView(u, viewer_u);
    VecDestroy(&u);
    DMDestroy(&DM_u);
    PetscViewerDestroy(&viewer_u);
    PetscFunctionReturn(0);
}

Vec parabolic_problem_x::get_U() { return U_up; }

void parabolic_problem_x::set_U(Vec const & U) { VecCopy(U, U_up); }

PetscErrorCode parabolic_problem_x::assemble_lhs() 
{
    PetscFunctionBegin;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icux[3], icux_right[3], iux_left;
    PetscReal const Ret = Re/dt;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;
    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = D_x /N[0];
    PetscReal const hy = D_y /N[1];
    PetscReal const hz = D_z/ N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    Vec penalization_local;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalization_local);
    DMGlobalToLocalBegin(dmGrid, mask_U, INSERT_VALUES, penalization_local);
    DMGlobalToLocalEnd(dmGrid, mask_U, INSERT_VALUES, penalization_local);
    DMStagVecGetArrayRead(dmGrid, penalization_local, &arrPenalization);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal val_penalization = arrPenalization[ez][ey][ex][iux_left];
                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    const PetscReal valA = 1.0;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                }      
                if (ex == 0) {
                    DMStagStencil row;
                    const PetscReal valA = 1.0;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = LEFT;
                    row.c = 0;

                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                } else {
                    DMStagStencil row, col[7];
                    PetscReal valA[7];
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
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
                            col[1].i = ex;
                            col[1].j = ey - 1;
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
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
                            col[1].i = ex;
                            col[1].j = ey - 1;
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
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
                            col[1].i = ex;
                            col[1].j = ey - 1;
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
                        }
                    } else if (ez == 0) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        col[5].k = ez + 1;
                        col[5].loc = LEFT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);                      
                    } else if (ez == N[2] - 1) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                    } else {
                        nEntries = 7;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                    }
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }
                
            }
        }
    }        
    
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, penalization_local, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalization_local);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    PetscFunctionReturn(0); 
}

PetscErrorCode parabolic_problem_x::solve_step(PetscReal const & theta)

{
    PetscFunctionBegin;   
    KSP       ksp;
    PC        pc;
    this->assemble_rhs(theta);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetTolerances(ksp, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, U_up);

    if(monitor_convergence) {
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        PetscInt iterations;
        KSPGetIterationNumber(ksp, &iterations);
        PetscReal residual_norm;
        KSPGetResidualNorm(ksp, &residual_norm);

        if (reason < 0) {
            PetscPrintf(PETSC_COMM_WORLD, "u-comp KSP did not converge. Reason: %s\n", KSPConvergedReasons[reason]);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, 
                        "u-comp KSP converged in %d iterations with a final residual norm of %g. Reason: %s\n", 
                        iterations, residual_norm, KSPConvergedReasons[reason]);
        }
    }

    KSPDestroy(&ksp);   
    PetscFunctionReturn(0);
}

PetscErrorCode parabolic_problem_x::solve_step(PetscReal const & theta, Vec const & U_up)
{
    PetscFunctionBegin;   
    KSP       ksp;
    PC        pc;
    this->assemble_rhs(theta, U_up);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetTolerances(ksp, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, U_up);
    if(monitor_convergence) {

        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        PetscInt iterations;
        KSPGetIterationNumber(ksp, &iterations);
        PetscReal residual_norm;
        KSPGetResidualNorm(ksp, &residual_norm);

        if (reason < 0) {
            PetscPrintf(PETSC_COMM_WORLD, "u-comp KSP did not converge. Reason: %s\n", KSPConvergedReasons[reason]);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, 
                        "u-comp KSP converged in %d iterations with a final residual norm of %g. Reason: %s\n", 
                        iterations, residual_norm, KSPConvergedReasons[reason]);
        }
    }


    KSPDestroy(&ksp);   
    PetscFunctionReturn(0);
}

PetscErrorCode parabolic_problem_x::solve()
{
    PetscFunctionBegin;
    assemble_lhs();
    for(size_t i = 0; i < iter; i++){            
        theta = (1/Re)*(i+1)*dt*pi*pi/3;
        this->solve_step(theta);
        std::cout << "Iteration: " << i << std::endl;
        Vec bench;
        DMCreateGlobalVector(dmGrid, &bench);
        CreateAnalyticalU(dmGrid, bench, theta);
        CheckSolution(U_up, bench, "U");
        VecDestroy(&bench);
        exodus(i);              
    }
    PetscFunctionReturn(0);
}


const PetscErrorCode parabolic_problem_y::assemble_rhs(PetscReal const & theta)
{
    PetscFunctionBegin;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up[3], iuy, iuy_up;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    PetscReal const hx = D_x /N[0];
    PetscReal const hz = D_z/ N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);        
    }

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);

    VecScale(V_up, -Re/dt);

    Vec local;
    PetscReal ****arrVec;
    DMCreateLocalVector(dmGrid, &local);
    DMGlobalToLocalBegin(dmGrid, V_up, INSERT_VALUES, local);
    DMGlobalToLocalEnd(dmGrid, V_up, INSERT_VALUES, local);
    DMStagVecGetArrayRead(dmGrid, local, &arrVec);

    Vec local_rhs;
    PetscReal ****arrRhs;
    DMGetLocalVector(dmGrid, &local_rhs);
    VecSet(local_rhs, 0.0);
    DMStagVecGetArray(dmGrid, local_rhs, &arrRhs);   

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid, mask_V, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid, mask_V, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid, penalizationLocal, &arrPenalization); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {  
                if (ey == N[1] - 1) {
                    if(arrPenalization[ez][ey][ex][iuy_up] == 0.0) {
                        arrRhs[ez][ey][ex][iuy_up] = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iuy_up] = 0.0;
                    }
                }                
                if (ey == 0) {
                    if(arrPenalization[ez][ey][ex][iuy] == 0.0){
                        arrRhs[ez][ey][ex][iuy] = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iuy] = 0.0;
                    }
                } else {
                    PetscReal valRhs;
                    if (ex == 0) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz, theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz, theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        } else {
                            PetscReal bc_1;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        }
                    } else if (ex == N[0] - 1) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz, theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz, theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        } else {
                            PetscReal bc_1;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        }
                    } else if (ez == 0) {
                        PetscReal bc_2;
                        bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz, theta);
                        valRhs = arrVec[ez][ey][ex][iuy] - bc_2/(hz*hz);
                        arrRhs[ez][ey][ex][iuy] = valRhs;

                    } else if (ez == N[2] - 1) {
                        PetscReal bc_1;
                        bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz, theta);
                        valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hz*hz);
                        arrRhs[ez][ey][ex][iuy] = valRhs;
                    } else {
                        valRhs = arrVec[ez][ey][ex][iuy];
                        arrRhs[ez][ey][ex][iuy] = valRhs;
                    }
                }                
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, local, &arrVec);
    DMRestoreLocalVector(dmGrid, &local);
    DMLocalToGlobal(dmGrid, local_rhs, INSERT_VALUES, rhs);
    DMStagVecRestoreArray(dmGrid, local_rhs, &arrRhs);
    DMRestoreLocalVector(dmGrid, &local_rhs);
    DMStagVecRestoreArrayRead(dmGrid, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalizationLocal);

    PetscFunctionReturn(0);  

}

const PetscErrorCode parabolic_problem_y::assemble_rhs(PetscReal const & theta, Vec const & V_up)
{
    PetscFunctionBegin;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up[3], iuy, iuy_up;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    PetscReal const hx = D_x /N[0];
    PetscReal const hz = D_z/ N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);        
    }

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);

    VecScale(V_up, -Re/dt);

    Vec local;
    PetscReal ****arrVec;
    DMCreateLocalVector(dmGrid, &local);
    DMGlobalToLocalBegin(dmGrid, V_up, INSERT_VALUES, local);
    DMGlobalToLocalEnd(dmGrid, V_up, INSERT_VALUES, local);
    DMStagVecGetArrayRead(dmGrid, local, &arrVec);

    Vec local_rhs;
    PetscReal ****arrRhs;
    DMGetLocalVector(dmGrid, &local_rhs);
    VecSet(local_rhs, 0.0);
    DMStagVecGetArray(dmGrid, local_rhs, &arrRhs);   

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid, mask_V, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid, mask_V, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid, penalizationLocal, &arrPenalization); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {  
                if (ey == N[1] - 1) {
                    if(arrPenalization[ez][ey][ex][iuy_up] == 0.0) {
                        arrRhs[ez][ey][ex][iuy_up] = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iuy_up] = 0.0;
                    }
                }                
                if (ey == 0) {
                    if(arrPenalization[ez][ey][ex][iuy] == 0.0){
                        arrRhs[ez][ey][ex][iuy] = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iuy] = 0.0;
                    }
                } else {
                    PetscReal valRhs;
                    if (ex == 0) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz, theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz, theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        } else {
                            PetscReal bc_1;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]-hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        }
                    } else if (ex == N[0] - 1) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz, theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz, theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx) - bc_2/(hz*hz);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        } else {
                            PetscReal bc_1;
                            bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]]+hx, arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuy] = valRhs;
                        }
                    } else if (ez == 0) {
                        PetscReal bc_2;
                        bc_2 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]-hz, theta);
                        valRhs = arrVec[ez][ey][ex][iuy] - bc_2/(hz*hz);
                        arrRhs[ez][ey][ex][iuy] = valRhs;

                    } else if (ez == N[2] - 1) {
                        PetscReal bc_1;
                        bc_1 = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]]+hz, theta);
                        valRhs = arrVec[ez][ey][ex][iuy] - bc_1/(hz*hz);
                        arrRhs[ez][ey][ex][iuy] = valRhs;
                    } else {
                        valRhs = arrVec[ez][ey][ex][iuy];
                        arrRhs[ez][ey][ex][iuy] = valRhs;
                    }
                }                
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, local, &arrVec);
    DMRestoreLocalVector(dmGrid, &local);
    DMLocalToGlobal(dmGrid, local_rhs, INSERT_VALUES, rhs);
    DMStagVecRestoreArray(dmGrid, local_rhs, &arrRhs);
    DMRestoreLocalVector(dmGrid, &local_rhs);
    DMStagVecRestoreArrayRead(dmGrid, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalizationLocal);

    PetscFunctionReturn(0);  

}

Vec parabolic_problem_y::get_V() { return V_up; }

void parabolic_problem_y::set_V(Vec const & V) { VecCopy(V, V_up);}

PetscErrorCode parabolic_problem_y::assemble_lhs()
{
    PetscFunctionBegin;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuy[3], icuy_up[3], iuy_down;
    PetscReal const Ret = Re/dt;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;
    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal const hx = D_x /N[0];
    PetscReal const hy = D_y /N[1];
    PetscReal const hz = D_z/ N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    Vec penalization_local;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalization_local);
    DMGlobalToLocalBegin(dmGrid, mask_V, INSERT_VALUES, penalization_local);
    DMGlobalToLocalEnd(dmGrid, mask_V, INSERT_VALUES, penalization_local);
    DMStagVecGetArrayRead(dmGrid, penalization_local, &arrPenalization);

    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);        
    }


    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) { 
                PetscReal val_penalization = arrPenalization[ez][ey][ex][iuy_down];
                    if (ey == N[1] - 1) {
                    DMStagStencil     row;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = UP;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                }                
                if (ey == 0) {
                    DMStagStencil     row;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = DOWN;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                } else {
                    DMStagStencil row, col[7];
                    PetscReal   valA[7];
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
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        }
                    } else if (ex == N[0] - 1) {
                        if (ez == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);                        
                        } else if (ez == N[2] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        }
                    } else if (ez == 0) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = DOWN;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        col[5].k   = ez + 1;
                        col[5].loc = DOWN;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = DOWN;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                    } else {
                        nEntries   = 7;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = DOWN;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                    }
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }                    
            }
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, penalization_local, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalization_local);
    
    PetscFunctionReturn(0);  
}

PetscErrorCode parabolic_problem_y::solve_step(PetscReal const & theta)
{
    PetscFunctionBegin;
    KSP       ksp;
    PC        pc;        
    this->assemble_rhs(theta);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetTolerances(ksp, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, V_up);

    if(monitor_convergence)
    {

        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        PetscInt iterations;
        KSPGetIterationNumber(ksp, &iterations);
        PetscReal residual_norm;
        KSPGetResidualNorm(ksp, &residual_norm);

        if (reason < 0) {
            PetscPrintf(PETSC_COMM_WORLD, "v-comp KSP did not converge. Reason: %s\n", KSPConvergedReasons[reason]);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, 
                        "v-comp KSP converged in %d iterations with a final residual norm of %g. Reason: %s\n", 
                        iterations, residual_norm, KSPConvergedReasons[reason]);
        }
    }
    

    KSPDestroy(&ksp);        
    PetscFunctionReturn(0);
}

PetscErrorCode parabolic_problem_y::solve_step(PetscReal const & theta, Vec const & V_up)
{
    PetscFunctionBegin;
    KSP       ksp;
    PC        pc;        
    this->assemble_rhs(theta, V_up);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetTolerances(ksp, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, V_up);    

    if(monitor_convergence)
    {
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        PetscInt iterations;
        KSPGetIterationNumber(ksp, &iterations);
        PetscReal residual_norm;
        KSPGetResidualNorm(ksp, &residual_norm);

        if (reason < 0) {
            PetscPrintf(PETSC_COMM_WORLD, "v-comp KSP did not converge. Reason: %s\n", KSPConvergedReasons[reason]);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, 
                        "v-comp KSP converged in %d iterations with a final residual norm of %g. Reason: %s\n", 
                        iterations, residual_norm, KSPConvergedReasons[reason]);
        }
    }
    KSPDestroy(&ksp);        
    PetscFunctionReturn(0);
}

PetscErrorCode parabolic_problem_y::solve()
{
    PetscFunctionBegin;
    assemble_lhs();

    for(size_t i = 0; i < iter; i++){               
        theta = (1/Re)*(i+1)*dt*pi*pi/3;
        this->solve_step(theta);
        Vec bench;
        DMCreateGlobalVector(dmGrid, &bench);
        CreateAnalyticalV(dmGrid, bench, theta);
        CheckSolution(V_up, bench, "V");
        VecDestroy(&bench);
        exodus(i);                   
    }
    PetscFunctionReturn(0);
}

PetscErrorCode parabolic_problem_y::exodus(size_t const & i)
{
    PetscFunctionBegin;
    PetscViewer viewer_v;
    DM DM_v;
    //DMStagCreateCompatibleDMStag(dmGrid_Staggered_y, 0, 0, 1, 0, &DM_v);
    Vec v;
    DMStagVecSplitToDMDA(dmGrid, V_up, DOWN, 0, &DM_v, &v);
    PetscObjectSetName((PetscObject)v, "y_component");
    char filename_v[50];
    sprintf(filename_v, "%sy_component%03zu.vtr", base_path, i);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid), filename_v, FILE_MODE_WRITE, &viewer_v);
    VecView(v, viewer_v);
    VecDestroy(&v);
    DMDestroy(&DM_v);
    PetscViewerDestroy(&viewer_v);
    PetscFunctionReturn(0);
}


PetscErrorCode const parabolic_problem_z::assemble_rhs(PetscReal const & theta) 
{
    PetscFunctionBegin;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuz[3], icuz_front[3], iuz, iuz_front;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;
    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);

    PetscReal hx = D_x /N[0];
    PetscReal hy = D_y /N[1];
    

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front);

    VecScale(W_up, -Re/dt);

    Vec local;
    PetscReal ****arrVec;
    DMCreateLocalVector(dmGrid, &local);
    DMGlobalToLocalBegin(dmGrid, W_up, INSERT_VALUES, local);
    DMGlobalToLocalEnd(dmGrid, W_up, INSERT_VALUES, local);
    DMStagVecGetArrayRead(dmGrid, local, &arrVec);

    Vec local_rhs;
    PetscReal ****arrRhs;
    DMGetLocalVector(dmGrid, &local_rhs);
    VecSet(local_rhs, 0.0);
    DMStagVecGetArray(dmGrid, local_rhs, &arrRhs);  

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid, mask_W, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid, mask_W, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid, penalizationLocal, &arrPenalization); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {              
                if (ez == N[2] - 1) {
                    if(arrPenalization[ez][ey][ex][iuz_front] == 0.0){
                        arrRhs[ez][ey][ex][iuz_front] = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iuz_front] = 0.0;
                    }
                }                
                if (ez == 0) {
                    if(arrPenalization[ez][ey][ex][iuz] == 0.0){
                        arrRhs[ez][ey][ex][iuz] = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iuz] = 0.0;
                    }
                } else {
                    PetscReal valRhs;
                    if (ex == 0) {
                        if (ey == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy) - bc_1/(hx*hx);                       
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        } else if (ey == N[1] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] + hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy) - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuz] = valRhs;                          
                        } else {
                            PetscReal bc_1;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        }
                    } else if (ex == N[0] - 1) {
                        if (ey == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] + hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy) - bc_1/(hx*hx);  
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        } else if (ey == N[1] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] + hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] + hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy) - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        } else {
                            PetscReal bc_1;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] + hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        }
                    } else if (ey == 0) {
                        PetscReal bc_2;
                        bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                        valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy);
                        arrRhs[ez][ey][ex][iuz] = valRhs;                         
                    } else if (ey == N[1] - 1) {
                        PetscReal bc_2;
                        bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] + hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                        valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy);
                        arrRhs[ez][ey][ex][iuz] = valRhs;
                    } else {
                        valRhs = arrVec[ez][ey][ex][iuz];
                        arrRhs[ez][ey][ex][iuz] = valRhs;                        
                    }
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, local, &arrVec);
    DMRestoreLocalVector(dmGrid, &local);
    DMLocalToGlobal(dmGrid, local_rhs, INSERT_VALUES, rhs);
    DMStagVecRestoreArray(dmGrid, local_rhs, &arrRhs);
    DMRestoreLocalVector(dmGrid, &local_rhs);
    DMStagVecRestoreArrayRead(dmGrid, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalizationLocal);

    PetscFunctionReturn(0); 
}

PetscErrorCode const parabolic_problem_z::assemble_rhs(PetscReal const & theta, Vec const & W_up)
{
    PetscFunctionBegin;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuz[3], icuz_front[3], iuz, iuz_front;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;
    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    
    PetscReal hx = D_x /N[0];
    PetscReal hy = D_y /N[1];
    

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front);

    VecScale(W_up, -Re/dt);

    Vec local;
    PetscReal ****arrVec;
    DMCreateLocalVector(dmGrid, &local);
    DMGlobalToLocalBegin(dmGrid, W_up, INSERT_VALUES, local);
    DMGlobalToLocalEnd(dmGrid, W_up, INSERT_VALUES, local);
    DMStagVecGetArrayRead(dmGrid, local, &arrVec);

    Vec local_rhs;
    PetscReal ****arrRhs;
    DMGetLocalVector(dmGrid, &local_rhs);
    VecSet(local_rhs, 0.0);
    DMStagVecGetArray(dmGrid, local_rhs, &arrRhs);  

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid, mask_W, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid, mask_W, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid, penalizationLocal, &arrPenalization); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {              
                if (ez == N[2] - 1) {
                    if(arrPenalization[ez][ey][ex][iuz_front] == 0.0){
                        arrRhs[ez][ey][ex][iuz_front] = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iuz_front] = 0.0;
                    }
                }                
                if (ez == 0) {
                    if(arrPenalization[ez][ey][ex][iuz] == 0.0){
                        arrRhs[ez][ey][ex][iuz] = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                    } else {
                        arrRhs[ez][ey][ex][iuz] = 0.0;
                    }
                } else {
                    PetscReal valRhs;
                    if (ex == 0) {
                        if (ey == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy) - bc_1/(hx*hx);                       
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        } else if (ey == N[1] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] + hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy) - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuz] = valRhs;                          
                        } else {
                            PetscReal bc_1;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] - hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        }
                    } else if (ex == N[0] - 1) {
                        if (ey == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] + hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy) - bc_1/(hx*hx);  
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        } else if (ey == N[1] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] + hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] + hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy) - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        } else {
                            PetscReal bc_1;
                            bc_1 = uzRef(arrCoord[ez][ey][ex][icuz[0]] + hx, arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
                            valRhs = arrVec[ez][ey][ex][iuz] - bc_1/(hx*hx);
                            arrRhs[ez][ey][ex][iuz] = valRhs;
                        }
                    } else if (ey == 0) {
                        PetscReal bc_2;
                        bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] - hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                        valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy);
                        arrRhs[ez][ey][ex][iuz] = valRhs;                         
                    } else if (ey == N[1] - 1) {
                        PetscReal bc_2;
                        bc_2 = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]] + hy, arrCoord[ez][ey][ex][icuz[2]], theta);
                        valRhs = arrVec[ez][ey][ex][iuz] - bc_2/(hy*hy);
                        arrRhs[ez][ey][ex][iuz] = valRhs;
                    } else {
                        valRhs = arrVec[ez][ey][ex][iuz];
                        arrRhs[ez][ey][ex][iuz] = valRhs;                        
                    }
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, local, &arrVec);
    DMRestoreLocalVector(dmGrid, &local);
    DMLocalToGlobal(dmGrid, local_rhs, INSERT_VALUES, rhs);
    DMStagVecRestoreArray(dmGrid, local_rhs, &arrRhs);
    DMRestoreLocalVector(dmGrid, &local_rhs);
    DMStagVecRestoreArrayRead(dmGrid, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalizationLocal);

    PetscFunctionReturn(0); 
}

Vec parabolic_problem_z::get_W() { return W_up; }

void parabolic_problem_z::set_W(Vec const & W) { VecCopy(W, W_up); }

PetscErrorCode parabolic_problem_z::assemble_lhs() 
{
    PetscFunctionBegin;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    PetscInt icuz[3], icuz_front[3], iuz_back;
    PetscReal const Ret = Re/dt;
    Vec coordLocal;
    DM dmCoord;
    PetscReal ****arrCoord;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    PetscReal hx = D_x /N[0];
    PetscReal hy = D_y /N[1];
    PetscReal hz = D_z/ N[2];

    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinatesLocal(dmGrid, &coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    Vec penalization_local;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid, &penalization_local);
    DMGlobalToLocalBegin(dmGrid, mask_W, INSERT_VALUES, penalization_local);
    DMGlobalToLocalEnd(dmGrid, mask_W, INSERT_VALUES, penalization_local);
    DMStagVecGetArrayRead(dmGrid, penalization_local, &arrPenalization);

    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }
    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {  

                PetscReal val_penalization = arrPenalization[ez][ey][ex][iuz_back];            
                if (ez == N[2] - 1) {
                    DMStagStencil     row;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = FRONT;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                }                
                if (ez == 0) {
                    DMStagStencil     row;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = BACK;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, 1, &row, &valA, INSERT_VALUES);
                } else {
                    DMStagStencil row, col[7];
                    PetscReal   valA[7];
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
                            valA[0]    = -2.0 / (hx * hx) - 2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
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
                        } else if (ey == N[1] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
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
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        }
                    } else if (ex == N[0] - 1) {
                        if (ey == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        } else if (ey == N[1] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
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
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                        }
                    } else if (ey == 0) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                    } else if (ey == N[1] - 1) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) - 2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
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
                    } else {
                        nEntries   = 7;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret - Re*val_penalization;
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
                    }
                    DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid, penalization_local, &arrPenalization);
    DMRestoreLocalVector(dmGrid, &penalization_local);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    PetscFunctionReturn(0); 
}

PetscErrorCode parabolic_problem_z::solve_step(PetscReal const & theta)
{
    PetscFunctionBegin;  
    KSP       ksp;
    PC        pc;
    this->assemble_rhs(theta);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetTolerances(ksp, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, W_up);

    if(monitor_convergence)
    {
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        PetscInt iterations;
        KSPGetIterationNumber(ksp, &iterations);
        PetscReal residual_norm;
        KSPGetResidualNorm(ksp, &residual_norm);

        if (reason < 0) {
            PetscPrintf(PETSC_COMM_WORLD, "w-comp KSP did not converge. Reason: %s\n", KSPConvergedReasons[reason]);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, 
                        "w-comp KSP converged in %d iterations with a final residual norm of %g. Reason: %s\n", 
                        iterations, residual_norm, KSPConvergedReasons[reason]);
        }
    }

    KSPDestroy(&ksp);  
    PetscFunctionReturn(0);
}

PetscErrorCode parabolic_problem_z::solve_step(PetscReal const & theta, Vec const & W_up)
{
    PetscFunctionBegin;  
    KSP       ksp;
    PC        pc;
    this->assemble_rhs(theta, W_up);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetTolerances(ksp, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs, W_up);

    if(monitor_convergence)
    {
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        PetscInt iterations;
        KSPGetIterationNumber(ksp, &iterations);
        PetscReal residual_norm;
        KSPGetResidualNorm(ksp, &residual_norm);

        if (reason < 0) {
            PetscPrintf(PETSC_COMM_WORLD, "w-comp KSP did not converge. Reason: %s\n", KSPConvergedReasons[reason]);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, 
                        "w-comp KSP converged in %d iterations with a final residual norm of %g. Reason: %s\n", 
                        iterations, residual_norm, KSPConvergedReasons[reason]);
        }
    }

    KSPDestroy(&ksp);  
    PetscFunctionReturn(0);
}

PetscErrorCode parabolic_problem_z::solve()
{
    PetscFunctionBegin;
    assemble_lhs();
    for(size_t i = 0; i < iter; i++){
            
            theta = (1/Re)*(i+1)*dt*pi*pi/3;
            this->solve_step(theta);

            Vec bench;
            DMCreateGlobalVector(dmGrid, &bench);
            CreateAnalyticalW(dmGrid, bench, theta);
            CheckSolution(W_up, bench, "W");
            VecDestroy(&bench);
            exodus(i);                 
    }
    PetscFunctionReturn(0);
}

PetscErrorCode parabolic_problem_z::exodus(size_t const & i)
{
    PetscFunctionBegin;
    PetscViewer viewer_w;
    DM DM_w;
    //DMStagCreateCompatibleDMStag(dmGrid_Staggered_z, 0, 0, 1, 0, &DM_w);
    Vec w;
    DMStagVecSplitToDMDA(dmGrid, W_up, BACK, 0, &DM_w, &w);
    PetscObjectSetName((PetscObject)w, "z_component");
    char filename_w[50];
    sprintf(filename_w, "%sz_component%03zu.vtr", base_path, i);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid), filename_w, FILE_MODE_WRITE, &viewer_w);
    VecView(w, viewer_w);
    VecDestroy(&w);
    DMDestroy(&DM_w);
    PetscViewerDestroy(&viewer_w);
    PetscFunctionReturn(0);
}
