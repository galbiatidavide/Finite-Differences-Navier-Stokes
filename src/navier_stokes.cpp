#include "navier_stokes.hpp"

/*
PetscErrorCode const navier_stokes_problem::update_bc_U(PetscReal const & theta)
{
    PetscInt icux_left[3], icux_right[3], iux_left, iux_right;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec; 

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_staggered_x, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_staggered_x, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid_staggered_x, &dmCoord);

    DMGetCoordinates(dmGrid_staggered_x, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_staggered_x, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid_staggered_x, RIGHT, 0, &iux_right);

    DMCreateLocalVector(dmGrid_staggered_x, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_x, U_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_x, U_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArray(dmGrid_staggered_x, vecLocal, &arrVec);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid_staggered_x, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_x, mask_U, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_x, mask_U, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_x, penalizationLocal, &arrPenalization);


    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ex == N[0] - 1) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iux_right] == 0.0) {
                        val = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iux_right] = val;                    
                } else if(ex == 0) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iux_left] == 0.0) {
                        val = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]], arrCoord[ez][ey][ex][icux_left[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iux_left] = val;
                }

            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid_staggered_x, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid_staggered_x, vecLocal, INSERT_VALUES, U_up);
    DMRestoreLocalVector(dmGrid_staggered_x, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArrayRead(dmGrid_staggered_x, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid_staggered_x, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const navier_stokes_problem::update_bc_V(PetscReal const & theta) 
{

    PetscInt icuy_down[3], icuy_up[3], iuy_down, iuy_up;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_staggered_y, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_staggered_y, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid_staggered_y, &dmCoord);

    DMGetCoordinates(dmGrid_staggered_y, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy_down[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
    } 
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_staggered_y, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid_staggered_y, UP, 0, &iuy_up);

    DMCreateLocalVector(dmGrid_staggered_y, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_y, V_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_y, V_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_y, vecLocal, &arrVec);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid_staggered_y, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_y, mask_V, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_y, mask_V, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_y, penalizationLocal, &arrPenalization);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey == N[1] - 1) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iuy_up] == 0.0) {
                        val = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iuy_up] = val;
                } else if(ey == 0) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iuy_down] == 0.0) {
                        val = uyRef(arrCoord[ez][ey][ex][icuy_down[0]], arrCoord[ez][ey][ex][icuy_down[1]], arrCoord[ez][ey][ex][icuy_down[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iuy_down] = val;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid_staggered_y, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid_staggered_y, vecLocal, INSERT_VALUES, V_up);
    DMRestoreLocalVector(dmGrid_staggered_y, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArrayRead(dmGrid_staggered_y, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid_staggered_y, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const navier_stokes_problem::update_bc_W(PetscReal const & theta) 
{

    PetscInt icuz_back[3], icuz_front[3], iuz_back, iuz_front;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_staggered_z, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_staggered_z, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid_staggered_z, &dmCoord);

    DMGetCoordinates(dmGrid_staggered_z, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz_back[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_staggered_z, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid_staggered_z, FRONT, 0, &iuz_front); 
    
    DMCreateLocalVector(dmGrid_staggered_z, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_z, W_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_z, W_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_z, vecLocal, &arrVec);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid_staggered_z, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_z, mask_W, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_z, mask_W, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_z, penalizationLocal, &arrPenalization);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez == N[2] - 1) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iuz_front] == 0.0) {
                        val = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iuz_front] = val;
                } else if(ez == 0) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iuz_back] == 0.0) {
                        val = uzRef(arrCoord[ez][ey][ex][icuz_back[0]], arrCoord[ez][ey][ex][icuz_back[1]], arrCoord[ez][ey][ex][icuz_back[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iuz_back] = val;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid_staggered_z, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid_staggered_z, vecLocal, INSERT_VALUES, W_up);
    DMRestoreLocalVector(dmGrid_staggered_z, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArrayRead(dmGrid_staggered_z, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid_staggered_z, &penalizationLocal);

    PetscFunctionReturn(0);
}
*/


navier_stokes_problem::navier_stokes_problem(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, DM const & dmGrid_centered, DM const & dmGrid_cent_rich, DM const & dmGrid_stag_transp, DM const dmGrid_shift_transp, Vec const & U_up, Vec const & V_up, Vec const W_up)
    : dmGrid_staggered_x(dmGrid_staggered_x), dmGrid_staggered_y(dmGrid_staggered_y), dmGrid_staggered_z(dmGrid_staggered_z), dmGrid_centered(dmGrid_centered), dmGrid_cent_rich(dmGrid_cent_rich), dmGrid_stag_transp(dmGrid_stag_transp), dmGrid_shift_transp(dmGrid_shift_transp), U_up(U_up), V_up(V_up), W_up(W_up)

{
    DMCreateGlobalVector(dmGrid_centered, &P);
    DMCreateGlobalVector(dmGrid_staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_staggered_z, &P_z);
    DMCreateGlobalVector(dmGrid_centered, &Magnitude);

    DMCreateGlobalVector(dmGrid_staggered_x, &mask_U);
    DMCreateGlobalVector(dmGrid_staggered_y, &mask_V);
    DMCreateGlobalVector(dmGrid_staggered_z, &mask_W);

    if(brinkman)
    {
        createMaskU(dmGrid_staggered_x, mask_U, vertices, faces);
        createMaskV(dmGrid_staggered_y, mask_V, vertices, faces);
        createMaskW(dmGrid_staggered_z, mask_W, vertices, faces);
    }
    else {
        VecSet(mask_U, 0.0);
        VecSet(mask_V, 0.0);
        VecSet(mask_W, 0.0);
    }
}

navier_stokes_problem::navier_stokes_problem()
{
    //Allocate the grids
    CreateGrid(&dmGrid_staggered_x, 0, 1, 0);
    CreateGrid(&dmGrid_staggered_y, 0, 1, 0);
    CreateGrid(&dmGrid_staggered_z, 0, 1, 0);
    CreateGrid(&dmGrid_centered, 0, 0, 1);
    CreateGrid(&dmGrid_cent_rich, 0, 1, 1);
    CreateGrid(&dmGrid_shift_transp, 1, 1, 0);
    CreateGrid(&dmGrid_stag_transp, 0, 1, 0);

    //Create parallel vectors
    DMCreateGlobalVector(dmGrid_staggered_x, &U_up);
    DMCreateGlobalVector(dmGrid_staggered_y, &V_up);
    DMCreateGlobalVector(dmGrid_staggered_z, &W_up);
    CreateAnalyticalU(dmGrid_staggered_x, U_up, 0);
    CreateAnalyticalV(dmGrid_staggered_y, V_up, 0);
    CreateAnalyticalW(dmGrid_staggered_z, W_up, 0);

    DMCreateGlobalVector(dmGrid_centered, &P);
    DMCreateGlobalVector(dmGrid_staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_staggered_z, &P_z);
    DMCreateGlobalVector(dmGrid_centered, &Magnitude);

    /*DMCreateGlobalVector(dmGrid_staggered_x, &U_prova);
    DMCreateGlobalVector(dmGrid_staggered_y, &V_prova);
    DMCreateGlobalVector(dmGrid_staggered_z, &W_prova);*/

    DMCreateGlobalVector(dmGrid_staggered_x, &mask_U);
    DMCreateGlobalVector(dmGrid_staggered_y, &mask_V);
    DMCreateGlobalVector(dmGrid_staggered_z, &mask_W);

    if(brinkman)
    {
        createMaskU(dmGrid_staggered_x, mask_U, vertices, faces);
        createMaskV(dmGrid_staggered_y, mask_V, vertices, faces);
        createMaskW(dmGrid_staggered_z, mask_W, vertices, faces);
    }
    else {
        VecSet(mask_U, 0.0);
        VecSet(mask_V, 0.0);
        VecSet(mask_W, 0.0);
    }

};

PetscErrorCode const navier_stokes_problem::update_velocity(PetscReal const & theta)
{
    PetscFunctionBegin;

    VecAXPY(U_up, -dt, P_x);
    VecAXPY(V_up, -dt, P_y);
    VecAXPY(W_up, -dt, P_z);

    /*VecCopy(U_pre, U_up);
    VecCopy(V_pre, V_up);
    VecCopy(W_pre, W_up);*/

    /*update_bc_U(theta);
    update_bc_V(theta);
    update_bc_W(theta);*/


    PetscFunctionReturn(0);
}

PetscErrorCode const navier_stokes_problem::assemble_magnitude(Vec & Magnitude_Shifted, Vec const & U, Vec const & V, Vec const & W) 
{
    PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
    PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMGetCoordinateDM(dmGrid_cent_rich, &dmCoord);

    DMGetCoordinates(dmGrid_cent_rich, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_cent_rich, LEFT, 0, &iu_left);
    DMStagGetLocationSlot(dmGrid_cent_rich, RIGHT, 0, &iu_right);
    DMStagGetLocationSlot(dmGrid_cent_rich, UP, 0, &iu_up);
    DMStagGetLocationSlot(dmGrid_cent_rich, DOWN, 0, &iu_down);
    DMStagGetLocationSlot(dmGrid_cent_rich, FRONT, 0, &iu_front);
    DMStagGetLocationSlot(dmGrid_cent_rich, BACK, 0, &iu_back);
    DMStagGetLocationSlot(dmGrid_cent_rich, ELEMENT, 0, &iu_element);

    DMCreateLocalVector(dmGrid_cent_rich, &vecULocal);
    DMGlobalToLocalBegin(dmGrid_cent_rich, U, INSERT_VALUES, vecULocal);
    DMGlobalToLocalEnd(dmGrid_cent_rich, U, INSERT_VALUES, vecULocal);
    DMStagVecGetArrayRead(dmGrid_cent_rich, vecULocal, &arrU);

    DMCreateLocalVector(dmGrid_cent_rich, &vecVLocal);
    DMGlobalToLocalBegin(dmGrid_cent_rich, V, INSERT_VALUES, vecVLocal);
    DMGlobalToLocalEnd(dmGrid_cent_rich, V, INSERT_VALUES, vecVLocal);
    DMStagVecGetArrayRead(dmGrid_cent_rich, vecVLocal, &arrV);

    DMCreateLocalVector(dmGrid_cent_rich, &vecWLocal);
    DMGlobalToLocalBegin(dmGrid_cent_rich, W, INSERT_VALUES, vecWLocal);
    DMGlobalToLocalEnd(dmGrid_cent_rich, W, INSERT_VALUES, vecWLocal);
    DMStagVecGetArrayRead(dmGrid_cent_rich, vecWLocal, &arrW);    

    DMGetLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMStagVecGetArray(dmGrid_cent_rich, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal inter, left, right, up, down, front, back;
                left = arrU[ez][ey][ex][iu_left];
                right = arrU[ez][ey][ex][iu_right];
                up = arrV[ez][ey][ex][iu_up];
                down = arrV[ez][ey][ex][iu_down];
                front = arrW[ez][ey][ex][iu_front];
                back = arrW[ez][ey][ex][iu_back];
                inter = sqrt(((right + left)*(right + left)) / 4 + ((up + down)*(up + down)) / 4 + ((front + back)*(front + back)) / 4);
                arrOut[ez][ey][ex][iu_element] = inter;
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecULocal, &arrU);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecVLocal, &arrV);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecWLocal, &arrW);
    DMStagVecRestoreArray(dmGrid_cent_rich, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_cent_rich, vecOutLocal, INSERT_VALUES, Magnitude_Shifted);    
    DMRestoreLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecULocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecVLocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecWLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);    
}

PetscErrorCode const navier_stokes_problem::compute_magnitude()
{    
    PetscFunctionBegin;
    Vec U_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &U_shifted);
    DMStagMigrateVec(dmGrid_staggered_x, U_up, dmGrid_cent_rich, U_shifted);
    Vec V_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &V_shifted);
    DMStagMigrateVec(dmGrid_staggered_y, V_up, dmGrid_cent_rich, V_shifted);
    Vec W_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &W_shifted);
    DMStagMigrateVec(dmGrid_staggered_z, W_up, dmGrid_cent_rich, W_shifted);

    Vec magnitude_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &magnitude_shifted);
    assemble_magnitude(magnitude_shifted, U_shifted, V_shifted, W_shifted);
    DMStagMigrateVec(dmGrid_cent_rich, magnitude_shifted, dmGrid_centered, Magnitude);

    VecDestroy(&magnitude_shifted);
    VecDestroy(&U_shifted);  
    VecDestroy(&V_shifted);
    VecDestroy(&W_shifted);            
    PetscFunctionReturn(0); 
}

PetscErrorCode navier_stokes_problem::exodus(size_t i){

        PetscFunctionBegin;
        PetscViewer viewer_magnitude;
        DM DM_magnitude;
        //DMStagCreateCompatibleDMStag(dmGrid_staggered_x, 0, 0, 1, 0, &DM_u);
        Vec magnitude;
        DMStagVecSplitToDMDA(dmGrid_centered, Magnitude, ELEMENT, 0, &DM_magnitude, &magnitude);
        PetscObjectSetName((PetscObject)magnitude, "magnitude");
        char filename_magnitude[50]; 
        sprintf(filename_magnitude, "%smagnitude%03zu.vtr", base_path, i);

        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_magnitude, FILE_MODE_WRITE, &viewer_magnitude);
        VecView(magnitude, viewer_magnitude);
        VecDestroy(&magnitude);
        DMDestroy(&DM_magnitude);
        PetscViewerDestroy(&viewer_magnitude);

        PetscViewer viewer_u;
        DM DM_u;
        //DMStagCreateCompatibleDMStag(dmGrid_staggered_x, 0, 0, 1, 0, &DM_u);
        Vec u;
        DMStagVecSplitToDMDA(dmGrid_staggered_x, U_up, LEFT, 0, &DM_u, &u);
        PetscObjectSetName((PetscObject)u, "x_component");
        char filename_u[50]; 
        sprintf(filename_u, "%sx_component%03zu.vtr", base_path, i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
        VecView(u, viewer_u);
        VecDestroy(&u);
        DMDestroy(&DM_u);
        PetscViewerDestroy(&viewer_u); 

        PetscViewer viewer_v;
        DM DM_v;
        //DMStagCreateCompatibleDMStag(dmGrid_staggered_y, 0, 0, 1, 0, &DM_v);
        Vec v;
        DMStagVecSplitToDMDA(dmGrid_staggered_y, V_up, DOWN, 0, &DM_v, &v);
        PetscObjectSetName((PetscObject)v, "y_component");
        char filename_v[50];
        sprintf(filename_v, "%sy_component%03zu.vtr", base_path, i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_staggered_y), filename_v, FILE_MODE_WRITE, &viewer_v);
        VecView(v, viewer_v);
        VecDestroy(&v);
        DMDestroy(&DM_v);
        PetscViewerDestroy(&viewer_v);

        PetscViewer viewer_w;
        DM DM_w;
        //DMStagCreateCompatibleDMStag(dmGrid_staggered_z, 0, 0, 1, 0, &DM_w);
        Vec w;
        DMStagVecSplitToDMDA(dmGrid_staggered_z, W_up, BACK, 0, &DM_w, &w);
        PetscObjectSetName((PetscObject)w, "z_component");
        char filename_w[50];
        sprintf(filename_w, "%sz_component%03zu.vtr", base_path, i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_staggered_z), filename_w, FILE_MODE_WRITE, &viewer_w);
        VecView(w, viewer_w);
        VecDestroy(&w);
        DMDestroy(&DM_w);
        PetscViewerDestroy(&viewer_w);

        PetscViewer viewer_p;
        DM DM_p;
        //DMStagCreateCompatibleDMStag(dmGrid_centered, 0, 0, 0, 1, &DM_p);
        Vec p;
        DMStagVecSplitToDMDA(dmGrid_centered, P, ELEMENT, 0, &DM_p, &p);
        PetscObjectSetName((PetscObject)p, "p");
        char filename_p[50];
        sprintf(filename_p, "%sp%03zu.vtr", base_path, i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_p, FILE_MODE_WRITE, &viewer_p);
        VecView(p, viewer_p);
        VecDestroy(&p);
        DMDestroy(&DM_p);
        PetscViewerDestroy(&viewer_p);
        PetscFunctionReturn(0);
}

PetscErrorCode const navier_stokes_problem::solve()
{
    PetscFunctionBegin;

    PrintSimulationParameters();

    transport_problem_x transport_x(dmGrid_shift_transp, dmGrid_stag_transp, dmGrid_cent_rich);
    transport_problem_y transport_y(dmGrid_shift_transp, dmGrid_stag_transp, dmGrid_cent_rich);
    transport_problem_z transport_z(dmGrid_shift_transp, dmGrid_stag_transp, dmGrid_cent_rich);

    //this->assemble_lhs();

    parabolic_problem_x parabolic_x(dmGrid_staggered_x);
    parabolic_problem_y parabolic_y(dmGrid_staggered_y); 
    parabolic_problem_z parabolic_z(dmGrid_staggered_z);


    poisson_problem poisson(dmGrid_staggered_x, dmGrid_staggered_y, dmGrid_staggered_z, dmGrid_centered, dmGrid_cent_rich);
    
    for(size_t i = 0; i < iter; ++i){

        theta = i * dt;

        transport_x.solve_step_x(theta, U_up, V_up, W_up);
        transport_y.solve_step_y(theta, U_up, V_up, W_up);
        transport_z.solve_step_z(theta, U_up, V_up, W_up);

        parabolic_x.solve_step(theta, U_up);
        parabolic_y.solve_step(theta, V_up);
        parabolic_z.solve_step(theta, W_up);

        poisson.manage_pressure(U_up, V_up, W_up, P);
        poisson.manage_pressure_x(P, P_x);
        poisson.manage_pressure_y(P, P_y);
        poisson.manage_pressure_z(P, P_z);



        //this->manage_pressure();
        /*Vec pressure_bench;
        DMCreateGlobalVector(dmGrid_centered, &pressure_bench);
        CreateReferencePressure(dmGrid_centered, pressure_bench, theta);
        CheckSolution(P, pressure_bench);
        VecDestroy(&pressure_bench);*/
        /*this->manage_pressure_x();
        this->manage_pressure_y();
        this->manage_pressure_z();*/
        //theta = (i+1)*dt;

        this->update_velocity(theta);

        if(check_convergence)
        {
            Vec solution_x;
            DMCreateGlobalVector(dmGrid_staggered_x, &solution_x);
            CreateAnalyticalU(dmGrid_staggered_x, solution_x, theta);
            CheckSolution(U_up, solution_x, "U");
            VecDestroy(&solution_x);
            Vec solution_y;
            DMCreateGlobalVector(dmGrid_staggered_y, &solution_y);
            CreateAnalyticalV(dmGrid_staggered_y, solution_y, theta);
            CheckSolution(V_up, solution_y, "V");
            VecDestroy(&solution_y);
            Vec solution_z;
            DMCreateGlobalVector(dmGrid_staggered_z, &solution_z);
            CreateAnalyticalW(dmGrid_staggered_z, solution_z, theta);
            CheckSolution(W_up, solution_z, "W");
            VecDestroy(&solution_z);

            /*Vec solution_p;
            DMCreateGlobalVector(dmGrid_centered, &solution_p);
            double mean_analytical;
            VecSum(solution_p, &mean_analytical);
            int size;
            VecGetSize(solution_p, &size);
            mean_analytical /= size;
            double mean_numerical;
            VecSum(P, &mean_numerical);
            mean_numerical /= size;
            VecShift(P, -mean_numerical);
            CreateAnalyticalP(dmGrid_centered, solution_p, theta);
            CheckSolution(P, solution_p, "P");
            VecDestroy(&solution_p);
            */




        }

        compute_magnitude();

        exodus(i);

        std::cout << "Iteration " << i << " completed." << std::endl; 




    

    }




    PetscFunctionReturn(0);
}

navier_stokes_problem::~navier_stokes_problem()
{
    VecDestroy(&P);
    VecDestroy(&P_x);
    VecDestroy(&P_y);
    VecDestroy(&P_z);
    VecDestroy(&Magnitude);
    VecDestroy(&U_up);
    VecDestroy(&V_up);
    VecDestroy(&W_up);
    DMDestroy(&dmGrid_staggered_x);
    DMDestroy(&dmGrid_staggered_y);
    DMDestroy(&dmGrid_staggered_z);
    DMDestroy(&dmGrid_centered);
    DMDestroy(&dmGrid_cent_rich);
    DMDestroy(&dmGrid_shift_transp);
    DMDestroy(&dmGrid_stag_transp);
    VecDestroy(&mask_U);
    VecDestroy(&mask_V);
    VecDestroy(&mask_W);

    std::cout << "Navier-Stokes Destructor Called" << std::endl;

}