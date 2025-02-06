#include "stokes.hpp"


PetscErrorCode const stokes_problem::update_velocity(PetscReal const & theta)
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

PetscErrorCode const stokes_problem::assemble_magnitude(Vec & Magnitude_Shifted, Vec const & U, Vec const & V, Vec const & W) 
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

PetscErrorCode const stokes_problem::compute_magnitude()
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

PetscErrorCode stokes_problem::exodus(size_t i){

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

PetscErrorCode const stokes_problem::solve()
{
    PetscFunctionBegin;

    PrintSimulationParameters();

    //this->assemble_lhs();

    parabolic_problem_x parabolic_x(dmGrid_staggered_x);
    parabolic_problem_y parabolic_y(dmGrid_staggered_y); 
    parabolic_problem_z parabolic_z(dmGrid_staggered_z);

    parabolic_x.assemble_lhs();
    parabolic_y.assemble_lhs();
    parabolic_z.assemble_lhs();  

    poisson_problem poisson(dmGrid_staggered_x, dmGrid_staggered_y, dmGrid_staggered_z, dmGrid_centered, dmGrid_cent_rich);
    
    for(size_t i = 0; i < iter; ++i){

        theta = i * dt;

        parabolic_x.solve_step(theta, U_up);
        parabolic_y.solve_step(theta, V_up);
        parabolic_z.solve_step(theta, W_up);

        poisson.manage_pressure(U_up, V_up, W_up, P);
        poisson.manage_pressure_x(P, P_x);
        poisson.manage_pressure_y(P, P_y);
        poisson.manage_pressure_z(P, P_z);



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
