//
// Created by dave on 26/10/23.
//

#include <iostream>
#include <chrono>

#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPoints.h>
#include <vector>
#include <array>
#include <fstream>
#include <unordered_map>
#include <cmath>

#include "Advection.cpp"
#include "Diffusion.cpp"
#include "Poisson.cpp"
#include "Solver.cpp"
#include "Boundary.cpp"
#include "PostProcessing.cpp"
#include "Reader.cpp"

//#include "Setting.cpp"

int main(int argc, char **argv)
{    
    auto start = std::chrono::high_resolution_clock::now();

    PetscInitialize(&argc, &argv, (char*)0, (char*)0);

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    if(rank == 0){
        std::cout << "Refinement level: " << nx << "x" << ny << "x" << nz << std::endl;
        std::cout << "Grid size: " << Lx << "x" << Ly << "x" << Lz << std::endl;
        std::cout << "Time step: " << dt << std::endl;
        std::cout << "Number of iterations: " << iter << std::endl;
        std::cout << "Reynolds number: " << Re << std::endl;
        std::cout << "Reference velocity: " << vRef << std::endl;
        std::cout << "Reference length: " << LRef << std::endl;
        std::cout << "Reference time: " << LRef/vRef << std::endl;
        std::cout << "Reference viscosity: " << vRef*LRef/Re << std::endl;
    }

    // Create necessary grids
    DM dmGrid_Shifted, dmGrid_Centered, dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Staggered; //need to declare to due to solving linear system laplacian messing things up
    {
        CreateGrid(&dmGrid_Shifted, 0, 1, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Centered, 0, 0, 1, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        CreateGrid(&dmGrid_Staggered_x, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
        DMClone(dmGrid_Staggered_x, &dmGrid_Staggered_y);
        DMClone(dmGrid_Staggered_x, &dmGrid_Staggered_z);
    }

    /******************* */

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <surface_cut.stl>" << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> faces;
    std::string filename = argv[1];
    reader(filename, vertices, faces);

    /*DM dmGrid, dm_cent;
    CreateGrid(&dm_cent, 0, 0, 1, nx, nx, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid, 0, 1, 0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);*/

    /*Vec vec_stag, vec_cent;
    DMCreateGlobalVector(dm_cent, &vec_cent);
    DMCreateGlobalVector(dmGrid, &vec_stag);*/

    Vec penalization_U, penalization_V, penalization_W;
    DMCreateGlobalVector(dmGrid_Staggered_x, &penalization_U);
    DMCreateGlobalVector(dmGrid_Staggered_y, &penalization_V);
    DMCreateGlobalVector(dmGrid_Staggered_z, &penalization_W);
    createMaskU(dmGrid_Staggered_x, penalization_U, vertices, faces);
    std::cout << "Mask U created" << std::endl;
    createMaskV(dmGrid_Staggered_y, penalization_V, vertices, faces);
    std::cout << "Mask V created" << std::endl;
    createMaskW(dmGrid_Staggered_z, penalization_W, vertices, faces);
    std::cout << "Mask W created" << std::endl;

    //staggered


    PetscViewer viewer_stag;
    DM dmGrid;
    //DMStagCreateCompatibleDMStag(dmSol_Staggered_y, 0, 0, 1, 0, &DM_v);
    Vec rd;
    DMStagVecSplitToDMDA(dmGrid_Staggered_x, penalization_V, DOWN, 0, &dmGrid, &rd);
    PetscObjectSetName((PetscObject)rd, "r_values_down");
    char filename_rd[50];
    sprintf(filename_rd, "results/r_values_down%03zu.vtr");
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid), filename_rd, FILE_MODE_WRITE, &viewer_stag);
    VecView(rd, viewer_stag);
    VecDestroy(&rd);
    DMDestroy(&dmGrid);
    PetscViewerDestroy(&viewer_stag);

    /********************************* */


    Vec U_0, V_0, W_0;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_0);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_0);





    DMCreateGlobalVector(dmGrid_Staggered_z, &W_0);
    /*PetscViewer viewer_Pi;
    DM da_solution_P;
    DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &da_solution_P);
    Vec P_grid;
    DMStagVecSplitToDMDA(dmGrid_Staggered_x, U_0, LEFT, 0, &da_solution_P, &P_grid);
    PetscObjectSetName((PetscObject)P_grid, "bravoo");
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "bravoo.vtr", FILE_MODE_WRITE, &viewer_Pi);
    VecView(P_grid, viewer_Pi);
    PetscViewerDestroy(&viewer_Pi);*/
    CreateAnalyticalU(dmGrid_Staggered_x, U_0, penalization_U, 0);
    CreateAnalyticalV(dmGrid_Staggered_y, V_0, penalization_V, 0);
    CreateAnalyticalW(dmGrid_Staggered_z, W_0, penalization_W, 0);

   
    Vec U_int, V_int, W_int;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_int);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_int);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_int);
    Vec U_pre, V_pre, W_pre;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_pre);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_pre);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_pre);
    Vec P, P_x, P_y, P_z;
    DMCreateGlobalVector(dmGrid_Centered, &P);
    DMCreateGlobalVector(dmGrid_Staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_Staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_Staggered_z, &P_z);
    Vec U_up, V_up, W_up;
    DMCreateGlobalVector(dmGrid_Staggered_x, &U_up);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_up);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_up);
    Vec Magnitude;
    DMCreateGlobalVector(dmGrid_Centered, &Magnitude);

    for(size_t i = 0; i < iter; ++i){
        if (i == 0){

            Vec U_0_second, V_0_second, W_0_second;
            Vec U_0_third, V_0_third, W_0_third;
            DMCreateGlobalVector(dmGrid_Staggered_x, &U_0_second);
            DMCreateGlobalVector(dmGrid_Staggered_y, &V_0_second);
            DMCreateGlobalVector(dmGrid_Staggered_z, &W_0_second);
            DMCreateGlobalVector(dmGrid_Staggered_x, &U_0_third);
            DMCreateGlobalVector(dmGrid_Staggered_y, &V_0_third);
            DMCreateGlobalVector(dmGrid_Staggered_z, &W_0_third);
            VecCopy(U_0, U_0_second);
            VecCopy(V_0, V_0_second);
            VecCopy(W_0, W_0_second);
            VecCopy(U_0, U_0_third);
            VecCopy(V_0, V_0_third);
            VecCopy(W_0, W_0_third);
            
            //theta = d*d*(i)*dt;
            //theta = (vRef/Re)*k*k*dt*(i);
            theta = dt*(i)/0.2;

            ManageAdvection_x(dt, U_int, U_0, V_0, W_0, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, penalization_U, penalization_V, penalization_W);
            ManageAdvection_y(dt, V_int, U_0_second, V_0_second, W_0_second, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, penalization_U, penalization_V, penalization_W);
            ManageAdvection_z(dt, W_int, U_0_third, V_0_third, W_0_third, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, penalization_U, penalization_V, penalization_W);



            /*Vec bench_x;
            DMCreateGlobalVector(dmGrid_Staggered_z, &bench_x);
            CreateReferenceSolutionTry(dmGrid_Staggered_z, &bench_x, theta);
            CheckSolution(W_int, bench_x);

            //std::cout << "Spirit of Kitty Hawk completed: advection done." << std::endl;
            Vec bench;
            DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
            CreateReferenceSolutionTry(dmGrid_Staggered_x, &bench);
            PetscViewer viewer_Pi;
            DM da_solution_P;
            DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &da_solution_P);
            Vec P_grid;
            DMStagVecSplitToDMDA(dmGrid_Staggered_x, bench, LEFT, 0, &da_solution_P, &P_grid);
            PetscObjectSetName((PetscObject)P_grid, "bravo");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "bravo.vtr", FILE_MODE_WRITE, &viewer_Pi);
            VecView(P_grid, viewer_Pi);
            PetscViewerDestroy(&viewer_Pi);
            CheckSolution(W_int, bench);
            theta = d*d*(i+1)*dt;*/

            ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_pre, V_pre, W_pre, U_int, V_int, W_int, penalization_U, penalization_V, penalization_W, theta);

            /*PetscViewer viewer_Pi_one;
            DM da_solution_P_one;
            DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &da_solution_P_one);
            Vec P_grid_one;
            DMStagVecSplitToDMDA(dmGrid_Staggered_z, U_pre, LEFT, 0, &da_solution_P_one, &P_grid_one);
            PetscObjectSetName((PetscObject)P_grid_one, "bravoU");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P_one), "bravoU.vtr", FILE_MODE_WRITE, &viewer_Pi_one);
            VecView(P_grid_one, viewer_Pi_one);
            PetscViewerDestroy(&viewer_Pi_one);

            PetscViewer viewer_Pi_two;
            DM da_solution_P_two;
            DMStagCreateCompatibleDMStag(dmGrid_Staggered_y, 0, 0, 1, 0, &da_solution_P_two);
            Vec P_grid_two;
            DMStagVecSplitToDMDA(dmGrid_Staggered_y, V_pre, DOWN, 0, &da_solution_P_two, &P_grid_two);
            PetscObjectSetName((PetscObject)P_grid_two, "bravoV");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P_two), "bravoV.vtr", FILE_MODE_WRITE, &viewer_Pi_two);
            VecView(P_grid_two, viewer_Pi_two);
            PetscViewerDestroy(&viewer_Pi_two);

            PetscViewer viewer_Pi_three;
            DM da_solution_P_three;
            DMStagCreateCompatibleDMStag(dmGrid_Staggered_z, 0, 0, 1, 0, &da_solution_P_three);
            Vec P_grid_three;
            DMStagVecSplitToDMDA(dmGrid_Staggered_z, W_pre, BACK, 0, &da_solution_P_three, &P_grid_three);
            PetscObjectSetName((PetscObject)P_grid_three, "bravoW");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P_three), "bravoW.vtr", FILE_MODE_WRITE, &viewer_Pi_three);
            VecView(P_grid_three, viewer_Pi_three);
            PetscViewerDestroy(&viewer_Pi_three);
            std::cout<<"DONEEEEE    "<<i<<std::endl;*/

            //std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
            ManagePressure(dmGrid_Centered, dmGrid_Shifted, dmGrid_Staggered, dt, P, U_pre, V_pre, W_pre);
            ManagePressure_x(dmGrid_Staggered_x, dmGrid_Centered, dmGrid_Shifted, P_x, P);
            ManagePressure_y(dmGrid_Staggered_y, dmGrid_Centered, dmGrid_Shifted, P_y, P);
            ManagePressure_z(dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, P_z, P);
            
            //std::cout << "Spirit of Oklahoma completed:   pressure done." << std::endl;
            //theta = (vRef/Re)*k*k*dt*(i+1);
            theta = dt*(i+1)/0.2;

            UpdateVelocity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, U_up, V_up, W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, theta, penalization_U, penalization_V, penalization_W);
            std::cout << "Spirit of California comdmGrid_centpleted: iteration " << i << " done" << std::endl;

            Vec bench;
            DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
            CreateReferenceSolutionTry(dmGrid_Staggered_x, bench, theta);
            /*PetscViewer viewer_Pi;
            DM da_solution_P;
            DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &da_solution_P);
            Vec P_grid;
            DMStagVecSplitToDMDA(dmGrid_Staggered_x, bench, LEFT, 0, &da_solution_P, &P_grid);
            PetscObjectSetName((PetscObject)P_grid, "bravo");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "bravo.vtr", FILE_MODE_WRITE, &viewer_Pi);
            VecView(P_grid, viewer_Pi);
            PetscViewerDestroy(&viewer_Pi);*/
            CheckSolution(U_up, bench);
            VecDestroy(&bench);

            ComputeMagnitude(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, Magnitude, U_up, V_up, W_up);

            PetscViewer viewer_magnitude;
            DM DM_magnitude;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec magnitude;
            DMStagVecSplitToDMDA(dmGrid_Centered, Magnitude, ELEMENT, 0, &DM_magnitude, &magnitude);
            PetscObjectSetName((PetscObject)magnitude, "magnitude");
            char filename_magnitude[50]; 
            sprintf(filename_magnitude, "results/magnitude%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_magnitude, FILE_MODE_WRITE, &viewer_magnitude);
            VecView(magnitude, viewer_magnitude);
            VecDestroy(&magnitude);
            DMDestroy(&DM_magnitude);
            PetscViewerDestroy(&viewer_magnitude);

            PetscViewer viewer_u;
            DM DM_u;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmGrid_Staggered_x, U_up, LEFT, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "x_component");
            char filename_u[50]; 
            sprintf(filename_u, "results/x_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u); 

            PetscViewer viewer_v;
            DM DM_v;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_y, 0, 0, 1, 0, &DM_v);
            Vec v;
            DMStagVecSplitToDMDA(dmGrid_Staggered_y, V_up, DOWN, 0, &DM_v, &v);
            PetscObjectSetName((PetscObject)v, "y_component");
            char filename_v[50];
            sprintf(filename_v, "results/y_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_y), filename_v, FILE_MODE_WRITE, &viewer_v);
            VecView(v, viewer_v);
            VecDestroy(&v);
            DMDestroy(&DM_v);
            PetscViewerDestroy(&viewer_v);

            PetscViewer viewer_w;
            DM DM_w;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_z, 0, 0, 1, 0, &DM_w);
            Vec w;
            DMStagVecSplitToDMDA(dmGrid_Staggered_z, W_up, BACK, 0, &DM_w, &w);
            PetscObjectSetName((PetscObject)w, "z_component");
            char filename_w[50];
            sprintf(filename_w, "results/z_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_z), filename_w, FILE_MODE_WRITE, &viewer_w);
            VecView(w, viewer_w);
            VecDestroy(&w);
            DMDestroy(&DM_w);
            PetscViewerDestroy(&viewer_w);

            PetscViewer viewer_p;
            DM DM_p;
            //DMStagCreateCompatibleDMStag(dmGrid_Centered, 0, 0, 0, 1, &DM_p);
            Vec p;
            DMStagVecSplitToDMDA(dmGrid_Centered, P, ELEMENT, 0, &DM_p, &p);
            PetscObjectSetName((PetscObject)p, "p");
            char filename_p[50];
            sprintf(filename_p, "results/p%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_p, FILE_MODE_WRITE, &viewer_p);
            VecView(p, viewer_p);
            VecDestroy(&p);
            DMDestroy(&DM_p);
            PetscViewerDestroy(&viewer_p);
            std::cout << "---------------------------------------------------------" << std::endl;

            VecDestroy(&U_0_second);
            VecDestroy(&V_0_second);
            VecDestroy(&W_0_second);
            VecDestroy(&U_0_third);
            VecDestroy(&V_0_third);
            VecDestroy(&W_0_third);

        } else {

            Vec U_up_second, V_up_second, W_up_second;
            Vec U_up_third, V_up_third, W_up_third;
            DMCreateGlobalVector(dmGrid_Staggered_x, &U_up_second);
            DMCreateGlobalVector(dmGrid_Staggered_y, &V_up_second);
            DMCreateGlobalVector(dmGrid_Staggered_z, &W_up_second);
            DMCreateGlobalVector(dmGrid_Staggered_x, &U_up_third);
            DMCreateGlobalVector(dmGrid_Staggered_y, &V_up_third);
            DMCreateGlobalVector(dmGrid_Staggered_z, &W_up_third);
            VecCopy(U_up, U_up_second);
            VecCopy(V_up, V_up_second);
            VecCopy(W_up, W_up_second);
            VecCopy(U_up, U_up_third);
            VecCopy(V_up, V_up_third);
            VecCopy(W_up, W_up_third);                    

            //theta = d*d*(i)*dt;
            //theta = (vRef/Re)*k*k*dt*(i);
            theta = dt*(i)/0.2;


            ManageAdvection_x(dt, U_int, U_up, V_up, W_up, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, penalization_U, penalization_V, penalization_W);
            ManageAdvection_y(dt, V_int, U_up_second, V_up_second, W_up_second, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, penalization_U, penalization_V, penalization_W);
            ManageAdvection_z(dt, W_int, U_up_third, V_up_third, W_up_third, nx, ny, nz, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz, theta, penalization_U, penalization_V, penalization_W);
            //std::cout << "Spirit of Kitty Hawk completed: advection done." << std::endl;

            //theta = d*d*(i+1)*dt;
            ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_pre, V_pre, W_pre, U_int, V_int, W_int, penalization_U, penalization_V, penalization_W, theta);

            //std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
            ManagePressure(dmGrid_Centered, dmGrid_Shifted, dmGrid_Staggered, dt, P, U_pre, V_pre, W_pre);
            ManagePressure_x(dmGrid_Staggered_x, dmGrid_Centered, dmGrid_Shifted, P_x, P);
            ManagePressure_y(dmGrid_Staggered_y, dmGrid_Centered, dmGrid_Shifted, P_y, P);
            ManagePressure_z(dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, P_z, P);
            //std::cout << "Spirit of Oklahoma completed:   pressure done." << std::endl;
            //theta = (vRef/Re)*k*k*dt*(i+1);
            theta = dt*(i+1)/0.2;

            UpdateVelocity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, U_up, V_up, W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, theta, penalization_U, penalization_V, penalization_W);
            //std::cout << "Spirit of California completed: iteration " << i << " done" << std::endl;
            Vec bench;
            DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
            CreateReferenceSolutionTry(dmGrid_Staggered_x, bench, theta);
            /*PetscViewer viewer_Pi;
            DM da_solution_P;
            DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &da_solution_P);
            Vec P_grid;
            DMStagVecSplitToDMDA(dmGrid_Staggered_x, bench, LEFT, 0, &da_solution_P, &P_grid);
            PetscObjectSetName((PetscObject)P_grid, "bravo");
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)da_solution_P), "bravo.vtr", FILE_MODE_WRITE, &viewer_Pi);
            VecView(P_grid, viewer_Pi);
            PetscViewerDestroy(&viewer_Pi);*/
            CheckSolution(U_up, bench);
            VecDestroy(&bench);


            ComputeMagnitude(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, Magnitude, U_up, V_up, W_up);

            PetscViewer viewer_magnitude;
            DM DM_magnitude;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec magnitude;
            DMStagVecSplitToDMDA(dmGrid_Centered, Magnitude, ELEMENT, 0, &DM_magnitude, &magnitude);
            PetscObjectSetName((PetscObject)magnitude, "magnitude");
            char filename_magnitude[50]; 
            sprintf(filename_magnitude, "results/magnitude%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_magnitude, FILE_MODE_WRITE, &viewer_magnitude);
            VecView(magnitude, viewer_magnitude);
            VecDestroy(&magnitude);
            DMDestroy(&DM_magnitude);
            PetscViewerDestroy(&viewer_magnitude);

            PetscViewer viewer_u;
            DM DM_u;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
            Vec u;
            DMStagVecSplitToDMDA(dmGrid_Staggered_x, U_up, LEFT, 0, &DM_u, &u);
            PetscObjectSetName((PetscObject)u, "x_component");
            char filename_u[50]; 
            sprintf(filename_u, "results/x_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
            VecView(u, viewer_u);
            VecDestroy(&u);
            DMDestroy(&DM_u);
            PetscViewerDestroy(&viewer_u); 

            PetscViewer viewer_v;
            DM DM_v;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_y, 0, 0, 1, 0, &DM_v);
            Vec v;
            DMStagVecSplitToDMDA(dmGrid_Staggered_y, V_up, DOWN, 0, &DM_v, &v);
            PetscObjectSetName((PetscObject)v, "y_component");
            char filename_v[50];
            sprintf(filename_v, "results/y_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_y), filename_v, FILE_MODE_WRITE, &viewer_v);
            VecView(v, viewer_v);
            VecDestroy(&v);
            DMDestroy(&DM_v);
            PetscViewerDestroy(&viewer_v);

            PetscViewer viewer_w;
            DM DM_w;
            //DMStagCreateCompatibleDMStag(dmGrid_Staggered_z, 0, 0, 1, 0, &DM_w);
            Vec w;
            DMStagVecSplitToDMDA(dmGrid_Staggered_z, W_up, BACK, 0, &DM_w, &w);
            PetscObjectSetName((PetscObject)w, "z_component");
            char filename_w[50];
            sprintf(filename_w, "results/z_component%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_z), filename_w, FILE_MODE_WRITE, &viewer_w);
            VecView(w, viewer_w);
            VecDestroy(&w);
            DMDestroy(&DM_w);
            PetscViewerDestroy(&viewer_w);

            PetscViewer viewer_p;
            DM DM_p;
            //DMStagCreateCompatibleDMStag(dmGrid_Centered, 0, 0, 0, 1, &DM_p);
            Vec p;
            DMStagVecSplitToDMDA(dmGrid_Centered, P, ELEMENT, 0, &DM_p, &p);
            PetscObjectSetName((PetscObject)p, "p");
            char filename_p[50];
            sprintf(filename_p, "results/p%03zu.vtr", i);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_p, FILE_MODE_WRITE, &viewer_p);
            VecView(p, viewer_p);
            VecDestroy(&p);
            DMDestroy(&DM_p);
            PetscViewerDestroy(&viewer_p);
            std::cout << "Iteration " << i << " completed." << std::endl;     
            std::cout << "---------------------------------------------------------" << std::endl;

            VecDestroy(&U_up_second);
            VecDestroy(&V_up_second);
            VecDestroy(&W_up_second);
            VecDestroy(&U_up_third);
            VecDestroy(&V_up_third);
            VecDestroy(&W_up_third);                
   
        }
    }


    VecDestroy(&U_0);
    VecDestroy(&V_0);
    VecDestroy(&W_0);
    VecDestroy(&U_int);
    VecDestroy(&V_int);
    VecDestroy(&W_int);
    VecDestroy(&U_pre);
    VecDestroy(&V_pre);
    VecDestroy(&W_pre);
    VecDestroy(&P);
    VecDestroy(&P_x);
    VecDestroy(&P_y);
    VecDestroy(&P_z);
    VecDestroy(&U_up);
    VecDestroy(&V_up);
    VecDestroy(&W_up);
    VecDestroy(&Magnitude); 
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_x);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_y);
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered_z);
    PetscObjectDestroy((PetscObject*)&dmGrid_Centered);
    PetscObjectDestroy((PetscObject*)&dmGrid_Shifted);  
    PetscObjectDestroy((PetscObject*)&dmGrid_Staggered);  

    PetscFinalize();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Trinity test successfully completed. Ad maiora!" << std::endl;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    PetscFunctionReturn(0); 
}




