#include <petsc.h>
#include <array>
#include <memory>
#include <cmath>
#include <limits>
#include <petscmat.h>
#include <petscdm.h>
#include <petscksp.h>
#include "grid.hpp"

// TO DO:
// 1. Parsare BC
// 2. Parsare Exact Solution 
// 3. Setter per Reynolds
// 4. capire pointer functions o functions in bc e exact solution


class Parabolic {
private:
    std::unique_ptr<Grid> grid;
    PetscReal dt; //Time step
    PetscReal T; //Final Time
    PetscReal time; //Current Time
    std::vector<Mat> lhs_comp;
    std::vector<Vec> rhs_comp;
    std::function<double(double, double, double, double)>  exactSolution[3];
    std::vector<Component> n_components; 
    unsigned int iter;

    //parameters 
    PetscReal Re = 10;
    PetscReal Ret = Re/dt;
    PetscReal hx, hy, hz;

public:

    // Costruttore e distruttore
    Parabolic(Params input): dt(input.dt), T(input.T), 
                            hx((input.intervals[0][1]-input.intervals[0][0])/input.n_discr[0]), 
                            hy((input.intervals[1][1]-input.intervals[1][0])/input.n_discr[1]), 
                            hz((input.intervals[2][1]-input.intervals[2][0])/input.n_discr[2]), 
                            iter(input.iter) {
        
        grid = std::make_unique<StaggeredGrid>(input);
        grid->bc_setUp("parabolic");
        //grid->save_grid();

        lhs_comp.resize(grid->components.size());
        rhs_comp.resize(grid->components.size());

        for(auto &lhs : lhs_comp){
            DMCreateMatrix(grid->dmGrid, &lhs);
        }

        for(auto &rhs : rhs_comp){
            DMCreateGlobalVector(grid->dmGrid, &rhs);
        }

         for (int i = 0; i < 3; ++i) {
            exactSolution[i] = grid->bc.bcFunctions[i];
        }


    };

    PetscErrorCode init_n(){

        PetscFunctionBegin;

        n_components.resize(3);
        Vec U_n, V_n, W_n;
        DMCreateGlobalVector(grid->dmGrid, &U_n);
        DMCreateGlobalVector(grid->dmGrid, &V_n);
        DMCreateGlobalVector(grid->dmGrid, &W_n);

        n_components[0] = {U_n, {grid->components[0].location[0], grid->components[0].location[1]}, grid->components[0].name + "_n"};
        n_components[1] = {V_n, {grid->components[1].location[0], grid->components[1].location[1]}, grid->components[1].name + "_n"};
        n_components[2] = {W_n, {grid->components[2].location[0], grid->components[2].location[1]}, grid->components[2].name + "_n"};  

        VecCopy(grid->components[0].variable, n_components[0].variable);
        VecCopy(grid->components[1].variable, n_components[1].variable);
        VecCopy(grid->components[2].variable, n_components[2].variable);  

        //grid->save_grid();
        output(0); 
        
        PetscFunctionReturn(0);
    }


    PetscErrorCode assemble_matrices();
    PetscErrorCode assemble_rhs(const unsigned int &time_step);

    PetscErrorCode solveTimeStep(const unsigned int &time_step){

    PetscFunctionBegin;

    for (size_t i = 0; i < grid->components.size(); ++i) {  
    KSP       ksp;
    PC        pc;
    
    if(time_step == 0){
        init_n();
        assemble_matrices(); 
    }
  
        assemble_rhs(time_step);
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetType(ksp, KSPCG);
        KSPSetOperators(ksp, lhs_comp[i], lhs_comp[i]);
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCFIELDSPLIT);
        PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, rhs_comp[i], grid->components[i].variable);
        KSPDestroy(&ksp);  
        VecCopy(grid->components[i].variable, n_components[i].variable); 

    }

    PetscFunctionReturn(0);

};

const PetscErrorCode Solve()
    {
        PetscFunctionBegin;

        for(size_t i = 0; i < iter; i++){            
            
            solveTimeStep(i);

            Vec bench;
            DMCreateGlobalVector(grid->dmGrid, &bench);
            grid->CreateReferenceSolutionTry(0, i);
            CheckSolution(grid->components[0].variable);
            output(i);
              
        }
        PetscFunctionReturn(0);
    }

    
PetscErrorCode saveMatrices();


protected:

    // Assemble the mass and stiffness matrices.
    

    // Assemble the right-hand side of the problem.
    //PetscErrorCode assemble_rhs(const double &time);

    // Solve the problem for one time step.
    //PetscErrorCode solveTimeStep();

    PetscErrorCode CheckSolution(Vec const & sol)
    {
    Vec       diff;
    PetscReal normsolRef, errAbs, errRel;
    PetscFunctionBegin;

    PetscFunctionBegin;
    VecDuplicate(sol, &diff);
    VecCopy(sol, diff);
    VecAXPY(diff, -1.0, grid->globalVec);
    VecNorm(diff, NORM_2, &errAbs);
    VecNorm(grid->globalVec, NORM_2, &normsolRef);
    errRel = errAbs / normsolRef;
    PetscPrintf(PETSC_COMM_WORLD, "Error (abs): %g\nError (rel): %g\n", (double)errAbs, (double)errRel);
    VecDestroy(&diff);
    PetscFunctionReturn(0);
    }

    PetscErrorCode output(const unsigned int &time_step) {
        
        PetscFunctionBegin
        for(unsigned int i = 0; i < grid->components.size(); i++){
            PetscViewer viewer;
            Vec r;
            DM pda;
            DMStagVecSplitToDMDA(grid->dmGrid, n_components.at(i).variable, n_components.at(i).location[0],  DM_BOUNDARY_NONE, &pda, &r);
            PetscObjectSetName((PetscObject)r, n_components.at(i).name.c_str());
            char filename_r[50];
            sprintf(filename_r, "results/%s_000%i.vtr", n_components.at(i).name.c_str(), time_step);
            PetscViewerVTKOpen(PetscObjectComm((PetscObject)pda), filename_r, FILE_MODE_WRITE, &viewer);
            VecView(r, viewer);
            VecDestroy(&r);
            DMDestroy(&pda);
            PetscViewerDestroy(&viewer);
        }

        PetscFunctionReturn(0);
    }

PetscErrorCode output_rhs(const unsigned int &time_step) {
        
        PetscFunctionBegin
        for(unsigned int i = 0; i < grid->components.size(); i++){
            // PetscViewer viewer;
            // Vec r;
            // DM pda;
            // DMStagVecSplitToDMDA(grid->dmGrid, rhs_comp[i], grid->components[i].location[0],  DM_BOUNDARY_NONE, &pda, &r);
            // PetscObjectSetName((PetscObject)r, grid->components.at(i).name.c_str());
            // char filename_r[50];
            // sprintf(filename_r, "results/rhs%i_00%i.vtr", i, time_step);
            // PetscViewerVTKOpen(PetscObjectComm((PetscObject)pda), filename_r, FILE_MODE_WRITE, &viewer);
            // VecView(r, viewer);
            // VecDestroy(&r);
            // DMDestroy(&pda);
            // PetscViewerDestroy(&viewer);
            PetscViewer viewer_2;
            Vec r;
            DM pda;

            DMStagVecSplitToDMDA(grid->dmGrid, rhs_comp[i], grid->components[i].location[0], DM_BOUNDARY_NONE, &pda, &r);
            PetscObjectSetName((PetscObject)r, "rhs");  // Set name of vector

            char filename_r[50];
            sprintf(filename_r, "results/rhs_00%i.txt", i);  // Change extension to .txt
            PetscViewerASCIIOpen(PetscObjectComm((PetscObject)pda), filename_r, &viewer_2); // Use ASCII viewer

            VecView(r, viewer_2);  // View the vector contents in text format

            // Cleanup
            VecDestroy(&r);
            DMDestroy(&pda);
            PetscViewerDestroy(&viewer_2);
        }

        PetscFunctionReturn(0);
}

public:

//destructors
~Parabolic() {

    for(auto c : n_components){
        VecDestroy(&c.variable);
    }

    for(auto lhs : lhs_comp){
        MatDestroy(&lhs);
    }

    for(auto rhs : rhs_comp){
        VecDestroy(&rhs);
    }

};

};


