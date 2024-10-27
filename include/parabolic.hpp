#include <petsc.h>
#include <array>
#include <memory>
#include <petscmat.h>
#include "Grid.hpp"

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

    //parameters 
    PetscReal Re = 1.0;
    PetscReal Ret = Re/dt;
    PetscReal hx, hy, hz;

public:

    // Costruttore e distruttore
    Parabolic(Params input): dt(input.dt), T(input.T), 
                            hx((input.intervals[0][1]-input.intervals[0][0])/input.n_discr[0]), 
                            hy((input.intervals[1][1]-input.intervals[1][0])/input.n_discr[1]), 
                            hz((input.intervals[2][1]-input.intervals[2][0])/input.n_discr[2]) {
        

        grid = std::make_unique<StaggeredGrid>(input);
        grid->bc_setUp("parabolic");

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

        init_n();

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
        
        PetscFunctionReturn(0);
    }

    // Solve the problem.
    void solve();

    //PetscErrorCode assemble_matrices();
    // Compute the error.
    //double compute_error(const VectorTools::NormType &norm_type);

    // Metodi per ottenere i risultati finali
    // double compute_error(){};
    // Vec GetSolutionX(){};
    // Vec GetSolutionY(){};
    // Vec GetSolutionZ(){};
    void SaveSolution(){
        grid->save_grid();
    };

    PetscErrorCode assemble_matrices();
    PetscErrorCode assemble_rhs(const double &time);



PetscErrorCode saveMatrices(const std::string& filename_prefix);


protected:

    // Assemble the mass and stiffness matrices.
    

    // Assemble the right-hand side of the problem.
    //PetscErrorCode assemble_rhs(const double &time);


    // Solve the problem for one time step.
    void solve_time_step(){};

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


