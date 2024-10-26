#include <petsc.h>
#include <array>
#include <memory>
#include <petscmat.h>
#include "Grid.hpp"

// TO DO:
// 1. Parsare BC
// 2. Parsare Exact Solution 
// 3. Setter per Reynolds


class Parabolic {
private:
    std::unique_ptr<Grid> grid;
    PetscReal dt; //Time step
    PetscReal T; //Final Time
    PetscReal time; //Current Time
    std::vector<Mat> lhs_comp;
    std::vector<Vec> rhs_comp;
    std::function<double(double, double, double, double)>*  exactSolution[3];
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
            exactSolution[i] = &grid->bc.bcFunctions[i];
        }

        


    };

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



PetscErrorCode saveMatrices(const std::string& filename_prefix);


protected:

    // Assemble the mass and stiffness matrices.
    

    // Assemble the right-hand side of the problem.
    PetscErrorCode assemble_rhs(const double &time);


    // Solve the problem for one time step.
    void solve_time_step(){};

public:

//destructors
~Parabolic() {};

};


