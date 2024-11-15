//#include "MeshHandler.hpp"
#include "Parabolic.hpp"


int main(int argc, char **argv) {
    
    PetscInitialize(&argc, &argv, (char*)0, (char*)0);

    Params input;
    input.n_discr = {4, 4, 4};
    input.intervals = {{{-0.5, 0.5}, {-0.5, 0.5}, {-0.5, 0.5}}};
    input.T = 1.0;
    input.dt = 0.5;
    input.Re = 10;

{

    // StaggeredGrid grid(input);
    // grid.bc_setUp("stokes");
    // grid.save_grid();
    // MeshHandler<StaggeredGrid> staggeredMeshHandler(input, "surface_cut.stl");
    // staggeredMeshHandler.bc_setUp("stationary");
    // staggeredMeshHandler.reader();
    // staggeredMeshHandler.give_penalty();
    // staggeredMeshHandler.save_grid();


    Parabolic pb(input);
    //pb.assemble_matrices();
    //pb.saveMatrices();
    pb.Solve();
}

    PetscFinalize();

    return 0;

}