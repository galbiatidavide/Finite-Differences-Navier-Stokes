#include "MeshHandler.hpp"


int main(int argc, char **argv) {
    
    PetscInitialize(&argc, &argv, (char*)0, (char*)0);

    Params input;
    input.dofs = {8, 2, 1, 0};
    input.n_discr = {10, 10, 10, 10};
    input.intervals = {{{0, 1}, {0, 1}, {0, 1}}};

    MeshHandler<StaggeredGrid> staggeredMeshHandler(input, "surface_cut.stl");
    staggeredMeshHandler.reader();
    staggeredMeshHandler.save_vertices();



    PetscFinalize();
    return 0;

}