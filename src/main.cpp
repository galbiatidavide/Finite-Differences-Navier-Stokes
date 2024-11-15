//#include "MeshHandler.hpp"
#include "parabolic.hpp"


int main(int argc, char **argv) {
    
    PetscInitialize(&argc, &argv, (char*)0, (char*)0);

    Params input;
    input.n_discr = {4, 4, 4};
    input.intervals = {{{-0.5, 0.5}, {-0.5, 0.5}, {-0.5, 0.5}}};
    input.T = 1.0;
    input.dt = 0.5;
    input.Re = 10;
    input.iter = 2;

{

    Parabolic pb(input);
    pb.Solve();

}

    PetscFinalize();

    return 0;

}