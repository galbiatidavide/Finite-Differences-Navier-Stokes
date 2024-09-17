#include <iostream>
#include <chrono>
#include <cmath>
#include <limits>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscvec.h>
#include <array>
#include <petscvec.h>


#define BACK_DOWN_LEFT   DMSTAG_BACK_DOWN_LEFT
#define BACK_DOWN        DMSTAG_BACK_DOWN
#define BACK_DOWN_RIGHT  DMSTAG_BACK_DOWN_RIGHT
#define BACK_LEFT        DMSTAG_BACK_LEFT
#define BACK             DMSTAG_BACK
#define BACK_RIGHT       DMSTAG_BACK_RIGHT
#define BACK_UP_LEFT     DMSTAG_BACK_UP_LEFT
#define BACK_UP          DMSTAG_BACK_UP
#define BACK_UP_RIGHT    DMSTAG_BACK_UP_RIGHT
#define DOWN_LEFT        DMSTAG_DOWN_LEFT
#define DOWN             DMSTAG_DOWN
#define DOWN_RIGHT       DMSTAG_DOWN_RIGHT
#define LEFT             DMSTAG_LEFT
#define ELEMENT          DMSTAG_ELEMENT
#define RIGHT            DMSTAG_RIGHT
#define UP_LEFT          DMSTAG_UP_LEFT
#define UP               DMSTAG_UP
#define UP_RIGHT         DMSTAG_UP_RIGHT
#define FRONT_DOWN_LEFT  DMSTAG_FRONT_DOWN_LEFT
#define FRONT_DOWN       DMSTAG_FRONT_DOWN
#define FRONT_DOWN_RIGHT DMSTAG_FRONT_DOWN_RIGHT
#define FRONT_LEFT       DMSTAG_FRONT_LEFT
#define FRONT            DMSTAG_FRONT
#define FRONT_RIGHT      DMSTAG_FRONT_RIGHT
#define FRONT_UP_LEFT    DMSTAG_FRONT_UP_LEFT
#define FRONT_UP         DMSTAG_FRONT_UP
#define FRONT_UP_RIGHT   DMSTAG_FRONT_UP_RIGHT


// Define the Params structure
struct Params {
    std::array<PetscScalar, 4> n_discr;
    std::array<PetscScalar, 4> dofs;
    std::array<std::array<PetscScalar, 2>, 3> intervals;
    const PetscInt stencilWidth = 1;
};



// Template class for Grid
template <typename GridType>
class Grid {
protected:
    DM dmGrid;
    Vec vecLocal;
    Params input;

public:
    // Constructor
    Grid(Params given_input) : input(given_input) {
        static_cast<GridType*>(this)->setDofs(input); // Call the derived class method to set dofs
        CreateGrid(&dmGrid, input);
    }

    // Create grid 
    PetscErrorCode CreateGrid(DM* dmGrid, Params input){
        PetscFunctionBegin;
        DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, input.n_discr[0], input.n_discr[1], input.n_discr[2], 
                        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 
                        input.dofs[0], input.dofs[1], input.dofs[2], input.dofs[3], 
                        DMSTAG_STENCIL_BOX, input.stencilWidth, NULL, NULL, NULL, dmGrid);

        DMSetFromOptions(*dmGrid);
        DMSetUp(*dmGrid);
        DMStagSetUniformCoordinatesExplicit(*dmGrid, input.intervals[0][0], input.intervals[0][1], 
                                                     input.intervals[1][0], input.intervals[1][1], 
                                                     input.intervals[2][0], input.intervals[2][1]);

        PetscFunctionReturn(0);
    };


    void print_input() {
        std::cout << input.dofs[0] << " " << input.dofs[1] << " " << input.dofs[2] << " " << input.dofs[3] << std::endl;
    }

    // Destructor
    virtual ~Grid() {};
};

// Derived template specialization for staggered grid
class StaggeredGrid : public Grid<StaggeredGrid> {
public:
    StaggeredGrid(Params given_input) : Grid(given_input) {}

    // Method to set dofs for staggered grid
    void setDofs(Params& input) {
        input.dofs = {0, 0, 1, 0}; // Set staggered grid dofs
    }

    ~StaggeredGrid() {};
};

class CenteredGrid : public Grid<CenteredGrid> {
public:
    CenteredGrid(Params given_input) : Grid(given_input) {}

    // Method to set dofs for centered grid
    void setDofs(Params& input) {
        input.dofs = {0, 0, 1, 1}; // Set centered grid dofs
    }

    ~CenteredGrid() {};
};

class ShiftedGrid : public Grid<ShiftedGrid> {
public:
    ShiftedGrid(Params given_input) : Grid(given_input) {}

    // Method to set dofs for shifted grid
    void setDofs(Params& input) {
        input.dofs = {0, 1, 1, 0}; // Set shifted grid dofs
    }

    ~ShiftedGrid() {};
};


















