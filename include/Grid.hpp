#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <limits>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscvec.h>
#include <array>
#include <vector>
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
    std::array<PetscInt, 4> n_discr;
    std::array<PetscInt, 4> dofs;
    std::array<std::array<PetscScalar, 2>, 3> intervals;
    const PetscInt stencilWidth = 1;
};



// Template class for Grid
template <typename GridType>
class Grid {
protected:
    DM dmGrid;
    Vec globalVec;
    Params input;
    PetscInt gridSize = input.n_discr[0] * input.n_discr[1] * input.n_discr[2];
    std::vector<std::array<double, 3>> coordinates;
    std::vector<DMStagStencilLocation> boundaryTypes;

public:
    // Constructor
    Grid(Params given_input) : input(given_input) {
        static_cast<GridType*>(this)->setDofs(input); // Call the derived class method to set dofs
        CreateGrid(&dmGrid, input);
        DMCreateGlobalVector(dmGrid, &globalVec);
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


    void get_coordinates(){
        static_cast<GridType*>(this)->setTypes();
        std::vector<PetscInt[3]> icux(boundaryTypes.size()); 
        
        PetscInt startx, starty, startz, N[3], ex, ey, ez, d;
        DM dmCoord;
        Vec coord, coordLocal, vecLocal;
        PetscReal ****arrCoord, ****arrVec;   

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &input.n_discr[0], &input.n_discr[1], &input.n_discr[2], NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);


        for (auto b : boundaryTypes) {
            for (auto ic : icux){
                for (d = 0; d < 3; ++d) {
                    DMStagGetLocationSlot(dmCoord, b, d, &ic[d]);
                }
            }
        }


        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
        DMCreateLocalVector(dmGrid, &vecLocal);
        DMGlobalToLocalBegin(dmGrid, globalVec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid, globalVec, INSERT_VALUES, vecLocal);
        DMStagVecGetArray(dmGrid, vecLocal, &arrVec);  

        for (ez = startz; ez < startz + input.n_discr[2]; ++ez) {
            for (ey = starty; ey < starty + input.n_discr[1]; ++ey) {
                for (ex = startx; ex < startx + input.n_discr[0]; ++ex) {
                    for(auto i : icux) {
                    coordinates.push_back({arrCoord[ez][ey][ex][i[0]], arrCoord[ez][ey][ex][i[1]], arrCoord[ez][ey][ex][i[2]]});
                    coordFile << arrCoord[ez][ey][ex][i[0]] << " " << arrCoord[ez][ey][ex][i[1]] << " " << arrCoord[ez][ey][ex][i[2]] << std::endl;
                    }
                }

            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMRestoreLocalVector(dmCoord, &coordLocal);
        DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
        DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, globalVec);
        DMRestoreLocalVector(dmGrid, &vecLocal);
    
    }

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

    void setTypes() {
        boundaryTypes = {RIGHT, LEFT, UP, DOWN, BACK, FRONT};
    }

    ~StaggeredGrid() {};

};


class CenteredGrid : public Grid<CenteredGrid> {
public:
    CenteredGrid(Params given_input) : Grid(given_input) {}

    // Method to set dofs for centered grid
    void setDofs(Params& input) {
        input.dofs = {0, 0, 0, 1}; // Set centered grid dofs
    }

    void setTypes() {
        boundaryTypes = {ELEMENT};
    }

    ~CenteredGrid() {};
};

class ShiftedGrid : public Grid<ShiftedGrid> {
public:
    ShiftedGrid(Params given_input) : Grid(given_input) {}

    // Method to set dofs for shifted grid
    void setDofs(Params& input) {
        input.dofs = {0, 1, 1, 0}; // Set shifted grid dofs lati e facce per termini misti del non lineare
    }

    void setTypes() {
        boundaryTypes = {RIGHT, LEFT, UP, DOWN, BACK, FRONT};
    }

    ~ShiftedGrid() {};
};


















