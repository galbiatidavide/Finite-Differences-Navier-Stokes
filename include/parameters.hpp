#include <iostream>
#include <chrono>
#include <cmath>
#include <limits>
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscdm.h>


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

PetscReal constexpr pi = 3.14159265358979323846;
PetscReal constexpr eps=std::numeric_limits<float>::max();

//Define flow parameters
PetscReal constexpr A{4*sqrt(2)/(3*sqrt(3))};
PetscReal constexpr vRef{4};
PetscReal constexpr k{2*pi};
//PetscReal constexpr c{(5/6)*pi};
//PetscReal constexpr d{(1/6)*pi};
//PetscReal constexpr vNorm{1e-2*(1/A)};
//PetscReal constexpr vRef{vNorm*A};*/
//PetscReal constexpr vRef = 1;
PetscReal constexpr a = pi/4;
PetscReal constexpr d = 1.5*pi;
//PetscReal constexpr A = 1;
//PetscReal constexpr B = 1;
//PetscReal constexpr C = 1;
//Define the domain & grid
PetscInt constexpr nx{32};
PetscInt constexpr ny{32};
PetscInt constexpr nz{32};
PetscReal constexpr Lx_0{-0.5};
PetscReal constexpr Ly_0{-0.5};
PetscReal constexpr Lz_0{-0.5};
PetscReal constexpr Lx{0.5};
PetscReal constexpr Ly{0.5};
PetscReal constexpr Lz{0.5};


//Define the time
PetscReal constexpr dt{0.001}; 
PetscReal constexpr iter{1000};

PetscReal constexpr Re{1};
PetscReal constexpr LRef{100};
PetscReal constexpr nu{vRef*LRef/Re};
PetscReal theta{0};




constexpr PetscReal uxRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    //return vRef*A*(sin(k*x - c)*cos(k*y - d)*sin(k*z) - cos(k*z - c)*sin(k*x - d)*sin(k*y))*exp(-theta);
    //return vRef*(A*sin(k*z) + C*cos(k*y))*exp(-theta);
    return -a*(exp(a*x)*sin(a*y + d*z) + exp(a*z)*cos(a*x + d*y))*exp(-theta); //navier stokes
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z; //parabolic
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return (-k*cos(k*x)*cos(k*z)*sin(k*y) - k*sin(k*x)*sin(k*y)*sin(k*z))*exp(-theta); //stokes
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return z*z;

    //return -0.1*y/(2*pi*sqrt((x)*(x) + (y)*(y))) + 0.1*cos(10*pi*z);
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
}

constexpr PetscReal uyRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    //return vRef*A*(sin(k*y - c)*cos(k*z - d)*sin(k*x) - cos(k*x - c)*sin(k*y - d)*sin(k*z))*exp(-theta);
    //return vRef*(B*sin(k*x) + A*cos(k*z))*exp(-theta);
    return -a*(exp(a*y)*sin(a*z + d*x) + exp(a*x)*cos(a*y + d*z))*exp(-theta);
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z;
    //return (-k*cos(k*x)*cos(k*y)*sin(k*z) - k*sin(k*x)*sin(k*y)*sin(k*z))*exp(-theta);
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return 0.1*x/(2*pi*sqrt((x)*(x) + (y)*(y))) + + 0.1*cos(10*pi*z);
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
}

constexpr PetscReal uzRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    //return vRef*A*(sin(k*z - c)*cos(k*x - d)*sin(k*y) - cos(k*y - c)*sin(k*z - d)*sin(k*x))*exp(-theta);
    //return vRef*(C*sin(k*y) + B*cos(k*x))*exp(-theta);
    return -a*(exp(a*z)*sin(a*x + d*y) + exp(a*y)*cos(a*z + d*x))*exp(-theta);
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z;
    //return (-k*cos(k*y)*cos(k*z)*sin(k*x) - k*sin(k*x)*sin(k*y)*sin(k*z))*exp(-theta);
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return 0;
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);


}

constexpr PetscReal solution(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    //return vRef*A*(sin(k*x - c)*cos(k*y - d)*sin(k*z) - cos(k*z - c)*sin(k*x - d)*sin(k*y))*exp(-theta);
    //return vRef*(A*sin(k*z) + C*cos(k*y))*exp(-theta);
    return -a*(exp(a*x)*sin(a*y + d*z) + exp(a*z)*cos(a*x + d*y))*exp(-theta);
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z; //parabolic
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) + 2*pi*0.1*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)*sin(2*pi*x) + 2*pi*0.1*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)*sin(2*pi*y) + 2*pi*0.1*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)*sin(2*pi*z);

    //return (-k*cos(k*x)*cos(k*z)*sin(k*y) - k*sin(k*x)*sin(k*y)*sin(k*z))*exp(-theta); //stokes
    //return -0.1*y/(2*pi*sqrt((x)*(x) + (y)*(y)));
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) + 0.05*4*pi*cos(2*pi*x)*sin(2*pi*x)*cos(2*pi*y)*cos(2*pi*y)*cos(2*pi*z)*cos(2*pi*z) +0.05*4*pi*cos(2*pi*x)*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*y)*cos(2*pi*z)*sin(2*pi*z)+0.05*4*pi*cos(2*pi*x)*cos(2*pi*x)*cos(2*pi*y)*sin(2*pi*y)*cos(2*pi*z)*cos(2*pi*z);
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);



}

constexpr PetscReal pRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    //return vRef*A*(sin(k*x - c)*cos(k*y - d)*sin(k*z) - cos(k*z - c)*sin(k*x - d)*sin(k*y))*exp(-theta);
    //return vRef*(A*sin(k*z) + C*cos(k*y))*exp(-theta);
    //return -vRef*a*(exp(a*x)*sin(a*y + d*z) + exp(a*z)*cos(a*x + d*y))*exp(-theta);
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z;
    //return (-k*cos(k*x)*cos(k*z)*sin(k*y) - k*sin(k*x)*sin(k*y)*sin(k*z))*exp(-theta); //stokes
    //return -0.1*y/(2*pi*sqrt((x)*(x) + (y)*(y)));
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) + 0.05*4*pi*cos(2*pi*x)*sin(2*pi*x)*cos(2*pi*y)*cos(2*pi*y)*cos(2*pi*z)*cos(2*pi*z) +0.05*4*pi*cos(2*pi*x)*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*y)*cos(2*pi*z)*sin(2*pi*z)+0.05*4*pi*cos(2*pi*x)*cos(2*pi*x)*cos(2*pi*y)*sin(2*pi*y)*cos(2*pi*z)*cos(2*pi*z);
    //return 2;
    return -0.5*a*a*(exp(2*a*x) + exp(2*a*y) + exp(2*a*z) + 2*sin(a*x + d*y)*cos(a*z + d*x)*exp(a*(y + z)) +  2*sin(a*y + d*z)*cos(a*x + d*y)*exp(a*(x + z)) +  2*sin(a*z + d*x)*cos(a*y + d*z)*exp(a*(y + x)))*exp(-2*theta);
    


}


PetscErrorCode CheckSolution(Vec const & sol, Vec const & solRef)
{
    Vec       diff;
    PetscReal normsolRef, errAbs, errRel;
    PetscFunctionBegin;

    PetscFunctionBegin;
    VecDuplicate(sol, &diff);
    VecCopy(sol, diff);
    VecAXPY(diff, -1.0, solRef);
    VecNorm(diff, NORM_2, &errAbs);
    VecNorm(solRef, NORM_2, &normsolRef);
    errRel = errAbs / normsolRef;
    PetscPrintf(PETSC_COMM_WORLD, "Error (abs): %g\nError (rel): %g\n", (double)errAbs, (double)errRel);
    VecDestroy(&diff);
    PetscFunctionReturn(0);
}


PetscErrorCode CreateReferenceSolutionTry(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iux, icux[3], iuy, icuy[3], iuz, icuz[3], iue, icue[3];
    DM              dmCoord;
    Vec             vecLocal, coord, coordLocal;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, LEFT, 0, &icux[0]);
    DMStagGetLocationSlot(dmCoord, LEFT, 1, &icux[1]);
    DMStagGetLocationSlot(dmCoord, LEFT, 2, &icux[2]); 
    DMStagGetLocationSlot(dmCoord, DOWN, 0, &icuy[0]);
    DMStagGetLocationSlot(dmCoord, DOWN, 1, &icuy[1]);
    DMStagGetLocationSlot(dmCoord, DOWN, 2, &icuy[2]);
    DMStagGetLocationSlot(dmCoord, BACK, 0, &icuz[0]);
    DMStagGetLocationSlot(dmCoord, BACK, 1, &icuz[1]);
    DMStagGetLocationSlot(dmCoord, BACK, 2, &icuz[2]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 0, &icue[0]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 1, &icue[1]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 2, &icue[2]);     
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iue);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iux] = solution(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode CreateReferencePressure(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iux, icux[3], iuy, icuy[3], iuz, icuz[3], iue, icue[3];
    DM              dmCoord;
    Vec             vecLocal, coord, coordLocal;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, LEFT, 0, &icux[0]);
    DMStagGetLocationSlot(dmCoord, LEFT, 1, &icux[1]);
    DMStagGetLocationSlot(dmCoord, LEFT, 2, &icux[2]); 
    DMStagGetLocationSlot(dmCoord, DOWN, 0, &icuy[0]);
    DMStagGetLocationSlot(dmCoord, DOWN, 1, &icuy[1]);
    DMStagGetLocationSlot(dmCoord, DOWN, 2, &icuy[2]);
    DMStagGetLocationSlot(dmCoord, BACK, 0, &icuz[0]);
    DMStagGetLocationSlot(dmCoord, BACK, 1, &icuz[1]);
    DMStagGetLocationSlot(dmCoord, BACK, 2, &icuz[2]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 0, &icue[0]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 1, &icue[1]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 2, &icue[2]);     
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iue);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iue] = pRef(arrCoord[ez][ey][ex][icue[0]], arrCoord[ez][ey][ex][icue[1]], arrCoord[ez][ey][ex][icue[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}


/*
PetscErrorCode CreateReferenceSolutionTryForce(DM const & dmGrid, Vec & vec, PetscReal const & theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, iux, icux[3], iuy, icuy[3], iuz, icuz[3], iue, icue[3];
    DM              dmCoord;
    Vec             vecLocal, coord, coordLocal;
    PetscReal ****arrVec, ****arrCoord;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmCoord, LEFT, 0, &icux[0]);
    DMStagGetLocationSlot(dmCoord, LEFT, 1, &icux[1]);
    DMStagGetLocationSlot(dmCoord, LEFT, 2, &icux[2]); 
    DMStagGetLocationSlot(dmCoord, DOWN, 0, &icuy[0]);
    DMStagGetLocationSlot(dmCoord, DOWN, 1, &icuy[1]);
    DMStagGetLocationSlot(dmCoord, DOWN, 2, &icuy[2]);
    DMStagGetLocationSlot(dmCoord, BACK, 0, &icuz[0]);
    DMStagGetLocationSlot(dmCoord, BACK, 1, &icuz[1]);
    DMStagGetLocationSlot(dmCoord, BACK, 2, &icuz[2]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 0, &icue[0]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 1, &icue[1]);
    DMStagGetLocationSlot(dmCoord, ELEMENT, 2, &icue[2]);     
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iue);
    DMGetLocalVector(dmGrid, &vecLocal);
    DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                arrVec[ez][ey][ex][iux] = force(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, vec);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &vecLocal);

    PetscFunctionReturn(0);
}


*/