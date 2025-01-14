#include <chrono>
#include <iostream>
#include <limits>
//PETSc
#include <petscdmstag.h>
#include <petscksp.h>
#include <petscmat.h>
//VTK
#include <vtkSTLReader.h>

#ifndef PROBLEM_SETTING_HPP
#define PROBLEM_SETTING_HPP

namespace problem_setting
{
    //Define macros for positions in staggered grid
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

    PetscReal constexpr a = pi/4;
    PetscReal constexpr d = 1.5*pi;
    PetscInt constexpr nx{32};
    PetscInt constexpr ny{32};
    PetscInt constexpr nz{32};
    /*PetscReal constexpr Lx_0{253};
    PetscReal constexpr Ly_0{235};
    PetscReal constexpr Lz_0{190};
    PetscReal constexpr Lx{292};
    PetscReal constexpr Ly{310};
    PetscReal constexpr Lz{232};*/
    PetscReal constexpr Lx_0{-0.5};
    PetscReal constexpr Ly_0{-0.5};
    PetscReal constexpr Lz_0{-0.5};
    PetscReal constexpr Lx{0.5};
    PetscReal constexpr Ly{0.5};
    PetscReal constexpr Lz{0.5};

    PetscReal constexpr dt{0.00625}; 
    PetscReal constexpr iter{15};

    PetscReal constexpr Re{1.0};

    inline bool brinkmann{false};

    inline PetscReal theta;
    //mask
    inline std::vector<std::array<double, 3>> vertices;
    inline std::vector<std::array<int, 3>> faces;
    inline std::string filename;

    constexpr const char *base_path = "results/";


    constexpr PetscReal uxRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
    {
        return -a*(exp(a*x)*sin(a*y + d*z) + exp(a*z)*cos(a*x + d*y))*exp(-d*d*theta); //navier stokes
        //quando fai test ricorda di cambiare il valore di theta e degli altri parametri e di sistemarlo
        //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z; //parabolic
        //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z); //transport
        //return 0;
    }

    constexpr PetscReal uyRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
    {
        return -a*(exp(a*y)*sin(a*z + d*x) + exp(a*x)*cos(a*y + d*z))*exp(-d*d*theta); // navier stokes
        //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z; //parabolic
        //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z); //transport
        //return -1.0;
    }

    constexpr PetscReal uzRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
    {
        return -a*(exp(a*z)*sin(a*x + d*y) + exp(a*y)*cos(a*z + d*x))*exp(-d*d*theta); //navier stokes
        //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z; // parabolic
        //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z); // transport
        //return 0;
    }

    constexpr PetscReal pRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
    {
        return -0.5*a*a*(exp(2*a*x) + exp(2*a*y) + exp(2*a*z) + 2*sin(a*x + d*y)*cos(a*z + d*x)*exp(a*(y + z)) +  2*sin(a*y + d*z)*cos(a*x + d*y)*exp(a*(x + z)) +  2*sin(a*z + d*x)*cos(a*y + d*z)*exp(a*(y + x)))*exp(-2*d*d*theta);
    }

    PetscErrorCode CheckSolution(Vec const & sol, Vec const & solRef, std::string const & comp);

    PetscErrorCode CreateAnalyticalU(DM const & dmGrid, Vec & vec, PetscReal const & theta);
    
    PetscErrorCode CreateAnalyticalV(DM const & dmGrid, Vec & vec, PetscReal const & theta);

    PetscErrorCode CreateAnalyticalW(DM const & dmGrid, Vec & vec, PetscReal const & theta);

    PetscErrorCode CreateAnalyticalP(DM const & dmGrid, Vec & vec, PetscReal const & theta);

    PetscErrorCode CreateGrid(DM * const dmGrid, PetscInt const & dof1, PetscInt const & dof2, PetscInt const & dof3, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz);

    void PrintSimulationParameters();

    bool rayIntersectsTriangle(const std::array<double, 3>& rayOrigin, const std::array<double, 3>& rayVector, const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2);

    bool isPointInsideMesh(const std::array<double, 3>& point, const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<int, 3>>& faces);

    void reader(const std::string& filename, std::vector<std::array<double, 3>>& vertices, std::vector<std::array<int, 3>>& faces);

    PetscErrorCode createMaskU(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces);

    PetscErrorCode createMaskV(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces);

    PetscErrorCode createMaskW(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces);


}

#endif // PROBLEM_SETTING_HPP