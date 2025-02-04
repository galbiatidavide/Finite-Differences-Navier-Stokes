
#include "macros.hpp"

#ifndef PROBLEM_SETTING_HPP
#define PROBLEM_SETTING_HPP

namespace problem_setting
{

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
    PetscReal constexpr D_x{Lx - Lx_0};
    PetscReal constexpr D_y{Ly - Ly_0};
    PetscReal constexpr D_z{Lz - Lz_0};

    PetscReal constexpr dt{0.00625/2};
    PetscReal constexpr iter{32};

    PetscReal constexpr Re{1};

    inline bool brinkmann{false};
    inline bool monitor_convergence{false};
    inline bool check_convergence{false};


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
        /*if(y == 310)
        {
            return -1.0;
        }
        else if (y == 235)
        {
            if((x - 265.8)*(x - 265.8) + (z - 206.4)*(z - 206.4) <= 3.5)
            {
                return -1.18;
            }
            else if((x - 268.58)*(x - 268.58) + (z - 215.03)*(z - 215.03) <= 3.5)
            {
                return -1.20;
            }
            else {
                return 0;
            }
        }
        else 
        {
            return 0;
        }*/


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
}

#endif // PROBLEM_SETTING_HPP