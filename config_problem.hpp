/// \file

#include "macros.hpp"

/**
 * @file config_problem.hpp
 * @brief Configuration file for setting up problem parameters.
 *
 * This file defines the parameters for simulations. 
 * Before starting any simulation, careful consideration of the parameters is required.
 * It is recommended to work in a dimensionless framework.
 */

#ifndef PROBLEM_SETTING_HPP
#define PROBLEM_SETTING_HPP

namespace problem_setting
{


/**
 * @brief Set this as base path to store results. Default is "results/"
 */
constexpr const char *base_path = "results/";

/**
 * @brief Set flag to (dis)able brinkmann penalization. Required for complex geometry handling.
 */
constexpr bool brinkmann{false};
/**
 * @brief Set flag to (dis)able convergence error of KSP solver.
 */
constexpr bool monitor_convergence{true};
/**
 * @brief Set flag to (dis)able convergence error of variable. Only implemented for NAvier-Stokes and parabolic problems.
 */
constexpr bool check_convergence{true};


/**
 * @brief Set the number of elements in x-direction.
 */
constexpr PetscInt nx{32};
/**
 * @brief Set the number of elements in y-direction.
 */
constexpr PetscInt ny{32};
/**
 * @brief Set the number of elements in z-direction.
 */
constexpr PetscInt nz{32};

/*constexpr PetscReal Lx_0{253};
constexpr PetscReal Ly_0{235};
constexpr PetscReal Lz_0{190};
constexpr PetscReal Lx{292};
constexpr PetscReal Ly{310};
constexpr PetscReal Lz{232};*/

/**
 * @brief Set the domain limits in x, y and z directions.
 */
constexpr PetscReal Lx_0{-0.5};
constexpr PetscReal Ly_0{-0.5};
constexpr PetscReal Lz_0{-0.5};
constexpr PetscReal Lx{0.5};
constexpr PetscReal Ly{0.5};
constexpr PetscReal Lz{0.5};

/**
 * @brief Set time-step.
 */
constexpr PetscReal dt{0.00625/2};
/**
 * @brief Set final time.
 */
constexpr PetscReal iter{32};
/**
 * @brief Set starting time.
 */
inline PetscReal theta;

/**
 * @brief Set Reynolds number. For non advective problems, set mu = 1/Re (adimensional framework)
 */
constexpr PetscReal Re{1};

/**
 * @brief Set flow parameters. You can declare as many constexpr variables as you need.
 */
constexpr PetscReal a = pi / 4;
constexpr PetscReal d = 1.5 * pi;

/**
 * @brief Computes the reference solution for the problem in the x-direction.
 *
 * This function provides an analytical reference solution for the given 
 * problem in the x-direction. It is time-dependent and useful for benchmarking 
 * numerical methods.
 *
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @param z Z-coordinate.
 * @param theta Time-dependent parameter.
 * @return Computed reference solution in the x-direction.
 */
constexpr PetscReal uxRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    return -a*(exp(a*x)*sin(a*y + d*z) + exp(a*z)*cos(a*x + d*y))*exp(-d*d*theta); //navier stokes
    //quando fai test ricorda di cambiare il valore di theta e degli altri parametri e di sistemarlo
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z; //parabolic
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z); //transport
    //return 0;
}

/**
 * @brief Computes the reference solution for the problem in the y-direction.
 *
 * This function provides an analytical reference solution for the given 
 * problem in the y-direction. It is used for validation and testing of 
 * numerical simulations.
 *
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @param z Z-coordinate.
 * @param theta Time-dependent parameter.
 * @return Computed reference solution in the y-direction.
 */
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

/**
 * @brief Computes the reference solution for the problem in the z-direction.
 *
 * This function provides an analytical reference solution for the given 
 * problem in the z-direction. It is useful for testing the consistency of 
 * numerical solvers.
 *
 * @param x X-coordinate.
 * @param y Y-coordinate.
 * @param z Z-coordinate.
 * @param theta Time-dependent parameter.
 * @return Computed reference solution in the z-direction.
 */
constexpr PetscReal uzRef(PetscReal const & x, PetscReal const & y, PetscReal const & z, PetscReal const & theta)
{
    return -a*(exp(a*z)*sin(a*x + d*y) + exp(a*y)*cos(a*z + d*x))*exp(-d*d*theta); //navier stokes
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z; // parabolic
    //return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z); // transport
    //return 0;
}




}

#endif // PROBLEM_SETTING_HPP