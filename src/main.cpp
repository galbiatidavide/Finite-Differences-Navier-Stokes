/******************************************************************************
 *                                                                            *
 *  Project:    Brinkman - Navier-Stokes multiple solver                      *                  
 *  Author:     Dave & Ale                                                    *
 *  Created:    26 October 2023                                               *
 *                                                                            *
 *  Copyright © 2023 Dave & Ale. All rights reserved.                         *
 *                                                                            *
 *  Redistribution and use in source and binary forms, with or without        *
 *  modification, are permitted provided that the following conditions are    *
 *  met:                                                                      *
 *                                                                            *
 *  1. Redistributions of source code must retain the above copyright notice, *
 *     this list of conditions, and the following disclaimer.                 *
 *  2. Redistributions in binary form must reproduce the above copyright      *
 *     notice, this list of conditions, and the following disclaimer in the   *
 *     documentation and/or other materials provided with the distribution.   *
 *                                                                            *
 *  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY EXPRESS OR IMPLIED          *
 *  WARRANTIES, INCLUDING BUT NOT LIMITED TO MERCHANTABILITY OR FITNESS FOR   *
 *  A PARTICULAR PURPOSE. THE AUTHOR SHALL NOT BE HELD LIABLE FOR ANY CLAIMS  *
 *  OR DAMAGES ARISING FROM THE USE OF THIS SOFTWARE.                         *
 *                                                                            *
 ******************************************************************************/

#include "navier_stokes.hpp"
#include "inviscid_euler.hpp"
#include "stokes.hpp"
#include "advection_diffusion.hpp"



int main(int argc, char **argv)
{   

    using namespace problem_setting;
    PetscInitialize(&argc, &argv, (char*)0, (char*)0);
    //std::freopen("output.txt", "w", stdout);


    /**
     * @brief Filename for the geometry.
     */
    /*if(argc == 2 and std::string(argv[1]) != "-objects_dump" and brinkman)
    {
        filename = argv[1];
        std::cout << "Reading geometry from " << argv[1] << std::endl;
        reader(filename, vertices, faces);        
    }
    else if(argc == 2 and std::string(argv[1]) != "-objects_dump" and !brinkman)
    {
        std::cout<<"Error!. Set brinkman flag = TRUE"<<std::endl;
        exit(1);
    }*/
    
    

    #ifdef COMPILE_PARABOLIC
    #ifndef COMPILE_NAVIER_STOKES
    #ifndef COMPILE_STOKES
    #ifndef COMPILE_ADVECTION_DIFFUSION
    {
        parabolic_problem_x parabolic_x;
        parabolic_x.solve();
    }
    #endif
    #endif
    #endif
    #endif

    /*
    #ifdef COMPILE_PARABOLIC
    #ifndef COMPILE_NAVIER_STOKES
    #ifndef COMPILE_STOKES
    #ifndef COMPILE_ADVECTION_DIFFUSION
    {
        parabolic_problem_y parabolic_y;
        parabolic_y.solve();
    }
    #endif
    #endif
    #endif
    #endif

    #ifdef COMPILE_PARABOLIC
    #ifndef COMPILE_NAVIER_STOKES
    #ifndef COMPILE_STOKES
    #ifndef COMPILE_ADVECTION_DIFFUSION
    {
        parabolic_problem_z parabolic_z;
        parabolic_z.solve();
    }
    #endif
    #endif
    #endif
    #endif
    */

    #ifdef COMPILE_NAVIER_STOKES
    {
        //poisson_problem poisson;
        //poisson.manage_pressure();
        navier_stokes_problem navier_stokes;
        navier_stokes.solve();

    }
    #endif

    #ifdef COMPILE_POISSON
    {
        auto start = std::chrono::high_resolution_clock::now();
        poisson_problem poisson;
        poisson.manage_pressure();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    }
    #endif

    #ifdef COMPILE_EULER
    {
        euler_problem euler;
        euler.solve();
    }
    #endif

    #ifdef COMPILE_STOKES
    {
        stokes_problem stokes;
        stokes.solve();
    }
    #endif

    #ifdef COMPILE_ADVECTION_DIFFUSION
    {
        advection_diffusion_problem advection_diffusion;
        advection_diffusion.solve();
    }
    #endif


    PetscFinalize();
    std::cout << "Test successfully completed. Ad maiora!" << std::endl;
    //std::fclose(stdout);

    PetscFunctionReturn(0); 
}
