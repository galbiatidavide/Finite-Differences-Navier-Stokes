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

#include <chrono>
#include <iostream>
#include <limits>
//PETSc
#include <petscdmstag.h>
#include <petscksp.h>
//VTK
#include <vtkSTLReader.h>

#include "parameters.hpp"
#include "navier_stokes.hpp"
#include "reader.cpp"


int main(int argc, char **argv)
{   

    auto start = std::chrono::high_resolution_clock::now();

    PetscInitialize(&argc, &argv, (char*)0, (char*)0);

    if(argc == 2)
    {
        /*std::vector<std::array<double, 3>> vertices;
        std::vector<std::array<int, 3>> faces;
        std::string filename = argv[1];
        reader(filename, vertices, faces);*/
    }

    {
        navier_stokes_problem navier_stokes;
        navier_stokes.solve();
    }

    PetscFinalize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    std::cout << "Test successfully completed. Ad maiora!" << std::endl;
    PetscFunctionReturn(0); 
}



