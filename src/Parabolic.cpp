#include "Parabolic.hpp"

PetscErrorCode Parabolic::assemble_matrices(){

PetscFunctionBegin;

for(unsigned int i = 0; i < grid->components.size(); i++){

            PetscInt        start[3], n[3], N[3], ex, ey, ez;

            DMStagGetCorners(grid->dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], NULL, NULL, NULL);
            DMStagGetGlobalSizes(grid->dmGrid, &N[0], &N[1], &N[2]);


            if(grid -> components[i].name == "x_component"){
            // std::cout << "----------------------------" << std::endl;
            // std::cout << "lhs assemble for x_component" << std::endl;
            for (ez = start[2]; ez < start[2] + n[2]; ++ez) {
                for (ey = start[1]; ey < start[1] + n[1]; ++ey) {
                    for (ex = start[0]; ex < start[0] + n[0]; ++ex) {

                /* Right Boundary velocity Dirichlet */
                if (ex == N[0] - 1) {
                    DMStagStencil row;
                    const PetscReal valA = 1.0;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = RIGHT;
                    row.c = 0;
                    DMStagMatSetValuesStencil(grid->dmGrid, lhs_comp[i], 1, &row, 1, &row, &valA, INSERT_VALUES);
                }      
                /* Equation on left face of this element */
                if (ex == 0) {
                    DMStagStencil row;
                    const PetscReal valA = 1.0;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = LEFT;
                    row.c = 0;
                    DMStagMatSetValuesStencil(grid->dmGrid, lhs_comp[i], 1, &row, 1, &row, &valA, INSERT_VALUES);
                    
                    } else {
                    DMStagStencil row, col[7];
                    PetscReal valA[7];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = LEFT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);

                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = LEFT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);                          
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);                      
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = LEFT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = LEFT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = LEFT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = LEFT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = LEFT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);                          
                        }
                    } else if (ez == 0) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = LEFT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = LEFT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = LEFT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = LEFT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = LEFT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);                     
                    } else if (ez == N[2] - 1) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = LEFT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = LEFT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = LEFT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = LEFT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = LEFT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);                    
                    } else {
                        nEntries = 7;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = LEFT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = LEFT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = LEFT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = LEFT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = LEFT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = LEFT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        col[6].i = ex;
                        col[6].j = ey;
                        col[6].k = ez + 1;
                        col[6].loc = LEFT;
                        col[6].c = 0;
                        valA[6] = 1.0 / (hz * hz);
                    }

                    DMStagMatSetValuesStencil(grid->dmGrid, lhs_comp[i], 1, &row, nEntries, col, valA, INSERT_VALUES);

                } //else
            } //for ex
        } //for ey
    } //for ez

    MatAssemblyBegin(lhs_comp[i], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(lhs_comp[i], MAT_FINAL_ASSEMBLY);
    // std::cout << "lhs assembled for x_component" << std::endl;
    // std::cout << "-----------------------------" << std::endl;

    } //if x_component

    if(grid->components[i].name == "y_component"){

    // std::cout << "lhs assemble for y_component" << std::endl;

        for (ez = start[2]; ez < start[2] + n[2]; ++ez) {
            for (ey = start[1]; ey < start[1] + n[1]; ++ey) {
                for (ex = start[0]; ex < start[0] + n[0]; ++ex) {

         if (ey == N[1] - 1) {
                    /* Top boundary velocity Dirichlet */
                    DMStagStencil     row;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = UP;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(grid -> dmGrid, lhs_comp[i], 1, &row, 1, &row, &valA, INSERT_VALUES);
                }                
                /* Equation on bottom face of this element */
                if (ey == 0) {
                    /* Bottom boundary velocity Dirichlet */
                    DMStagStencil     row;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = DOWN;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(grid -> dmGrid, lhs_comp[i], 1, &row, 1, &row, &valA, INSERT_VALUES);
                } else {
                    /* Y-momentum equation, (v_xx + v_yy + v_zz) - p_y = f^y */
                    DMStagStencil row, col[7];
                    PetscReal   valA[7];
                    PetscInt      nEntries;
                    row.i   = ex;
                    row.j   = ey;
                    row.k   = ez;
                    row.loc = DOWN;
                    row.c   = 0;
                    if (ex == 0) {
                        if (ez == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = DOWN;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz);
                        }
                    } else if (ex == N[0] - 1) {
                        if (ez == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex - 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);                        
                        } else if (ez == N[2] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex - 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = DOWN;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = DOWN;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = DOWN;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex - 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = DOWN;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = DOWN;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = DOWN;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = DOWN;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = DOWN;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = DOWN;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = DOWN;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = DOWN;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = DOWN;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = DOWN;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = DOWN;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = DOWN;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = DOWN;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = DOWN;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez - 1;
                        col[5].loc = DOWN;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);                  
                    } else {
                        nEntries   = 7;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = DOWN;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = DOWN;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = DOWN;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = DOWN;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = DOWN;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez - 1;
                        col[5].loc = DOWN;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        col[6].i   = ex;
                        col[6].j   = ey;
                        col[6].k   = ez + 1;
                        col[6].loc = DOWN;
                        col[6].c   = 0;
                        valA[6]    = 1.0 / (hz * hz);                       
                    }

                    DMStagMatSetValuesStencil(grid->dmGrid, lhs_comp[i], 1, &row, nEntries, col, valA, INSERT_VALUES);

                } //else
                
            } //for ez
        } //for ey
    } //for ex

    MatAssemblyBegin(lhs_comp[i], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(lhs_comp[i], MAT_FINAL_ASSEMBLY);

    // std::cout << "lhs assembled for y_component" << std::endl;
    // std::cout << "-----------------------------" << std::endl;

} //if y_component

if(grid->components[i].name == "z_component"){

// std::cout << "lhs assemble for z_component" << std::endl;

for (ez = start[2]; ez < start[2] + n[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0]; ++ex) {              
                if (ez == N[2] - 1) {
                    /* Front boundary velocity Dirichlet */
                    DMStagStencil     row;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = FRONT;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(grid->dmGrid, lhs_comp[i], 1, &row, 1, &row, &valA, INSERT_VALUES);
                }                
                /* Equation on back face of this element */
                if (ez == 0) {
                    /* Back boundary velocity Dirichlet */
                    DMStagStencil     row;
                    const PetscReal valA = 1.0;
                    row.i                  = ex;
                    row.j                  = ey;
                    row.k                  = ez;
                    row.loc                = BACK;
                    row.c                  = 0;
                    DMStagMatSetValuesStencil(grid->dmGrid, lhs_comp[i], 1, &row, 1, &row, &valA, INSERT_VALUES);
                } else {
                    /* Z-momentum equation, (w_xx + w_yy + w_zz) - p_z = f^z */
                    DMStagStencil row, col[7];
                    PetscReal   valA[7];
                    PetscInt      nEntries;
                    row.i   = ex;
                    row.j   = ey;
                    row.k   = ez;
                    row.loc = BACK;
                    row.c   = 0;
                    if (ex == 0) {
                        if (ey == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) - 2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz); 
                        } else if (ey == N[1] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex + 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);                            
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex + 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = BACK;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz);
                        }
                    } else if (ex == N[0] - 1) {
                        if (ey == 0) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey + 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz); 
                        } else if (ey == N[1] - 1) {
                            nEntries   = 5;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex - 1;
                            col[2].j   = ey;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hx * hx);
                            col[3].i   = ex;
                            col[3].j   = ey;
                            col[3].k   = ez - 1;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hz * hz);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez + 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                        } else {
                            nEntries   = 6;
                            col[0].i   = ex;
                            col[0].j   = ey;
                            col[0].k   = ez;
                            col[0].loc = BACK;
                            col[0].c   = 0;
                            valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                            col[1].i   = ex;
                            col[1].j   = ey - 1;
                            col[1].k   = ez;
                            col[1].loc = BACK;
                            col[1].c   = 0;
                            valA[1]    = 1.0 / (hy * hy);
                            col[2].i   = ex;
                            col[2].j   = ey + 1;
                            col[2].k   = ez;
                            col[2].loc = BACK;
                            col[2].c   = 0;
                            valA[2]    = 1.0 / (hy * hy);
                            col[3].i   = ex - 1;
                            col[3].j   = ey;
                            col[3].k   = ez;
                            col[3].loc = BACK;
                            col[3].c   = 0;
                            valA[3]    = 1.0 / (hx * hx);
                            col[4].i   = ex;
                            col[4].j   = ey;
                            col[4].k   = ez - 1;
                            col[4].loc = BACK;
                            col[4].c   = 0;
                            valA[4]    = 1.0 / (hz * hz);
                            col[5].i   = ex;
                            col[5].j   = ey;
                            col[5].k   = ez + 1;
                            col[5].loc = BACK;
                            col[5].c   = 0;
                            valA[5]    = 1.0 / (hz * hz); 
                        }
                    } else if (ey == 0) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey + 1;
                        col[1].k   = ez;
                        col[1].loc = BACK;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex - 1;
                        col[2].j   = ey;
                        col[2].k   = ez;
                        col[2].loc = BACK;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hx * hx);
                        col[3].i   = ex + 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = BACK;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez - 1;
                        col[4].loc = BACK;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = BACK;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);                       
                    } else if (ey == N[1] - 1) {
                        nEntries   = 6;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) - 2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = BACK;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex - 1;
                        col[2].j   = ey;
                        col[2].k   = ez;
                        col[2].loc = BACK;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hx * hx);
                        col[3].i   = ex + 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = BACK;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex;
                        col[4].j   = ey;
                        col[4].k   = ez - 1;
                        col[4].loc = BACK;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hz * hz);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez + 1;
                        col[5].loc = BACK;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                    } else {
                        nEntries   = 7;
                        col[0].i   = ex;
                        col[0].j   = ey;
                        col[0].k   = ez;
                        col[0].loc = BACK;
                        col[0].c   = 0;
                        valA[0]    = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz) - Ret;
                        col[1].i   = ex;
                        col[1].j   = ey - 1;
                        col[1].k   = ez;
                        col[1].loc = BACK;
                        col[1].c   = 0;
                        valA[1]    = 1.0 / (hy * hy);
                        col[2].i   = ex;
                        col[2].j   = ey + 1;
                        col[2].k   = ez;
                        col[2].loc = BACK;
                        col[2].c   = 0;
                        valA[2]    = 1.0 / (hy * hy);
                        col[3].i   = ex - 1;
                        col[3].j   = ey;
                        col[3].k   = ez;
                        col[3].loc = BACK;
                        col[3].c   = 0;
                        valA[3]    = 1.0 / (hx * hx);
                        col[4].i   = ex + 1;
                        col[4].j   = ey;
                        col[4].k   = ez;
                        col[4].loc = BACK;
                        col[4].c   = 0;
                        valA[4]    = 1.0 / (hx * hx);
                        col[5].i   = ex;
                        col[5].j   = ey;
                        col[5].k   = ez - 1;
                        col[5].loc = BACK;
                        col[5].c   = 0;
                        valA[5]    = 1.0 / (hz * hz);
                        col[6].i   = ex;
                        col[6].j   = ey;
                        col[6].k   = ez + 1;
                        col[6].loc = BACK;
                        col[6].c   = 0;
                        valA[6]    = 1.0 / (hz * hz);                         
                    }

                    DMStagMatSetValuesStencil(grid->dmGrid, lhs_comp[i], 1, &row, nEntries, col, valA, INSERT_VALUES);
                } //else
            } //for ex
        } //for ey
    } //for ez

    MatAssemblyBegin(lhs_comp[i], MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(lhs_comp[i], MAT_FINAL_ASSEMBLY);

    // std::cout << "lhs assembled for z_component" << std::endl;
    // std::cout << "-----------------------------" << std::endl;

} //if z_component
} //for componenti 

std::cout << "lhs assembled for all components" << std::endl;
std::cout << "---------------------------------" << std::endl;

PetscFunctionReturn(0); 

}


PetscErrorCode Parabolic::assemble_rhs(const unsigned int &time_step){

    PetscReal theta = grid->bc.get_time(time_step);
    std::cout << "theta = " << theta << " at time " << time << std::endl;

    PetscFunctionBegin;
    
    for(unsigned int i = 0; i < grid->components.size(); i++){

        PetscInt        start[3], n[3], nExtra[3], N[3], ex, ey, ez, iuxStart, iuxEnd, icuxStart[3], icuxEnd[3];
        DM              dmCoord;
        Vec             coord, coordLocal;
        PetscReal       ****arrCoord, ****arrVecLocal_n, ****arrRhs;


        DMStagGetCorners(grid->dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
        DMStagGetGlobalSizes(grid->dmGrid, &N[0], &N[1], &N[2]);

        DMGetCoordinateDM(grid->dmGrid, &dmCoord);
        DMGetCoordinates(grid->dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);


        for (int d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, grid->components[i].location[0], d, &icuxStart[d]);
            DMStagGetLocationSlot(dmCoord, grid->components[i].location[1], d, &icuxEnd[d]);
        }

        DMStagGetLocationSlot(grid->dmGrid, grid->components[i].location[0], 0, &iuxStart);
        DMStagGetLocationSlot(grid->dmGrid, grid->components[i].location[1], 0, &iuxEnd);

        Vec local_n;
        DMCreateLocalVector(grid->dmGrid, &local_n);
        DMGlobalToLocalBegin(grid->dmGrid, n_components[i].variable, INSERT_VALUES, local_n);
        DMGlobalToLocalEnd(grid->dmGrid, n_components[i].variable, INSERT_VALUES, local_n);
        DMStagVecGetArrayRead(grid->dmGrid, local_n, &arrVecLocal_n);

        Vec vecLocalRhs;
        DMGetLocalVector(grid->dmGrid, &vecLocalRhs);
        DMStagVecGetArray(grid->dmGrid, vecLocalRhs, &arrRhs);

       //Set rhs of x_component
        if (grid->components[i].name == "x_component") {

            for (ez = start[2]; ez < start[2] + n[2]; ++ez) {
                for (ey = start[1]; ey < start[1] + n[1]; ++ey) {
                    for(ex = start[0]; ex < start[0] +n[0]; ++ex){
                

                if (ex == N[0] - 1) {
                    arrRhs[ez][ey][ex][iuxEnd] = exactSolution[i](arrCoord[ez][ey][ex][icuxEnd[0]], arrCoord[ez][ey][ex][icuxEnd[1]], arrCoord[ez][ey][ex][icuxEnd[2]], theta);}   
                if (ex == 0) {
                    arrRhs[ez][ey][ex][iuxStart] = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);} 
                else {
                    if (ey == 0) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]]-hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]-hz, theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hy*hy) - bc_2/(hz*hz);                             
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]] - hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]+hz, theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hy*hy) - bc_2/(hz*hz);
                        } else {
                            PetscReal bc_2;
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]] - hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_2/(hy*hy);                                                     
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]]+hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]-hz, theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hy*hy) - bc_2/(hz*hz);                                                  
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]] + hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]] + hz, theta);                            
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hy*hy) - bc_2/(hz*hz);       
                        } else {
                            PetscReal bc_2;
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]]+hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_2/(hy*hy);                            
                        }
                    } else if (ez == 0) {
                        PetscReal bc_1;
                        bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]-hz, theta);
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hz*hz);                       
                    } else if (ez == N[2] - 1) {
                        PetscReal bc_1;
                        bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]+hz, theta);
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hz*hz);                        
                    } else {
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart];
                    }

                } //else

            } //for ex
        } //for ey
    } //for ez

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(grid->dmGrid, local_n, &arrVecLocal_n);
    DMStagVecRestoreArray(grid->dmGrid, vecLocalRhs, &arrRhs);
    DMLocalToGlobal(grid->dmGrid, vecLocalRhs, INSERT_VALUES, rhs_comp[i]);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(grid->dmGrid, &local_n);
    DMRestoreLocalVector(grid->dmGrid, &vecLocalRhs);

        // PetscViewer viewer_2;
        // Vec r;
        // DM pda;

        // DMStagVecSplitToDMDA(grid->dmGrid, rhs_comp[i], grid->components[i].location[0], DM_BOUNDARY_NONE, &pda, &r);
        // PetscObjectSetName((PetscObject)r, "rhs");  // Set name of vector

        // char filename_r[50];
        // sprintf(filename_r, "results/rhs_00%i.txt", i);  // Change extension to .txt
        // PetscViewerASCIIOpen(PetscObjectComm((PetscObject)pda), filename_r, &viewer_2); // Use ASCII viewer

        // VecView(r, viewer_2);  // View the vector contents in text format

        // // Cleanup
        // VecDestroy(&r);
        // DMDestroy(&pda);
        // PetscViewerDestroy(&viewer_2);

    } //if x_component

    if(grid->components[i].name == "y_component"){

        for (ez = start[2]; ez < start[2] + n[2]; ++ez) {
                for (ey = start[1]; ey < start[1] + n[1]; ++ey) {
                    for(ex = start[0]; ex < start[0] +n[0]; ++ex){
                
                if (ey == N[1] - 1) {
                    arrRhs[ez][ey][ex][iuxEnd] = exactSolution[i](arrCoord[ez][ey][ex][icuxEnd[0]], arrCoord[ez][ey][ex][icuxEnd[1]], arrCoord[ez][ey][ex][icuxEnd[2]], theta);}   
                if (ey == 0) {
                    arrRhs[ez][ey][ex][iuxStart] = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);} 
                else {
                    if (ex == 0) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]-hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]-hz, theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx) - bc_2/(hz*hz);                             
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]-hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]+hz, theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx) - bc_2/(hz*hz);
                        } else {
                            PetscReal bc_1;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]-hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx);                                                     
                        }
                    } else if (ex == N[0] - 1) {
                        if (ez == 0) {
                            PetscReal bc_1, bc_2;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]+hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]-hz, theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx) - bc_2/(hz*hz);                                                  
                        } else if (ez == N[2] - 1) {
                            PetscReal bc_1, bc_2;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]+hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]+hz, theta);                           
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx) - bc_2/(hz*hz);       
                        } else {
                            PetscReal bc_1;
                            bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]+hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                            arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx);                            
                        }
                    } else if (ez == 0) {
                        PetscReal bc_2;
                        bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]-hz, theta);
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_2/(hz*hz);                       
                    } else if (ez == N[2] - 1) {
                        PetscReal bc_2;
                        bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]]+hz, theta);
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_2/(hz*hz);                        
                    } else {
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart];
                    }

                } //else

            } //for ex
        } //for ey
    } //for ez

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(grid->dmGrid, local_n, &arrVecLocal_n);
    DMStagVecRestoreArray(grid->dmGrid, vecLocalRhs, &arrRhs);
    DMLocalToGlobal(grid->dmGrid, vecLocalRhs, INSERT_VALUES, rhs_comp[i]);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(grid->dmGrid, &local_n);
    DMRestoreLocalVector(grid->dmGrid, &vecLocalRhs);

    // PetscViewer viewer_2;
    // Vec r;
    // DM pda;

    // DMStagVecSplitToDMDA(grid->dmGrid, rhs_comp[i], grid->components[i].location[0], DM_BOUNDARY_NONE, &pda, &r);
    // PetscObjectSetName((PetscObject)r, "rhs");  // Set name of vector

    // char filename_r[50];
    // sprintf(filename_r, "results/rhs_00%i.txt", i);  // Change extension to .txt
    // PetscViewerASCIIOpen(PetscObjectComm((PetscObject)pda), filename_r, &viewer_2); // Use ASCII viewer

    // VecView(r, viewer_2);  // View the vector contents in text format

    // // Cleanup
    // VecDestroy(&r);
    // DMDestroy(&pda);
    // PetscViewerDestroy(&viewer_2);

    } //if y_component


    if(grid->components[i].name == "z_component"){

    for (ez = start[2]; ez < start[2] + n[2]; ++ez) {
            for (ey = start[1]; ey < start[1] + n[1]; ++ey) {
                for(ex = start[0]; ex < start[0] +n[0]; ++ex){
            
            if (ez == N[2] - 1) {
                arrRhs[ez][ey][ex][iuxEnd] = exactSolution[i](arrCoord[ez][ey][ex][icuxEnd[0]], arrCoord[ez][ey][ex][icuxEnd[1]], arrCoord[ez][ey][ex][icuxEnd[2]], theta);}   
            if (ez == 0) {
                arrRhs[ez][ey][ex][iuxStart] = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);} 
            else {
                if (ex == 0) {
                    if (ey == 0) {
                        PetscReal bc_1, bc_2;
                        bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]-hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                        bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]]-hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx) - bc_2/(hy*hy);                             
                    } else if (ey == N[1] - 1) {
                        PetscReal bc_1, bc_2;
                        bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]-hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                        bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]]+hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx) - bc_2/(hy*hy);
                    } else {
                        PetscReal bc_1;
                        bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]-hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx);                                                     
                    }
                } else if (ex == N[0] - 1) {
                    if (ey == 0) {
                        PetscReal bc_1, bc_2;
                        bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]+hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                        bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]]-hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx) - bc_2/(hy*hy);                                                  
                    } else if (ey == N[1] - 1) {
                        PetscReal bc_1, bc_2;
                        bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]+hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                        bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]]+hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);                           
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx) - bc_2/(hy*hy);       
                    } else {
                        PetscReal bc_1;
                        bc_1 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]]+hx, arrCoord[ez][ey][ex][icuxStart[1]], arrCoord[ez][ey][ex][icuxStart[2]], theta);
                        arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_1/(hx*hx);                            
                    }
                } else if (ey == 0) {
                    PetscReal bc_2;
                    bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]]-hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                    arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_2/(hy*hy);                       
                } else if (ey == N[1] - 1) {
                    PetscReal bc_2;
                    bc_2 = exactSolution[i](arrCoord[ez][ey][ex][icuxStart[0]], arrCoord[ez][ey][ex][icuxStart[1]]+hy, arrCoord[ez][ey][ex][icuxStart[2]], theta);
                    arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart] - bc_2/(hy*hy);                        
                } else {
                    arrRhs[ez][ey][ex][iuxStart] = -Ret*arrVecLocal_n[ez][ey][ex][iuxStart];
                }

            } //else

        } //for ex
    } //for ey
} //for ez

DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
DMStagVecRestoreArrayRead(grid->dmGrid, local_n, &arrVecLocal_n);
DMStagVecRestoreArray(grid->dmGrid, vecLocalRhs, &arrRhs);
DMLocalToGlobal(grid->dmGrid, vecLocalRhs, INSERT_VALUES, rhs_comp[i]);
DMRestoreLocalVector(dmCoord, &coordLocal);
DMRestoreLocalVector(grid->dmGrid, &local_n);
DMRestoreLocalVector(grid->dmGrid, &vecLocalRhs);

    // PetscViewer viewer_2;
    // Vec r;
    // DM pda;

    // DMStagVecSplitToDMDA(grid->dmGrid, rhs_comp[i], grid->components[i].location[0], DM_BOUNDARY_NONE, &pda, &r);
    // PetscObjectSetName((PetscObject)r, "rhs");  // Set name of vector

    // char filename_r[50];
    // sprintf(filename_r, "results/rhs_00%i.txt", i);  // Change extension to .txt
    // PetscViewerASCIIOpen(PetscObjectComm((PetscObject)pda), filename_r, &viewer_2); // Use ASCII viewer

    // VecView(r, viewer_2);  // View the vector contents in text format

    // // Cleanup
    // VecDestroy(&r);
    // DMDestroy(&pda);
    // PetscViewerDestroy(&viewer_2);

} //if z_component

} //for componenti

std::cout << "rhs assembled for all components at time " << theta << std::endl;
std::cout << "-----------------------------------------" << std::endl;


PetscFunctionReturn(0);

}


PetscErrorCode Parabolic::solveTimeStep(const unsigned int &time_step){

    PetscFunctionBegin;

    for (size_t i = 0; i < grid->components.size(); ++i) {  
    KSP       ksp;
    PC        pc;
    
    assemble_rhs(time_step);
    output_rhs(time_step);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    saveMatrices();
    KSPSetOperators(ksp, lhs_comp[i], lhs_comp[i]);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCFIELDSPLIT);
    PCFieldSplitSetDetectSaddlePoint(pc, PETSC_TRUE);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, rhs_comp[i], grid->components[i].variable);
    KSPDestroy(&ksp);
    VecCopy(grid->components[i].variable, n_components[i].variable);
    output(time_step);
    }

    PetscFunctionReturn(0);
}


PetscErrorCode Parabolic::Solve(){
    
        PetscFunctionBegin;

        assemble_matrices();
        //saveMatrices();

        PetscInt timeStep = 0;
    
        while(timeStep*dt < T-0.5){
            solveTimeStep(timeStep);
            timeStep++;
        }
    
        PetscFunctionReturn(0);
}



PetscErrorCode Parabolic::saveMatrices() {
    PetscFunctionBegin;

    for (size_t i = 0; i < lhs_comp.size(); ++i) {

        std::string filename =  "dio_" + std::to_string(i) + ".txt";
        PetscViewer viewer;

        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename.c_str(), &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
        
        MatView(lhs_comp[i], viewer);

        PetscViewerPopFormat(viewer);
        PetscViewerDestroy(&viewer); 

    }

    PetscFunctionReturn(0);
}



