#include "poisson.hpp"


poisson_problem::poisson_problem(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, DM const & dmGrid_centered, DM const & dmGrid_cent_rich)
    : dmGrid_staggered_x(dmGrid_staggered_x), dmGrid_staggered_y(dmGrid_staggered_y), dmGrid_staggered_z(dmGrid_staggered_z), dmGrid_centered(dmGrid_centered), dmGrid_cent_rich(dmGrid_cent_rich)

{
    DMCreateGlobalVector(dmGrid_centered, &P);
    DMCreateGlobalVector(dmGrid_staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_staggered_z, &P_z);
    DMCreateGlobalVector(dmGrid_staggered_x, &U_up);
    DMCreateGlobalVector(dmGrid_staggered_y, &V_up);
    DMCreateGlobalVector(dmGrid_staggered_z, &W_up);

    DMCreateMatrix(dmGrid_centered, &A);
    this->assemble_lhs();
}

poisson_problem::poisson_problem()
{
    //Allocate the grids
    CreateGrid(&dmGrid_staggered_x, 0, 1, 0);
    CreateGrid(&dmGrid_staggered_y, 0, 1, 0);
    CreateGrid(&dmGrid_staggered_z, 0, 1, 0);
    CreateGrid(&dmGrid_centered, 0, 0, 1);
    CreateGrid(&dmGrid_cent_rich, 0, 1, 1);

    //Create parallel vectors
    DMCreateGlobalVector(dmGrid_staggered_x, &U_up);
    DMCreateGlobalVector(dmGrid_staggered_y, &V_up);
    DMCreateGlobalVector(dmGrid_staggered_z, &W_up);
    CreateAnalyticalU(dmGrid_staggered_x, U_up, 0);
    CreateAnalyticalV(dmGrid_staggered_y, V_up, 0);
    CreateAnalyticalW(dmGrid_staggered_z, W_up, 0);

    DMCreateGlobalVector(dmGrid_centered, &P);
    DMCreateGlobalVector(dmGrid_staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_staggered_z, &P_z);
    DMCreateMatrix(dmGrid_centered, &A);

    this->assemble_lhs();

};

PetscErrorCode const poisson_problem::assemble_lhs() 
{
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;

    PetscFunctionBegin;
    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    DMStagGetCorners(dmGrid_centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_centered, &N[0], &N[1], &N[2]);
    PetscReal const hx = D_x/ N[0];
    PetscReal const hy = D_y/ N[1];
    PetscReal const hz = D_z/ N[2];

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ex == N[0] - 1) {

                    DMStagStencil row, col[6];
                    PetscReal valA[6];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {                  
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx); 
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);                            
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                } 
                else if (ex == 0) {
                    DMStagStencil row, col[6];
                    PetscReal valA[6];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez + 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 4;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex + 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex;
                            col[3].j = ey;
                            col[3].k = ez - 1;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hz * hz);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez + 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries = 5;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                    } else {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -1.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex + 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex;
                        col[4].j = ey;
                        col[4].k = ez - 1;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                } else {
                    DMStagStencil row, col[7];
                    PetscReal valA[7];
                    PetscInt nEntries;
                    row.i = ex;
                    row.j = ey;
                    row.k = ez;
                    row.loc = ELEMENT;
                    row.c = 0;
                    if (ey == 0) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey + 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                        }
                    } else if (ey == N[1] - 1) {
                        if (ez == 0) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez + 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else if (ez == N[2] - 1) {
                            nEntries = 5;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                        } else {
                            nEntries = 6;
                            col[0].i = ex;
                            col[0].j = ey;
                            col[0].k = ez;
                            col[0].loc = ELEMENT;
                            col[0].c = 0;
                            valA[0] = -2.0 / (hx * hx) + -1.0 / (hy * hy) - 2.0 / (hz * hz);
                            col[1].i = ex;
                            col[1].j = ey - 1;
                            col[1].k = ez;
                            col[1].loc = ELEMENT;
                            col[1].c = 0;
                            valA[1] = 1.0 / (hy * hy);
                            col[2].i = ex - 1;
                            col[2].j = ey;
                            col[2].k = ez;
                            col[2].loc = ELEMENT;
                            col[2].c = 0;
                            valA[2] = 1.0 / (hx * hx);
                            col[3].i = ex + 1;
                            col[3].j = ey;
                            col[3].k = ez;
                            col[3].loc = ELEMENT;
                            col[3].c = 0;
                            valA[3] = 1.0 / (hx * hx);
                            col[4].i = ex;
                            col[4].j = ey;
                            col[4].k = ez - 1;
                            col[4].loc = ELEMENT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                        }
                    } else if (ez == 0) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    } else if (ez == N[2] - 1) {
                        nEntries = 6;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 1.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    } else {
                        nEntries = 7;
                        col[0].i = ex;
                        col[0].j = ey;
                        col[0].k = ez;
                        col[0].loc = ELEMENT;
                        col[0].c = 0;
                        valA[0] = -2.0 / (hx * hx) + -2.0 / (hy * hy) - 2.0 / (hz * hz);
                        col[1].i = ex;
                        col[1].j = ey - 1;
                        col[1].k = ez;
                        col[1].loc = ELEMENT;
                        col[1].c = 0;
                        valA[1] = 1.0 / (hy * hy);
                        col[2].i = ex;
                        col[2].j = ey + 1;
                        col[2].k = ez;
                        col[2].loc = ELEMENT;
                        col[2].c = 0;
                        valA[2] = 1.0 / (hy * hy);
                        col[3].i = ex - 1;
                        col[3].j = ey;
                        col[3].k = ez;
                        col[3].loc = ELEMENT;
                        col[3].c = 0;
                        valA[3] = 1.0 / (hx * hx);
                        col[4].i = ex + 1;
                        col[4].j = ey;
                        col[4].k = ez;
                        col[4].loc = ELEMENT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hx * hx);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez - 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                        col[6].i = ex;
                        col[6].j = ey;
                        col[6].k = ez + 1;
                        col[6].loc = ELEMENT;
                        col[6].c = 0;
                        valA[6] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }                
            }
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    PetscFunctionReturn(0); 
}

PetscErrorCode const poisson_problem::assemble_divergence(Vec & div, Vec const & U, Vec const &  V, Vec const & W) 
{
    PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_cent_rich, &N[0], &N[1], &N[2]);
    PetscReal const hx = D_x/ N[0];
    PetscReal const hy = D_y/ N[1];
    PetscReal const hz = D_z/ N[2];
    DMGetCoordinateDM(dmGrid_cent_rich, &dmCoord);

    DMGetCoordinates(dmGrid_cent_rich, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_cent_rich, LEFT, 0, &iu_left);
    DMStagGetLocationSlot(dmGrid_cent_rich, RIGHT, 0, &iu_right);
    DMStagGetLocationSlot(dmGrid_cent_rich, UP, 0, &iu_up);
    DMStagGetLocationSlot(dmGrid_cent_rich, DOWN, 0, &iu_down);
    DMStagGetLocationSlot(dmGrid_cent_rich, FRONT, 0, &iu_front);
    DMStagGetLocationSlot(dmGrid_cent_rich, BACK, 0, &iu_back);
    DMStagGetLocationSlot(dmGrid_cent_rich, ELEMENT, 0, &iu_element);

    DMCreateLocalVector(dmGrid_cent_rich, &vecULocal);
    DMGlobalToLocalBegin(dmGrid_cent_rich, U, INSERT_VALUES, vecULocal);
    DMGlobalToLocalEnd(dmGrid_cent_rich, U, INSERT_VALUES, vecULocal);
    DMStagVecGetArrayRead(dmGrid_cent_rich, vecULocal, &arrU);

    DMCreateLocalVector(dmGrid_cent_rich, &vecVLocal);
    DMGlobalToLocalBegin(dmGrid_cent_rich, V, INSERT_VALUES, vecVLocal);
    DMGlobalToLocalEnd(dmGrid_cent_rich, V, INSERT_VALUES, vecVLocal);
    DMStagVecGetArrayRead(dmGrid_cent_rich, vecVLocal, &arrV);

    DMCreateLocalVector(dmGrid_cent_rich, &vecWLocal);
    DMGlobalToLocalBegin(dmGrid_cent_rich, W, INSERT_VALUES, vecWLocal);
    DMGlobalToLocalEnd(dmGrid_cent_rich, W, INSERT_VALUES, vecWLocal);
    DMStagVecGetArrayRead(dmGrid_cent_rich, vecWLocal, &arrW);    

    DMGetLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMStagVecGetArray(dmGrid_cent_rich, vecOutLocal, &arrOut);

    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                PetscReal inter, left, right, up, down, front, back;
                left = arrU[ez][ey][ex][iu_left];
                right = arrU[ez][ey][ex][iu_right];
                up = arrV[ez][ey][ex][iu_up];
                down = arrV[ez][ey][ex][iu_down];
                front = arrW[ez][ey][ex][iu_front];
                back = arrW[ez][ey][ex][iu_back];
                inter = ((up - down) / hy + (right - left) / hx + (front - back) / hz)/dt;
                arrOut[ez][ey][ex][iu_element] = inter;
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecULocal, &arrU);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecVLocal, &arrV);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecWLocal, &arrW);
    DMStagVecRestoreArray(dmGrid_cent_rich, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_cent_rich, vecOutLocal, INSERT_VALUES, div);    
    DMRestoreLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecULocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecVLocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecWLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const poisson_problem::compute_divergence(Vec & div, Vec const & U_n, Vec const & V_n, Vec const & W_n) 
{
        PetscFunctionBegin;

        Vec U_shifted, V_shifted, W_shifted;
        DMCreateGlobalVector(dmGrid_cent_rich, &U_shifted);
        DMCreateGlobalVector(dmGrid_cent_rich, &V_shifted);
        DMCreateGlobalVector(dmGrid_cent_rich, &W_shifted);
        DMStagMigrateVec(dmGrid_staggered_x, U_n, dmGrid_cent_rich, U_shifted);
        DMStagMigrateVec(dmGrid_staggered_y, V_n, dmGrid_cent_rich, V_shifted);
        DMStagMigrateVec(dmGrid_staggered_z, W_n, dmGrid_cent_rich, W_shifted);

        Vec div_shifted;
        DMCreateGlobalVector(dmGrid_cent_rich, &div_shifted);
        assemble_divergence(div_shifted, U_shifted, V_shifted, W_shifted);
        DMStagMigrateVec(dmGrid_cent_rich, div_shifted, dmGrid_centered, div);

        VecDestroy(&U_shifted);
        VecDestroy(&V_shifted);
        VecDestroy(&W_shifted);
        VecDestroy(&div_shifted);

        PetscFunctionReturn(0); 
}

PetscErrorCode const poisson_problem::derive_x_P(Vec & P_x_shifted, Vec const & vec)
{
    PetscInt iux_left, iux_right, iux_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_cent_rich, &N[0], &N[1], &N[2]);
    PetscReal const hx = D_x/ N[0];
    DMGetCoordinateDM(dmGrid_cent_rich, &dmCoord);

    DMGetCoordinates(dmGrid_cent_rich, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_cent_rich, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid_cent_rich, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid_cent_rich, ELEMENT, 0, &iux_element);

    DMCreateLocalVector(dmGrid_cent_rich, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_cent_rich, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_cent_rich, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_cent_rich, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMStagVecGetArray(dmGrid_cent_rich, vecOutLocal, &arrOut);
  
    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
               
                if (ex != 0) {
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey][ex - 1][iux_element];
                    next = arrVec[ez][ey][ex][iux_element];
                    inter = (next - prev) / hx;
                    arrOut[ez][ey][ex][iux_left] = inter;      
                }
                if(ex == 0) {
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iux_element];
                    second = arrVec[ez][ey][ex + 1][iux_element];
                    third = arrVec[ez][ey][ex + 2][iux_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hx);
                    arrOut[ez][ey][ex][iux_left] = inter;
                }        
                if(ex == N[0] - 1){
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iux_element];
                    second = arrVec[ez][ey][ex - 1][iux_element];
                    third = arrVec[ez][ey][ex - 2][iux_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hx);
                    arrOut[ez][ey][ex][iux_right] = -inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid_cent_rich, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_cent_rich, vecOutLocal, INSERT_VALUES, P_x_shifted);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const poisson_problem::derive_y_P(Vec & P_y_shifted, Vec const & vec)
{
    PetscInt iuy_up, iuy_down, iuy_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_cent_rich, &N[0], &N[1], &N[2]);
    PetscReal const hy = D_y/ N[1];
    DMGetCoordinateDM(dmGrid_cent_rich, &dmCoord);

    DMGetCoordinates(dmGrid_cent_rich, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_cent_rich, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid_cent_rich, UP, 0, &iuy_up);
    DMStagGetLocationSlot(dmGrid_cent_rich, ELEMENT, 0, &iuy_element);

    DMCreateLocalVector(dmGrid_cent_rich, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_cent_rich, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_cent_rich, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_cent_rich, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMStagVecGetArray(dmGrid_cent_rich, vecOutLocal, &arrOut); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if(ey != 0){
                    PetscReal inter, prev, next;
                    prev = arrVec[ez][ey - 1][ex][iuy_element];
                    next = arrVec[ez][ey][ex][iuy_element];
                    inter = (next - prev) / hy;
                    arrOut[ez][ey][ex][iuy_down] = inter;
                }
                if(ey == 0) {
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iuy_element];
                    second = arrVec[ez][ey + 1][ex][iuy_element];
                    third = arrVec[ez][ey + 2][ex][iuy_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hy);
                    arrOut[ez][ey][ex][iuy_down] = inter;
                }
                if(ey == N[1] - 1){
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iuy_element];
                    second = arrVec[ez][ey - 1][ex][iuy_element];
                    third = arrVec[ez][ey - 2][ex][iuy_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hy);
                    arrOut[ez][ey][ex][iuy_up] = -inter;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid_cent_rich, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_cent_rich, vecOutLocal, INSERT_VALUES, P_y_shifted);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const poisson_problem::derive_z_P(Vec & P_z_shifted, Vec const & vec)
{
    PetscInt iuz_back, iuz_front, iuz_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_cent_rich, &N[0], &N[1], &N[2]);
    PetscReal const hz = D_z/ N[2];
    DMGetCoordinateDM(dmGrid_cent_rich, &dmCoord);

    DMGetCoordinates(dmGrid_cent_rich, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_cent_rich, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid_cent_rich, FRONT, 0, &iuz_front); 
    DMStagGetLocationSlot(dmGrid_cent_rich, ELEMENT, 0, &iuz_element);
    
    DMCreateLocalVector(dmGrid_cent_rich, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_cent_rich, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_cent_rich, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_cent_rich, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMStagVecGetArray(dmGrid_cent_rich, vecOutLocal, &arrOut); 

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                if (ez != 0) {
                    PetscReal inter, prev, next;
                    prev = arrVec[ez - 1][ey][ex][iuz_element];
                    next = arrVec[ez][ey][ex][iuz_element];
                    inter = (next - prev) / hz;
                    arrOut[ez][ey][ex][iuz_back] = inter;
                }
                if(ez == 0) {
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iuz_element];
                    second = arrVec[ez + 1][ey][ex][iuz_element];
                    third = arrVec[ez + 2][ey][ex][iuz_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hz);
                    arrOut[ez][ey][ex][iuz_back] = inter;
                }
                if(ez == N[2] - 1){
                    PetscReal first, second, third, inter;
                    first = arrVec[ez][ey][ex][iuz_element];
                    second = arrVec[ez - 1][ey][ex][iuz_element];
                    third = arrVec[ez - 2][ey][ex][iuz_element];
                    inter = (-2.0 * first + 3.0 * second - third) / (hz);
                    arrOut[ez][ey][ex][iuz_front] = -inter;
                }
                
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid_cent_rich, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_cent_rich, vecOutLocal, INSERT_VALUES, P_z_shifted);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0); 
}

PetscErrorCode const poisson_problem::manage_pressure_x(std::optional<std::reference_wrapper<Vec>> P_opt, std::optional<std::reference_wrapper<Vec>> P_x_opt)
{
    PetscFunctionBegin;

    Vec &P = P_opt ? P_opt.value().get() : this->P;
    Vec &P_x = P_x_opt ? P_x_opt.value().get() : this->P_x;

    Vec P_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &P_shifted);
    DMStagMigrateVec(dmGrid_centered, P, dmGrid_cent_rich, P_shifted);
    
    Vec P_x_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &P_x_shifted);
    derive_x_P(P_x_shifted, P_shifted);
    DMStagMigrateVec(dmGrid_cent_rich, P_x_shifted, dmGrid_staggered_x, P_x);

    VecDestroy(&P_x_shifted);
    VecDestroy(&P_shifted);
    /*DMDestroy(&dmGrid_centered);
    DMDestroy(&dmGrid_cent_rich);*/

    PetscFunctionReturn(0); 
}

PetscErrorCode const poisson_problem::manage_pressure_y(std::optional<std::reference_wrapper<Vec>> P_opt, std::optional<std::reference_wrapper<Vec>> P_y_opt)
{
    PetscFunctionBegin;
    Vec &P = P_opt ? P_opt.value().get() : this->P;
    Vec &P_y = P_y_opt ? P_y_opt.value().get() : this->P_y;

    Vec P_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &P_shifted);
    DMStagMigrateVec(dmGrid_centered, P, dmGrid_cent_rich, P_shifted);

    Vec P_y_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &P_y_shifted);
    derive_y_P(P_y_shifted, P_shifted);
    DMStagMigrateVec(dmGrid_cent_rich, P_y_shifted, dmGrid_staggered_y, P_y);

    VecDestroy(&P_y_shifted);
    VecDestroy(&P_shifted);

    PetscFunctionReturn(0); 
}

PetscErrorCode const poisson_problem::manage_pressure_z(std::optional<std::reference_wrapper<Vec>> P_opt, std::optional<std::reference_wrapper<Vec>> P_z_opt)
{
    PetscFunctionBegin;

    Vec &P = P_opt ? P_opt.value().get() : this->P;
    Vec &P_z = P_opt ? P_z_opt.value().get() : this->P_z;

    Vec P_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &P_shifted);
    DMStagMigrateVec(dmGrid_centered, P, dmGrid_cent_rich, P_shifted);

    Vec P_z_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &P_z_shifted);
    derive_z_P(P_z_shifted, P_shifted);
    DMStagMigrateVec(dmGrid_cent_rich, P_z_shifted, dmGrid_staggered_z, P_z);

    VecDestroy(&P_z_shifted);
    VecDestroy(&P_shifted);

    PetscFunctionReturn(0); 
}

PetscErrorCode const poisson_problem::manage_pressure(std::optional<std::reference_wrapper<Vec>> U_opt,
std::optional<std::reference_wrapper<Vec>> V_up_opt,
std::optional<std::reference_wrapper<Vec>> W_up_opt,  
std::optional<std::reference_wrapper<Vec>> P_opt)
{
    KSP ksp;
    PC  pc;

    PetscFunctionBegin;

    Vec &U_up = U_opt ? U_opt.value().get() : this->U_up;
    Vec &V_up = V_up_opt ? V_up_opt.value().get() : this->V_up;
    Vec &W_up = W_up_opt ? W_up_opt.value().get() : this->W_up;
    Vec &P = P_opt ? P_opt.value().get() : this->P;  

    Vec div;
    DMCreateGlobalVector(dmGrid_centered, &div);
    compute_divergence(div, U_up, V_up, W_up);     

    /*PetscReal mean;
    PetscInt size;
    VecSum(div, &mean);
    VecGetSize(div, &size);
    mean = mean / size;
    VecShift(div, -mean);*/
    AttachNullspace(dmGrid_centered, A);
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCGAMG);
    //PCHYPRESetType(pc, "euclid");
    KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, div, P);

    if(monitor_convergence) {
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        PetscInt iterations;
        KSPGetIterationNumber(ksp, &iterations);
        PetscReal residual_norm;
        KSPGetResidualNorm(ksp, &residual_norm);

        if (reason < 0) {
            PetscPrintf(PETSC_COMM_WORLD, "p-field KSP did not converge. Reason: %s\n", KSPConvergedReasons[reason]);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, 
                        "p-field KSP converged in %d iterations with a final residual norm of %g. Reason: %s\n", 
                        iterations, residual_norm, KSPConvergedReasons[reason]);
        }
    }


    VecDestroy(&div);
    KSPDestroy(&ksp);
       
    PetscFunctionReturn(0); 
}

PetscErrorCode const poisson_problem::exodus(size_t i){

        PetscFunctionBegin;

        PetscViewer viewer_p;
        DM DM_p;
        //DMStagCreateCompatibleDMStag(dmGrid_centered, 0, 0, 0, 1, &DM_p);
        Vec p;
        DMStagVecSplitToDMDA(dmGrid_centered, P, ELEMENT, 0, &DM_p, &p);
        PetscObjectSetName((PetscObject)p, "p");
        char filename_p[50];
        sprintf(filename_p, "%sp%03zu.vtr", base_path, i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_p, FILE_MODE_WRITE, &viewer_p);
        VecView(p, viewer_p);
        VecDestroy(&p);
        DMDestroy(&DM_p);
        PetscViewerDestroy(&viewer_p);
        std::cout << "Iteration " << i << " completed." << std::endl; 

        PetscFunctionReturn(0);
}


poisson_problem::~poisson_problem()
{
    VecDestroy(&P);
    VecDestroy(&P_x);
    VecDestroy(&P_y);
    VecDestroy(&P_z);
    VecDestroy(&U_up);
    VecDestroy(&V_up);
    VecDestroy(&W_up);
    MatDestroy(&A);
    DMDestroy(&dmGrid_staggered_x);
    DMDestroy(&dmGrid_staggered_y);
    DMDestroy(&dmGrid_staggered_z);
    DMDestroy(&dmGrid_centered);
    DMDestroy(&dmGrid_cent_rich);
    std::cout << "Poisson Destructor Called" << std::endl;
}

PetscErrorCode poisson_problem::AttachNullspace(DM dmSol, Mat A)
{
    DM           dmPressure;
    Vec          constantPressure, basis;
    PetscReal    nrm;
    MatNullSpace matNullSpace;

    PetscFunctionBeginUser;
    DMStagCreateCompatibleDMStag(dmSol, 0, 0, 1, 0, &dmPressure);
    DMGetGlobalVector(dmPressure, &constantPressure);
    VecSet(constantPressure, 1.0);
    VecNorm(constantPressure, NORM_2, &nrm);
    VecScale(constantPressure, 1.0 / nrm);
    DMCreateGlobalVector(dmSol, &basis);
    DMStagMigrateVec(dmPressure, constantPressure, dmSol, basis);
    MatNullSpaceCreate(PetscObjectComm((PetscObject)dmSol), PETSC_FALSE, 1, &basis, &matNullSpace);
    VecDestroy(&basis);
    VecDestroy(&constantPressure);
    MatSetNullSpace(A, matNullSpace);
    MatNullSpaceDestroy(&matNullSpace);
    PetscFunctionReturn(0);
}