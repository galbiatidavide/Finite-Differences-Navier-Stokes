#include "navier_stokes.hpp"

PetscErrorCode const navier_stokes_problem::assemble_lhs() 
{
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_centered, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];

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
                            //valA[0] = -1.0 / (hx * hx) + -1.0 / (hy * hy) - 1.0 / (hz * hz);
                            valA[0] = 1.0;
                            /*col[1].i = ex;
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
                            valA[3] = 1.0 / (hz * hz);*/
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

PetscErrorCode const navier_stokes_problem::assemble_divergence(Vec & div, Vec const & U, Vec const &  V, Vec const & W) 
{
    PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_cent_rich, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];
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

PetscErrorCode const navier_stokes_problem::compute_divergence(Vec & div, Vec const & U_n, Vec const & V_n, Vec const & W_n) 
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

PetscErrorCode const navier_stokes_problem::derive_x_P(Vec & P_x_shifted, Vec const & vec)
{
    PetscInt iux_left, iux_right, iux_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_cent_rich, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
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

PetscErrorCode const navier_stokes_problem::derive_y_P(Vec & P_y_shifted, Vec const & vec)
{
    PetscInt iuy_up, iuy_down, iuy_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_cent_rich, &N[0], &N[1], &N[2]);
    PetscReal const hy = 1.0 / N[1];
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

PetscErrorCode const navier_stokes_problem::derive_z_P(Vec & P_z_shifted, Vec const & vec)
{
    PetscInt iuz_back, iuz_front, iuz_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_cent_rich, &N[0], &N[1], &N[2]);
    PetscReal const hz = 1.0 / N[2];
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

PetscErrorCode const navier_stokes_problem::manage_pressure_x()
{
    PetscFunctionBegin;

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

PetscErrorCode const navier_stokes_problem::manage_pressure_y()
{
    PetscFunctionBegin

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

PetscErrorCode const navier_stokes_problem::manage_pressure_z()
{
    PetscFunctionBegin

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

PetscErrorCode const navier_stokes_problem::manage_pressure()
{
    KSP ksp;
    PC  pc;

    PetscFunctionBegin;

    //Stiamo usando direttamente U_pre, V_pre, W_pre come reference. Dovrebbe essere sicuro ma non mi fiderei. Pare funzionare comunque.
    /*Vec U_n, V_n, W_n;
    DMCreateGlobalVector(dmGrid_Staggered, &U_n);
    DMCreateGlobalVector(dmGrid_Staggered, &V_n);
    DMCreateGlobalVector(dmGrid_Staggered, &W_n);
    VecCopy(U_pre, U_n);
    VecCopy(V_pre, V_n);
    VecCopy(W_pre, W_n);*/

    Vec div;
    DMCreateGlobalVector(dmGrid_centered, &div);
    compute_divergence(div, U_up, V_up, W_up);     

    /*Vec force;
    DMCreateGlobalVector(dmGrid_centered, &force);
    CreateReferenceSolutionTryForce(dmGrid_centered, force, 0);
    CheckSolution(div, force);*/

    //assemble_lhs();
    /*PetscReal mean;
    PetscInt size;
    VecSum(div, &mean);
    VecGetSize(div, &size);
    mean = mean / size;
    VecShift(div, -mean);*/

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCBJACOBI);
    //PCHYPRESetType(pc, "euclid");
    KSPSetTolerances(ksp, PETSC_DEFAULT, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, div, P);

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


    VecDestroy(&div);
    KSPDestroy(&ksp);
       
    PetscFunctionReturn(0); 
}

PetscErrorCode const navier_stokes_problem::update_bc_U(PetscReal const & theta)
{
    PetscInt icux_left[3], icux_right[3], iux_left, iux_right;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec; 

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_staggered_x, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_staggered_x, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid_staggered_x, &dmCoord);

    DMGetCoordinates(dmGrid_staggered_x, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_staggered_x, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid_staggered_x, RIGHT, 0, &iux_right);

    DMCreateLocalVector(dmGrid_staggered_x, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_x, U_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_x, U_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArray(dmGrid_staggered_x, vecLocal, &arrVec);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid_staggered_x, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_x, mask_U, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_x, mask_U, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_x, penalizationLocal, &arrPenalization);


    for (ez = startz; ez < startz + nz; ++ez) { 
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ex == N[0] - 1) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iux_right] == 0.0) {
                        val = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iux_right] = val;                    
                } else if(ex == 0) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iux_left] == 0.0) {
                        val = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]], arrCoord[ez][ey][ex][icux_left[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iux_left] = val;
                }

            }
        }
    }


    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid_staggered_x, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid_staggered_x, vecLocal, INSERT_VALUES, U_up);
    DMRestoreLocalVector(dmGrid_staggered_x, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArrayRead(dmGrid_staggered_x, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid_staggered_x, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const navier_stokes_problem::update_bc_V(PetscReal const & theta) 
{

    PetscInt icuy_down[3], icuy_up[3], iuy_down, iuy_up;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_staggered_y, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_staggered_y, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid_staggered_y, &dmCoord);

    DMGetCoordinates(dmGrid_staggered_y, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy_down[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
    } 
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_staggered_y, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid_staggered_y, UP, 0, &iuy_up);

    DMCreateLocalVector(dmGrid_staggered_y, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_y, V_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_y, V_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_y, vecLocal, &arrVec);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid_staggered_y, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_y, mask_V, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_y, mask_V, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_y, penalizationLocal, &arrPenalization);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ey == N[1] - 1) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iuy_up] == 0.0) {
                        val = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iuy_up] = val;
                } else if(ey == 0) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iuy_down] == 0.0) {
                        val = uyRef(arrCoord[ez][ey][ex][icuy_down[0]], arrCoord[ez][ey][ex][icuy_down[1]], arrCoord[ez][ey][ex][icuy_down[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iuy_down] = val;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid_staggered_y, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid_staggered_y, vecLocal, INSERT_VALUES, V_up);
    DMRestoreLocalVector(dmGrid_staggered_y, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArrayRead(dmGrid_staggered_y, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid_staggered_y, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const navier_stokes_problem::update_bc_W(PetscReal const & theta) 
{

    PetscInt icuz_back[3], icuz_front[3], iuz_back, iuz_front;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_staggered_z, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_staggered_z, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid_staggered_z, &dmCoord);

    DMGetCoordinates(dmGrid_staggered_z, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz_back[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_staggered_z, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid_staggered_z, FRONT, 0, &iuz_front); 
    
    DMCreateLocalVector(dmGrid_staggered_z, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_z, W_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_z, W_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_z, vecLocal, &arrVec);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid_staggered_z, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid_staggered_z, mask_W, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid_staggered_z, mask_W, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid_staggered_z, penalizationLocal, &arrPenalization);

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {

                if (ez == N[2] - 1) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iuz_front] == 0.0) {
                        val = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iuz_front] = val;
                } else if(ez == 0) {
                    PetscReal val;
                    if(arrPenalization[ez][ey][ex][iuz_back] == 0.0) {
                        val = uzRef(arrCoord[ez][ey][ex][icuz_back[0]], arrCoord[ez][ey][ex][icuz_back[1]], arrCoord[ez][ey][ex][icuz_back[2]], theta);
                    } else {
                        val = 0.0;
                    }
                    arrVec[ez][ey][ex][iuz_back] = val;
                }
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid_staggered_z, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid_staggered_z, vecLocal, INSERT_VALUES, W_up);
    DMRestoreLocalVector(dmGrid_staggered_z, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArrayRead(dmGrid_staggered_z, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid_staggered_z, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const navier_stokes_problem::update_velocity(PetscReal const & theta)
{
    PetscFunctionBegin;

    VecAXPY(U_up, -dt, P_x);
    VecAXPY(V_up, -dt, P_y);
    VecAXPY(W_up, -dt, P_z);

    /*VecCopy(U_pre, U_up);
    VecCopy(V_pre, V_up);
    VecCopy(W_pre, W_up);*/

    update_bc_U(theta);
    update_bc_V(theta);
    update_bc_W(theta);


    PetscFunctionReturn(0);
}

PetscErrorCode const navier_stokes_problem::assemble_magnitude(Vec & Magnitude_Shifted, Vec const & U, Vec const & V, Vec const & W) 
{
    PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
    PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_cent_rich, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
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
                inter = sqrt(((right + left)*(right + left)) / 4 + ((up + down)*(up + down)) / 4 + ((front + back)*(front + back)) / 4);
                arrOut[ez][ey][ex][iu_element] = inter;
            }
        }
    }

    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecULocal, &arrU);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecVLocal, &arrV);
    DMStagVecRestoreArrayRead(dmGrid_cent_rich, vecWLocal, &arrW);
    DMStagVecRestoreArray(dmGrid_cent_rich, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_cent_rich, vecOutLocal, INSERT_VALUES, Magnitude_Shifted);    
    DMRestoreLocalVector(dmGrid_cent_rich, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecULocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecVLocal);
    DMRestoreLocalVector(dmGrid_cent_rich, &vecWLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);    
}

PetscErrorCode const navier_stokes_problem::compute_magnitude()
{    
    PetscFunctionBegin;
    Vec U_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &U_shifted);
    DMStagMigrateVec(dmGrid_staggered_x, U_up, dmGrid_cent_rich, U_shifted);
    Vec V_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &V_shifted);
    DMStagMigrateVec(dmGrid_staggered_y, V_up, dmGrid_cent_rich, V_shifted);
    Vec W_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &W_shifted);
    DMStagMigrateVec(dmGrid_staggered_z, W_up, dmGrid_cent_rich, W_shifted);

    Vec magnitude_shifted;
    DMCreateGlobalVector(dmGrid_cent_rich, &magnitude_shifted);
    assemble_magnitude(magnitude_shifted, U_shifted, V_shifted, W_shifted);
    DMStagMigrateVec(dmGrid_cent_rich, magnitude_shifted, dmGrid_centered, Magnitude);

    VecDestroy(&magnitude_shifted);
    VecDestroy(&U_shifted);  
    VecDestroy(&V_shifted);
    VecDestroy(&W_shifted);            
    PetscFunctionReturn(0); 
}

PetscErrorCode navier_stokes_problem::exodus(size_t i){

        PetscFunctionBegin;
        PetscViewer viewer_magnitude;
        DM DM_magnitude;
        //DMStagCreateCompatibleDMStag(dmGrid_staggered_x, 0, 0, 1, 0, &DM_u);
        Vec magnitude;
        DMStagVecSplitToDMDA(dmGrid_centered, Magnitude, ELEMENT, 0, &DM_magnitude, &magnitude);
        PetscObjectSetName((PetscObject)magnitude, "magnitude");
        char filename_magnitude[50]; 
        sprintf(filename_magnitude, "%smagnitude%03zu.vtr", base_path, i);

        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_centered), filename_magnitude, FILE_MODE_WRITE, &viewer_magnitude);
        VecView(magnitude, viewer_magnitude);
        VecDestroy(&magnitude);
        DMDestroy(&DM_magnitude);
        PetscViewerDestroy(&viewer_magnitude);

        PetscViewer viewer_u;
        DM DM_u;
        //DMStagCreateCompatibleDMStag(dmGrid_staggered_x, 0, 0, 1, 0, &DM_u);
        Vec u;
        DMStagVecSplitToDMDA(dmGrid_staggered_x, U_up, LEFT, 0, &DM_u, &u);
        PetscObjectSetName((PetscObject)u, "x_component");
        char filename_u[50]; 
        sprintf(filename_u, "%sx_component%03zu.vtr", base_path, i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
        VecView(u, viewer_u);
        VecDestroy(&u);
        DMDestroy(&DM_u);
        PetscViewerDestroy(&viewer_u); 

        PetscViewer viewer_v;
        DM DM_v;
        //DMStagCreateCompatibleDMStag(dmGrid_staggered_y, 0, 0, 1, 0, &DM_v);
        Vec v;
        DMStagVecSplitToDMDA(dmGrid_staggered_y, V_up, DOWN, 0, &DM_v, &v);
        PetscObjectSetName((PetscObject)v, "y_component");
        char filename_v[50];
        sprintf(filename_v, "%sy_component%03zu.vtr", base_path, i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_staggered_y), filename_v, FILE_MODE_WRITE, &viewer_v);
        VecView(v, viewer_v);
        VecDestroy(&v);
        DMDestroy(&DM_v);
        PetscViewerDestroy(&viewer_v);

        PetscViewer viewer_w;
        DM DM_w;
        //DMStagCreateCompatibleDMStag(dmGrid_staggered_z, 0, 0, 1, 0, &DM_w);
        Vec w;
        DMStagVecSplitToDMDA(dmGrid_staggered_z, W_up, BACK, 0, &DM_w, &w);
        PetscObjectSetName((PetscObject)w, "z_component");
        char filename_w[50];
        sprintf(filename_w, "%sz_component%03zu.vtr", base_path, i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_staggered_z), filename_w, FILE_MODE_WRITE, &viewer_w);
        VecView(w, viewer_w);
        VecDestroy(&w);
        DMDestroy(&DM_w);
        PetscViewerDestroy(&viewer_w);

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

PetscErrorCode const navier_stokes_problem::solve()
{
    PetscFunctionBegin;

    PrintSimulationParameters();

    transport_problem_x transport_x(dmGrid_shift_transp, dmGrid_stag_transp, dmGrid_cent_rich);
    transport_problem_y transport_y(dmGrid_shift_transp, dmGrid_stag_transp, dmGrid_cent_rich);
    transport_problem_z transport_z(dmGrid_shift_transp, dmGrid_stag_transp, dmGrid_cent_rich);

    this->assemble_lhs();

    parabolic_problem_x parabolic_x(dmGrid_staggered_x);
    parabolic_problem_y parabolic_y(dmGrid_staggered_y); 
    parabolic_problem_z parabolic_z(dmGrid_staggered_z);

    parabolic_x.assemble_lhs();
    parabolic_y.assemble_lhs();
    parabolic_z.assemble_lhs();  
    
    for(size_t i = 0; i < iter; ++i){

        theta = i * dt;

        transport_x.solve_step_x(theta, U_up, V_up, W_up);
        transport_y.solve_step_y(theta, U_up, V_up, W_up);
        transport_z.solve_step_z(theta, U_up, V_up, W_up);

        parabolic_x.solve_step(theta, U_up);
        parabolic_y.solve_step(theta, V_up);
        parabolic_z.solve_step(theta, W_up);

        this->manage_pressure();
        /*Vec pressure_bench;
        DMCreateGlobalVector(dmGrid_centered, &pressure_bench);
        CreateReferencePressure(dmGrid_centered, pressure_bench, theta);
        CheckSolution(P, pressure_bench);
        VecDestroy(&pressure_bench);*/
        this->manage_pressure_x();
        this->manage_pressure_y();
        this->manage_pressure_z();
        theta = (i+1)*dt;
        this->update_velocity(theta);
        Vec solution_x;
        DMCreateGlobalVector(dmGrid_staggered_x, &solution_x);
        CreateAnalyticalU(dmGrid_staggered_x, solution_x, theta);
        CheckSolution(U_up, solution_x, "U");
        VecDestroy(&solution_x);
        Vec solution_y;
        DMCreateGlobalVector(dmGrid_staggered_y, &solution_y);
        CreateAnalyticalV(dmGrid_staggered_y, solution_y, theta);
        CheckSolution(V_up, solution_y, "V");
        VecDestroy(&solution_y);
        Vec solution_z;
        DMCreateGlobalVector(dmGrid_staggered_z, &solution_z);
        CreateAnalyticalW(dmGrid_staggered_z, solution_z, theta);
        CheckSolution(W_up, solution_z, "W");
        VecDestroy(&solution_z);

        compute_magnitude();

        exodus(i);



    

    }




    PetscFunctionReturn(0);
}
