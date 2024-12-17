#include "parabolic.hpp"
#include "Setting.hpp"
#include "Transport.hpp"

#ifndef NAVIER_STOKES_PROBLEM_HPP
#define NAVIER_STOKES_PROBLEM_HPP

class navier_stokes_problem: public parabolic_problem_x, public parabolic_problem_y, public parabolic_problem_z
{
private:
DM const & dmGrid_Staggered_x;
DM const & dmGrid_Staggered_y;
DM const & dmGrid_Staggered_z;
DM const & dmGrid_Centered;
DM const & dmGrid_Shifted;

Vec P;
Vec P_x;
Vec P_y;
Vec P_z;
Vec Magnitude;

Mat A;

PetscReal const & dt;
PetscReal const & Re;
PetscInt const & iter;

Vec U_prova;
Vec V_prova;
Vec W_prova;

Vec penalization_u;
Vec penalization_v;
Vec penalization_w;

ProblemSetting<Transport> setting_transport;


PetscErrorCode const assemble_lhs() 
{
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_Centered, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_Centered, &N[0], &N[1], &N[2]);
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
                    DMStagMatSetValuesStencil(dmGrid_Centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
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
                        col[4].loc = LEFT;
                        col[4].c = 0;
                        valA[4] = 1.0 / (hz * hz);
                        col[5].i = ex;
                        col[5].j = ey;
                        col[5].k = ez + 1;
                        col[5].loc = ELEMENT;
                        col[5].c = 0;
                        valA[5] = 1.0 / (hz * hz);
                    }
                    DMStagMatSetValuesStencil(dmGrid_Centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
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
                    DMStagMatSetValuesStencil(dmGrid_Centered, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                }                
            }
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    PetscFunctionReturn(0); 
}

PetscErrorCode const assemble_divergence(Vec & div, Vec const & U, Vec const &  V, Vec const & W) 
{
    PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    PetscReal const hy = 1.0 / N[1];
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);

    DMGetCoordinates(dmGrid_Shifted, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &iu_left);
    DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &iu_right);
    DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &iu_up);
    DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &iu_down);
    DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &iu_front);
    DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &iu_back);
    DMStagGetLocationSlot(dmGrid_Shifted, ELEMENT, 0, &iu_element);

    DMCreateLocalVector(dmGrid_Shifted, &vecULocal);
    DMGlobalToLocalBegin(dmGrid_Shifted, U, INSERT_VALUES, vecULocal);
    DMGlobalToLocalEnd(dmGrid_Shifted, U, INSERT_VALUES, vecULocal);
    DMStagVecGetArrayRead(dmGrid_Shifted, vecULocal, &arrU);

    DMCreateLocalVector(dmGrid_Shifted, &vecVLocal);
    DMGlobalToLocalBegin(dmGrid_Shifted, V, INSERT_VALUES, vecVLocal);
    DMGlobalToLocalEnd(dmGrid_Shifted, V, INSERT_VALUES, vecVLocal);
    DMStagVecGetArrayRead(dmGrid_Shifted, vecVLocal, &arrV);

    DMCreateLocalVector(dmGrid_Shifted, &vecWLocal);
    DMGlobalToLocalBegin(dmGrid_Shifted, W, INSERT_VALUES, vecWLocal);
    DMGlobalToLocalEnd(dmGrid_Shifted, W, INSERT_VALUES, vecWLocal);
    DMStagVecGetArrayRead(dmGrid_Shifted, vecWLocal, &arrW);    

    DMGetLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMStagVecGetArray(dmGrid_Shifted, vecOutLocal, &arrOut);

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
    DMStagVecRestoreArrayRead(dmGrid_Shifted, vecULocal, &arrU);
    DMStagVecRestoreArrayRead(dmGrid_Shifted, vecVLocal, &arrV);
    DMStagVecRestoreArrayRead(dmGrid_Shifted, vecWLocal, &arrW);
    DMStagVecRestoreArray(dmGrid_Shifted, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_Shifted, vecOutLocal, INSERT_VALUES, div);    
    DMRestoreLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_Shifted, &vecULocal);
    DMRestoreLocalVector(dmGrid_Shifted, &vecVLocal);
    DMRestoreLocalVector(dmGrid_Shifted, &vecWLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const compute_divergence(Vec & div, Vec const & U_n, Vec const & V_n, Vec const & W_n) 
{
        PetscFunctionBegin;

        Vec U_shifted, V_shifted, W_shifted;
        DMCreateGlobalVector(dmGrid_Shifted, &U_shifted);
        DMCreateGlobalVector(dmGrid_Shifted, &V_shifted);
        DMCreateGlobalVector(dmGrid_Shifted, &W_shifted);
        DMStagMigrateVec(dmGrid_Staggered_x, U_n, dmGrid_Shifted, U_shifted);
        DMStagMigrateVec(dmGrid_Staggered_y, V_n, dmGrid_Shifted, V_shifted);
        DMStagMigrateVec(dmGrid_Staggered_z, W_n, dmGrid_Shifted, W_shifted);

        Vec div_shifted;
        DMCreateGlobalVector(dmGrid_Shifted, &div_shifted);
        assemble_divergence(div_shifted, U_shifted, V_shifted, W_shifted);
        DMStagMigrateVec(dmGrid_Shifted, div_shifted, dmGrid_Centered, div);

        VecDestroy(&U_shifted);
        VecDestroy(&V_shifted);
        VecDestroy(&W_shifted);
        VecDestroy(&div_shifted);

        PetscFunctionReturn(0); 
}

PetscErrorCode const derive_x_P(Vec & P_x_shifted, Vec const & vec)
{
    PetscInt iux_left, iux_right, iux_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
    PetscReal const hx = 1.0 / N[0];
    DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);

    DMGetCoordinates(dmGrid_Shifted, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid_Shifted, ELEMENT, 0, &iux_element);

    DMCreateLocalVector(dmGrid_Shifted, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_Shifted, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_Shifted, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_Shifted, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMStagVecGetArray(dmGrid_Shifted, vecOutLocal, &arrOut);
  
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
    DMStagVecRestoreArrayRead(dmGrid_Shifted, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid_Shifted, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_Shifted, vecOutLocal, INSERT_VALUES, P_x_shifted);
    DMRestoreLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_Shifted, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const derive_y_P(Vec & P_y_shifted, Vec const & vec)
{
    PetscInt iuy_up, iuy_down, iuy_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
    PetscReal const hy = 1.0 / N[1];
    DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);

    DMGetCoordinates(dmGrid_Shifted, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &iuy_up);
    DMStagGetLocationSlot(dmGrid_Shifted, ELEMENT, 0, &iuy_element);

    DMCreateLocalVector(dmGrid_Shifted, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_Shifted, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_Shifted, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_Shifted, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMStagVecGetArray(dmGrid_Shifted, vecOutLocal, &arrOut); 

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
    DMStagVecRestoreArrayRead(dmGrid_Shifted, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid_Shifted, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_Shifted, vecOutLocal, INSERT_VALUES, P_y_shifted);
    DMRestoreLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_Shifted, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const derive_z_P(Vec & P_z_shifted, Vec const & vec)
{
    PetscInt iuz_back, iuz_front, iuz_element;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_Shifted, &N[0], &N[1], &N[2]);
    PetscReal const hz = 1.0 / N[2];
    DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);

    DMGetCoordinates(dmGrid_Shifted, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &iuz_front); 
    DMStagGetLocationSlot(dmGrid_Shifted, ELEMENT, 0, &iuz_element);
    
    DMCreateLocalVector(dmGrid_Shifted, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_Shifted, vec, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_Shifted, vec, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_Shifted, vecLocal, &arrVec);

    DMGetLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMStagVecGetArray(dmGrid_Shifted, vecOutLocal, &arrOut); 

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
    DMStagVecRestoreArrayRead(dmGrid_Shifted, vecLocal, &arrVec);
    DMStagVecRestoreArray(dmGrid_Shifted, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_Shifted, vecOutLocal, INSERT_VALUES, P_z_shifted);
    DMRestoreLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_Shifted, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0); 
}

PetscErrorCode const manage_pressure_x()
{
    PetscFunctionBegin;

    Vec P_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &P_shifted);
    DMStagMigrateVec(dmGrid_Centered, P, dmGrid_Shifted, P_shifted);
    
    Vec P_x_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &P_x_shifted);
    derive_x_P(P_x_shifted, P_shifted);
    DMStagMigrateVec(dmGrid_Shifted, P_x_shifted, dmGrid_Staggered_x, P_x);

    VecDestroy(&P_x_shifted);
    VecDestroy(&P_shifted);
    /*DMDestroy(&dmGrid_centered);
    DMDestroy(&dmGrid_shifted);*/

    PetscFunctionReturn(0); 
}

PetscErrorCode const manage_pressure_y()
{
    PetscFunctionBegin

    Vec P_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &P_shifted);
    DMStagMigrateVec(dmGrid_Centered, P, dmGrid_Shifted, P_shifted);

    Vec P_y_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &P_y_shifted);
    derive_y_P(P_y_shifted, P_shifted);
    DMStagMigrateVec(dmGrid_Shifted, P_y_shifted, dmGrid_Staggered_y, P_y);

    VecDestroy(&P_y_shifted);
    VecDestroy(&P_shifted);

    PetscFunctionReturn(0); 
}

PetscErrorCode const manage_pressure_z()
{
    PetscFunctionBegin

    Vec P_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &P_shifted);
    DMStagMigrateVec(dmGrid_Centered, P, dmGrid_Shifted, P_shifted);

    Vec P_z_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &P_z_shifted);
    derive_z_P(P_z_shifted, P_shifted);
    DMStagMigrateVec(dmGrid_Shifted, P_z_shifted, dmGrid_Staggered_z, P_z);

    VecDestroy(&P_z_shifted);
    VecDestroy(&P_shifted);

    PetscFunctionReturn(0); 
}

PetscErrorCode const manage_pressure()
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
    DMCreateGlobalVector(dmGrid_Centered, &div);
    compute_divergence(div, U_up, V_up, W_up);     

    /*Vec force;
    DMCreateGlobalVector(dmGrid_centered, &force);
    CreateReferenceSolutionTryForce(dmGrid_centered, force, 0);
    CheckSolution(div, force);*/

    //assemble_lhs();
    PetscReal mean;
    PetscInt size;
    VecSum(div, &mean);
    VecGetSize(div, &size);
    mean = mean / size;
    VecShift(div, -mean);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);
    KSPSetOperators(ksp, A, A);
    KSPGetPC(ksp, &pc);
    //PCSetType(pc, PCNONE);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, div, P);

    VecDestroy(&div);
    KSPDestroy(&ksp);
       
    PetscFunctionReturn(0); 
}

PetscErrorCode const update_bc_U(PetscReal const & theta)
{
    PetscInt icux_left[3], icux_right[3], iux_left, iux_right;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec; 

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_Staggered_x, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_Staggered_x, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid_Staggered_x, &dmCoord);

    DMGetCoordinates(dmGrid_Staggered_x, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_Staggered_x, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid_Staggered_x, RIGHT, 0, &iux_right);

    DMCreateLocalVector(dmGrid_Staggered_x, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_Staggered_x, U_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_Staggered_x, U_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArray(dmGrid_Staggered_x, vecLocal, &arrVec);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid_Staggered_x, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid_Staggered_x, penalization_u, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid_Staggered_x, penalization_u, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid_Staggered_x, penalizationLocal, &arrPenalization);


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
    DMStagVecRestoreArray(dmGrid_Staggered_x, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid_Staggered_x, vecLocal, INSERT_VALUES, U_up);
    DMRestoreLocalVector(dmGrid_Staggered_x, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArrayRead(dmGrid_Staggered_x, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid_Staggered_x, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const update_bc_V(PetscReal const & theta) 
{

    PetscInt icuy_down[3], icuy_up[3], iuy_down, iuy_up;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_Staggered_y, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_Staggered_y, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid_Staggered_y, &dmCoord);

    DMGetCoordinates(dmGrid_Staggered_y, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy_down[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
    } 
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_Staggered_y, DOWN, 0, &iuy_down);
    DMStagGetLocationSlot(dmGrid_Staggered_y, UP, 0, &iuy_up);

    DMCreateLocalVector(dmGrid_Staggered_y, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_Staggered_y, V_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_Staggered_y, V_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_Staggered_y, vecLocal, &arrVec);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid_Staggered_y, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid_Staggered_y, penalization_v, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid_Staggered_y, penalization_v, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid_Staggered_y, penalizationLocal, &arrPenalization);

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
    DMStagVecRestoreArray(dmGrid_Staggered_y, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid_Staggered_y, vecLocal, INSERT_VALUES, V_up);
    DMRestoreLocalVector(dmGrid_Staggered_y, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArrayRead(dmGrid_Staggered_y, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid_Staggered_y, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const update_bc_W(PetscReal const & theta) 
{

    PetscInt icuz_back[3], icuz_front[3], iuz_back, iuz_front;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec vecLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrVec;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_Staggered_z, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid_Staggered_z, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid_Staggered_z, &dmCoord);

    DMGetCoordinates(dmGrid_Staggered_z, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz_back[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
    }  
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_Staggered_z, BACK, 0, &iuz_back);
    DMStagGetLocationSlot(dmGrid_Staggered_z, FRONT, 0, &iuz_front); 
    
    DMCreateLocalVector(dmGrid_Staggered_z, &vecLocal);
    DMGlobalToLocalBegin(dmGrid_Staggered_z, W_up, INSERT_VALUES, vecLocal);
    DMGlobalToLocalEnd(dmGrid_Staggered_z, W_up, INSERT_VALUES, vecLocal);
    DMStagVecGetArrayRead(dmGrid_Staggered_z, vecLocal, &arrVec);

    Vec penalizationLocal;
    PetscReal ****arrPenalization;
    DMCreateLocalVector(dmGrid_Staggered_z, &penalizationLocal);
    DMGlobalToLocalBegin(dmGrid_Staggered_z, penalization_w, INSERT_VALUES, penalizationLocal);
    DMGlobalToLocalEnd(dmGrid_Staggered_z, penalization_w, INSERT_VALUES, penalizationLocal);
    DMStagVecGetArrayRead(dmGrid_Staggered_z, penalizationLocal, &arrPenalization);

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
    DMStagVecRestoreArray(dmGrid_Staggered_z, vecLocal, &arrVec);
    DMLocalToGlobal(dmGrid_Staggered_z, vecLocal, INSERT_VALUES, W_up);
    DMRestoreLocalVector(dmGrid_Staggered_z, &vecLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArrayRead(dmGrid_Staggered_z, penalizationLocal, &arrPenalization);
    DMRestoreLocalVector(dmGrid_Staggered_z, &penalizationLocal);

    PetscFunctionReturn(0);
}

PetscErrorCode const update_velocity(PetscReal const & theta)
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

PetscErrorCode const assemble_magnitude(Vec & Magnitude_Shifted, Vec const & U, Vec const & V, Vec const & W) 
{
    PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
    PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez;
    DM dmCoord;
    Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
    PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid_Shifted, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMGetCoordinateDM(dmGrid_Shifted, &dmCoord);

    DMGetCoordinates(dmGrid_Shifted, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid_Shifted, LEFT, 0, &iu_left);
    DMStagGetLocationSlot(dmGrid_Shifted, RIGHT, 0, &iu_right);
    DMStagGetLocationSlot(dmGrid_Shifted, UP, 0, &iu_up);
    DMStagGetLocationSlot(dmGrid_Shifted, DOWN, 0, &iu_down);
    DMStagGetLocationSlot(dmGrid_Shifted, FRONT, 0, &iu_front);
    DMStagGetLocationSlot(dmGrid_Shifted, BACK, 0, &iu_back);
    DMStagGetLocationSlot(dmGrid_Shifted, ELEMENT, 0, &iu_element);

    DMCreateLocalVector(dmGrid_Shifted, &vecULocal);
    DMGlobalToLocalBegin(dmGrid_Shifted, U, INSERT_VALUES, vecULocal);
    DMGlobalToLocalEnd(dmGrid_Shifted, U, INSERT_VALUES, vecULocal);
    DMStagVecGetArrayRead(dmGrid_Shifted, vecULocal, &arrU);

    DMCreateLocalVector(dmGrid_Shifted, &vecVLocal);
    DMGlobalToLocalBegin(dmGrid_Shifted, V, INSERT_VALUES, vecVLocal);
    DMGlobalToLocalEnd(dmGrid_Shifted, V, INSERT_VALUES, vecVLocal);
    DMStagVecGetArrayRead(dmGrid_Shifted, vecVLocal, &arrV);

    DMCreateLocalVector(dmGrid_Shifted, &vecWLocal);
    DMGlobalToLocalBegin(dmGrid_Shifted, W, INSERT_VALUES, vecWLocal);
    DMGlobalToLocalEnd(dmGrid_Shifted, W, INSERT_VALUES, vecWLocal);
    DMStagVecGetArrayRead(dmGrid_Shifted, vecWLocal, &arrW);    

    DMGetLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMStagVecGetArray(dmGrid_Shifted, vecOutLocal, &arrOut);

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
    DMStagVecRestoreArrayRead(dmGrid_Shifted, vecULocal, &arrU);
    DMStagVecRestoreArrayRead(dmGrid_Shifted, vecVLocal, &arrV);
    DMStagVecRestoreArrayRead(dmGrid_Shifted, vecWLocal, &arrW);
    DMStagVecRestoreArray(dmGrid_Shifted, vecOutLocal, &arrOut);
    DMLocalToGlobal(dmGrid_Shifted, vecOutLocal, INSERT_VALUES, Magnitude_Shifted);    
    DMRestoreLocalVector(dmGrid_Shifted, &vecOutLocal);
    DMRestoreLocalVector(dmGrid_Shifted, &vecULocal);
    DMRestoreLocalVector(dmGrid_Shifted, &vecVLocal);
    DMRestoreLocalVector(dmGrid_Shifted, &vecWLocal);
    DMRestoreLocalVector(dmCoord, &coordLocal);

    PetscFunctionReturn(0);    
}

PetscErrorCode const compute_magnitude()
{    
    PetscFunctionBegin;
    Vec U_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &U_shifted);
    DMStagMigrateVec(dmGrid_Staggered_x, U_up, dmGrid_Shifted, U_shifted);
    Vec V_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &V_shifted);
    DMStagMigrateVec(dmGrid_Staggered_y, V_up, dmGrid_Shifted, V_shifted);
    Vec W_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &W_shifted);
    DMStagMigrateVec(dmGrid_Staggered_z, W_up, dmGrid_Shifted, W_shifted);

    Vec magnitude_shifted;
    DMCreateGlobalVector(dmGrid_Shifted, &magnitude_shifted);
    assemble_magnitude(magnitude_shifted, U_shifted, V_shifted, W_shifted);
    DMStagMigrateVec(dmGrid_Shifted, magnitude_shifted, dmGrid_Centered, Magnitude);

    VecDestroy(&magnitude_shifted);
    VecDestroy(&U_shifted);  
    VecDestroy(&V_shifted);
    VecDestroy(&W_shifted);            
    PetscFunctionReturn(0); 
}

public:

navier_stokes_problem(ProblemSetting<Transport> const & setting_transport, Vec const & U_0_, Vec const & V_0_, Vec const & W_0_, Vec const & penalization_u, Vec const & penalization_v, Vec const & penalization_w, DM const & dmGrid_Staggered_x, DM const & dmGrid_Staggered_y, DM const & dmGrid_Staggered_z, DM const & dmGrid_Centered, DM const & dmGrid_Shifted, PetscReal const & dt, PetscInt const & iter, PetscReal const & Re) :
parabolic_problem_x(U_0_, penalization_u, dmGrid_Staggered_x, Re, dt, iter),
parabolic_problem_y(V_0_, penalization_v, dmGrid_Staggered_y, Re, dt, iter),
parabolic_problem_z(W_0_, penalization_w, dmGrid_Staggered_z, Re, dt, iter),
dmGrid_Staggered_x(dmGrid_Staggered_x),
dmGrid_Staggered_y(dmGrid_Staggered_y),
dmGrid_Staggered_z(dmGrid_Staggered_z),
dmGrid_Centered(dmGrid_Centered),
dmGrid_Shifted(dmGrid_Shifted),
dt(dt), iter(iter), Re(Re),
setting_transport(setting_transport),
penalization_u(penalization_u),
penalization_v(penalization_v),
penalization_w(penalization_w)
{
    DMCreateGlobalVector(dmGrid_Centered, &P);
    DMCreateGlobalVector(dmGrid_Staggered_x, &P_x);
    DMCreateGlobalVector(dmGrid_Staggered_y, &P_y);
    DMCreateGlobalVector(dmGrid_Staggered_z, &P_z);
    DMCreateGlobalVector(dmGrid_Centered, &Magnitude);
    DMCreateMatrix(dmGrid_Centered, &A);
    /*DMCreateGlobalVector(dmGrid_Staggered_x, &U_prova);
    DMCreateGlobalVector(dmGrid_Staggered_y, &V_prova);
    DMCreateGlobalVector(dmGrid_Staggered_z, &W_prova);*/

};

PetscErrorCode const solve()
{
    PetscFunctionBegin;

    this->assemble_lhs();

    parabolic_problem_x::assemble_lhs();
    parabolic_problem_y::assemble_lhs();
    parabolic_problem_z::assemble_lhs();

    
    for(size_t i = 0; i < 100; ++i){

        theta = d*d*(i)*dt;

        transport_problem_x transport_x(setting_transport, U_up, V_up, W_up, penalization_u, penalization_v, penalization_w);
        transport_problem_y transport_y(setting_transport, U_up, V_up, W_up, penalization_u, penalization_v, penalization_w);
        transport_problem_z transport_z(setting_transport, U_up, V_up, W_up, penalization_u, penalization_v, penalization_w);
        transport_x.solve_step_x(theta);
        transport_y.solve_step_y(theta);
        transport_z.solve_step_z(theta);

        U_prova = transport_x.get_U();
        V_prova = transport_y.get_V();
        W_prova = transport_z.get_W();

        VecCopy(U_prova, U_up);
        VecCopy(V_prova, V_up);
        VecCopy(W_prova, W_up);

        parabolic_problem_x::solve_step(theta);
        parabolic_problem_y::solve_step(theta);
        parabolic_problem_z::solve_step(theta);

        this->manage_pressure();
        /*Vec pressure_bench;
        DMCreateGlobalVector(dmGrid_Centered, &pressure_bench);
        CreateReferencePressure(dmGrid_Centered, pressure_bench, theta);
        CheckSolution(P, pressure_bench);
        VecDestroy(&pressure_bench);*/
        this->manage_pressure_x();
        this->manage_pressure_y();
        this->manage_pressure_z();
        theta = d*d*(i+1)*dt;
        this->update_velocity(theta);
        Vec bench;
        DMCreateGlobalVector(dmGrid_Staggered_y, &bench);
        CreateReferenceSolutionTry(dmGrid_Staggered_y, bench, theta);
        CheckSolution(V_up, bench);
        VecDestroy(&bench);
        compute_magnitude();

        PetscViewer viewer_magnitude;
        DM DM_magnitude;
        //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
        Vec magnitude;
        DMStagVecSplitToDMDA(dmGrid_Centered, Magnitude, ELEMENT, 0, &DM_magnitude, &magnitude);
        PetscObjectSetName((PetscObject)magnitude, "magnitude");
        char filename_magnitude[50]; 
        sprintf(filename_magnitude, "results/magnitude%03zu.vtr", i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_magnitude, FILE_MODE_WRITE, &viewer_magnitude);
        VecView(magnitude, viewer_magnitude);
        VecDestroy(&magnitude);
        DMDestroy(&DM_magnitude);
        PetscViewerDestroy(&viewer_magnitude);

        PetscViewer viewer_u;
        DM DM_u;
        //DMStagCreateCompatibleDMStag(dmGrid_Staggered_x, 0, 0, 1, 0, &DM_u);
        Vec u;
        DMStagVecSplitToDMDA(dmGrid_Staggered_x, U_up, LEFT, 0, &DM_u, &u);
        PetscObjectSetName((PetscObject)u, "x_component");
        char filename_u[50]; 
        sprintf(filename_u, "results/x_component%03zu.vtr", i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_x), filename_u, FILE_MODE_WRITE, &viewer_u);
        VecView(u, viewer_u);
        VecDestroy(&u);
        DMDestroy(&DM_u);
        PetscViewerDestroy(&viewer_u); 

        PetscViewer viewer_v;
        DM DM_v;
        //DMStagCreateCompatibleDMStag(dmGrid_Staggered_y, 0, 0, 1, 0, &DM_v);
        Vec v;
        DMStagVecSplitToDMDA(dmGrid_Staggered_y, V_up, DOWN, 0, &DM_v, &v);
        PetscObjectSetName((PetscObject)v, "y_component");
        char filename_v[50];
        sprintf(filename_v, "results/y_component%03zu.vtr", i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_y), filename_v, FILE_MODE_WRITE, &viewer_v);
        VecView(v, viewer_v);
        VecDestroy(&v);
        DMDestroy(&DM_v);
        PetscViewerDestroy(&viewer_v);

        PetscViewer viewer_w;
        DM DM_w;
        //DMStagCreateCompatibleDMStag(dmGrid_Staggered_z, 0, 0, 1, 0, &DM_w);
        Vec w;
        DMStagVecSplitToDMDA(dmGrid_Staggered_z, W_up, BACK, 0, &DM_w, &w);
        PetscObjectSetName((PetscObject)w, "z_component");
        char filename_w[50];
        sprintf(filename_w, "results/z_component%03zu.vtr", i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Staggered_z), filename_w, FILE_MODE_WRITE, &viewer_w);
        VecView(w, viewer_w);
        VecDestroy(&w);
        DMDestroy(&DM_w);
        PetscViewerDestroy(&viewer_w);

        PetscViewer viewer_p;
        DM DM_p;
        //DMStagCreateCompatibleDMStag(dmGrid_Centered, 0, 0, 0, 1, &DM_p);
        Vec p;
        DMStagVecSplitToDMDA(dmGrid_Centered, P, ELEMENT, 0, &DM_p, &p);
        PetscObjectSetName((PetscObject)p, "p");
        char filename_p[50];
        sprintf(filename_p, "results/p%03zu.vtr", i);
        PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid_Centered), filename_p, FILE_MODE_WRITE, &viewer_p);
        VecView(p, viewer_p);
        VecDestroy(&p);
        DMDestroy(&DM_p);
        PetscViewerDestroy(&viewer_p);
        std::cout << "Iteration " << i << " completed." << std::endl;     
        std::cout << "---------------------------------------------------------" << std::endl;




    }




    PetscFunctionReturn(0);
}

virtual ~navier_stokes_problem()
{
    VecDestroy(&P);
    VecDestroy(&P_x);
    VecDestroy(&P_y);
    VecDestroy(&P_z);
    VecDestroy(&Magnitude);
    MatDestroy(&A);
    /*VecDestroy(&U_prova);
    VecDestroy(&V_prova);
    VecDestroy(&W_prova);*/
    // NO NEED TO CALL: HO STAMPATO QUANDO ENTRA NEL DISTRUTTORE DI PARABOLIC_PROBLEM E LO FA DOPO cout << "Destructor called" << std::endl;
    /*pb_x.~parabolic_problem_x();
    pb_y.~parabolic_problem_y();
    pb_z.~parabolic_problem_z();*/
    std::cout << "Destructor called" << std::endl;

}

};

#endif // NAVIER_STOKES_PROBLEM_HPP


