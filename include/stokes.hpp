
#ifndef STOKES_PROBLEM_HPP
#define STOKES_PROBLEM_HPP

class stokes_problem: public parabolic_problem_x, public parabolic_problem_y, public parabolic_problem_z {
private:

/*DM const & dmGrid_Staggered_x;
DM const & dmGrid_Staggered_y;
DM const & dmGrid_Staggered_z;*/
DM const & dmGrid_Staggered;
DM const & dmGrid_Centered;
DM const & dmGrid_Shifted;
Vec U_pre, V_pre, W_pre;
//Vec U_up, V_up, W_up;
Vec Magnitude;
Vec P, P_x, P_y, P_z;
//Vec U_0, V_0, W_0;




public:
    stokes_problem(DM const & dmGrid_Staggered_x, DM const & dmGrid_Staggered_y, DM const & dmGrid_Staggered_z, DM const & dmGrid_Staggered, DM const & dmGrid_Centered, DM const & dmGrid_Shifted, Vec const & U_init, Vec const & V_init, Vec const & W_init) : 
    parabolic_problem_x(dmGrid_Staggered_x, U_init), parabolic_problem_y(dmGrid_Staggered_y, V_init), parabolic_problem_z(dmGrid_Staggered_z, W_init),
    dmGrid_Staggered(dmGrid_Staggered), dmGrid_Centered(dmGrid_Centered), dmGrid_Shifted(dmGrid_Shifted)
    {

        DMCreateGlobalVector(dmGrid_Staggered_x, &U_pre);
        DMCreateGlobalVector(dmGrid_Staggered_y, &V_pre);
        DMCreateGlobalVector(dmGrid_Staggered_z, &W_pre);
        /*DMCreateGlobalVector(dmGrid_Staggered_x, &U_up);
        DMCreateGlobalVector(dmGrid_Staggered_y, &V_up);
        DMCreateGlobalVector(dmGrid_Staggered_z, &W_up);*/
        DMCreateGlobalVector(dmGrid_Centered, &Magnitude);
        DMCreateGlobalVector(dmGrid_Centered, &P);
        DMCreateGlobalVector(dmGrid_Staggered_x, &P_x);
        DMCreateGlobalVector(dmGrid_Staggered_y, &P_y);
        DMCreateGlobalVector(dmGrid_Staggered_z, &P_z);
        /*DMCreateGlobalVector(dmGrid_Staggered_x, &U_0);
        DMCreateGlobalVector(dmGrid_Staggered_y, &V_0);
        DMCreateGlobalVector(dmGrid_Staggered_z, &W_0);
        VecCopy(U_init, U_0);
        VecCopy(V_init, V_0);
        VecCopy(W_init, W_0);*/

    }

    /*
    annotazione importante: la condizione di compatibilita' per il toerema di Neumann in teoria non e' rispettata in Chorin Teman
    */

    PetscErrorCode AttachNullspace(DM const & dmGrid, Mat & A)
    {
    Vec          basis;
    MatNullSpace matNullSpace;
    PetscFunctionBegin;
    DMCreateGlobalVector(dmGrid, &basis);
    VecSet(basis, 1.0);
    MatNullSpaceCreate(PetscObjectComm((PetscObject)dmGrid), PETSC_FALSE, 1, &basis, &matNullSpace);
    VecDestroy(&basis);
    MatSetNullSpace(A, matNullSpace);
    MatNullSpaceDestroy(&matNullSpace);
    PetscFunctionReturn(0);
    }

    PetscErrorCode Assemble_P(DM const & dmGrid, Mat & A, Vec & rhs, Vec const & source) 
    {
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;

        PetscFunctionBegin;

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
        PetscReal const hx = 1.0 / N[0];
        PetscReal const hy = 1.0 / N[1];
        PetscReal const hz = 1.0 / N[2];

        Vec local;
        DMCreateLocalVector(dmGrid,&local);
        DMGlobalToLocalBegin(dmGrid,source,INSERT_VALUES,local);
        DMGlobalToLocalEnd(dmGrid,source,INSERT_VALUES,local);

        for (ez = startz; ez < startz + nz; ++ez) { /* With DMStag, always iterate x fastest, y second fastest, z slowest */
            for (ey = starty; ey < starty + ny; ++ey) {
                for (ex = startx; ex < startx + nx; ++ex) {
                    if (ex == N[0] - 1) {

                        DMStagStencil row, col[6];
                        PetscReal valA[6], valRhs;
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                           
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                          
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                      
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                            col[4].loc = LEFT;
                            col[4].c = 0;
                            valA[4] = 1.0 / (hz * hz);
                            col[5].i = ex;
                            col[5].j = ey;
                            col[5].k = ez + 1;
                            col[5].loc = ELEMENT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        }
                        DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                    } 
                    else if (ex == 0) {
                        DMStagStencil row, col[6];
                        PetscReal valA[6], valRhs;
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                               
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                          
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                      
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        }
                        DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);                              
                    } else {
                        /* X-momentum interior equation : (u_xx + u_yy + u_zz) - p_x = f^x */
                        DMStagStencil row, col[7];
                        PetscReal valA[7], valRhs;
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                               
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                          
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);                      
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                                /* Missing up term */
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
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
                                DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                                DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                            
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
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);                        
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
                            col[5].loc = LEFT;
                            col[5].c = 0;
                            valA[5] = 1.0 / (hz * hz);
                            col[6].i = ex;
                            col[6].j = ey;
                            col[6].k = ez + 1;
                            col[6].loc = ELEMENT;
                            col[6].c = 0;
                            valA[6] = 1.0 / (hz * hz);
                            DMStagVecGetValuesStencil(dmGrid, local, 1, &row, &valRhs);
                            DMStagVecSetValuesStencil(dmGrid, rhs, 1, &row, &valRhs, INSERT_VALUES);
                        }
                        DMStagMatSetValuesStencil(dmGrid, A, 1, &row, nEntries, col, valA, INSERT_VALUES);
                    }                
                }
            }
        }
        VecDestroy(&local);
        MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(rhs);
        VecAssemblyEnd(rhs);

        PetscFunctionReturn(0); 
    }

    PetscErrorCode AssembleDivergence(DM const & dmGrid, Vec & div, Vec const & U, Vec const &  V, Vec const & W, PetscReal const & dt) 
    {
        PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
        DM dmCoord;
        Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

        PetscFunctionBegin;

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
        PetscReal const hx = 1.0 / N[0];
        PetscReal const hy = 1.0 / N[1];
        PetscReal const hz = 1.0 / N[2];
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid, LEFT, 0, &iu_left);
        DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iu_right);
        DMStagGetLocationSlot(dmGrid, UP, 0, &iu_up);
        DMStagGetLocationSlot(dmGrid, DOWN, 0, &iu_down);
        DMStagGetLocationSlot(dmGrid, FRONT, 0, &iu_front);
        DMStagGetLocationSlot(dmGrid, BACK, 0, &iu_back);
        DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iu_element);

        DMCreateLocalVector(dmGrid, &vecULocal);
        DMGlobalToLocalBegin(dmGrid, U, INSERT_VALUES, vecULocal);
        DMGlobalToLocalEnd(dmGrid, U, INSERT_VALUES, vecULocal);
        DMStagVecGetArrayRead(dmGrid, vecULocal, &arrU);

        DMCreateLocalVector(dmGrid, &vecVLocal);
        DMGlobalToLocalBegin(dmGrid, V, INSERT_VALUES, vecVLocal);
        DMGlobalToLocalEnd(dmGrid, V, INSERT_VALUES, vecVLocal);
        DMStagVecGetArrayRead(dmGrid, vecVLocal, &arrV);

        DMCreateLocalVector(dmGrid, &vecWLocal);
        DMGlobalToLocalBegin(dmGrid, W, INSERT_VALUES, vecWLocal);
        DMGlobalToLocalEnd(dmGrid, W, INSERT_VALUES, vecWLocal);
        DMStagVecGetArrayRead(dmGrid, vecWLocal, &arrW);    

        DMGetLocalVector(dmGrid, &vecOutLocal);
        DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

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
        DMStagVecRestoreArrayRead(dmGrid, vecULocal, &arrU);
        DMStagVecRestoreArrayRead(dmGrid, vecVLocal, &arrV);
        DMStagVecRestoreArrayRead(dmGrid, vecWLocal, &arrW);
        DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
        DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, div);    
        DMRestoreLocalVector(dmGrid, &vecOutLocal);
        DMRestoreLocalVector(dmGrid, &vecULocal);
        DMRestoreLocalVector(dmGrid, &vecVLocal);
        DMRestoreLocalVector(dmGrid, &vecWLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

    PetscErrorCode ComputeDivergence(DM const & dmGrid_centered, DM const & dmGrid_shifted, DM const & dmGrid_staggered, Vec & div, Vec const & U_n, Vec const & V_n, Vec const & W_n, PetscReal const & dt) 
    {
            PetscFunctionBegin;

            Vec U_shifted, V_shifted, W_shifted;
            DMCreateGlobalVector(dmGrid_shifted, &U_shifted);
            DMCreateGlobalVector(dmGrid_shifted, &V_shifted);
            DMCreateGlobalVector(dmGrid_shifted, &W_shifted);
            DMStagMigrateVec(dmGrid_staggered, U_n, dmGrid_shifted, U_shifted);
            DMStagMigrateVec(dmGrid_staggered, V_n, dmGrid_shifted, V_shifted);
            DMStagMigrateVec(dmGrid_staggered, W_n, dmGrid_shifted, W_shifted);

            Vec div_shifted;
            DMCreateGlobalVector(dmGrid_shifted, &div_shifted);
            AssembleDivergence(dmGrid_shifted, div_shifted, U_shifted, V_shifted, W_shifted, dt);
            DMStagMigrateVec(dmGrid_shifted, div_shifted, dmGrid_centered, div);

            VecDestroy(&U_shifted);
            VecDestroy(&V_shifted);
            VecDestroy(&W_shifted);
            VecDestroy(&div_shifted);

            PetscFunctionReturn(0); 
    }

    PetscErrorCode Derive_x_P(DM const & dmGrid, Vec & P_x, Vec const & vec)
    {
        PetscInt iux_left, iux_right, iux_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
        DM dmCoord;
        Vec vecLocal, vecOutLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec, ****arrOut;    

        PetscFunctionBegin;

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
        PetscReal const hx = 1.0 / N[0];
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
        DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);
        DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iux_element);

        DMCreateLocalVector(dmGrid, &vecLocal);
        DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

        DMGetLocalVector(dmGrid, &vecOutLocal);
        DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);
    
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
                        arrOut[ez][ey][ex][iux_left] = 0.0;
                    }        
                    if(ex == N[0] - 1){
                        arrOut[ez][ey][ex][iux_right] = 0.0;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
        DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
        DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, P_x);
        DMRestoreLocalVector(dmGrid, &vecOutLocal);
        DMRestoreLocalVector(dmGrid, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

    PetscErrorCode Derive_y_P(DM const & dmGrid, Vec & P_y, Vec const & vec)
    {
        PetscInt iuy_up, iuy_down, iuy_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
        DM dmCoord;
        Vec vecLocal, vecOutLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec, ****arrOut;    

        PetscFunctionBegin;

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
        PetscReal const hy = 1.0 / N[1];
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);
        DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);
        DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iuy_element);

        DMCreateLocalVector(dmGrid, &vecLocal);
        DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

        DMGetLocalVector(dmGrid, &vecOutLocal);
        DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut); 

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
                        arrOut[ez][ey][ex][iuy_down] = 0.0;
                    }
                    if(ey == N[1] - 1){
                        arrOut[ez][ey][ex][iuy_up] = 0.0;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
        DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
        DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, P_y);
        DMRestoreLocalVector(dmGrid, &vecOutLocal);
        DMRestoreLocalVector(dmGrid, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

    PetscErrorCode Derive_z_P(DM const & dmGrid, Vec & P_z, Vec const & vec)
    {
        PetscInt iuz_back, iuz_front, iuz_element;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez;
        DM dmCoord;
        Vec vecLocal, vecOutLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec, ****arrOut;    

        PetscFunctionBegin;

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
        PetscReal const hz = 1.0 / N[2];
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);
        DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front); 
        DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iuz_element);
        
        DMCreateLocalVector(dmGrid, &vecLocal);
        DMGlobalToLocalBegin(dmGrid, vec, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid, vec, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

        DMGetLocalVector(dmGrid, &vecOutLocal);
        DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut); 

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
                        arrOut[ez][ey][ex][iuz_back] = 0.0;
                    }
                    if(ez == N[2] - 1){
                        arrOut[ez][ey][ex][iuz_front] = 0.0;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArrayRead(dmGrid, vecLocal, &arrVec);
        DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
        DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, P_z);
        DMRestoreLocalVector(dmGrid, &vecOutLocal);
        DMRestoreLocalVector(dmGrid, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0); 
    }

    //capire se preconditioner va bene ref articolo islamico
    PetscErrorCode ManagePressure(DM const & dmGrid_centered, DM const & dmGrid_shifted, DM const & dmGrid_staggered, PetscReal const & dt, Vec & P, Vec const & U_n, Vec const & V_n, Vec const & W_n)
    {
        Mat A;
        Vec rhs;
        KSP ksp;
        PC  pc;

        PetscFunctionBegin;

        DMCreateGlobalVector(dmGrid_centered, &rhs);
        DMCreateMatrix(dmGrid_centered, &A);

        Vec div;
        DMCreateGlobalVector(dmGrid_centered, &div);
        ComputeDivergence(dmGrid_centered, dmGrid_shifted, dmGrid_staggered, div, U_n, V_n, W_n, dt); 

        Assemble_P(dmGrid_centered, A, rhs, div);
        VecDestroy(&div);

        AttachNullspace(dmGrid_centered, A);

        /*//questo sistema non e' precondizionato: il preconditioner e' nullo. Bisogna farlo ed e' importante: punto piu' lento del codice
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetType(ksp, KSPGMRES);
        KSPSetOperators(ksp, A, A);
        //KSPGetPC(ksp, &pc);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, rhs, P);*/

        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetType(ksp, KSPGMRES);
        KSPSetOperators(ksp, A, A);
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "pilut");
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, rhs, P);

        MatDestroy(&A);
        VecDestroy(&rhs);
        KSPDestroy(&ksp);
        
        PetscFunctionReturn(0); 
    }

    //sistemare il create compatible dmstag (evito di passare tutti quei parametri) ma non so perche non riesco a farlo
    PetscErrorCode ManagePressure_x(DM const & dmGrid_staggered, DM const & dmGrid_centered, DM const & dmGrid_shifted, Vec & P_x, Vec const & P)
    {
        PetscFunctionBegin;

        Vec P_shifted;
        DMCreateGlobalVector(dmGrid_shifted, &P_shifted);
        DMStagMigrateVec(dmGrid_centered, P, dmGrid_shifted, P_shifted);
        
        Vec P_x_shifted;
        DMCreateGlobalVector(dmGrid_shifted, &P_x_shifted);
        Derive_x_P(dmGrid_shifted, P_x_shifted, P_shifted);
        DMStagMigrateVec(dmGrid_shifted, P_x_shifted, dmGrid_staggered, P_x);

        VecDestroy(&P_x_shifted);
        VecDestroy(&P_shifted);
        /*DMDestroy(&dmGrid_centered);
        DMDestroy(&dmGrid_shifted);*/

        PetscFunctionReturn(0); 
    }

    PetscErrorCode ManagePressure_y(DM const & dmGrid_staggered, DM const & dmGrid_centered, DM const & dmGrid_shifted, Vec & P_y, Vec const & P)
    {
        PetscFunctionBegin

        Vec P_shifted;
        DMCreateGlobalVector(dmGrid_shifted, &P_shifted);
        DMStagMigrateVec(dmGrid_centered, P, dmGrid_shifted, P_shifted);

        Vec P_y_shifted;
        DMCreateGlobalVector(dmGrid_shifted, &P_y_shifted);
        Derive_y_P(dmGrid_shifted, P_y_shifted, P_shifted);
        DMStagMigrateVec(dmGrid_shifted, P_y_shifted, dmGrid_staggered, P_y);

        VecDestroy(&P_y_shifted);
        VecDestroy(&P_shifted);

        PetscFunctionReturn(0); 
    }

    PetscErrorCode ManagePressure_z(DM const & dmGrid_staggered, DM const & dmGrid_centered, DM const & dmGrid_shifted, Vec & P_z, Vec const & P)
    {
        PetscFunctionBegin

        Vec P_shifted;
        DMCreateGlobalVector(dmGrid_shifted, &P_shifted);
        DMStagMigrateVec(dmGrid_centered, P, dmGrid_shifted, P_shifted);

        Vec P_z_shifted;
        DMCreateGlobalVector(dmGrid_shifted, &P_z_shifted);
        Derive_z_P(dmGrid_shifted, P_z_shifted, P_shifted);
        DMStagMigrateVec(dmGrid_shifted, P_z_shifted, dmGrid_staggered, P_z);

        VecDestroy(&P_z_shifted);
        VecDestroy(&P_shifted);

        PetscFunctionReturn(0); 
    }


    PetscErrorCode UpdatebcU(DM const & dmGrid, Vec & U_up, PetscReal const & theta)
    {
        PetscInt icux_left[3], icux_right[3], iux_left, iux_right;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
        DM dmCoord;
        Vec vecLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec; 

        PetscFunctionBegin;

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
            DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        }  
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
        DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);

        DMCreateLocalVector(dmGrid, &vecLocal);
        DMGlobalToLocalBegin(dmGrid, U_up, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid, U_up, INSERT_VALUES, vecLocal);
        DMStagVecGetArray(dmGrid, vecLocal, &arrVec);

        for (ez = startz; ez < startz + nz; ++ez) { 
            for (ey = starty; ey < starty + ny; ++ey) {
                for (ex = startx; ex < startx + nx; ++ex) {

                    if (ex == N[0] - 1) {
                        PetscReal val;
                        val = uxRef(arrCoord[ez][ey][ex][icux_right[0]], arrCoord[ez][ey][ex][icux_right[1]], arrCoord[ez][ey][ex][icux_right[2]], theta);
                        arrVec[ez][ey][ex][iux_right] = val;
                        
                    } else if(ex == 0) {
                        PetscReal val;
                        val = uxRef(arrCoord[ez][ey][ex][icux_left[0]], arrCoord[ez][ey][ex][icux_left[1]], arrCoord[ez][ey][ex][icux_left[2]], theta);
                        arrVec[ez][ey][ex][iux_left] = val;
                    }

                }
            }
        }


        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
        DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, U_up);
        DMRestoreLocalVector(dmGrid, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

    PetscErrorCode UpdatebcV(DM const & dmGrid, Vec & V_up, PetscReal const & theta) 
    {

        PetscInt icuy_down[3], icuy_up[3], iuy_down, iuy_up;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
        DM dmCoord;
        Vec vecLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec;   

        PetscFunctionBegin;

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy_down[d]);
            DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        } 
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);
        DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);

        DMCreateLocalVector(dmGrid, &vecLocal);
        DMGlobalToLocalBegin(dmGrid, V_up, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid, V_up, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

        for (ez = startz; ez < startz + nz; ++ez) {
            for (ey = starty; ey < starty + ny; ++ey) {
                for (ex = startx; ex < startx + nx; ++ex) {

                    if (ey == N[1] - 1) {
                        PetscReal val;
                        val = uyRef(arrCoord[ez][ey][ex][icuy_up[0]], arrCoord[ez][ey][ex][icuy_up[1]], arrCoord[ez][ey][ex][icuy_up[2]], theta);
                        arrVec[ez][ey][ex][iuy_up] = val;
                    } else if(ey == 0) {
                        PetscReal val;
                        val = uyRef(arrCoord[ez][ey][ex][icuy_down[0]], arrCoord[ez][ey][ex][icuy_down[1]], arrCoord[ez][ey][ex][icuy_down[2]], theta);
                        arrVec[ez][ey][ex][iuy_down] = val;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
        DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, V_up);
        DMRestoreLocalVector(dmGrid, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }

    PetscErrorCode UpdatebcW(DM const & dmGrid, Vec & W_up, PetscReal const & theta) 
    {

        PetscInt icuz_back[3], icuz_front[3], iuz_back, iuz_front;
        PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
        DM dmCoord;
        Vec vecLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrVec;   

        PetscFunctionBegin;

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

        for (d = 0; d < 3; ++d) {
            DMStagGetLocationSlot(dmCoord, BACK, d, &icuz_back[d]);
            DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
        }  
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);
        DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front); 
        
        DMCreateLocalVector(dmGrid, &vecLocal);
        DMGlobalToLocalBegin(dmGrid, W_up, INSERT_VALUES, vecLocal);
        DMGlobalToLocalEnd(dmGrid, W_up, INSERT_VALUES, vecLocal);
        DMStagVecGetArrayRead(dmGrid, vecLocal, &arrVec);

        for (ez = startz; ez < startz + nz; ++ez) {
            for (ey = starty; ey < starty + ny; ++ey) {
                for (ex = startx; ex < startx + nx; ++ex) {

                    if (ez == N[2] - 1) {
                        PetscReal val;
                        val = uzRef(arrCoord[ez][ey][ex][icuz_front[0]], arrCoord[ez][ey][ex][icuz_front[1]], arrCoord[ez][ey][ex][icuz_front[2]], theta);
                        arrVec[ez][ey][ex][iuz_front] = val;
                    } else if(ez == 0) {
                        PetscReal val;
                        val = uzRef(arrCoord[ez][ey][ex][icuz_back[0]], arrCoord[ez][ey][ex][icuz_back[1]], arrCoord[ez][ey][ex][icuz_back[2]], theta);
                        arrVec[ez][ey][ex][iuz_back] = val;
                    }
                }
            }
        }

        DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
        DMStagVecRestoreArray(dmGrid, vecLocal, &arrVec);
        DMLocalToGlobal(dmGrid, vecLocal, INSERT_VALUES, W_up);
        DMRestoreLocalVector(dmGrid, &vecLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);
    }


    PetscErrorCode UpdateVelocity(DM const & dmGrid_staggered_x, DM const & dmGrid_staggered_y, DM const & dmGrid_staggered_z, PetscReal const & dt, Vec & U_up, Vec & V_up, Vec & W_up, Vec const & P_x, Vec const & P_y, Vec const & P_z, Vec const & U_pre, Vec const & V_pre, Vec const & W_pre, PetscReal const & theta)
    {
        PetscFunctionBegin;

        VecAXPY(U_pre, -dt, P_x);
        VecAXPY(V_pre, -dt, P_y);
        VecAXPY(W_pre, -dt, P_z);

        VecCopy(U_pre, U_up);
        VecCopy(V_pre, V_up);
        VecCopy(W_pre, W_up);

        UpdatebcU(dmGrid_staggered_x, U_up, theta);
        UpdatebcV(dmGrid_staggered_y, V_up, theta);
        UpdatebcW(dmGrid_staggered_z, W_up, theta);


        PetscFunctionReturn(0);
    }



    PetscErrorCode AssembleMagnitude(DM const & dmGrid, Vec & magnitude, Vec const & U, Vec const & V, Vec const & W) 
    {
        PetscInt iu_left, iu_right, iu_up, iu_down, iu_front, iu_back, iu_element;
        PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez;
        DM dmCoord;
        Vec vecULocal, vecVLocal, vecWLocal, vecOutLocal, coord, coordLocal;
        PetscReal ****arrCoord, ****arrU, ****arrV, ****arrW, ****arrOut;    

        PetscFunctionBegin;

        DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
        DMGetCoordinateDM(dmGrid, &dmCoord);

        DMGetCoordinates(dmGrid, &coord);
        DMGetLocalVector(dmCoord, &coordLocal);
        DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
        DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

        DMStagGetLocationSlot(dmGrid, LEFT, 0, &iu_left);
        DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iu_right);
        DMStagGetLocationSlot(dmGrid, UP, 0, &iu_up);
        DMStagGetLocationSlot(dmGrid, DOWN, 0, &iu_down);
        DMStagGetLocationSlot(dmGrid, FRONT, 0, &iu_front);
        DMStagGetLocationSlot(dmGrid, BACK, 0, &iu_back);
        DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iu_element);

        DMCreateLocalVector(dmGrid, &vecULocal);
        DMGlobalToLocalBegin(dmGrid, U, INSERT_VALUES, vecULocal);
        DMGlobalToLocalEnd(dmGrid, U, INSERT_VALUES, vecULocal);
        DMStagVecGetArrayRead(dmGrid, vecULocal, &arrU);

        DMCreateLocalVector(dmGrid, &vecVLocal);
        DMGlobalToLocalBegin(dmGrid, V, INSERT_VALUES, vecVLocal);
        DMGlobalToLocalEnd(dmGrid, V, INSERT_VALUES, vecVLocal);
        DMStagVecGetArrayRead(dmGrid, vecVLocal, &arrV);

        DMCreateLocalVector(dmGrid, &vecWLocal);
        DMGlobalToLocalBegin(dmGrid, W, INSERT_VALUES, vecWLocal);
        DMGlobalToLocalEnd(dmGrid, W, INSERT_VALUES, vecWLocal);
        DMStagVecGetArrayRead(dmGrid, vecWLocal, &arrW);    

        DMGetLocalVector(dmGrid, &vecOutLocal);
        DMStagVecGetArray(dmGrid, vecOutLocal, &arrOut);

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
        DMStagVecRestoreArrayRead(dmGrid, vecULocal, &arrU);
        DMStagVecRestoreArrayRead(dmGrid, vecVLocal, &arrV);
        DMStagVecRestoreArrayRead(dmGrid, vecWLocal, &arrW);
        DMStagVecRestoreArray(dmGrid, vecOutLocal, &arrOut);
        DMLocalToGlobal(dmGrid, vecOutLocal, INSERT_VALUES, magnitude);    
        DMRestoreLocalVector(dmGrid, &vecOutLocal);
        DMRestoreLocalVector(dmGrid, &vecULocal);
        DMRestoreLocalVector(dmGrid, &vecVLocal);
        DMRestoreLocalVector(dmGrid, &vecWLocal);
        DMRestoreLocalVector(dmCoord, &coordLocal);

        PetscFunctionReturn(0);    
    }

    PetscErrorCode ComputeMagnitude(DM const & dmGrid_Staggered_x, DM const & dmGrid_Staggered_y, DM const & dmGrid_Staggered_z, DM const & dmGrid_Centered, DM const & dmGrid_Shifted, Vec & magnitude, Vec const & U, Vec const & V, Vec const & W)
    {
        
        PetscFunctionBegin;

        Vec U_shifted;
        DMCreateGlobalVector(dmGrid_Shifted, &U_shifted);
        DMStagMigrateVec(dmGrid_Staggered_x, U, dmGrid_Shifted, U_shifted);
        Vec V_shifted;
        DMCreateGlobalVector(dmGrid_Shifted, &V_shifted);
        DMStagMigrateVec(dmGrid_Staggered_y, V, dmGrid_Shifted, V_shifted);
        Vec W_shifted;
        DMCreateGlobalVector(dmGrid_Shifted, &W_shifted);
        DMStagMigrateVec(dmGrid_Staggered_z, W, dmGrid_Shifted, W_shifted);


        Vec magnitude_shifted;
        DMCreateGlobalVector(dmGrid_Shifted, &magnitude_shifted);
        AssembleMagnitude(dmGrid_Shifted, magnitude_shifted, U_shifted, V_shifted, W_shifted);
        DMStagMigrateVec(dmGrid_Shifted, magnitude_shifted, dmGrid_Centered, magnitude);

        VecDestroy(&magnitude_shifted);
        VecDestroy(&U_shifted);  
        VecDestroy(&V_shifted);
        VecDestroy(&W_shifted);
                
        PetscFunctionReturn(0); 

    }

   
    PetscErrorCode solve()
    {
        PetscFunctionBegin;
    
        for(size_t i = 0; i < iter; ++i){
            if (i == 0){

                theta = 3*(vRef/Re)*k*k*dt*(i);

                parabolic_problem_x::ManageViscosity(dmGrid_Staggered_x, dt, Re, U_pre, U, theta);
                parabolic_problem_y::ManageViscosity(dmGrid_Staggered_y, dt, Re, V_pre, V, theta);
                parabolic_problem_z::ManageViscosity(dmGrid_Staggered_z, dt, Re, W_pre, W, theta);
                //ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_pre, V_pre, W_pre, U, V, W, theta);

                //std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
                ManagePressure(dmGrid_Centered, dmGrid_Shifted, dmGrid_Staggered, dt, P, U_pre, V_pre, W_pre);
                ManagePressure_x(dmGrid_Staggered_x, dmGrid_Centered, dmGrid_Shifted, P_x, P);
                ManagePressure_y(dmGrid_Staggered_y, dmGrid_Centered, dmGrid_Shifted, P_y, P);
                ManagePressure_z(dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, P_z, P);
                
                theta = 3*(vRef/Re)*k*k*dt*(i+1);

                UpdateVelocity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, U_up, V_up, W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, theta);

                Vec bench;
                DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
                CreateReferenceSolutionTry(dmGrid_Staggered_x, bench, theta);
                CheckSolution(U_up, bench);
                VecDestroy(&bench);

                ComputeMagnitude(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, Magnitude, U_up, V_up, W_up);

                dump(i);
            } else {

                
                theta = 3*(vRef/Re)*k*k*dt*(i);

                parabolic_problem_x::ManageViscosity(dmGrid_Staggered_x, dt, Re, U_pre, U_up, theta);
                parabolic_problem_y::ManageViscosity(dmGrid_Staggered_y, dt, Re, V_pre, V_up, theta);
                parabolic_problem_z::ManageViscosity(dmGrid_Staggered_z, dt, Re, W_pre, W_up, theta);
                //ManageViscosity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, Re, U_pre, V_pre, W_pre, U_up, V_up, W_up, theta);

                //std::cout << "Spirit of Nebraska completed:   diffusion done." << std::endl;
                ManagePressure(dmGrid_Centered, dmGrid_Shifted, dmGrid_Staggered, dt, P, U_pre, V_pre, W_pre);
                ManagePressure_x(dmGrid_Staggered_x, dmGrid_Centered, dmGrid_Shifted, P_x, P);
                ManagePressure_y(dmGrid_Staggered_y, dmGrid_Centered, dmGrid_Shifted, P_y, P);
                ManagePressure_z(dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, P_z, P);
                theta = 3*(vRef/Re)*k*k*dt*(i+1);
                UpdateVelocity(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dt, U_up, V_up, W_up, P_x, P_y, P_z, U_pre, V_pre, W_pre, theta);
                //std::cout << "Spirit of California completed: iteration " << i << " done" << std::endl;
                Vec bench;
                DMCreateGlobalVector(dmGrid_Staggered_x, &bench);
                CreateReferenceSolutionTry(dmGrid_Staggered_x, bench, theta);
                CheckSolution(U_up, bench);
                VecDestroy(&bench);


                ComputeMagnitude(dmGrid_Staggered_x, dmGrid_Staggered_y, dmGrid_Staggered_z, dmGrid_Centered, dmGrid_Shifted, Magnitude, U_up, V_up, W_up);
                dump(i);

            }
        }

        PetscFunctionReturn(0);
    };

    PetscErrorCode dump(size_t const & i) override
    {
        PetscFunctionBegin;

        parabolic_problem_x::dump(i);
        parabolic_problem_y::dump(i);
        parabolic_problem_z::dump(i);
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
        PetscFunctionReturn(0); 
    
    }

    ~stokes_problem()
    {
        /*VecDestroy(&U_0);
        VecDestroy(&V_0);
        VecDestroy(&W_0);*/
        VecDestroy(&U_pre);
        VecDestroy(&V_pre);
        VecDestroy(&W_pre);
        /*VecDestroy(&U_up);
        VecDestroy(&V_up);
        VecDestroy(&W_up);*/
        VecDestroy(&P);
        VecDestroy(&P_x);
        VecDestroy(&P_y);
        VecDestroy(&P_z);
        VecDestroy(&Magnitude);
    }
};

#endif // STOKES_PROBLEM_HPP
