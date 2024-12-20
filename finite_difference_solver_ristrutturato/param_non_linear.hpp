#define pi 3.1415926535



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

static PetscScalar a = pi/4;
static PetscScalar d = pi/2;
static PetscScalar par = 1;

PetscReal theta = 0.0;
//static PetscScalar A = 4*sqrt(2)/(3*sqrt(3));
/*static PetscScalar k = 2*pi;
static PetscScalar c = (5/6)*pi;
static PetscScalar d = (1/6)*pi;
static PetscScalar A = 1;
static PetscScalar B = 1;
static PetscScalar C = 1;*/


static PetscScalar uxRef(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar theta)

{
    //return  cos(pi*x)*sin(2*pi*y);
    //return  -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y))*exp(-theta)*par;
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    //return (sin(k*x - c)*cos(k*y - d)*sin(k*z) - cos(k*z - c)*sin(k*x - d)*sin(k*y))*exp(-theta)*1e-4;
    //return (A*sin(k*z) + C*cos(k*y))*exp(-theta);
    return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z;

}

static PetscScalar uyRef(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar theta)
{
    //return  -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z))*exp(-theta)*par;
    //return (sin(k*y - c)*cos(k*z - d)*sin(k*x) - cos(k*x - c)*sin(k*y - d)*sin(k*z))*exp(-theta)*1e-4;
    return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return -cos(2*pi*x)*sin(2*pi*y)*exp(-theta);
    //return (B*sin(k*x) + A*cos(k*z))*exp(-theta);
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);


}

static PetscScalar uzRef(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar theta)
{
    return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z);
    //return (sin(k*z - c)*cos(k*x - d)*sin(k*y) - cos(k*y - c)*sin(k*z - d)*sin(k*x))*exp(-theta)*1e-4;
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    //return 0.0*x + 0.0*y + 0.0*z;
    //return (C*sin(k*y) + B*cos(k*x))*exp(-theta);
    //return cos(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);

}


static PetscScalar solution(PetscScalar x, PetscScalar y, PetscScalar z, PetscScalar theta)
{
    //return (sin(k*x - c)*cos(k*y - d)*sin(k*z) - cos(k*z - c)*sin(k*x - d)*sin(k*y))*exp(-theta)*1e-4;
    //return -84*pi*pi*sin(2*pi*x)*cos(4*pi*y)*cos(8*pi*z);
    //return 0.0*x + 0.0*y + 0.0*z;
    //return sin((pi/3)*(x+y+z))*exp(-theta) + x*y*z;
    //return (A*sin(k*z) + C*cos(k*y))*exp(-theta);
    return cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) + 0.01*4*pi*cos(2*pi*x)*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)*cos(2*pi*z)*sin(2*pi*y) + 0.01*4*pi*cos(2*pi*x)*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*y)*cos(2*pi*z)*sin(2*pi*z) + 0.01*4*pi*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*y)*cos(2*pi*z)*cos(2*pi*z)*sin(2*pi*x);

}



static PetscErrorCode CreateReferenceSolutionFirst(DM dmGrid, Vec *pSolRef, PetscScalar theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmGrid, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmGrid, &solRefLocal);
    DMStagVecGetArray(dmGrid, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[1] + n[1] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iux] = uxRef(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, solRefLocal, &arrSol);
    DMLocalToGlobal(dmGrid, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &solRefLocal);

    return 0;
}

static PetscErrorCode CreateReferenceSolutionSecond(DM dmGrid, Vec *pSolRef, PetscScalar theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmGrid, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmGrid, &solRefLocal);
    DMStagVecGetArray(dmGrid, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[0] + n[0] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iuy] = uyRef(arrCoord[ez][ey][ex][icuy[0]], arrCoord[ez][ey][ex][icuy[1]], arrCoord[ez][ey][ex][icuy[2]], theta);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, solRefLocal, &arrSol);
    DMLocalToGlobal(dmGrid, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &solRefLocal);

    return 0;
}

static PetscErrorCode CreateReferenceSolutionThird(DM dmGrid, Vec *pSolRef, PetscScalar theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmGrid, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmGrid, &solRefLocal);
    DMStagVecGetArray(dmGrid, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[0] + n[0] && ey < start[1] + n[1]) 
                arrSol[ez][ey][ex][iuz] = uzRef(arrCoord[ez][ey][ex][icuz[0]], arrCoord[ez][ey][ex][icuz[1]], arrCoord[ez][ey][ex][icuz[2]], theta);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, solRefLocal, &arrSol);
    DMLocalToGlobal(dmGrid, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &solRefLocal);

    return 0;
}

static PetscErrorCode CreateReferenceSolutionTry(DM dmGrid, Vec *pSolRef, PetscScalar theta)
{
    PetscInt        start[3], n[3], nExtra[3], ex, ey, ez, d;
    PetscInt        iux, iuy, iuz, icux[3], icuy[3], icuz[3], iue, icue[3];
    Vec             solRef, solRefLocal, coord, coordLocal;
    DM              dmCoord;
    PetscScalar ****arrSol, ****arrCoord;

    PetscFunctionBeginUser;
    DMCreateGlobalVector(dmGrid, pSolRef);
    solRef = *pSolRef;
    DMStagGetCorners(dmGrid, &start[0], &start[1], &start[2], &n[0], &n[1], &n[2], &nExtra[0], &nExtra[1], &nExtra[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);
    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz);
    DMStagGetLocationSlot(dmGrid, ELEMENT, 0, &iue);
    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz[d]);
        DMStagGetLocationSlot(dmCoord, ELEMENT, d, &icue[d]);
    }
    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);
    DMGetLocalVector(dmGrid, &solRefLocal);
    DMStagVecGetArray(dmGrid, solRefLocal, &arrSol);
    for (ez = start[2]; ez < start[2] + n[2] + nExtra[2]; ++ez) {
        for (ey = start[1]; ey < start[1] + n[1] + nExtra[1]; ++ey) {
            for (ex = start[0]; ex < start[0] + n[0] + nExtra[0]; ++ex) {
                //if (ex < start[1] + n[1] && ey < start[2] + n[2]) 
                arrSol[ez][ey][ex][iux] = solution(arrCoord[ez][ey][ex][icux[0]], arrCoord[ez][ey][ex][icux[1]], arrCoord[ez][ey][ex][icux[2]], theta);
            }
        }
    }
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMStagVecRestoreArray(dmGrid, solRefLocal, &arrSol);
    DMLocalToGlobal(dmGrid, solRefLocal, INSERT_VALUES, solRef);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMRestoreLocalVector(dmGrid, &solRefLocal);

    return 0;
}

// Create Domain routines
static PetscErrorCode CreateGrid(DM* pdmGrid, PetscInt dof1, PetscInt dof2, PetscInt dof3, PetscScalar nx, PetscScalar ny, PetscScalar nz, PetscReal Lx_0, PetscReal Lx, PetscReal Ly_0, PetscReal Ly, PetscReal Lz_0, PetscReal Lz)
{
    DM dmGrid;
    const PetscInt dof0 = 0;
    const PetscInt stencilWidth = 1;
    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, nx, ny, nz, PETSC_DECIDE,
                   PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL,
                   NULL, pdmGrid);
    dmGrid = *pdmGrid;
    DMSetFromOptions(dmGrid);
    DMSetUp(dmGrid);
    DMStagSetUniformCoordinatesExplicit(dmGrid, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    //PetscObjectDestroy((PetscObject*)&dmGrid);

    return 0;
}


PetscErrorCode CheckSolution(Vec const & sol, Vec const & solRef)
{
    Vec       diff;
    PetscReal normsolRef, errAbs, errRel;
    PetscFunctionBegin;

    PetscFunctionBegin;
    VecDuplicate(sol, &diff);
    VecCopy(sol, diff);
    VecAXPY(diff, -1.0, solRef);
    VecNorm(diff, NORM_2, &errAbs);
    VecNorm(solRef, NORM_2, &normsolRef);
    errRel = errAbs / normsolRef;
    PetscPrintf(PETSC_COMM_WORLD, "Error (abs): %g\nError (rel): %g\n", (double)errAbs, (double)errRel);
    VecDestroy(&diff);
    PetscFunctionReturn(0);
}


