#include <petsc.h>
#include <iostream>


static char help[] = "DMStag vector operations example in 3D.\n\n";

PetscErrorCode ComputeDerivative(DM dm, Vec x, Vec dy) {
    PetscErrorCode ierr;
    Vec            local_x, local_dy;
    PetscInt       startx, starty, startz, nx, ny, nz;
    PetscScalar    ****arr_x, ****arr_dy;

    // Create local vectors for input and output
    ierr = DMCreateLocalVector(dm, &local_x); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(dm, &local_dy); CHKERRQ(ierr);

    // Scatter global to local
    ierr = DMGlobalToLocalBegin(dm, x, INSERT_VALUES, local_x); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(dm, x, INSERT_VALUES, local_x); CHKERRQ(ierr);

    ierr = DMStagGetCorners(dm, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL); CHKERRQ(ierr);

    // Access arrays for direct manipulation
    ierr = DMDAVecGetArrayDOF(dm, local_x, &arr_x); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(dm, local_dy, &arr_dy); CHKERRQ(ierr);

    for (PetscInt k = startz; k < startz + nz; ++k) {
        for (PetscInt j = starty; j < starty + ny; ++j) {
            for (PetscInt i = startx; i < startx + nx; ++i) {
                PetscScalar center_val = 0.0, below_val = 0.0;

                // Get the center value
                ierr = DMStagVecGetArrayDOF(dm, local_x, k, j, i, &center_val); CHKERRQ(ierr);

                // Get the below value (boundary condition: value outside grid is 0)
                if (j > starty) {
                    ierr = DMStagVecGetArrayDOF(dm, local_x, k, j - 1, i, &below_val); CHKERRQ(ierr);
                }

                // Calculate the derivative and store it in DOWN position
                PetscScalar dy_val = center_val - below_val;

                // Store in DOWN position
                ierr = DMStagVecSetArrayDOF(dm, local_dy, k, j, i, dy_val); CHKERRQ(ierr);

                // Store in UP position for the cell below
                if (j > starty) {
                    ierr = DMStagVecSetArrayDOF(dm, local_dy, k, j - 1, i, dy_val); CHKERRQ(ierr);
                }
            }
        }
    }

    ierr = DMDAVecRestoreArrayDOF(dm, local_x, &arr_x); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(dm, local_dy, &arr_dy); CHKERRQ(ierr);

    // Scatter local to global
    ierr = DMLocalToGlobalBegin(dm, local_dy, INSERT_VALUES, dy); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(dm, local_dy, INSERT_VALUES, dy); CHKERRQ(ierr);

    // Destroy local vectors
    ierr = VecDestroy(&local_x); CHKERRQ(ierr);
    ierr = VecDestroy(&local_dy); CHKERRQ(ierr);

    return 0;
}

int main(int argc, char **argv) {
    PetscErrorCode ierr;
    DM             dm;
    Vec            x, dy;
    DMBoundaryType bx = DM_BOUNDARY_NONE, by = DM_BOUNDARY_NONE, bz = DM_BOUNDARY_NONE;

    ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;

    ierr = DMStagCreate3d(PETSC_COMM_WORLD, bx, by, bz, 10, 10, 10, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 0, 1, 1, 1, DMSTAG_STENCIL_BOX, 1, NULL, NULL, NULL, &dm); CHKERRQ(ierr);
    ierr = DMSetUp(dm); CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(dm, &x); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dm, &dy); CHKERRQ(ierr);

    // Set the value 1 at the center of each cell
    PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez;
    DMStagStencil point;

    ierr = DMStagGetCorners(dm, &startx, &starty, &startz, &nx, &ny, &nz, &ex, &ey, &ez); CHKERRQ(ierr);

    for (PetscInt k = startz; k < startz + nz; ++k) {
        for (PetscInt j = starty; j < starty + ny; ++j) {
            for (PetscInt i = startx; i < startx + nx; ++i) {
                point.i = i;
                point.j = j;
                point.k = k;
                point.loc = DMSTAG_ELEMENT;
                PetscScalar value = 1.0;
                ierr = DMStagVecSetValuesStencil(dm, x, 1, &point, &value, INSERT_VALUES); CHKERRQ(ierr);
            }
        }
    }

    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

    // Compute the derivative
    ierr = ComputeDerivative(dm, x, dy); CHKERRQ(ierr);

    /*// Destroy PETSc objects and finalize
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&dy); CHKERRQ(ierr);
    ierr = DMDestroy(&dm); CHKERRQ(ierr);*/
    ierr = PetscFinalize();

    return ierr;
}
