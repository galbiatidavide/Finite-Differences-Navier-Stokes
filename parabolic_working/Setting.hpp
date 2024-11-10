#ifndef SETTING_HPP
#define SETTING_HPP

#include <functional>
#include <petscdm.h>
#include <petscdmstag.h>

class Setting {
public:
    Setting(PetscInt nx_, PetscInt ny_, PetscInt nz_,
            PetscReal Lx_0_, PetscReal Ly_0_, PetscReal Lz_0_,
            PetscReal Lx_,   PetscReal Ly_,   PetscReal Lz_,
            PetscReal dt_,   PetscReal iter_)
        : nx(nx_), ny(ny_), nz(nz_),
          Lx_0(Lx_0_), Ly_0(Ly_0_), Lz_0(Lz_0_),
          Lx(Lx_),     Ly(Ly_),     Lz(Lz_),
          dt(dt_),     iter(iter_) {}  

    const PetscInt nx, ny, nz;
    const PetscReal Lx_0, Ly_0, Lz_0, Lx, Ly, Lz;
    const PetscReal dt, iter;
};

template <typename ProblemType>
class ProblemSetting : public Setting {
public:
    using Setting::Setting;
};

template <>
class ProblemSetting<struct Stokes> : public Setting {
public:
    using Setting::Setting;
};

template <>
class ProblemSetting<struct NavierStokes> : public Setting {
public:
    using Setting::Setting;
};

template <>
class ProblemSetting<struct Euler> : public Setting {
public:
    using Setting::Setting;
};

template <>
class ProblemSetting<struct Parabolic> : public Setting {
public:
    ProblemSetting(PetscInt nx_, PetscInt ny_, PetscInt nz_,
                   PetscReal Lx_0_, PetscReal Ly_0_, PetscReal Lz_0_,
                   PetscReal Lx_, PetscReal Ly_, PetscReal Lz_,
                   PetscReal dt_, PetscReal iter_, PetscReal mu_)
        : Setting(nx_, ny_, nz_, Lx_0_, Ly_0_, Lz_0_, Lx_, Ly_, Lz_, dt_, iter_),
          mu(mu_) {
        DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                       nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                       0, 0, 1, 0, DMSTAG_STENCIL_BOX, 1, NULL, NULL, NULL, &dmGrid);
        DMSetFromOptions(dmGrid);
        DMSetUp(dmGrid);
        DMStagSetUniformCoordinatesExplicit(dmGrid, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    }

    const PetscReal mu;
    DM dmGrid;

    ~ProblemSetting() {
        DMDestroy(&dmGrid);
    }
};

template <>
class ProblemSetting<struct Poisson> : public Setting {
public:
    using Setting::Setting;
};

#endif // SETTING_HPP
