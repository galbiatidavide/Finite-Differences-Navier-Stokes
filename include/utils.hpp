#include "macros.hpp"
#include "config_problem.hpp"

using namespace problem_setting;

#ifndef UTILS_HPP
#define UTILS_HPP

PetscErrorCode CheckSolution(Vec const & sol, Vec const & solRef, std::string const & comp);

PetscErrorCode CreateAnalyticalU(DM const & dmGrid, Vec & vec, PetscReal const & theta);

PetscErrorCode CreateAnalyticalV(DM const & dmGrid, Vec & vec, PetscReal const & theta);

PetscErrorCode CreateAnalyticalW(DM const & dmGrid, Vec & vec, PetscReal const & theta);

//PetscErrorCode CreateAnalyticalP(DM const & dmGrid, Vec & vec, PetscReal const & theta);

PetscErrorCode CreateGrid(DM * const dmGrid, PetscInt const & dof1, PetscInt const & dof2, PetscInt const & dof3, PetscInt const & nx, PetscInt const & ny, PetscInt const & nz, PetscReal const & Lx_0, PetscReal const & Lx, PetscReal const & Ly_0, PetscReal const & Ly, PetscReal const & Lz_0, PetscReal const & Lz);

void PrintSimulationParameters();

bool rayIntersectsTriangle(const std::array<double, 3>& rayOrigin, const std::array<double, 3>& rayVector, const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2);

bool isPointInsideMesh(const std::array<double, 3>& point, const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<int, 3>>& faces);

void reader(const std::string& filename, std::vector<std::array<double, 3>>& vertices, std::vector<std::array<int, 3>>& faces);

PetscErrorCode createMaskU(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces);

PetscErrorCode createMaskV(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces);

PetscErrorCode createMaskW(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces);

PetscReal constexpr D_x{Lx - Lx_0};
PetscReal constexpr D_y{Ly - Ly_0};
PetscReal constexpr D_z{Lz - Lz_0};

//mask
inline std::vector<std::array<double, 3>> vertices;
inline std::vector<std::array<int, 3>> faces;
inline std::string filename;

#endif // UTILS_HPP