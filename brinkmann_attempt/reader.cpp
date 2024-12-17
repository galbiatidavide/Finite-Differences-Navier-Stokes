#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPoints.h>
#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <petscdm.h>
#include <petscdmstag.h>

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

struct GridCell {
    std::vector<std::array<double, 3>> vertices;
};


bool rayIntersectsTriangle(const std::array<double, 3>& rayOrigin, const std::array<double, 3>& rayVector, const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2) {
    const double EPSILON = 1e-8;
    std::array<double, 3> edge1, edge2, h, s, q;
    double a, f, u, v;

    // Calculate edges
    for (int i = 0; i < 3; ++i) {
        edge1[i] = v1[i] - v0[i];
        edge2[i] = v2[i] - v0[i];
    }

    // Calculate determinant
    h[0] = rayVector[1] * edge2[2] - rayVector[2] * edge2[1];
    h[1] = rayVector[2] * edge2[0] - rayVector[0] * edge2[2];
    h[2] = rayVector[0] * edge2[1] - rayVector[1] * edge2[0];

    a = edge1[0] * h[0] + edge1[1] * h[1] + edge1[2] * h[2];

    if (a > -EPSILON && a < EPSILON)
        return false; // This means the ray is parallel to the triangle.

    f = 1.0 / a;

    // Calculate u parameter and test bound
    for (int i = 0; i < 3; ++i) {
        s[i] = rayOrigin[i] - v0[i];
    }

    u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
    if (u < 0.0 || u > 1.0)
        return false;

    // Calculate v parameter and test bound
    q[0] = s[1] * edge1[2] - s[2] * edge1[1];
    q[1] = s[2] * edge1[0] - s[0] * edge1[2];
    q[2] = s[0] * edge1[1] - s[1] * edge1[0];

    v = f * (rayVector[0] * q[0] + rayVector[1] * q[1] + rayVector[2] * q[2]);
    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage, we can compute t to find out where the intersection point is on the line.
    double t = f * (edge2[0] * q[0] + edge2[1] * q[1] + edge2[2] * q[2]);
    if (t > EPSILON) // ray intersection
        return true;

    return false; // No intersection
}

bool isPointInsideMesh(const std::array<double, 3>& point, const std::vector<std::array<double, 3>>& vertices, const std::vector<std::array<int, 3>>& faces) {
    // Define a ray direction. We'll use the positive x-direction.
    std::array<double, 3> rayDirection = {1.0, 0.0, 0.0};

    int intersectionCount = 0;

    // Iterate over all faces (triangles)
    for (const auto& face : faces) {
        const std::array<double, 3>& v0 = vertices[face[0]];
        const std::array<double, 3>& v1 = vertices[face[1]];
        const std::array<double, 3>& v2 = vertices[face[2]];

        if (rayIntersectsTriangle(point, rayDirection, v0, v1, v2)) {
            intersectionCount++;
        }
    }

    // If the intersection count is odd, the point is inside; otherwise, it's outside.
    return (intersectionCount % 2 == 1);
}

void createGridSTL(const std::vector<std::array<double, 3>>& vertices, std::unordered_map<std::string, GridCell>& grid, double cellSize, const std::array<double, 3>& gridOrigin)
{
    
    for (const auto& vertex : vertices) {
        int x_idx = static_cast<int>((vertex[0] - gridOrigin[0]) / cellSize);
        int y_idx = static_cast<int>((vertex[1] - gridOrigin[1]) / cellSize);
        int z_idx = static_cast<int>((vertex[2] - gridOrigin[2]) / cellSize);

        std::string key = std::to_string(x_idx) + "_" + std::to_string(y_idx) + "_" + std::to_string(z_idx);
        grid[key].vertices.push_back(vertex);
    }
}

void saveGridToFile(const std::unordered_map<std::string, GridCell>& grid, const std::string& filename) {
    std::ofstream gridFile(filename);
    if (gridFile.is_open()) {
        for (const auto& cell : grid) {
            gridFile << "Cell " << cell.first << ":\n";
            for (const auto& v : cell.second.vertices) {
                gridFile << "    " << v[0] << " " << v[1] << " " << v[2] << std::endl;
            }
            gridFile << std::endl;
        }
        gridFile.close();
        std::cout << "Grid information saved to " << filename << std::endl;
    } else {
        std::cerr << "Unable to open file " << filename << " for writing" << std::endl;
    }
}

void reader(const std::string& filename, std::vector<std::array<double, 3>>& vertices, std::vector<std::array<int, 3>>& faces) {
    // Create a reader for the STL file
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // Get the polydata object
    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

    // Extract vertices
    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
    vertices.resize(points->GetNumberOfPoints());
    //std::vector<std::array<double, 3>> vertices(points->GetNumberOfPoints());

    

    double x_min = std::numeric_limits<double>::max();
    double y_min = std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();

    double x_max = std::numeric_limits<double>::lowest();
    double y_max = std::numeric_limits<double>::lowest();
    double z_max = std::numeric_limits<double>::lowest();

    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        double p[3];
        points->GetPoint(i, p);
        vertices[i] = {{ p[0], p[1], p[2] }};

        x_min = std::min(x_min, p[0]);
        y_min = std::min(y_min, p[1]);
        z_min = std::min(z_min, p[2]);

        x_max = std::max(x_max, p[0]);
        y_max = std::max(y_max, p[1]);
        z_max = std::max(z_max, p[2]);
    }

    // Calculate the maximum dimension of the bounding box
    double max_dimension = std::max({x_max - x_min, y_max - y_min, z_max - z_min});

    // Normalize vertices using the maximum dimension
    for (auto& v : vertices) {
        v[0] = (v[0] - x_min) / max_dimension;
        v[1] = (v[1] - y_min) / max_dimension;
        v[2] = (v[2] - z_min) / max_dimension;
    }


    // Define the grid resolution (e.g., 10x10x10 grid)
    int gridResolution = 10;
    double cellSize = 1.0 / gridResolution; // Since vertices are normalized, the grid spans [0, 1]
    std::array<double, 3> gridOrigin = {0.0, 0.0, 0.0}; // Grid starts at the origin

    // Create a grid and store vertices in the corresponding cells
    std::unordered_map<std::string, GridCell> grid;
    createGridSTL(vertices, grid, cellSize, gridOrigin);

    // Save grid information to a file
    saveGridToFile(grid, "grid_info.txt");

    // Save vertices to a file (optional)
    std::ofstream verticesFile("vertices.txt");
    if (verticesFile.is_open()) {
        for (const auto& v : vertices) {
            verticesFile << v[0] << " " << v[1] << " " << v[2] << std::endl;
        }
        verticesFile.close();
        std::cout << "Vertices saved to vertices.txt" << std::endl;
    } else {
        std::cerr << "Unable to open file vertices.txt for writing" << std::endl;
    }

    // Save faces to a file (optional)
    vtkSmartPointer<vtkCellArray> triangles = polyData->GetPolys();
    //std::vector<std::array<vtkIdType, 3>> faces;

    triangles->InitTraversal();
    vtkIdType npts;
    const vtkIdType *pts;
    while (triangles->GetNextCell(npts, pts)) {
        if (npts == 3) {
            faces.push_back({{ pts[0], pts[1], pts[2] }});
        }
    }


    std::ofstream facesFile("faces.txt");
    if (facesFile.is_open()) {
        for (const auto& f : faces) {
            facesFile << f[0] << " " << f[1] << " " << f[2] << std::endl;
        }
        facesFile.close();
        std::cout << "Faces saved to faces.txt" << std::endl;
    } else {
        std::cerr << "Unable to open file faces.txt for writing" << std::endl;
    }

}

PetscErrorCode createMaskU(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces) {

    PetscInt icux_right[3], icux_left[3], iux_right, iux_left;
    PetscInt startx, starty, startz, N[3], ex, ey, ez, nx, ny, nz, d;
    DM dmCoord;
    Vec coord, coordLocal, vec_stag_local;
    PetscReal ****arrCoord, ****arrVecStag;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
    }  

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);

    DMCreateLocalVector(dmGrid, &vec_stag_local);
    DMGlobalToLocalBegin(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMGlobalToLocalEnd(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMStagVecGetArray(dmGrid, vec_stag_local, &arrVecStag);

    std::vector<std::vector<PetscInt>> icux;
    icux.push_back({icux_right[0], icux_right[1], icux_right[2]});
    icux.push_back({icux_left[0], icux_left[1], icux_left[2]});

    std::vector<PetscInt> iux = {iux_right, iux_left};

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                for(auto i : icux) {
                    for (auto j : iux){
                    std::array<double, 3> point = {arrCoord[ez][ey][ex][i[0]], arrCoord[ez][ey][ex][i[1]], arrCoord[ez][ey][ex][i[2]]};
                    if (isPointInsideMesh(point, vertices, faces)) {
                        arrVecStag[ez][ey][ex][j] = 0.0;
                    }
                    else {
                        arrVecStag[ez][ey][ex][j] = eps;
                        }
                    }
                }

            }
        }
    }
    
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArray(dmGrid, vec_stag_local, &arrVecStag);
    DMLocalToGlobal(dmGrid, vec_stag_local, INSERT_VALUES, vec_stag);
    DMRestoreLocalVector(dmGrid, &vec_stag_local);

    PetscFunctionReturn(0);
}

PetscErrorCode createMaskV(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces) {

    PetscInt icuy_up[3], icuy_down[3], iuy_up, iuy_down;
    PetscInt startx, starty, startz, N[3], ex, ey, ez, nx, ny, nz, d;
    DM dmCoord;
    Vec coord, coordLocal, vec_stag_local;
    PetscReal ****arrCoord, ****arrVecStag;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, UP, d, &icuy_up[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icuy_down[d]);
    }  

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, UP, 0, &iuy_up);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iuy_down);

    DMCreateLocalVector(dmGrid, &vec_stag_local);
    DMGlobalToLocalBegin(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMGlobalToLocalEnd(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMStagVecGetArray(dmGrid, vec_stag_local, &arrVecStag);

    std::vector<std::vector<PetscInt>> icuy;
    icuy.push_back({icuy_up[0], icuy_up[1], icuy_up[2]});
    icuy.push_back({icuy_down[0], icuy_down[1], icuy_down[2]});

    std::vector<PetscInt> iuy = {iuy_up, iuy_down};

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                for(auto i : icuy) {
                    for (auto j : iuy){
                    std::array<double, 3> point = {arrCoord[ez][ey][ex][i[0]], arrCoord[ez][ey][ex][i[1]], arrCoord[ez][ey][ex][i[2]]};
                    if (isPointInsideMesh(point, vertices, faces)) {
                        arrVecStag[ez][ey][ex][j] = 0.0;
                    }
                    else {
                        arrVecStag[ez][ey][ex][j] = eps;
                        }
                    }
                }

            }
        }
    }
    
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArray(dmGrid, vec_stag_local, &arrVecStag);
    DMLocalToGlobal(dmGrid, vec_stag_local, INSERT_VALUES, vec_stag);
    DMRestoreLocalVector(dmGrid, &vec_stag_local);

    PetscFunctionReturn(0);
}

PetscErrorCode createMaskW(DM const & dmGrid, Vec & vec_stag, std::vector<std::array<double, 3>> const & vertices, std::vector<std::array<int, 3>> const & faces) {

    PetscInt icuz_front[3], icuz_back[3], iuz_front, iuz_back;
    PetscInt startx, starty, startz, N[3], ex, ey, ez, nx, ny, nz, d;
    DM dmCoord;
    Vec coord, coordLocal, vec_stag_local;
    PetscReal ****arrCoord, ****arrVecStag;   

    PetscFunctionBegin;

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);

    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icuz_front[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icuz_back[d]);
    }  

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iuz_front);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iuz_back);

    DMCreateLocalVector(dmGrid, &vec_stag_local);
    DMGlobalToLocalBegin(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMGlobalToLocalEnd(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMStagVecGetArray(dmGrid, vec_stag_local, &arrVecStag);

    std::vector<std::vector<PetscInt>> icuz;
    icuz.push_back({icuz_front[0], icuz_front[1], icuz_front[2]});
    icuz.push_back({icuz_back[0], icuz_back[1], icuz_back[2]});

    std::vector<PetscInt> iuz = {iuz_front, iuz_back};

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                for(auto i : icuz) {
                    for (auto j : iuz){
                    std::array<double, 3> point = {arrCoord[ez][ey][ex][i[0]], arrCoord[ez][ey][ex][i[1]], arrCoord[ez][ey][ex][i[2]]};
                    if (isPointInsideMesh(point, vertices, faces)) {
                        arrVecStag[ez][ey][ex][j] = 0.0;
                    }
                    else {
                        arrVecStag[ez][ey][ex][j] = eps;
                        }
                    }
                }

            }
        }
    }
    
    DMStagVecRestoreArrayRead(dmCoord, coordLocal, &arrCoord);
    DMRestoreLocalVector(dmCoord, &coordLocal);
    DMStagVecRestoreArray(dmGrid, vec_stag_local, &arrVecStag);
    DMLocalToGlobal(dmGrid, vec_stag_local, INSERT_VALUES, vec_stag);
    DMRestoreLocalVector(dmGrid, &vec_stag_local);

    PetscFunctionReturn(0);
}





    


    
