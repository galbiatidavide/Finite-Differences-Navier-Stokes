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
#include <petscksp.h>


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


bool rayIntersectsTriangle(const std::array<double, 3>& rayOrigin,
                           const std::array<double, 3>& rayVector,
                           const std::array<double, 3>& v0,
                           const std::array<double, 3>& v1,
                           const std::array<double, 3>& v2) {
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


bool isPointInsideMesh(const std::array<double, 3>& point,
                       const std::vector<std::array<double, 3>>& vertices,
                       const std::vector<std::array<int, 3>>& faces) {
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


static PetscErrorCode CreateGrid(DM* pdmSol, PetscInt dof1, PetscInt dof2, PetscInt dof3, PetscScalar nx, PetscScalar ny, PetscScalar nz, PetscReal Lx_0, PetscReal Lx, PetscReal Ly_0, PetscReal Ly, PetscReal Lz_0, PetscReal Lz)
{
    DM dmSol;
    const PetscInt dof0 = 0;
    const PetscInt stencilWidth = 1;
    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, nx, ny, nz, PETSC_DECIDE,
                   PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL,
                   NULL, pdmSol);
    dmSol = *pdmSol;
    DMSetFromOptions(dmSol);
    DMSetUp(dmSol);
    DMStagSetUniformCoordinatesExplicit(dmSol, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    return 0;
}

void createGridSTL(const std::vector<std::array<double, 3>>& vertices, 
                std::unordered_map<std::string, GridCell>& grid, 
                double cellSize, 
                const std::array<double, 3>& gridOrigin) {
    
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


int main(int argc, char* argv[]) {

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <surface_cut.stl>" << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> faces;
    std::string filename = argv[1];
    reader(filename, vertices, faces);

    PetscInitialize(&argc, &argv, (char*)0, "Help text for this program");

    static PetscInt nx_  = 40;
    static PetscInt ny_  = 40;
    static PetscInt nz_  = 40;

    static PetscReal Lx_0   = 0;
    static PetscReal Ly_0  = 0.2;
    static PetscReal Lz_0  = 0;
    static PetscReal Lx = 0.9;
    static PetscReal Ly = 0.9;
    static PetscReal Lz = 0.9;

    DM dmGrid, dm_cent;
    CreateGrid(&dm_cent, 0, 0, 1, nx_, ny_, nz_, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);
    CreateGrid(&dmGrid, 0, 1, 0, nx_, ny_, nz_, Lx_0, Lx, Ly_0, Ly, Lz_0, Lz);

    Vec vec_stag, vec_cent;
    DMCreateGlobalVector(dm_cent, &vec_cent);
    DMCreateGlobalVector(dmGrid, &vec_stag);

    //staggered
    PetscInt icux_right[3], icux_left[3], icux_up[3], icux_down[3], icux_back[3], icux_front[3], iux_right, iux_left, iux_up, iux_down, iux_back, iux_front;
    PetscInt startx, starty, startz, N[3], nx, ny, nz, ex, ey, ez, d;
    DM dmCoord;
    Vec coord, coordLocal, vec_stag_local;
    PetscReal ****arrCoord, ****arrVecStag;   

    DMStagGetCorners(dmGrid, &startx, &starty, &startz, &nx, &ny, &nz, NULL, NULL, NULL);
    DMStagGetGlobalSizes(dmGrid, &N[0], &N[1], &N[2]);
    DMGetCoordinateDM(dmGrid, &dmCoord);

    DMGetCoordinates(dmGrid, &coord);
    DMGetLocalVector(dmCoord, &coordLocal);
    DMGlobalToLocal(dmCoord, coord, INSERT_VALUES, coordLocal);


    for (d = 0; d < 3; ++d) {
        DMStagGetLocationSlot(dmCoord, RIGHT, d, &icux_right[d]);
        DMStagGetLocationSlot(dmCoord, LEFT, d, &icux_left[d]);
        DMStagGetLocationSlot(dmCoord, UP, d, &icux_up[d]);
        DMStagGetLocationSlot(dmCoord, DOWN, d, &icux_down[d]);
        DMStagGetLocationSlot(dmCoord, BACK, d, &icux_back[d]);
        DMStagGetLocationSlot(dmCoord, FRONT, d, &icux_front[d]);
    }  

    DMStagVecGetArrayRead(dmCoord, coordLocal, &arrCoord);

    DMStagGetLocationSlot(dmGrid, RIGHT, 0, &iux_right);
    DMStagGetLocationSlot(dmGrid, LEFT, 0, &iux_left);
    DMStagGetLocationSlot(dmGrid, UP, 0, &iux_up);
    DMStagGetLocationSlot(dmGrid, DOWN, 0, &iux_down);
    DMStagGetLocationSlot(dmGrid, BACK, 0, &iux_back);
    DMStagGetLocationSlot(dmGrid, FRONT, 0, &iux_front);

    DMCreateLocalVector(dmGrid, &vec_stag_local);
    DMGlobalToLocalBegin(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMGlobalToLocalEnd(dmGrid, vec_stag, INSERT_VALUES, vec_stag_local);
    DMStagVecGetArray(dmGrid, vec_stag_local, &arrVecStag);

    //LEFT, DOWN, BACK se sono alla fine anche destra sopra e davanti 

    std::vector<std::vector<PetscInt>> icux;
    icux.push_back({icux_right[0], icux_right[1], icux_right[2]});
    icux.push_back({icux_left[0], icux_left[1], icux_left[2]});
    icux.push_back({icux_up[0], icux_up[1], icux_up[2]});
    icux.push_back({icux_down[0], icux_down[1], icux_down[2]});
    icux.push_back({icux_back[0], icux_back[1], icux_back[2]});
    icux.push_back({icux_front[0], icux_front[1], icux_front[2]});


    std::vector<PetscInt> iux = {iux_right, iux_left, iux_up, iux_down, iux_back, iux_front};

    for (ez = startz; ez < startz + nz; ++ez) {
        for (ey = starty; ey < starty + ny; ++ey) {
            for (ex = startx; ex < startx + nx; ++ex) {
                for(auto i : icux) {
                    for (auto j : iux){
                    std::array<double, 3> point = {arrCoord[ez][ey][ex][i[0]], arrCoord[ez][ey][ex][i[1]], arrCoord[ez][ey][ex][i[2]]};
                    if (isPointInsideMesh(point, vertices, faces)) {
                        arrVecStag[ez][ey][ex][j] = 1.0;
                    }
                    else {
                        arrVecStag[ez][ey][ex][j] = 1000.0;
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

    PetscViewer viewer_stag;
    DM dmGrid;
    //DMStagCreateCompatibleDMStag(dmSol_Staggered_y, 0, 0, 1, 0, &DM_v);
    Vec rd;
    DMStagVecSplitToDMDA(dmGrid, vec_stag, DOWN, 0, &dmGrid, &rd);
    PetscObjectSetName((PetscObject)rd, "r_values_down");
    char filename_rd[50];
    sprintf(filename_rd, "results/r_values_down%03zu.vtr");
    PetscViewerVTKOpen(PetscObjectComm((PetscObject)dmGrid), filename_rd, FILE_MODE_WRITE, &viewer_stag);
    VecView(rd, viewer_stag);
    VecDestroy(&rd);
    DMDestroy(&dmGrid);
    PetscViewerDestroy(&viewer_stag);
    

    PetscFinalize();

    return EXIT_SUCCESS;
}